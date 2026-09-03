import { describe, it, expect } from 'vitest';
import { IntrLine2SegmentMesh2FI } from '../src/IntrLine2SegmentMesh2.js';
import {
    IntrRay2SegmentMesh2FI,
    defaultIntrRay2SegmentMesh2FIResult
} from '../src/IntrRay2SegmentMesh2.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { SegmentMesh } from '../src/SegmentMesh.js';
import { Vector, add, length, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

// A unit square as a closed contiguous mesh with vertices
// (0,0), (1,0), (1,1), (0,1).
function squareMesh(): SegmentMesh {
    return SegmentMesh.fromContiguous(
        [vec(0, 0), vec(1, 0), vec(1, 1), vec(0, 1)], false);
}

const fi = new IntrRay2SegmentMesh2FI();

describe('IntrRay2SegmentMesh2', () => {
    it('has an empty default result', () => {
        expect(defaultIntrRay2SegmentMesh2FIResult().intersections)
            .toHaveLength(0);
    });

    it('finds both crossings of a square with a ray outside it', () => {
        const result = fi.find(ray([-1, 0.5], [1, 0]), squareMesh());
        expect(result.intersections).toHaveLength(2);
        // Sorted by the line parameter, so the near crossing is first.
        expect(result.intersections[0].rayParameter).toBeCloseTo(1, 12);
        expect(result.intersections[1].rayParameter).toBeCloseTo(2, 12);
        expect(result.intersections[0].point.values[0]).toBeCloseTo(0, 12);
        expect(result.intersections[1].point.values[0]).toBeCloseTo(1, 12);
        // The mesh segment parameter interpolates the named endpoints.
        for (const hit of result.intersections) {
            const mesh = squareMesh();
            const V0 = mesh.getVertices()[hit.indexPair[0]];
            const V1 = mesh.getVertices()[hit.indexPair[1]];
            const P = add(mul(1 - hit.meshSegmentParameter, V0),
                mul(hit.meshSegmentParameter, V1));
            expect(length(sub(P, hit.point))).toBeCloseTo(0, 12);
        }
    });

    it('discards the crossings behind the ray origin', () => {
        // The ray starts inside the square and points in +x, so only the
        // x = 1 edge is hit.
        const result = fi.find(ray([0.5, 0.5], [1, 0]), squareMesh());
        expect(result.intersections).toHaveLength(1);
        expect(result.intersections[0].rayParameter).toBeCloseTo(0.5, 12);
        expect(result.intersections[0].point.values[0]).toBeCloseTo(1, 12);
        // The same ray reversed hits the x = 0 edge.
        const reversed = fi.find(ray([0.5, 0.5], [-1, 0]), squareMesh());
        expect(reversed.intersections).toHaveLength(1);
        expect(reversed.intersections[0].point.values[0]).toBeCloseTo(0, 12);
    });

    it('keeps a crossing exactly at the ray origin', () => {
        // The ray origin is on the x = 0 edge; the line parameter there is 0,
        // which is not discarded.
        const result = fi.find(ray([0, 0.5], [1, 0]), squareMesh());
        expect(result.intersections).toHaveLength(2);
        expect(result.intersections[0].rayParameter).toBeCloseTo(0, 12);
    });

    it('reports nothing for a ray that misses the mesh', () => {
        expect(fi.find(ray([-1, 5], [1, 0]), squareMesh()).intersections)
            .toHaveLength(0);
        expect(fi.find(ray([-1, 0.5], [-1, 0]), squareMesh()).intersections)
            .toHaveLength(0);
    });

    it('matches the line query filtered by parameter on random inputs', () => {
        let state = 60660;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        const lsQuery = new IntrLine2SegmentMesh2FI();
        // A zigzag open polyline.
        const vertices: Vector[] = [];
        for (let i = 0; i < 10; ++i) {
            vertices.push(vec(i - 5, (i % 2 === 0 ? -1 : 1) * 1.5));
        }
        const mesh = SegmentMesh.fromContiguous(vertices, true);

        let numHits = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const R = ray([rand() * 6, rand() * 3], [rand(), rand() + 0.001]);
            const result = fi.find(R, mesh);
            const lineResult = lsQuery.find(
                Line.fromOriginDirection(R.origin, R.direction), mesh);
            const expected = lineResult.intersections.filter(
                (object) => object.lineParameter >= 0);
            expect(result.intersections).toHaveLength(expected.length);
            numHits += result.intersections.length;

            for (let i = 0; i < expected.length; ++i) {
                const actual = result.intersections[i];
                expect(actual.indexPair).toEqual(expected[i].indexPair);
                expect(actual.rayParameter)
                    .toBeCloseTo(expected[i].lineParameter, 12);
                expect(actual.meshSegmentParameter)
                    .toBeCloseTo(expected[i].meshSegmentParameter, 12);
                // The point is on the ray and on the named mesh segment.
                expect(actual.rayParameter).toBeGreaterThanOrEqual(0);
                const onRay = add(R.origin,
                    mul(actual.rayParameter, R.direction));
                expect(length(sub(actual.point, onRay))).toBeCloseTo(0, 10);
                const V0 = mesh.getVertices()[actual.indexPair[0]];
                const V1 = mesh.getVertices()[actual.indexPair[1]];
                const onSegment = add(mul(1 - actual.meshSegmentParameter, V0),
                    mul(actual.meshSegmentParameter, V1));
                expect(length(sub(actual.point, onSegment)))
                    .toBeCloseTo(0, 10);
            }
        }
        expect(numHits).toBeGreaterThan(100);
    });
});
