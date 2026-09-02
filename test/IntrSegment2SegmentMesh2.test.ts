import { describe, it, expect } from 'vitest';
import { IntrLine2SegmentMesh2FI } from '../src/IntrLine2SegmentMesh2';
import {
    IntrSegment2SegmentMesh2FI,
    defaultIntrSegment2SegmentMesh2FIResult
} from '../src/IntrSegment2SegmentMesh2';
import { Line } from '../src/Line';
import { Segment } from '../src/Segment';
import { SegmentMesh } from '../src/SegmentMesh';
import { Vector, add, length, mul, sub } from '../src/Vector';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

// A unit square as a closed contiguous mesh.
function squareMesh(): SegmentMesh {
    return SegmentMesh.fromContiguous(
        [vec(0, 0), vec(1, 0), vec(1, 1), vec(0, 1)], false);
}

const fi = new IntrSegment2SegmentMesh2FI();

describe('IntrSegment2SegmentMesh2', () => {
    it('has an empty default result', () => {
        expect(defaultIntrSegment2SegmentMesh2FIResult().intersections)
            .toHaveLength(0);
    });

    it('finds both crossings of a square and reports [0,1] parameters', () => {
        // The segment spans x in [-1, 2] at y = 0.5, so the crossings at
        // x = 0 and x = 1 are at parameters 1/3 and 2/3.
        const result = fi.find(seg([-1, 0.5], [2, 0.5]), squareMesh());
        expect(result.intersections).toHaveLength(2);
        expect(result.intersections[0].segmentParameter).toBeCloseTo(1 / 3, 12);
        expect(result.intersections[1].segmentParameter).toBeCloseTo(2 / 3, 12);
        expect(result.intersections[0].point.values[0]).toBeCloseTo(0, 12);
        expect(result.intersections[1].point.values[0]).toBeCloseTo(1, 12);
    });

    it('discards crossings beyond the segment endpoints', () => {
        // The segment stops inside the square, so only the x = 0 edge is hit.
        const result = fi.find(seg([-1, 0.5], [0.5, 0.5]), squareMesh());
        expect(result.intersections).toHaveLength(1);
        expect(result.intersections[0].segmentParameter).toBeCloseTo(2 / 3, 12);
        // A segment fully inside the square hits nothing.
        expect(fi.find(seg([0.2, 0.5], [0.8, 0.5]), squareMesh()).intersections)
            .toHaveLength(0);
        // A segment fully outside on the line of the square hits nothing.
        expect(fi.find(seg([-3, 0.5], [-2, 0.5]), squareMesh()).intersections)
            .toHaveLength(0);
    });

    it('keeps crossings exactly at the segment endpoints', () => {
        const result = fi.find(seg([0, 0.5], [1, 0.5]), squareMesh());
        expect(result.intersections).toHaveLength(2);
        expect(result.intersections[0].segmentParameter).toBeCloseTo(0, 12);
        expect(result.intersections[1].segmentParameter).toBeCloseTo(1, 12);
    });

    it('reports nothing for a segment that misses the mesh', () => {
        expect(fi.find(seg([-1, 5], [2, 5]), squareMesh()).intersections)
            .toHaveLength(0);
    });

    it('matches the line query filtered by parameter on random inputs', () => {
        let state = 909090;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        const lsQuery = new IntrLine2SegmentMesh2FI();
        const vertices: Vector[] = [];
        for (let i = 0; i < 10; ++i) {
            vertices.push(vec(i - 5, (i % 2 === 0 ? -1 : 1) * 1.5));
        }
        const mesh = SegmentMesh.fromContiguous(vertices, true);

        let numHits = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const p0 = [rand() * 6, rand() * 3];
            const p1 = [rand() * 6, rand() * 3];
            if (p0[0] === p1[0] && p0[1] === p1[1]) {
                continue;
            }
            const S = seg(p0, p1);
            const result = fi.find(S, mesh);

            const direction = sub(S.p[1], S.p[0]);
            const lineResult = lsQuery.find(
                Line.fromOriginDirection(S.p[0], direction), mesh);
            const expected = lineResult.intersections.filter(
                (object) => object.lineParameter >= 0
                    && object.lineParameter <= 1);
            expect(result.intersections).toHaveLength(expected.length);
            numHits += result.intersections.length;

            for (let i = 0; i < expected.length; ++i) {
                const actual = result.intersections[i];
                expect(actual.indexPair).toEqual(expected[i].indexPair);
                expect(actual.segmentParameter)
                    .toBeCloseTo(expected[i].lineParameter, 12);
                expect(actual.segmentParameter).toBeGreaterThanOrEqual(0);
                expect(actual.segmentParameter).toBeLessThanOrEqual(1);
                // The point is on the segment and on the named mesh segment.
                const onSegment = add(S.p[0],
                    mul(actual.segmentParameter, direction));
                expect(length(sub(actual.point, onSegment)))
                    .toBeCloseTo(0, 10);
                const V0 = mesh.getVertices()[actual.indexPair[0]];
                const V1 = mesh.getVertices()[actual.indexPair[1]];
                const onMeshSegment = add(
                    mul(1 - actual.meshSegmentParameter, V0),
                    mul(actual.meshSegmentParameter, V1));
                expect(length(sub(actual.point, onMeshSegment)))
                    .toBeCloseTo(0, 10);
            }
        }
        expect(numHits).toBeGreaterThan(50);
    });
});
