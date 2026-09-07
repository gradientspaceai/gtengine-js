import { describe, it, expect } from 'vitest';
import { IntrLine2SegmentMesh2FI } from '../src/IntrLine2SegmentMesh2.js';
import {
    IntrRay2SegmentMesh2FI,
    defaultIntrRay2SegmentMesh2FIResult
} from '../src/IntrRay2SegmentMesh2.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { SegmentMesh } from '../src/SegmentMesh.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, ray as arbRay,
    rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrRay2SegmentMesh2.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrRay2SegmentMesh2 verification', () => {
    const fiQ = new IntrRay2SegmentMesh2FI();
    const lsQ = new IntrLine2SegmentMesh2FI();

    // A closed polygonal mesh whose vertices are on a jittered circle, so
    // that segments are well separated and rays hit several of them.
    const arbMesh = fc.tuple(fc.integer({ min: 3, max: 9 }),
        fc.array(fc.double({ min: 0.5, max: 2, noNaN: true }),
            { minLength: 9, maxLength: 9 }))
        .map(([n, radii]) => {
            const vs: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                const a = (2 * Math.PI * i) / n;
                const r = radii[i];
                vs.push(Vector.fromArray([r * Math.cos(a), r * Math.sin(a)]));
            }
            return SegmentMesh.fromContiguous(vs, false);
        });

    it('is the line-mesh query filtered to nonnegative parameters', () => {
        check(fc.tuple(arbRay(2), arbMesh), ([R, M]) => {
            const line = Line.fromOriginDirection(R.origin, R.direction);
            const ls = lsQ.find(line, M);
            const kept = ls.intersections.filter(o => o.lineParameter >= 0);
            const r = fiQ.find(R, M);
            expect(r.intersections.length).toBe(kept.length);
            for (let i = 0; i < kept.length; ++i) {
                expect(r.intersections[i].rayParameter)
                    .toBe(kept[i].lineParameter);
                expect(r.intersections[i].meshSegmentParameter)
                    .toBe(kept[i].meshSegmentParameter);
                expect(r.intersections[i].indexPair)
                    .toEqual(kept[i].indexPair);
                expect(r.intersections[i].point.values)
                    .toEqual(kept[i].point.values);
            }
        });
    });

    it('the reported points are on the ray and on the mesh segment', () => {
        check(fc.tuple(arbRay(2), arbMesh), ([R, M]) => {
            const vertices = M.getVertices();
            const r = fiQ.find(R, M);
            for (const o of r.intersections) {
                expect(o.rayParameter).toBeGreaterThanOrEqual(0);
                expect(o.meshSegmentParameter).toBeGreaterThanOrEqual(0);
                expect(o.meshSegmentParameter).toBeLessThanOrEqual(1);
                const scale = 1 + Math.abs(o.rayParameter)
                    + length(R.origin);
                expectVectorClose(o.point,
                    add(R.origin, mul(o.rayParameter, R.direction)),
                    1e-9 * scale, 1e-9);
                const V0 = vertices[o.indexPair[0]];
                const V1 = vertices[o.indexPair[1]];
                const t = o.meshSegmentParameter;
                expectVectorClose(o.point,
                    add(mul(1 - t, V0), mul(t, V1)), 1e-9 * scale, 1e-9);
            }
        });
    });

    it('reports the intersections in increasing ray parameter', () => {
        check(fc.tuple(arbRay(2), arbMesh), ([R, M]) => {
            const r = fiQ.find(R, M);
            for (let i = 1; i < r.intersections.length; ++i) {
                expect(r.intersections[i].rayParameter)
                    .toBeGreaterThanOrEqual(
                        r.intersections[i - 1].rayParameter);
            }
        });
    });

    it('a ray from inside a closed mesh crosses it an odd number of times', () => {
        check(fc.tuple(arbMesh, unitVector(2)), ([M, d]) => {
            // The meshes above are star shaped about the origin, so a ray
            // from the origin leaves the polygon exactly once. Skip the rare
            // draws that pass through a vertex.
            const R = Ray.fromOriginDirection(Vector.zero(2), d);
            const r = fiQ.find(R, M);
            for (const o of r.intersections) {
                if (o.meshSegmentParameter < 1e-9
                    || o.meshSegmentParameter > 1 - 1e-9) {
                    return;
                }
            }
            expect(r.intersections.length).toBe(1);
        });
    });

    it('does not alias the input mesh or the line-query result', () => {
        const M = SegmentMesh.fromContiguous([
            Vector.fromArray([-1, -1]), Vector.fromArray([1, -1]),
            Vector.fromArray([1, 1]), Vector.fromArray([-1, 1])], false);
        const R = Ray.fromOriginDirection(Vector.fromArray([-5, 0]),
            Vector.fromArray([1, 0]));
        const r = fiQ.find(R, M);
        expect(r.intersections.length).toBe(2);
        const before = M.getVertices().map(v => v.values.slice());
        r.intersections[0].point.values[0] = 999;
        r.intersections[0].indexPair[0] = 7;
        expect(M.getVertices().map(v => v.values.slice())).toEqual(before);
        // A second call is unaffected by the mutation of the first result.
        const r2 = fiQ.find(R, M);
        expect(r2.intersections[0].point.values[0]).toBeCloseTo(-1, 12);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbRay(2), arbMesh, rotationFrame(2),
            wellScaledVector(2)),
            ([R, M, Rot, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(Rot[0], v), dot(Rot[1], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const r0 = fiQ.find(R, M);
                // Skip near-vertex and near-origin hits, where the strict
                // filters can flip under a rigid motion.
                for (const o of r0.intersections) {
                    if (o.meshSegmentParameter < 1e-6
                        || o.meshSegmentParameter > 1 - 1e-6
                        || o.rayParameter < 1e-6) {
                        return;
                    }
                }
                const M2 = SegmentMesh.fromContiguous(
                    M.getVertices().map(map), false);
                const R2 = Ray.fromOriginDirection(map(R.origin),
                    rot(R.direction));
                const r1 = fiQ.find(R2, M2);
                expect(r1.intersections.length).toBe(r0.intersections.length);
                for (let i = 0; i < r0.intersections.length; ++i) {
                    expectClose(r1.intersections[i].rayParameter,
                        r0.intersections[i].rayParameter, 1e-7, 1e-7);
                }
            });
    });

    it('reports nothing for a ray pointing away from the mesh', () => {
        check(fc.tuple(arbMesh, unitVector(2)), ([M, d]) => {
            const R = Ray.fromOriginDirection(mul(20, d), d);
            expect(fiQ.find(R, M).intersections.length).toBe(0);
        });
    });

    it('has an empty default result that is not shared', () => {
        const a = defaultIntrRay2SegmentMesh2FIResult();
        const b = defaultIntrRay2SegmentMesh2FIResult();
        expect(a.intersections).toEqual([]);
        expect(a.intersections).not.toBe(b.intersections);
    });
});
