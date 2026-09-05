import { describe, it, expect } from 'vitest';
import { IntrLine2SegmentMesh2FI } from '../src/IntrLine2SegmentMesh2.js';
import { Line } from '../src/Line.js';
import { SegmentMesh } from '../src/SegmentMesh.js';
import { Vector, add, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(origin: number[], direction: number[]): Line {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Line.fromOriginDirection(Vector.fromArray(origin), d);
}

// The unit square as a closed polyline.
function unitSquare(): SegmentMesh {
    return SegmentMesh.fromContiguous(
        [vec(0, 0), vec(1, 0), vec(1, 1), vec(0, 1)], false);
}

describe('IntrLine2SegmentMesh2FI', () => {
    const query = new IntrLine2SegmentMesh2FI();

    it('finds the two crossings of a horizontal line through a square', () => {
        const mesh = unitSquare();
        const result = query.find(line([-5, 0.5], [1, 0]), mesh);
        expect(result.intersections.length).toBe(2);

        // Sorted by line parameter, so the left crossing comes first.
        expect(result.intersections[0].point.values[0]).toBeCloseTo(0, 12);
        expect(result.intersections[0].point.values[1]).toBeCloseTo(0.5, 12);
        expect(result.intersections[1].point.values[0]).toBeCloseTo(1, 12);
        expect(result.intersections[1].point.values[1]).toBeCloseTo(0.5, 12);
        expect(result.intersections[0].lineParameter)
            .toBeLessThan(result.intersections[1].lineParameter);

        // The left crossing is on the segment <V[3],V[0]>, the right one is
        // on <V[1],V[2]>.
        expect(result.intersections[0].indexPair).toEqual([3, 0]);
        expect(result.intersections[1].indexPair).toEqual([1, 2]);
    });

    it('reports the mesh-segment parameter of each crossing', () => {
        const mesh = SegmentMesh.fromDisjoint([vec(0, -1), vec(0, 3)]);
        const result = query.find(line([-1, 0], [1, 0]), mesh);
        expect(result.intersections.length).toBe(1);
        // The crossing is at (0,0), which is 1/4 of the way along the
        // segment from (0,-1) to (0,3).
        expect(result.intersections[0].meshSegmentParameter)
            .toBeCloseTo(0.25, 12);
        expect(result.intersections[0].lineParameter).toBeCloseTo(1, 12);
    });

    it('reports both endpoints when the line and a mesh segment are coincident', () => {
        const mesh = SegmentMesh.fromDisjoint([vec(2, 0), vec(5, 0)]);
        const result = query.find(line([0, 0], [1, 0]), mesh);
        expect(result.intersections.length).toBe(2);
        expect(result.intersections[0].point.values[0]).toBeCloseTo(2, 12);
        expect(result.intersections[1].point.values[0]).toBeCloseTo(5, 12);
        expect(result.intersections[0].lineParameter)
            .toBeLessThan(result.intersections[1].lineParameter);
    });

    it('returns an empty list when the line misses every segment', () => {
        const mesh = unitSquare();
        const result = query.find(line([-5, 5], [1, 0]), mesh);
        expect(result.intersections.length).toBe(0);
    });

    it('sorts every intersection by line parameter', () => {
        // A zig-zag polyline crossed by a diagonal line.
        const mesh = SegmentMesh.fromContiguous([
            vec(-3, -1), vec(-2, 1), vec(-1, -1), vec(0, 1), vec(1, -1),
            vec(2, 1), vec(3, -1)
        ], true);
        const result = query.find(line([-4, 0], [1, 0]), mesh);
        expect(result.intersections.length).toBe(6);
        for (let i = 1; i < result.intersections.length; ++i) {
            expect(result.intersections[i - 1].lineParameter)
                .toBeLessThanOrEqual(result.intersections[i].lineParameter);
        }
        for (const record of result.intersections) {
            expect(record.point.values[1]).toBeCloseTo(0, 12);
        }
    });

    it('reports points consistent with both parameterizations on random lines', () => {
        let seed = 13579;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const vertices: Vector[] = [];
        for (let i = 0; i < 12; ++i) {
            const theta = (2 * Math.PI * i) / 12;
            vertices.push(vec(2 * Math.cos(theta), 1.3 * Math.sin(theta)));
        }
        const mesh = SegmentMesh.fromContiguous(vertices, false);
        const meshVertices = mesh.getVertices();

        for (let trial = 0; trial < 150; ++trial) {
            const l = line([rand() * 6 - 3, rand() * 6 - 3],
                [rand() * 2 - 1, rand() * 2 - 1]);
            const result = query.find(l, mesh);
            for (const record of result.intersections) {
                const onLine = add(l.origin,
                    mul(record.lineParameter, l.direction));
                expect(sub(onLine, record.point).values[0]).toBeCloseTo(0, 8);
                expect(sub(onLine, record.point).values[1]).toBeCloseTo(0, 8);

                const p0 = meshVertices[record.indexPair[0]];
                const p1 = meshVertices[record.indexPair[1]];
                const onSegment = add(p0,
                    mul(record.meshSegmentParameter, sub(p1, p0)));
                expect(sub(onSegment, record.point).values[0])
                    .toBeCloseTo(0, 8);
                expect(sub(onSegment, record.point).values[1])
                    .toBeCloseTo(0, 8);
                expect(record.meshSegmentParameter).toBeGreaterThanOrEqual(0);
                expect(record.meshSegmentParameter).toBeLessThanOrEqual(1);
            }
            for (let i = 1; i < result.intersections.length; ++i) {
                expect(result.intersections[i - 1].lineParameter)
                    .toBeLessThanOrEqual(
                        result.intersections[i].lineParameter);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrLine2SegmentMesh2.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, wellScaled
} from './helpers/arbitraries.js';
import { IntrLine2Segment2FI } from '../src/IntrLine2Segment2.js';
import { Segment } from '../src/Segment.js';
import { dot, length } from '../src/Vector.js';
import { IntrRay2SegmentMesh2FI } from '../src/IntrRay2SegmentMesh2.js';
import { Ray } from '../src/Ray.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so a direction built from
// (cos a, sin a) never has a subnormal component; a subnormal component makes
// the exact DotPerp tests inside the line-line query underflow.
const angle2 = () => wellScaled(-Math.PI, Math.PI);

// A mesh whose vertices are on an integer lattice, so the collinearity tests
// that upstream performs with exact DotPerp comparisons are reproducible.
const latticeMesh = fc.array(
    fc.tuple(fc.integer({ min: -6, max: 6 }), fc.integer({ min: -6, max: 6 })),
    { minLength: 2, maxLength: 9 })
    .map(pts => SegmentMesh.fromContiguous(
        pts.map(([x, y]) => vec(x, y)), true));

const latticeLine = fc.tuple(
    fc.integer({ min: -6, max: 6 }), fc.integer({ min: -6, max: 6 }),
    fc.integer({ min: -4, max: 4 }), fc.integer({ min: -4, max: 4 }))
    .filter(([, , dx, dy]) => dx !== 0 || dy !== 0)
    .map(([ox, oy, dx, dy]) =>
        Line.fromOriginDirection(vec(ox, oy), vec(dx, dy)));

const generalMesh = fc.array(
    fc.tuple(wellScaled(-6, 6), wellScaled(-6, 6)),
    { minLength: 2, maxLength: 8 })
    .map(pts => SegmentMesh.fromContiguous(
        pts.map(([x, y]) => vec(x, y)), true));

const generalLine = fc.tuple(wellScaled(-5, 5), wellScaled(-5, 5), angle2())
    .map(([ox, oy, a]) =>
        Line.fromOriginDirection(vec(ox, oy),
            vec(Math.cos(a), Math.sin(a))));

describe('IntrLine2SegmentMesh2 verification', () => {
    const query = new IntrLine2SegmentMesh2FI();
    const lsq = new IntrLine2Segment2FI();

    // The independent reference: run the line-segment query per mesh segment
    // and collect the records the header documents.
    function reference(l: Line, mesh: SegmentMesh):
        { indexPair: [number, number], lineParameter: number,
            meshSegmentParameter: number, point: Vector }[] {
        const out: { indexPair: [number, number], lineParameter: number,
            meshSegmentParameter: number, point: Vector }[] = [];
        const vertices = mesh.getVertices();
        const indices = mesh.getIndices();
        for (let i = 0; i < indices.length; ++i) {
            const seg = Segment.fromEndpoints(vertices[indices[i][0]],
                vertices[indices[i][1]]);
            const r = lsq.find(l, seg);
            if (!r.intersect) { continue; }
            if (r.numIntersections === 1) {
                out.push({
                    indexPair: [indices[i][0], indices[i][1]],
                    lineParameter: r.lineParameter[0],
                    meshSegmentParameter: r.segmentParameter[0],
                    point: r.point.clone()
                });
            }
            else {
                for (let j = 0; j < 2; ++j) {
                    out.push({
                        indexPair: [indices[i][0], indices[i][1]],
                        lineParameter: r.lineParameter[j],
                        meshSegmentParameter: r.segmentParameter[j],
                        point: seg.p[j].clone()
                    });
                }
            }
        }
        return out;
    }

    it('the hits are the union of the per-segment line-segment hits', () => {
        check(fc.tuple(latticeLine, latticeMesh), ([l, mesh]) => {
            const got = query.find(l, mesh).intersections;
            const want = reference(l, mesh);
            expect(got.length).toBe(want.length);
            // Both lists carry the same multiset of records; the query sorts
            // by line parameter, so compare sorted copies.
            const key = (r: { indexPair: [number, number],
                lineParameter: number, meshSegmentParameter: number }) =>
                `${r.indexPair[0]},${r.indexPair[1]},${r.lineParameter},`
                + `${r.meshSegmentParameter}`;
            const a = got.map(key).sort();
            const b = want.map(key).sort();
            expect(a).toEqual(b);
        });
    });

    it('the records are sorted by line parameter, as upstream documents',
        () => {
            check(fc.tuple(generalLine, generalMesh), ([l, mesh]) => {
                const got = query.find(l, mesh).intersections;
                for (let i = 1; i < got.length; ++i) {
                    expect(got[i - 1].lineParameter)
                        .toBeLessThanOrEqual(got[i].lineParameter);
                }
            });
        });

    it('every record is consistent with its line and mesh segment', () => {
        check(fc.tuple(generalLine, generalMesh), ([l, mesh]) => {
            const vertices = mesh.getVertices();
            for (const r of query.find(l, mesh).intersections) {
                expect(r.indexPair[0]).toBeGreaterThanOrEqual(0);
                expect(r.indexPair[1]).toBeLessThan(vertices.length);
                expect(Number.isFinite(r.lineParameter)).toBe(true);
                expect(Number.isFinite(r.meshSegmentParameter)).toBe(true);
                expect(r.meshSegmentParameter)
                    .toBeGreaterThanOrEqual(-1e-12);
                expect(r.meshSegmentParameter).toBeLessThanOrEqual(1 + 1e-12);

                // The point is on the mesh segment at the reported parameter,
                // which the header states is for (1-t)*p0 + t*p1.
                const p0 = vertices[r.indexPair[0]];
                const p1 = vertices[r.indexPair[1]];
                expectVectorClose(r.point,
                    add(p0, mul(r.meshSegmentParameter, sub(p1, p0))),
                    1e-8, 1e-9);

                // The point is on the line at the reported parameter --
                // EXCEPT for a mesh segment collinear with the line, where
                // upstream copies IntrLine2Segment2's sentinel parameters
                // +-max() into the record while storing the true endpoint as
                // the point. See the dedicated test below.
                if (Math.abs(r.lineParameter) >= 1e100) {
                    continue;
                }
                expectVectorClose(r.point,
                    add(l.origin, mul(r.lineParameter, l.direction)),
                    1e-8, 1e-9);
            }
        });
    });

    it('a collinear mesh segment reports the +-max() line-parameter sentinel',
        () => {
            // Upstream bug (preserved): IntrLine2Segment2's collinear result
            // sets lineParameter = { -max(), +max() } with the point marked
            // invalid, and IntrLine2SegmentMesh2 copies those sentinels into
            // the record while replacing the point with the real endpoint.
            // The record is then internally inconsistent, and the sort by
            // line parameter puts both endpoints at the extremes of the list
            // rather than at their geometric positions.
            const mesh = SegmentMesh.fromDisjoint([vec(2, 0), vec(5, 0)]);
            const got = query.find(line([0, 0], [1, 0]), mesh).intersections;
            expect(got.length).toBe(2);
            expect(got[0].lineParameter).toBe(-Number.MAX_VALUE);
            expect(got[1].lineParameter).toBe(Number.MAX_VALUE);
            expectVectorClose(got[0].point, vec(2, 0), 0, 0);
            expectVectorClose(got[1].point, vec(5, 0), 0, 0);
            expect(got[0].meshSegmentParameter).toBe(0);
            expect(got[1].meshSegmentParameter).toBe(1);

            // The downstream consequence in IntrRay2SegmentMesh2, which keeps
            // records with lineParameter >= 0: the ray from (0,0) along +x
            // meets the collinear segment [(2,0),(5,0)] at both endpoints,
            // but only the +max() record survives the filter, so (2,0) is
            // dropped. The mirrored ray keeps a point that is behind it.
            const ray2 = new IntrRay2SegmentMesh2FI();
            const forward = ray2.find(
                Ray.fromOriginDirection(vec(0, 0), vec(1, 0)), mesh)
                .intersections;
            expect(forward.length).toBe(1);
            expectVectorClose(forward[0].point, vec(5, 0), 0, 0);

            const backward = ray2.find(
                Ray.fromOriginDirection(vec(0, 0), vec(-1, 0)), mesh)
                .intersections;
            expect(backward.length).toBe(1);
            // (5,0) is behind this ray, yet it is reported.
            expectVectorClose(backward[0].point, vec(5, 0), 0, 0);
        });

    it('a line through a sampled mesh point reports that point', () => {
        check(fc.tuple(latticeMesh,
            fc.nat({ max: 7 }),
            fc.double({ min: 0.15, max: 0.85, noNaN: true,
                noDefaultInfinity: true }),
            angle2()),
        ([mesh, which, u, ang]) => {
            const vertices = mesh.getVertices();
            const indices = mesh.getIndices();
            const idx = indices[which % indices.length];
            const p0 = vertices[idx[0]], p1 = vertices[idx[1]];
            if (length(sub(p1, p0)) < 1e-6) { return; }
            const target = add(p0, mul(u, sub(p1, p0)));
            const d = vec(Math.cos(ang), Math.sin(ang));
            // Skip a line that is (nearly) parallel to the segment; the hit
            // is then the coincident case with a different record shape.
            const e = sub(p1, p0);
            if (Math.abs(d.get(0) * e.get(1) - d.get(1) * e.get(0))
                < 1e-6 * length(e)) {
                return;
            }
            const l = Line.fromOriginDirection(target, d);
            const got = query.find(l, mesh).intersections;
            let found = false;
            for (const r of got) {
                if (length(sub(r.point, target)) < 1e-8) { found = true; }
            }
            expect(found).toBe(true);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(generalLine, generalMesh, angle2(),
            wellScaled(-3, 3), wellScaled(-3, 3)),
        ([l, mesh, ang, tx, ty]) => {
            const ca = Math.cos(ang), sa = Math.sin(ang);
            const rot = (v: Vector): Vector => vec(
                ca * v.get(0) - sa * v.get(1),
                sa * v.get(0) + ca * v.get(1));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty));
            const l2 = Line.fromOriginDirection(xf(l.origin), rot(l.direction));
            const mesh2 = SegmentMesh.fromContiguous(
                mesh.getVertices().map(xf), true);
            const g0 = query.find(l, mesh).intersections;
            const g1 = query.find(l2, mesh2).intersections;
            if (g0.length !== g1.length) {
                return;   // a grazing vertex may flip the count
            }
            // Skip collinear mesh segments: their records carry the +-max()
            // line-parameter sentinel (see the dedicated test above), which
            // no rigid motion can map onto another.
            for (const r of g0.concat(g1)) {
                if (Math.abs(r.lineParameter) >= 1e100) { return; }
            }
            for (let i = 0; i < g0.length; ++i) {
                expectClose(g1[i].lineParameter, g0[i].lineParameter,
                    1e-7, 1e-8);
                expectVectorClose(g1[i].point, xf(g0[i].point), 1e-7, 1e-8);
            }
        });
    });

    it('a coincident mesh segment contributes both of its endpoints', () => {
        const mesh = SegmentMesh.fromContiguous(
            [vec(0, 0), vec(2, 0), vec(2, 2)], true);
        const got = query.find(line([-3, 0], [1, 0]), mesh).intersections;
        // The first segment lies on the line and reports both endpoints; the
        // second segment meets the line only at (2,0).
        expect(got.length).toBe(3);
        expect(got[0].lineParameter).toBeLessThanOrEqual(got[1].lineParameter);
        expect(got[1].lineParameter).toBeLessThanOrEqual(got[2].lineParameter);
        expectVectorClose(got[0].point, vec(0, 0), 1e-12, 1e-12);
        for (const r of got.slice(1)) {
            expectVectorClose(r.point, vec(2, 0), 1e-12, 1e-12);
        }
    });

    it('an empty-topology mesh yields no intersections', () => {
        const mesh = new SegmentMesh();
        const got = query.find(line([0, 0], [1, 0]), mesh).intersections;
        expect(got.length).toBe(0);
    });

    it('a line that misses every segment yields no intersections', () => {
        const mesh = unitSquare();
        expect(query.find(line([-5, 9], [1, 0]), mesh).intersections.length)
            .toBe(0);
    });

    it('records do not alias the mesh vertices', () => {
        const mesh = SegmentMesh.fromContiguous(
            [vec(0, 0), vec(2, 0)], true);
        const got = query.find(line([-3, 0], [1, 0]), mesh).intersections;
        expect(got.length).toBe(2);
        got[0].point.set(0, 1234);
        expect(mesh.getVertices()[0].get(0)).toBe(0);
        expect(dot(mesh.getVertices()[0], mesh.getVertices()[0])).toBe(0);
    });
});
