import { describe, expect, it } from 'vitest';
import {
    IntrTetrahedron3Tetrahedron3TI,
    defaultIntrTetrahedron3Tetrahedron3TIResult,
    intrTetrahedron3Tetrahedron3InvalidIndex
} from '../src/IntrTetrahedron3Tetrahedron3.js';
import { Tetrahedron3 } from '../src/Tetrahedron3.js';
import { Vector, add, dot, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

const V3 = (x: number, y: number, z: number) => Vector.fromArray([x, y, z]);
const INVALID = intrTetrahedron3Tetrahedron3InvalidIndex;

function translate(tetra: Tetrahedron3, t: Vector): Tetrahedron3 {
    return Tetrahedron3.fromArray(tetra.v.map(v => add(v, t)));
}

// The query (like upstream) assumes the tetrahedron faces are
// counterclockwise ordered when viewed from outside, which is equivalent to
// a positive signed volume. Random vertices give either orientation, so swap
// two vertices when the volume is negative.
function orient(t: Tetrahedron3): Tetrahedron3 {
    const volume = dot(sub(t.v[1], t.v[0]),
        cross(sub(t.v[2], t.v[0]), sub(t.v[3], t.v[0])));
    return volume >= 0 ? t
        : Tetrahedron3.fromVertices(t.v[0], t.v[2], t.v[1], t.v[3]);
}

function unitTetra(): Tetrahedron3 {
    return Tetrahedron3.fromVertices(V3(0, 0, 0), V3(1, 0, 0), V3(0, 1, 0),
        V3(0, 0, 1));
}

// An independent separating-axis reference: project both tetrahedra onto
// each candidate axis and test the resulting intervals for overlap. The
// candidate axes are the 8 face normals and the 36 edge-pair cross products,
// which is a complete set for two convex polytopes.
function referenceIntersect(t0: Tetrahedron3, t1: Tetrahedron3): boolean {
    const project = (t: Tetrahedron3, axis: Vector) => {
        let min = Number.MAX_VALUE, max = -Number.MAX_VALUE;
        for (const v of t.v) {
            const d = dot(axis, v);
            min = Math.min(min, d);
            max = Math.max(max, d);
        }
        return { min, max };
    };

    const axes: Vector[] = [];
    for (const t of [t0, t1]) {
        for (let f = 0; f < 4; ++f) {
            const idx = Tetrahedron3.getFaceIndices(f);
            axes.push(cross(sub(t.v[idx[1]], t.v[idx[0]]),
                sub(t.v[idx[2]], t.v[idx[0]])));
        }
    }
    for (let i0 = 0; i0 < 6; ++i0) {
        const a = Tetrahedron3.getEdgeIndices(i0);
        const e0 = sub(t0.v[a[1]], t0.v[a[0]]);
        for (let i1 = 0; i1 < 6; ++i1) {
            const b = Tetrahedron3.getEdgeIndices(i1);
            const e1 = sub(t1.v[b[1]], t1.v[b[0]]);
            axes.push(cross(e0, e1));
        }
    }

    for (const axis of axes) {
        const lenSqr = dot(axis, axis);
        if (lenSqr <= 1e-24) {
            continue;
        }
        const p0 = project(t0, axis);
        const p1 = project(t1, axis);
        if (p0.max < p1.min || p1.max < p0.min) {
            return false;
        }
    }
    return true;
}

describe('IntrTetrahedron3Tetrahedron3TI default result', () => {
    it('matches the upstream default constructor', () => {
        const r = defaultIntrTetrahedron3Tetrahedron3TIResult();
        expect(r.intersect).toBe(false);
        expect(r.separating).toEqual([INVALID, INVALID]);
    });
});

describe('IntrTetrahedron3Tetrahedron3TI known configurations', () => {
    const ti = new IntrTetrahedron3Tetrahedron3TI();

    it('reports intersection of a tetrahedron with itself', () => {
        const t = unitTetra();
        const r = ti.test(t, t);
        expect(r.intersect).toBe(true);
        expect(r.separating).toEqual([INVALID, INVALID]);
    });

    it('reports intersection when one tetrahedron contains the other', () => {
        const big = Tetrahedron3.fromVertices(V3(-4, -4, -4), V3(8, -4, -4),
            V3(-4, 8, -4), V3(-4, -4, 8));
        const small = translate(unitTetra(), V3(-0.5, -0.5, -0.5));
        expect(ti.test(big, small).intersect).toBe(true);
        expect(ti.test(small, big).intersect).toBe(true);
    });

    it('separates tetrahedra that are far apart, by a face normal', () => {
        const t0 = unitTetra();
        const t1 = translate(unitTetra(), V3(10, 0, 0));
        const r = ti.test(t0, t1);
        expect(r.intersect).toBe(false);
        // Exactly one of the two indices is a valid face index.
        const valid0 = r.separating[0] !== INVALID;
        const valid1 = r.separating[1] !== INVALID;
        expect(valid0 !== valid1).toBe(true);
        const i = valid0 ? r.separating[0] : r.separating[1];
        expect(i).toBeGreaterThanOrEqual(0);
        expect(i).toBeLessThan(4);
    });

    it('is symmetric in its arguments', () => {
        const t0 = unitTetra();
        for (const d of [0.1, 0.5, 0.9, 1.0, 1.5, 3.0]) {
            const t1 = translate(unitTetra(), V3(d, d * 0.25, -d * 0.5));
            expect(ti.test(t0, t1).intersect)
                .toBe(ti.test(t1, t0).intersect);
        }
    });

    it('reports no intersection for face-touching tetrahedra', () => {
        // Two tetrahedra sharing the plane z = 0, one above and one below.
        // The face-normal axis has all of the second tetrahedron's projected
        // heights >= 0 with one strictly positive, so the query reports
        // separation (a touching configuration is not an overlap).
        const t0 = unitTetra();
        const t1 = Tetrahedron3.fromVertices(V3(0, 0, 0), V3(1, 0, 0),
            V3(0, 1, 0), V3(0, 0, -1));
        expect(ti.test(t0, t1).intersect).toBe(false);
        expect(ti.test(t1, t0).intersect).toBe(false);
    });

    it('reports overlap when the two tetrahedra share interior points', () => {
        const t0 = unitTetra();
        const t1 = translate(unitTetra(), V3(0.1, 0.1, 0.1));
        expect(ti.test(t0, t1).intersect).toBe(true);
    });

    it('does not throw for a degenerate (flat) tetrahedron', () => {
        const flat = Tetrahedron3.fromVertices(V3(0, 0, 0), V3(1, 0, 0),
            V3(0, 1, 0), V3(1, 1, 0));
        expect(() => new IntrTetrahedron3Tetrahedron3TI()
            .test(flat, unitTetra())).not.toThrow();
    });

    it('clamps the epsilon argument to [0,1]', () => {
        const t0 = unitTetra();
        const t1 = translate(unitTetra(), V3(1.4, 0, 0));
        // epsilon = 1 skips every edge-edge axis; the face normals alone
        // already separate this configuration.
        expect(ti.test(t0, t1, 1).intersect)
            .toBe(ti.test(t0, t1, 5).intersect);
        expect(ti.test(t0, t1, 0).intersect)
            .toBe(ti.test(t0, t1, -5).intersect);
    });
});

describe('IntrTetrahedron3Tetrahedron3TI edge-edge separation', () => {
    const ti = new IntrTetrahedron3Tetrahedron3TI();

    it('detects separation that only an edge-edge axis reveals', () => {
        // A regression case found by comparing the port with a faithful
        // transcription of the upstream (unnormalized) parallelism test:
        // upstream skips the edge-edge axes because |Dot(E0,E1)| >= 1 for
        // these edge lengths, so upstream reports a false intersection.
        const t0 = Tetrahedron3.fromVertices(
            V3(-1.1178684230511395, 0.8056851638507494, -0.4620938293924992),
            V3(1.8546957987615356, -1.2099614139692678, 1.4172073618589005),
            V3(-0.6184301369909337, -0.8468503946656596, -1.3403615948466405),
            V3(1.9352771062102527, -0.019677161248250474, -0.4333620064115906));
        const t1 = Tetrahedron3.fromVertices(
            V3(1.0989275011694652, -0.9575566029909797, 0.05797815604972589),
            V3(4.566261773261363, -1.9749088268610224, 0.005514266903285936),
            V3(3.1081347507928196, 0.3338751822402122, -1.7254567126396376),
            V3(3.419419889067961, -0.19753169743229249, 1.5311250721715042));

        const r = ti.test(t0, t1);
        expect(r.intersect).toBe(false);
        expect(r.separating[0]).toBe(4);
        expect(r.separating[1]).toBe(2);

        // Verify the reported axis independently: the projections of the two
        // tetrahedra onto it must be disjoint intervals.
        const a = Tetrahedron3.getEdgeIndices(r.separating[0]);
        const b = Tetrahedron3.getEdgeIndices(r.separating[1]);
        const axis = cross(sub(t0.v[a[1]], t0.v[a[0]]),
            sub(t1.v[b[1]], t1.v[b[0]]));
        const proj = (t: Tetrahedron3) => {
            const ds = t.v.map(v => dot(axis, v));
            return { min: Math.min(...ds), max: Math.max(...ds) };
        };
        const p0 = proj(t0), p1 = proj(t1);
        expect(p0.max < p1.min || p1.max < p0.min).toBe(true);
    });
});

describe('IntrTetrahedron3Tetrahedron3TI randomized cross-check', () => {
    it('agrees with an independent interval-based separating-axis test',
        () => {
            const ti = new IntrTetrahedron3Tetrahedron3TI();
            let seed = 20260101;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) & 0x7fffffff;
                return seed / 0x7fffffff;
            };
            const rv = (s: number) => V3((rand() - 0.5) * s,
                (rand() - 0.5) * s, (rand() - 0.5) * s);

            let hits = 0, misses = 0;
            for (let k = 0; k < 2500; ++k) {
                const t0 = orient(Tetrahedron3.fromVertices(rv(4), rv(4),
                    rv(4), rv(4)));
                const offset = rv(8);
                const t1 = translate(orient(
                    Tetrahedron3.fromVertices(rv(4), rv(4), rv(4), rv(4))),
                    offset);

                const expected = referenceIntersect(t0, t1);
                expect(ti.test(t0, t1).intersect).toBe(expected);
                // Symmetry of the query.
                expect(ti.test(t1, t0).intersect).toBe(expected);
                if (expected) {
                    ++hits;
                } else {
                    ++misses;
                }
            }
            expect(hits).toBeGreaterThan(100);
            expect(misses).toBeGreaterThan(1000);
        });

    it('reports intersection whenever the two tetrahedra share a point',
        () => {
            const ti = new IntrTetrahedron3Tetrahedron3TI();
            let seed = 555111;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) & 0x7fffffff;
                return seed / 0x7fffffff;
            };
            const rv = (s: number) => V3((rand() - 0.5) * s,
                (rand() - 0.5) * s, (rand() - 0.5) * s);

            // Build a strictly interior point of a tetrahedron with positive
            // barycentric weights, then translate the tetrahedron so that
            // this point coincides with a shared point Q.
            const shift = (t: Tetrahedron3, q: Vector) => {
                const w = [rand() + 0.2, rand() + 0.2, rand() + 0.2,
                    rand() + 0.2];
                const sum = w[0] + w[1] + w[2] + w[3];
                let p = V3(0, 0, 0);
                for (let i = 0; i < 4; ++i) {
                    p = add(p, Vector.fromArray(
                        t.v[i].values.map(x => (x * w[i]) / sum)));
                }
                return translate(t, sub(q, p));
            };

            for (let k = 0; k < 1000; ++k) {
                const q = rv(10);
                const t0 = shift(orient(
                    Tetrahedron3.fromVertices(rv(3), rv(3), rv(3), rv(3))), q);
                const t1 = shift(orient(
                    Tetrahedron3.fromVertices(rv(3), rv(3), rv(3), rv(3))), q);
                expect(ti.test(t0, t1).intersect).toBe(true);
            }
        });
});
