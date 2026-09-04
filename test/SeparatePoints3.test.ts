import { describe, expect, it } from 'vitest';
import { SeparatePoints3 } from '../src/SeparatePoints3.js';
import { Vector, dot } from '../src/Vector.js';
import { check, fc, latticeVector } from './helpers/arbitraries.js';
import { exactDyadic } from './helpers/exact.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The eight corners of the axis-aligned box [x0,x1] x [y0,y1] x [z0,z1].
function box(x0: number, x1: number, y0: number, y1: number,
    z0: number, z1: number): Vector[] {
    const points: Vector[] = [];
    for (const x of [x0, x1]) {
        for (const y of [y0, y1]) {
            for (const z of [z0, z1]) {
                points.push(v3(x, y, z));
            }
        }
    }
    return points;
}

// Rotate a point about the axis (1,1,1)/sqrt(3) by the given angle, using
// Rodrigues' formula. This produces point sets that are not axis aligned, so
// the plane normals and constants are inexact floating-point numbers.
function rotate(p: Vector, angle: number): Vector {
    const s = 1 / Math.sqrt(3);
    const k = [s, s, s];
    const c = Math.cos(angle), sn = Math.sin(angle);
    const v = p.values;
    const kd = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
    const kx = [
        k[1] * v[2] - k[2] * v[1],
        k[2] * v[0] - k[0] * v[2],
        k[0] * v[1] - k[1] * v[0]
    ];
    return v3(
        v[0] * c + kx[0] * sn + k[0] * kd * (1 - c),
        v[1] * c + kx[1] * sn + k[1] * kd * (1 - c),
        v[2] * c + kx[2] * sn + k[2] * kd * (1 - c));
}

function translate(points: readonly Vector[], t: Vector): Vector[] {
    return points.map(p => v3(p.values[0] + t.values[0],
        p.values[1] + t.values[1], p.values[2] + t.values[2]));
}

// Verify that the reported plane really separates the two point sets: every
// point of one set is on or below the plane and every point of the other set
// is on or above it, with at least one point strictly off the plane in each
// set (otherwise the "separation" is vacuous).
function verifySeparation(plane: { normal: Vector; constant: number },
    points0: readonly Vector[], points1: readonly Vector[]): void {
    const n = plane.normal;
    // The normal must be unit length.
    expect(dot(n, n)).toBeCloseTo(1, 12);

    const signed = (points: readonly Vector[]): number[] =>
        points.map(p => dot(n, p) - plane.constant);
    const s0 = signed(points0);
    const s1 = signed(points1);
    const tol = 1e-12 * (1 + Math.max(
        ...s0.map(Math.abs), ...s1.map(Math.abs)));

    const zeroBelowOneAbove =
        s0.every(s => s <= tol) && s1.every(s => s >= -tol);
    const zeroAboveOneBelow =
        s0.every(s => s >= -tol) && s1.every(s => s <= tol);
    expect(zeroBelowOneAbove || zeroAboveOneBelow).toBe(true);

    // The plane is not degenerate: the two sets are not both on the plane.
    expect(Math.max(...s0.map(Math.abs)) > tol
        || Math.max(...s1.map(Math.abs)) > tol).toBe(true);
}

// A deterministic pseudorandom generator so the tests are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function randomCloud(rand: () => number, count: number, cx: number,
    cy: number, cz: number, radius: number): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        points.push(v3(
            cx + radius * (2 * rand() - 1),
            cy + radius * (2 * rand() - 1),
            cz + radius * (2 * rand() - 1)));
    }
    return points;
}

describe('SeparatePoints3', () => {
    const query = new SeparatePoints3();

    it('separates two disjoint axis-aligned boxes', () => {
        const points0 = box(0, 1, 0, 1, 0, 1);
        const points1 = box(3, 4, 0, 1, 0, 1);
        const result = query.compute(points0, points1);
        expect(result.separated).toBe(true);
        verifySeparation(result.separatingPlane, points0, points1);
    });

    it('separates boxes that are disjoint along a diagonal direction', () => {
        const points0 = box(0, 1, 0, 1, 0, 1);
        const points1 = box(2, 3, 2, 3, 2, 3);
        const result = query.compute(points0, points1);
        expect(result.separated).toBe(true);
        verifySeparation(result.separatingPlane, points0, points1);
    });

    it('separates boxes that touch face to face', () => {
        // The boxes share the plane x = 1. The supporting plane of the face
        // x = 1 of box 0 is a separating plane in the sense of the algorithm
        // (the hulls have disjoint interiors).
        const points0 = box(0, 1, 0, 1, 0, 1);
        const points1 = box(1, 2, 0, 1, 0, 1);
        const result = query.compute(points0, points1);
        expect(result.separated).toBe(true);
        verifySeparation(result.separatingPlane, points0, points1);
    });

    it('separates boxes that nearly touch', () => {
        const gap = 1e-8;
        const points0 = box(0, 1, 0, 1, 0, 1);
        const points1 = box(1 + gap, 2 + gap, 0, 1, 0, 1);
        const result = query.compute(points0, points1);
        expect(result.separated).toBe(true);
        verifySeparation(result.separatingPlane, points0, points1);
    });

    it('does not separate boxes that overlap slightly', () => {
        const overlap = 1e-8;
        const points0 = box(0, 1, 0, 1, 0, 1);
        const points1 = box(1 - overlap, 2, 0, 1, 0, 1);
        expect(query.compute(points0, points1).separated).toBe(false);
    });

    it('does not separate boxes that overlap heavily', () => {
        const points0 = box(0, 2, 0, 2, 0, 2);
        const points1 = box(1, 3, 1, 3, 1, 3);
        expect(query.compute(points0, points1).separated).toBe(false);
    });

    it('does not separate a box contained in another box', () => {
        const points0 = box(0, 10, 0, 10, 0, 10);
        const points1 = box(4, 6, 4, 6, 4, 6);
        expect(query.compute(points0, points1).separated).toBe(false);
    });

    it('does not separate rotated overlapping boxes (the round-off case)', () => {
        // Rotating the boxes makes every candidate plane normal and constant
        // inexact. Upstream's floating-point side tests report a spurious
        // separation for configurations like this; the exact orientation
        // predicate of the port does not.
        const base = box(-1, 1, -1, 1, -1, 1);
        const points0 = base.map(p => rotate(p, 0.3));
        const points1 = translate(base.map(p => rotate(p, 1.1)),
            v3(0.5, 0.25, -0.125));
        expect(query.compute(points0, points1).separated).toBe(false);
    });

    it('does not separate the overlapping clouds upstream misclassifies', () => {
        // A concrete configuration for which the upstream floating-point side
        // tests report a separating plane. The reported plane is the
        // supporting plane of a face of hull 0, and both point sets lie on
        // its negative side (upstream's WhichSide returns +1 for hull 0 only
        // because a vertex of the face itself evaluates to a tiny positive
        // number). See the port notes in src/SeparatePoints3.ts.
        const points0 = [
            v3(0.9347284217365086, -0.7016730923205614, 0.0681410632096231),
            v3(0.9753749435767531, 0.45009304070845246, -0.40927881747484207),
            v3(-0.35152136581018567, 0.3707107463851571, -0.22273730067536235),
            v3(0.6667292881757021, 0.04053660621866584, 0.6666020704433322),
            v3(0.2834406332112849, 0.49213195964694023, 0.42226726887747645),
            v3(0.8978642271831632, -0.07511199964210391, 0.6709316726773977),
            v3(-0.9854007088579237, 0.35722421016544104, -0.8994384235702455),
            v3(0.7301426827907562, -0.7787917708046734, 0.09983729664236307)
        ];
        const points1 = [
            v3(0.2809420191682874, 0.9698736068792642, 0.42090822020545593),
            v3(-0.615513819642365, -0.18522269446402784, -0.25008522290736435),
            v3(0.02358881300315252, 0.12438158569857483, -0.18564756670966742),
            v3(0.1132712844759225, 0.840209798142314, -0.23032965175807474),
            v3(-0.3343439157120883, 0.14911684924736623, 0.27891098530963054),
            v3(-0.5629287933930754, 0.9055949041619897, -1.0867322510108353),
            v3(0.12913488326594225, 0.1969857557676733, -0.7294633005745709),
            v3(-0.771140193939209, -0.17589914873242374, -0.9751263245940208)
        ];
        expect(query.compute(points0, points1).separated).toBe(false);
        expect(query.compute(points1, points0).separated).toBe(false);
    });

    it('separates rotated boxes that are far apart', () => {
        const base = box(-1, 1, -1, 1, -1, 1);
        const points0 = base.map(p => rotate(p, 0.3));
        const points1 = translate(base.map(p => rotate(p, 1.1)),
            v3(6, 3, -2));
        const result = query.compute(points0, points1);
        expect(result.separated).toBe(true);
        verifySeparation(result.separatingPlane, points0, points1);
    });

    it('finds an edge-edge separating plane for crossed slabs', () => {
        // Two long thin boxes whose axes are skew: no face of either hull
        // separates them, so the cross-product planes are required.
        const points0 = [
            v3(-4, -0.25, -0.25), v3(-4, -0.25, 0.25),
            v3(-4, 0.25, -0.25), v3(-4, 0.25, 0.25),
            v3(4, -0.25, -0.25), v3(4, -0.25, 0.25),
            v3(4, 0.25, -0.25), v3(4, 0.25, 0.25)
        ];
        const points1 = [
            v3(-0.25, -4, 1.5), v3(0.25, -4, 1.5),
            v3(-0.25, -4, 2.0), v3(0.25, -4, 2.0),
            v3(-0.25, 4, 1.5), v3(0.25, 4, 1.5),
            v3(-0.25, 4, 2.0), v3(0.25, 4, 2.0)
        ];
        const result = query.compute(points0, points1);
        expect(result.separated).toBe(true);
        verifySeparation(result.separatingPlane, points0, points1);
    });

    it('reports no separation when a hull is degenerate', () => {
        // A planar point set has hull dimension 2, so the query bails out.
        const planar = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0), v3(2, 3, 0)
        ];
        const solid = box(10, 11, 10, 11, 10, 11);
        expect(query.compute(planar, solid).separated).toBe(false);
        expect(query.compute(solid, planar).separated).toBe(false);

        // A collinear point set has hull dimension 1.
        const collinear = [
            v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2), v3(3, 3, 3)
        ];
        expect(query.compute(collinear, solid).separated).toBe(false);

        // A single repeated point has hull dimension 0.
        const single = [v3(5, 5, 5), v3(5, 5, 5), v3(5, 5, 5), v3(5, 5, 5)];
        expect(query.compute(single, solid).separated).toBe(false);
    });

    it('is symmetric in its two arguments', () => {
        const cases: [Vector[], Vector[]][] = [
            [box(0, 1, 0, 1, 0, 1), box(3, 4, 0, 1, 0, 1)],
            [box(0, 2, 0, 2, 0, 2), box(1, 3, 1, 3, 1, 3)],
            [box(0, 10, 0, 10, 0, 10), box(4, 6, 4, 6, 4, 6)]
        ];
        for (const [p0, p1] of cases) {
            const forward = query.compute(p0, p1);
            const reverse = query.compute(p1, p0);
            expect(reverse.separated).toBe(forward.separated);
            if (forward.separated) {
                verifySeparation(forward.separatingPlane, p0, p1);
                verifySeparation(reverse.separatingPlane, p1, p0);
            }
        }
    });

    it('handles randomized separated and overlapping clouds', () => {
        const rand = makeRandom(20260902);

        // Clouds in disjoint boxes are always separated and the reported
        // plane is verified point by point.
        for (let trial = 0; trial < 15; ++trial) {
            const points0 = randomCloud(rand, 20, 0, 0, 0, 1);
            const points1 = randomCloud(rand, 20, 5, 0.5, -0.5, 1);
            const result = query.compute(points0, points1);
            expect(result.separated).toBe(true);
            verifySeparation(result.separatingPlane, points0, points1);
        }

        // Clouds about a common center overlap, so no separation exists.
        for (let trial = 0; trial < 15; ++trial) {
            const points0 = randomCloud(rand, 20, 0, 0, 0, 1);
            const points1 = randomCloud(rand, 20, 0.1, -0.05, 0.05, 1);
            expect(query.compute(points0, points1).separated).toBe(false);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V13): property-based cross-checks against an exact oracle.
// ---------------------------------------------------------------------------

// An exact separating-plane search that is independent of the port. Instead of
// the convex hulls it enumerates every candidate plane the separating-axis
// theorem can require directly from the input points: the plane through each
// triple of points of one set, and the plane through a point of set 0 whose
// normal is the cross product of a point pair from set 0 with a point pair from
// set 1. Every hull face normal and every hull edge-edge cross product is in
// this (much larger) family, so for point sets whose hulls are three
// dimensional the search succeeds exactly when the two convex hulls have
// disjoint interiors -- the specification of the query. All arithmetic is
// bigint on the exact dyadic representation of the coordinates, so the
// decision has no round-off.
function exactSeparated(points0: readonly Vector[],
    points1: readonly Vector[]): boolean {
    const all: number[] = [];
    for (const p of points0) { all.push(...p.values); }
    for (const p of points1) { all.push(...p.values); }
    const s = exactDyadic(all);
    const pack = (offset: number, count: number): bigint[][] => {
        const out: bigint[][] = [];
        for (let i = 0; i < count; ++i) {
            const k = 3 * (offset + i);
            out.push([s[k], s[k + 1], s[k + 2]]);
        }
        return out;
    };
    const P0 = pack(0, points0.length);
    const P1 = pack(points0.length, points1.length);

    const subv = (a: bigint[], b: bigint[]): bigint[] =>
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    const crossv = (a: bigint[], b: bigint[]): bigint[] => [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ];
    // -1, 0, +1 as the sign of the exact signed distance; 2 when every point
    // of the set lies on the plane and 0 when the plane splits the set.
    const side = (P: bigint[][], n: bigint[], o: bigint[]): number => {
        let pos = 0, neg = 0;
        for (const q of P) {
            const d = subv(q, o);
            const e = n[0] * d[0] + n[1] * d[1] + n[2] * d[2];
            if (e > 0n) { ++pos; } else if (e < 0n) { ++neg; }
            if (pos !== 0 && neg !== 0) { return 0; }
        }
        return pos !== 0 ? 1 : (neg !== 0 ? -1 : 2);
    };
    const separates = (n: bigint[], o: bigint[]): boolean => {
        if (n[0] === 0n && n[1] === 0n && n[2] === 0n) { return false; }
        const s0 = side(P0, n, o);
        const s1 = side(P1, n, o);
        if (s0 === 0 || s1 === 0) { return false; }
        if (s0 === 2 && s1 === 2) { return false; }   // both sets on the plane
        if (s0 === 2 || s1 === 2) { return true; }
        return s0 * s1 < 0;
    };

    for (const P of [P0, P1]) {
        for (let a = 0; a < P.length; ++a) {
            for (let b = a + 1; b < P.length; ++b) {
                for (let c = b + 1; c < P.length; ++c) {
                    if (separates(crossv(subv(P[b], P[a]), subv(P[c], P[a])),
                        P[a])) {
                        return true;
                    }
                }
            }
        }
    }
    for (let a = 0; a < P0.length; ++a) {
        for (let b = a + 1; b < P0.length; ++b) {
            const d0 = subv(P0[b], P0[a]);
            for (let c = 0; c < P1.length; ++c) {
                for (let d = c + 1; d < P1.length; ++d) {
                    if (separates(crossv(d0, subv(P1[d], P1[c])), P0[a])) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

// A lattice point set whose convex hull is guaranteed to be three dimensional:
// a unit tetrahedron at a lattice base plus extra lattice points near it. Small
// integers are exact in binary64 and stay exact through the hull predicates and
// the bigint oracle above.
const latticeSolid = (spread: number): fc.Arbitrary<Vector[]> =>
    fc.tuple(latticeVector(3, -spread, spread),
        fc.array(latticeVector(3, -2, 2), { minLength: 3, maxLength: 5 }))
        .map(([base, extra]) => {
            const at = (dx: number, dy: number, dz: number): Vector =>
                v3(base.get(0) + dx, base.get(1) + dy, base.get(2) + dz);
            const points = [at(0, 0, 0), at(1, 0, 0), at(0, 1, 0), at(0, 0, 1)];
            for (const e of extra) {
                points.push(at(e.get(0), e.get(1), e.get(2)));
            }
            return points;
        });

describe('SeparatePoints3 verification', () => {
    const query = new SeparatePoints3();

    // With this generator roughly half of the draws are separated and half
    // overlap, so both answers are exercised.
    it('agrees with an exact separating-plane search on lattice input', () => {
        check(fc.tuple(latticeSolid(1), latticeSolid(1)), ([p0, p1]) => {
            const actual = query.compute(p0, p1).separated;
            expect(actual).toBe(exactSeparated(p0, p1));
        }, 120);
    }, 30000);

    it('is symmetric under an argument swap', () => {
        check(fc.tuple(latticeSolid(1), latticeSolid(2)), ([p0, p1]) => {
            expect(query.compute(p1, p0).separated)
                .toBe(query.compute(p0, p1).separated);
        }, 120);
    }, 30000);

    it('reports a plane that really separates the sets', () => {
        check(fc.tuple(latticeSolid(4), latticeSolid(4)), ([p0, p1]) => {
            const result = query.compute(p0, p1);
            if (result.separated) {
                verifySeparation(result.separatingPlane, p0, p1);
            }
        }, 120);
    }, 30000);

    it('never separates sets that share a common interior point', () => {
        // Both sets contain the cube [-1,1]^3 shifted to a common center, so
        // their hulls have a common interior and no separation exists.
        check(fc.tuple(latticeVector(3, -5, 5),
            fc.array(latticeVector(3, -3, 3), { minLength: 4, maxLength: 6 }),
            fc.array(latticeVector(3, -3, 3), { minLength: 4, maxLength: 6 })),
        ([c, e0, e1]) => {
            const shell = (extra: Vector[]): Vector[] => {
                const points: Vector[] = [];
                for (const x of [-1, 1]) {
                    for (const y of [-1, 1]) {
                        for (const z of [-1, 1]) {
                            points.push(v3(c.get(0) + x, c.get(1) + y,
                                c.get(2) + z));
                        }
                    }
                }
                for (const p of extra) {
                    points.push(v3(c.get(0) + 2 * p.get(0),
                        c.get(1) + 2 * p.get(1), c.get(2) + 2 * p.get(2)));
                }
                return points;
            };
            expect(query.compute(shell(e0), shell(e1)).separated).toBe(false);
        }, 60);
    }, 30000);

    it('reports no separation for degenerate (non-3D) hulls', () => {
        check(fc.tuple(fc.array(latticeVector(3, -6, 6),
            { minLength: 4, maxLength: 8 }), latticeSolid(3)),
        ([flat, solid]) => {
            // Project the first set onto the plane z = 0: its hull has
            // dimension at most 2, so the query must decline.
            const planar = flat.map(p => v3(p.get(0), p.get(1), 0));
            expect(query.compute(planar, solid).separated).toBe(false);
            expect(query.compute(solid, planar).separated).toBe(false);
        }, 60);
    }, 30000);
});
