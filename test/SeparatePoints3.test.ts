import { describe, expect, it } from 'vitest';
import { SeparatePoints3 } from '../src/SeparatePoints3';
import { Vector, dot } from '../src/Vector';

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
