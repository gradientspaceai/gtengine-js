import { describe, expect, it } from 'vitest';
import { MinimumVolumeSphere3 } from '../src/MinimumVolumeSphere3';
import { Vector, sub, dot } from '../src/Vector';

function cross(u: Vector, v: Vector): Vector {
    return Vector.fromArray([
        u.get(1) * v.get(2) - u.get(2) * v.get(1),
        u.get(2) * v.get(0) - u.get(0) * v.get(2),
        u.get(0) * v.get(1) - u.get(1) * v.get(0)]);
}

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function distance(p: Vector, q: Vector): number {
    const d = sub(p, q);
    return Math.sqrt(dot(d, d));
}

// The circumcenter of a triangle in 3D, computed with the cross-product
// formula C = P0 + ((|a|^2 b - |b|^2 a) x (a x b)) / (2 |a x b|^2), which is
// a different derivation from the barycentric linear solve used by the port.
function circumcenter3(P0: Vector, P1: Vector, P2: Vector): Vector | null {
    const a = sub(P1, P0);
    const b = sub(P2, P0);
    const axb = cross(a, b);
    const denom = 2 * dot(axb, axb);
    if (denom < 1e-24) {
        return null;
    }
    const aa = dot(a, a);
    const bb = dot(b, b);
    const num = cross(
        Vector.fromArray([
            aa * b.get(0) - bb * a.get(0),
            aa * b.get(1) - bb * a.get(1),
            aa * b.get(2) - bb * a.get(2)]),
        axb);
    return Vector.fromArray([
        P0.get(0) + num.get(0) / denom,
        P0.get(1) + num.get(1) / denom,
        P0.get(2) + num.get(2) / denom]);
}

// The circumcenter of a tetrahedron via Cramer's rule on the linear system
// 2*(Pi - P0).C = |Pi|^2 - |P0|^2, i = 1,2,3.
function circumcenter4(P0: Vector, P1: Vector, P2: Vector, P3: Vector):
    Vector | null {
    const rows: number[][] = [];
    const rhs: number[] = [];
    for (const P of [P1, P2, P3]) {
        rows.push([
            2 * (P.get(0) - P0.get(0)),
            2 * (P.get(1) - P0.get(1)),
            2 * (P.get(2) - P0.get(2))]);
        rhs.push(dot(P, P) - dot(P0, P0));
    }
    const det3 = (m: number[][]): number =>
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    const d = det3(rows);
    if (Math.abs(d) < 1e-14) {
        return null;
    }
    const solution: number[] = [];
    for (let c = 0; c < 3; ++c) {
        const m = rows.map((row, r) =>
            row.map((value, col) => (col === c ? rhs[r] : value)));
        solution.push(det3(m) / d);
    }
    return Vector.fromArray(solution);
}

// An independent brute-force minimum-volume sphere: enumerate the spheres
// defined by every pair, triple and quadruple of points and keep the
// smallest one that contains every point.
function bruteForceMinimumSphere(points: readonly Vector[]):
    { center: Vector, radius: number } {
    const n = points.length;
    if (n === 1) {
        return { center: points[0].clone(), radius: 0 };
    }

    let best: { center: Vector, radius: number } | null = null;
    const consider = (center: Vector | null, radius: number): void => {
        if (center === null) {
            return;
        }
        const tol = 1e-9 * Math.max(1, radius);
        for (const p of points) {
            if (distance(p, center) > radius + tol) {
                return;
            }
        }
        if (best === null || radius < best.radius) {
            best = { center, radius };
        }
    };

    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            const center = Vector.fromArray([
                0.5 * (points[i].get(0) + points[j].get(0)),
                0.5 * (points[i].get(1) + points[j].get(1)),
                0.5 * (points[i].get(2) + points[j].get(2))]);
            consider(center, 0.5 * distance(points[i], points[j]));

            for (let k = j + 1; k < n; ++k) {
                const c3 = circumcenter3(points[i], points[j], points[k]);
                if (c3 !== null) {
                    consider(c3, distance(c3, points[i]));
                }
                for (let l = k + 1; l < n; ++l) {
                    const c4 = circumcenter4(points[i], points[j], points[k],
                        points[l]);
                    if (c4 !== null) {
                        consider(c4, distance(c4, points[i]));
                    }
                }
            }
        }
    }

    if (best === null) {
        throw new Error('brute force failed');
    }
    return best;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (1664525 * state + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('MinimumVolumeSphere3', () => {
    it('handles a single point', () => {
        const mvs = new MinimumVolumeSphere3();
        const { minimal, success } = mvs.compute([v3(1, 2, 3)]);
        expect(success).toBe(true);
        expect(minimal.radius).toBe(0);
        expect(minimal.center.get(2)).toBe(3);
        expect(mvs.numSupport).toBe(1);
    });

    it('computes the diametral sphere of two points', () => {
        const mvs = new MinimumVolumeSphere3();
        const { minimal, success } = mvs.compute([v3(0, 0, -3), v3(0, 0, 3)]);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(3, 12);
        expect(minimal.center.get(2)).toBeCloseTo(0, 12);
        expect(mvs.numSupport).toBe(2);
    });

    it('computes the circumsphere of a regular tetrahedron', () => {
        // Vertices of a regular tetrahedron inscribed in a sphere of radius
        // sqrt(3) centered at the origin.
        const points = [
            v3(1, 1, 1), v3(1, -1, -1), v3(-1, 1, -1), v3(-1, -1, 1)
        ];
        const mvs = new MinimumVolumeSphere3();
        const { minimal, success } = mvs.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(Math.sqrt(3), 10);
        expect(minimal.center.get(0)).toBeCloseTo(0, 10);
        expect(minimal.center.get(1)).toBeCloseTo(0, 10);
        expect(minimal.center.get(2)).toBeCloseTo(0, 10);
        expect(mvs.numSupport).toBe(4);
    });

    it('computes the sphere of the unit cube corners', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            points.push(v3(i & 1, (i >> 1) & 1, (i >> 2) & 1));
        }
        const { minimal, success } = new MinimumVolumeSphere3()
            .compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(0.5 * Math.sqrt(3), 10);
        for (let j = 0; j < 3; ++j) {
            expect(minimal.center.get(j)).toBeCloseTo(0.5, 10);
        }
    });

    it('uses only 3 support points for a flat obtuse configuration', () => {
        // Three points forming an acute triangle in the z = 0 plane, plus a
        // point well inside the circumsphere.
        const points = [
            v3(1, 0, 0),
            v3(-0.5, Math.sqrt(3) / 2, 0),
            v3(-0.5, -Math.sqrt(3) / 2, 0),
            v3(0, 0, 0.1)
        ];
        const mvs = new MinimumVolumeSphere3();
        const { minimal, success } = mvs.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(1, 10);
        expect(mvs.numSupport).toBe(3);
    });

    it('ignores duplicate points', () => {
        const points = [
            v3(0, 0, 0), v3(2, 0, 0), v3(0, 0, 0), v3(2, 0, 0),
            v3(1, 1, 0), v3(1, 1, 0)
        ];
        const mvs = new MinimumVolumeSphere3();
        const { minimal, success } = mvs.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(1, 10);
        expect(minimal.center.get(0)).toBeCloseTo(1, 10);
        expect(minimal.center.get(1)).toBeCloseTo(0, 10);
        for (let i = 0; i < mvs.numSupport; ++i) {
            const index = mvs.support[i];
            expect(index).toBeLessThan(points.length);
            expect(distance(points[index], minimal.center))
                .toBeCloseTo(minimal.radius, 9);
        }
    });

    it('handles a set of identical points', () => {
        const mvs = new MinimumVolumeSphere3();
        const { minimal, success } = mvs.compute(
            [v3(-1, 4, 7), v3(-1, 4, 7)]);
        expect(success).toBe(true);
        expect(minimal.radius).toBe(0);
        expect(mvs.numSupport).toBe(1);
    });

    it('handles collinear points', () => {
        const points = [
            v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2), v3(-2, -2, -2)
        ];
        const mvs = new MinimumVolumeSphere3();
        const { minimal, success } = mvs.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(2 * Math.sqrt(3), 9);
        expect(minimal.center.get(0)).toBeCloseTo(0, 9);
        expect(mvs.numSupport).toBe(2);
    });

    it('handles coplanar cocircular points', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 6; ++i) {
            const angle = 2 * Math.PI * i / 6;
            points.push(v3(3 * Math.cos(angle), 3 * Math.sin(angle), 2));
        }
        const { minimal, success } = new MinimumVolumeSphere3()
            .compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(3, 8);
        expect(minimal.center.get(2)).toBeCloseTo(2, 8);
    });

    it('throws for an empty point set', () => {
        expect(() => new MinimumVolumeSphere3().compute([]))
            .toThrow('Input must contain points.');
    });

    it('throws for non-3D points', () => {
        expect(() => new MinimumVolumeSphere3()
            .compute([Vector.fromArray([1, 2])])).toThrow();
    });

    it('matches brute force on random point sets', () => {
        const random = makeRandom(31337);
        for (let trial = 0; trial < 25; ++trial) {
            const numPoints = 4 + (trial % 6);
            const points: Vector[] = [];
            for (let i = 0; i < numPoints; ++i) {
                points.push(v3(4 * random() - 2, 4 * random() - 2,
                    4 * random() - 2));
            }

            const mvs = new MinimumVolumeSphere3();
            const { minimal, success } = mvs.compute(points);
            expect(success).toBe(true);

            const expected = bruteForceMinimumSphere(points);
            expect(minimal.radius).toBeCloseTo(expected.radius, 7);
            for (let j = 0; j < 3; ++j) {
                expect(minimal.center.get(j))
                    .toBeCloseTo(expected.center.get(j), 6);
            }

            for (const p of points) {
                expect(distance(p, minimal.center))
                    .toBeLessThanOrEqual(minimal.radius + 1e-9);
            }

            expect(mvs.numSupport).toBeGreaterThanOrEqual(2);
            for (let i = 0; i < mvs.numSupport; ++i) {
                expect(distance(points[mvs.support[i]], minimal.center))
                    .toBeCloseTo(minimal.radius, 6);
            }
        }
    });

    it('is minimal: moving the center leaves a point outside', () => {
        const random = makeRandom(9001);
        for (let trial = 0; trial < 15; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 10; ++i) {
                points.push(v3(random(), random(), random()));
            }
            const { minimal, success } = new MinimumVolumeSphere3()
                .compute(points);
            expect(success).toBe(true);

            const eps = 1e-4;
            for (let axis = 0; axis < 3; ++axis) {
                for (const sign of [-1, 1]) {
                    const c = minimal.center.clone();
                    c.set(axis, c.get(axis) + sign * eps);
                    let maxDist = 0;
                    for (const p of points) {
                        maxDist = Math.max(maxDist, distance(p, c));
                    }
                    expect(maxDist).toBeGreaterThan(minimal.radius - 1e-12);
                }
            }
        }
    });
});
