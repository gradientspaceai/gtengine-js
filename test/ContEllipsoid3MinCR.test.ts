import { describe, it, expect } from 'vitest';
import { getContainerEllipsoid3MinCR } from '../src/ContEllipsoid3MinCR.js';
import { Matrix } from '../src/Matrix.js';
import { Vector, sub } from '../src/Vector.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A rotation matrix whose columns are the ellipsoid axes: rotate about z by
// 'a' and then about x by 'b'.
function rot(a: number, b: number): Matrix {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const Rz = Matrix.fromArray(3, 3, [ca, -sa, 0, sa, ca, 0, 0, 0, 1]);
    const Rx = Matrix.fromArray(3, 3, [1, 0, 0, 0, cb, -sb, 0, sb, cb]);
    const M = new Matrix(3, 3);
    for (let r = 0; r < 3; ++r) {
        for (let c = 0; c < 3; ++c) {
            let s = 0;
            for (let k = 0; k < 3; ++k) { s += Rx.get(r, k) * Rz.get(k, c); }
            M.set(r, c, s);
        }
    }
    return M;
}

// The coordinates (u,v,w) = R^T*(X-C).
function local(X: Vector, C: Vector, R: Matrix): number[] {
    const d = sub(X, C);
    const u: number[] = [];
    for (let j = 0; j < 3; ++j) {
        let s = 0;
        for (let i = 0; i < 3; ++i) { s += R.get(i, j) * d.values[i]; }
        u.push(s);
    }
    return u;
}

// The quadratic form (X-C)^T R D R^T (X-C).
function form(X: Vector, C: Vector, R: Matrix, D: readonly number[]): number {
    const u = local(X, C, R);
    return D[0] * u[0] * u[0] + D[1] * u[1] * u[1] + D[2] * u[2] * u[2];
}

// The constraint coefficients A[i] = (u^2, v^2, w^2).
function constraints(points: readonly Vector[], C: Vector,
    R: Matrix): number[][] {
    return points.map(P => {
        const u = local(P, C, R);
        return [u[0] * u[0], u[1] * u[1], u[2] * u[2]];
    });
}

// An independent solver for the same problem: maximize D[0]*D[1]*D[2] subject
// to A[i].D <= 1 and D >= 0. The objective is strictly increasing in each
// variable, so the maximizer has one, two or three active constraints. This
// enumerates the candidates for each case - the tangency point of a single
// plane, the maximum of x*y*z along the line shared by a pair of planes, and
// the common point of a triple of planes - keeps the feasible ones and
// returns the largest product. It shares no algorithmic idea with the ported
// facet/edge walk.
function maxProductByEnumeration(A: readonly number[][]): number[] {
    const feasible = (D: readonly number[]): boolean => {
        for (let k = 0; k < 3; ++k) {
            if (!(D[k] >= 0) || !isFinite(D[k])) { return false; }
        }
        for (const a of A) {
            if (a[0] * D[0] + a[1] * D[1] + a[2] * D[2] > 1 + 1e-9) {
                return false;
            }
        }
        return true;
    };

    let best = [0, 0, 0];
    const consider = (D: number[]): void => {
        if (feasible(D) && D[0] * D[1] * D[2] > best[0] * best[1] * best[2]) {
            best = D;
        }
    };

    const n = A.length;
    for (let i = 0; i < n; ++i) {
        // The maximum of x*y*z on a single plane is at x = 1/(3*a0), etc.
        if (A[i][0] > 0 && A[i][1] > 0 && A[i][2] > 0) {
            consider([1 / (3 * A[i][0]), 1 / (3 * A[i][1]),
                1 / (3 * A[i][2])]);
        }

        for (let j = i + 1; j < n; ++j) {
            // Maximize x*y*z on the line A[i].D = A[j].D = 1. The line is
            // P + t*W with W = A[i] x A[j]; P is found by setting the
            // largest-magnitude component of W to zero and solving the
            // resulting 2x2 system. The derivative of the cubic in t is a
            // quadratic whose roots are the candidates.
            const W = [
                A[i][1] * A[j][2] - A[i][2] * A[j][1],
                A[i][2] * A[j][0] - A[i][0] * A[j][2],
                A[i][0] * A[j][1] - A[i][1] * A[j][0]
            ];
            let kFix = 0;
            for (let k = 1; k < 3; ++k) {
                if (Math.abs(W[k]) > Math.abs(W[kFix])) { kFix = k; }
            }
            if (W[kFix] === 0) { continue; }
            const k0 = (kFix + 1) % 3, k1 = (kFix + 2) % 3;
            // Solve A[i][k0]*p0 + A[i][k1]*p1 = 1 (and the same for j); the
            // determinant of that system is +-W[kFix].
            const det = A[i][k0] * A[j][k1] - A[i][k1] * A[j][k0];
            if (det === 0) { continue; }
            const P = [0, 0, 0];
            P[kFix] = 0;
            P[k0] = (A[j][k1] - A[i][k1]) / det;
            P[k1] = (A[i][k0] - A[j][k0]) / det;

            // f(t) = (P0+t*W0)(P1+t*W1)(P2+t*W2); f'(t) = c0 + c1*t + c2*t^2.
            const c2 = 3 * W[0] * W[1] * W[2];
            const c1 = 2 * (P[2] * W[0] * W[1] + P[1] * W[0] * W[2] +
                P[0] * W[1] * W[2]);
            const c0 = P[0] * P[1] * W[2] + P[0] * P[2] * W[1] +
                P[1] * P[2] * W[0];
            const roots: number[] = [];
            if (c2 !== 0) {
                const disc = c1 * c1 - 4 * c0 * c2;
                if (disc >= 0) {
                    const sq = Math.sqrt(disc);
                    roots.push((-c1 + sq) / (2 * c2), (-c1 - sq) / (2 * c2));
                }
            } else if (c1 !== 0) {
                roots.push(-c0 / c1);
            }
            for (const t of roots) {
                consider([P[0] + t * W[0], P[1] + t * W[1],
                    P[2] + t * W[2]]);
            }

            for (let k = j + 1; k < n; ++k) {
                // The common point of three planes.
                const D = solve3x3([A[i], A[j], A[k]]);
                if (D !== null) { consider(D); }
            }
        }
    }
    return best;
}

// Solve M*D = (1,1,1) by Cramer's rule; null if M is singular.
function solve3x3(M: readonly number[][]): number[] | null {
    const det =
        M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
        M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
        M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    if (det === 0) { return null; }
    const D: number[] = [];
    for (let c = 0; c < 3; ++c) {
        const N = M.map(row => row.slice());
        for (let r = 0; r < 3; ++r) { N[r][c] = 1; }
        const d =
            N[0][0] * (N[1][1] * N[2][2] - N[1][2] * N[2][1]) -
            N[0][1] * (N[1][0] * N[2][2] - N[1][2] * N[2][0]) +
            N[0][2] * (N[1][0] * N[2][1] - N[1][1] * N[2][0]);
        D.push(d / det);
    }
    return D;
}

function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1103515245) + 12345) >>> 0;
        return s / 4294967296;
    };
}

describe('getContainerEllipsoid3MinCR', () => {
    it('recovers the extents of an axis-aligned ellipsoid sample', () => {
        // Points on x^2/9 + y^2/4 + z^2 = 1, including the six vertices. The
        // vertex constraints bound D[0] <= 1/9, D[1] <= 1/4 and D[2] <= 1,
        // and that D satisfies every other constraint with equality, so it is
        // the maximizer of the product.
        const a = 3, b = 2, c = 1;
        const points: Vector[] = [];
        for (let i = 0; i <= 8; ++i) {
            const phi = Math.PI * i / 8;
            for (let j = 0; j < 16; ++j) {
                const theta = 2 * Math.PI * j / 16;
                points.push(v(
                    a * Math.sin(phi) * Math.cos(theta),
                    b * Math.sin(phi) * Math.sin(theta),
                    c * Math.cos(phi)));
            }
        }
        const C = v(0, 0, 0), R = rot(0, 0);
        const D = getContainerEllipsoid3MinCR(points, C, R);
        expect(D[0]).toBeCloseTo(1 / (a * a), 9);
        expect(D[1]).toBeCloseTo(1 / (b * b), 9);
        expect(D[2]).toBeCloseTo(1 / (c * c), 9);
        for (const P of points) {
            expect(form(P, C, R, D)).toBeLessThanOrEqual(1 + 1e-10);
        }
    });

    it('recovers the extents in a rotated, translated frame', () => {
        const a = 1.5, b = 4, c = 0.8;
        const C = v(2, -1, 0.5);
        const R = rot(0.6, -1.1);
        const points: Vector[] = [];
        const addLocal = (u: number, w: number, s: number): void => {
            const X = C.values.slice();
            for (let i = 0; i < 3; ++i) {
                X[i] += u * R.get(i, 0) + w * R.get(i, 1) + s * R.get(i, 2);
            }
            points.push(Vector.fromArray(X));
        };
        for (let i = 0; i <= 6; ++i) {
            const phi = Math.PI * i / 6;
            for (let j = 0; j < 12; ++j) {
                const theta = 2 * Math.PI * j / 12;
                addLocal(a * Math.sin(phi) * Math.cos(theta),
                    b * Math.sin(phi) * Math.sin(theta),
                    c * Math.cos(phi));
            }
        }
        const D = getContainerEllipsoid3MinCR(points, C, R);
        expect(D[0]).toBeCloseTo(1 / (a * a), 9);
        expect(D[1]).toBeCloseTo(1 / (b * b), 9);
        expect(D[2]).toBeCloseTo(1 / (c * c), 9);
        for (const P of points) {
            expect(form(P, C, R, D)).toBeLessThanOrEqual(1 + 1e-10);
        }
    });

    it('handles a single point (tangency at one third of each intercept)',
        () => {
            // One constraint a0*x+a1*y+a2*z = 1; the maximum of x*y*z on it
            // is at x = 1/(3*a0), y = 1/(3*a1), z = 1/(3*a2).
            const C = v(1, 2, 3), R = rot(0, 0);
            const P = v(4, 4, 5);  // (u,v,w) = (3,2,2), A = (9,4,4)
            const D = getContainerEllipsoid3MinCR([P], C, R);
            expect(D[0]).toBeCloseTo(1 / 27, 9);
            expect(D[1]).toBeCloseTo(1 / 12, 9);
            expect(D[2]).toBeCloseTo(1 / 12, 9);
            expect(form(P, C, R, D)).toBeCloseTo(1, 9);
        });

    it('is symmetric under reflections of the sample (box corners)', () => {
        // The eight corners (+-2, +-1, +-3) give a single distinct
        // constraint (4,1,9); the maximum is at 1/(3*a_k).
        const points: Vector[] = [];
        for (const sx of [-1, 1]) {
            for (const sy of [-1, 1]) {
                for (const sz of [-1, 1]) {
                    points.push(v(2 * sx, 1 * sy, 3 * sz));
                }
            }
        }
        const C = v(0, 0, 0), R = rot(0, 0);
        const D = getContainerEllipsoid3MinCR(points, C, R);
        expect(D[0]).toBeCloseTo(1 / 12, 9);
        expect(D[1]).toBeCloseTo(1 / 3, 9);
        expect(D[2]).toBeCloseTo(1 / 27, 9);
        for (const P of points) {
            expect(form(P, C, R, D)).toBeCloseTo(1, 9);
        }
    });

    it('ignores interior points', () => {
        const outer = [v(2, 0.5, 1), v(-1, 2, 0.7), v(0.4, -1, 2.5),
            v(-1.5, -1.5, -1.5), v(1, 1, -2)];
        const C = v(0, 0, 0), R = rot(0, 0);
        const D0 = getContainerEllipsoid3MinCR(outer, C, R);
        const inner = outer.concat([v(0.1, 0.1, 0.1), v(-0.2, 0.05, 0.1),
            v(0.3, -0.1, 0.2)]);
        const D1 = getContainerEllipsoid3MinCR(inner, C, R);
        for (let k = 0; k < 3; ++k) {
            expect(D1[k]).toBeCloseTo(D0[k], 9);
        }
    });

    it('is deterministic (the jitter uses a fixed-seed generator)', () => {
        const points = [v(2, 0.5, 1), v(-1, 2, 0.7), v(0.4, -1, 2.5),
            v(-1.5, -1.5, -1.5)];
        const C = v(0, 0, 0), R = rot(0.3, 0.2);
        const D0 = getContainerEllipsoid3MinCR(points, C, R);
        const D1 = getContainerEllipsoid3MinCR(points, C, R);
        expect(D1).toEqual(D0);
    });

    it('matches an independent exhaustive-candidate solver on random clouds',
        () => {
            const rnd = makeRandom(987654321);
            const numTrials = 60;
            let numOptimal = 0;
            for (let trial = 0; trial < numTrials; ++trial) {
                const R = rot(2 * Math.PI * rnd(), 2 * Math.PI * rnd());
                const C = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
                const numPoints = 4 + Math.floor(8 * rnd());
                const points: Vector[] = [];
                for (let i = 0; i < numPoints; ++i) {
                    points.push(v(C.values[0] + 4 * rnd() - 2,
                        C.values[1] + 4 * rnd() - 2,
                        C.values[2] + 4 * rnd() - 2));
                }
                // Keep the problem bounded in every direction.
                const bound = (sx: number, sy: number, sz: number): void => {
                    const X = C.values.slice();
                    for (let i = 0; i < 3; ++i) {
                        X[i] += 1.3 * sx * R.get(i, 0) +
                            1.1 * sy * R.get(i, 1) + 0.9 * sz * R.get(i, 2);
                    }
                    points.push(Vector.fromArray(X));
                };
                bound(1, 1, 1);
                bound(-1, 1, -1);

                const D = getContainerEllipsoid3MinCR(points, C, R);
                for (let k = 0; k < 3; ++k) {
                    expect(D[k]).toBeGreaterThan(0);
                    expect(isFinite(D[k])).toBe(true);
                }

                // Every point is inside or on the ellipsoid.
                let maxForm = 0;
                for (const P of points) {
                    const f = form(P, C, R, D);
                    expect(f).toBeLessThanOrEqual(1 + 1e-9);
                    maxForm = Math.max(maxForm, f);
                }

                // At least one constraint is active (the fit is tight).
                expect(maxForm).toBeGreaterThan(1 - 1e-8);

                // The product never exceeds the true maximum (the fit is
                // never smaller than the minimum-volume ellipsoid). The
                // ported algorithm jitters the constraints by up to 1e-12,
                // so the agreement is not exact.
                const E = maxProductByEnumeration(constraints(points, C, R));
                const pD = D[0] * D[1] * D[2];
                const pE = E[0] * E[1] * E[2];
                expect(pD).toBeLessThan(pE * (1 + 1e-6));
                if (pD > pE * (1 - 1e-6)) {
                    ++numOptimal;
                }
            }

            // Upstream's facet/edge walk can stall at a vertex of the
            // constraint polytope (its FindEdgeMax returns unchanged when
            // tMax is 0), so it is a heuristic: it always returns a
            // containing ellipsoid but not always the minimum-volume one.
            // See the port note in ContEllipsoid3MinCR.ts.
            expect(numOptimal).toBeGreaterThan(0.8 * numTrials);
        });

    it('degenerates gracefully when the points lie on a coordinate plane',
        () => {
            // Every point has w = 0, so nothing bounds D[2]. Upstream's
            // jitter of up to 1e-12 makes the third coefficient positive
            // anyway, so the fit is enormous along the third axis but still
            // contains the points.
            const points = [v(1, 0, 0), v(0, 1, 0), v(-1, -1, 0)];
            const C = v(0, 0, 0), R = rot(0, 0);
            const D = getContainerEllipsoid3MinCR(points, C, R);
            expect(D[2]).toBeGreaterThan(1e6);
            for (const P of points) {
                expect(form(P, C, R, D)).toBeLessThanOrEqual(1 + 1e-9);
            }
        });

    it('throws for an empty point set', () => {
        expect(() => getContainerEllipsoid3MinCR([], v(0, 0, 0), rot(0, 0)))
            .toThrow('no points');
    });

    it('validates the dimensions of its inputs', () => {
        expect(() => getContainerEllipsoid3MinCR([v(1, 1, 1)],
            Vector.fromArray([0, 0]), rot(0, 0))).toThrow('3D');
        expect(() => getContainerEllipsoid3MinCR(
            [Vector.fromArray([1, 1])], v(0, 0, 0), rot(0, 0)))
            .toThrow('3D');
        expect(() => getContainerEllipsoid3MinCR([v(1, 1, 1)], v(0, 0, 0),
            new Matrix(2, 2))).toThrow('3x3');
    });
});
