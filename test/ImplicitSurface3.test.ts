import { describe, it, expect } from 'vitest';
import { ImplicitSurface3 } from '../src/ImplicitSurface3.js';
import { Matrix } from '../src/Matrix.js';
import { Vector, add, dot, length, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, nonzero, rotationFrame,
    wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

// F(x,y,z) = x^2 + y^2 + z^2 - r^2.
class Sphere extends ImplicitSurface3 {
    constructor(private r: number) {
        super();
    }

    f(p: Vector): number {
        const [x, y, z] = p.values;
        return x * x + y * y + z * z - this.r * this.r;
    }

    fx(p: Vector): number { return 2 * p.values[0]; }
    fy(p: Vector): number { return 2 * p.values[1]; }
    fz(p: Vector): number { return 2 * p.values[2]; }
    fxx(_p: Vector): number { return 2; }
    fxy(_p: Vector): number { return 0; }
    fxz(_p: Vector): number { return 0; }
    fyy(_p: Vector): number { return 2; }
    fyz(_p: Vector): number { return 0; }
    fzz(_p: Vector): number { return 2; }
}

// F(x,y,z) = (x/a)^2 + (y/b)^2 + (z/c)^2 - 1.
class Ellipsoid extends ImplicitSurface3 {
    constructor(private a: number, private b: number, private c: number) {
        super();
    }

    f(p: Vector): number {
        const [x, y, z] = p.values;
        return x * x / (this.a * this.a) + y * y / (this.b * this.b) +
            z * z / (this.c * this.c) - 1;
    }

    fx(p: Vector): number { return 2 * p.values[0] / (this.a * this.a); }
    fy(p: Vector): number { return 2 * p.values[1] / (this.b * this.b); }
    fz(p: Vector): number { return 2 * p.values[2] / (this.c * this.c); }
    fxx(_p: Vector): number { return 2 / (this.a * this.a); }
    fxy(_p: Vector): number { return 0; }
    fxz(_p: Vector): number { return 0; }
    fyy(_p: Vector): number { return 2 / (this.b * this.b); }
    fyz(_p: Vector): number { return 0; }
    fzz(_p: Vector): number { return 2 / (this.c * this.c); }
}

// A surface with mixed second-order partials, F = x*y + y*z + z*x - 1, to
// exercise the off-diagonal Hessian entries.
class Mixed extends ImplicitSurface3 {
    f(p: Vector): number {
        const [x, y, z] = p.values;
        return x * y + y * z + z * x - 1;
    }

    fx(p: Vector): number { return p.values[1] + p.values[2]; }
    fy(p: Vector): number { return p.values[0] + p.values[2]; }
    fz(p: Vector): number { return p.values[0] + p.values[1]; }
    fxx(_p: Vector): number { return 0; }
    fxy(_p: Vector): number { return 1; }
    fxz(_p: Vector): number { return 1; }
    fyy(_p: Vector): number { return 0; }
    fyz(_p: Vector): number { return 1; }
    fzz(_p: Vector): number { return 0; }
}

function V3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Central-difference approximation of dF/dx_i.
function fdGradient(surface: ImplicitSurface3, p: Vector, h: number): Vector {
    const g = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        const pp = p.clone();
        const pm = p.clone();
        pp.values[i] += h;
        pm.values[i] -= h;
        g.values[i] = (surface.f(pp) - surface.f(pm)) / (2 * h);
    }
    return g;
}

// Central-difference approximation of d^2F/(dx_i dx_j).
function fdHessian(surface: ImplicitSurface3, p: Vector, h: number): number[][] {
    const H: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; ++i) {
        for (let j = 0; j < 3; ++j) {
            if (i === j) {
                const pp = p.clone();
                const pm = p.clone();
                pp.values[i] += h;
                pm.values[i] -= h;
                H[i][j] = (surface.f(pp) - 2 * surface.f(p) + surface.f(pm)) / (h * h);
            }
            else {
                const app = p.clone(); app.values[i] += h; app.values[j] += h;
                const apm = p.clone(); apm.values[i] += h; apm.values[j] -= h;
                const amp = p.clone(); amp.values[i] -= h; amp.values[j] += h;
                const amm = p.clone(); amm.values[i] -= h; amm.values[j] -= h;
                H[i][j] = (surface.f(app) - surface.f(apm) - surface.f(amp) +
                    surface.f(amm)) / (4 * h * h);
            }
        }
    }
    return H;
}

describe('ImplicitSurface3', () => {
    it('reports whether a point is on the surface', () => {
        const sphere = new Sphere(2);
        expect(sphere.isOnSurface(V3(2, 0, 0), 1e-12)).toBe(true);
        expect(sphere.isOnSurface(V3(0, 0, 2), 1e-12)).toBe(true);
        expect(sphere.isOnSurface(V3(0, 0, 2.5), 1e-12)).toBe(false);
        // The tolerance is on |F|, not on the distance to the surface.
        expect(sphere.isOnSurface(V3(0, 0, 2.5), 3)).toBe(true);
    });

    it('assembles the gradient from the first-order partials', () => {
        const ellipsoid = new Ellipsoid(1, 2, 3);
        const p = V3(0.3, -0.5, 1.1);
        const g = ellipsoid.getGradient(p);
        expect(g.size).toBe(3);
        expect(g.values[0]).toBeCloseTo(2 * 0.3 / 1, 14);
        expect(g.values[1]).toBeCloseTo(2 * -0.5 / 4, 14);
        expect(g.values[2]).toBeCloseTo(2 * 1.1 / 9, 14);
    });

    it('assembles a symmetric Hessian from the second-order partials', () => {
        const mixed = new Mixed();
        const H = mixed.getHessian(V3(1, 2, 3));
        expect(H.numRows).toBe(3);
        expect(H.numCols).toBe(3);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(H.get(r, c)).toBe(r === c ? 0 : 1);
                expect(H.get(r, c)).toBe(H.get(c, r));
            }
        }
    });

    it('matches finite differences for the gradient and the Hessian', () => {
        const surfaces: ImplicitSurface3[] = [
            new Sphere(2), new Ellipsoid(1, 2, 3), new Mixed()
        ];
        const points = [
            V3(0.7, -1.3, 0.5), V3(-0.2, 0.9, 1.7), V3(1.1, 1.1, -0.4)
        ];
        const h = 1e-4;
        for (const surface of surfaces) {
            for (const p of points) {
                const g = surface.getGradient(p);
                const gFD = fdGradient(surface, p, h);
                for (let i = 0; i < 3; ++i) {
                    expect(g.values[i]).toBeCloseTo(gFD.values[i], 6);
                }

                const H = surface.getHessian(p);
                const HFD = fdHessian(surface, p, h);
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        expect(H.get(r, c)).toBeCloseTo(HFD[r][c], 4);
                    }
                }
            }
        }
    });

    it('computes a right-handed orthonormal frame whose normal is the gradient', () => {
        const surfaces: Array<{ surface: ImplicitSurface3, points: Vector[] }> = [
            {
                surface: new Sphere(2),
                points: [V3(2, 0, 0), V3(0, 0, 2), V3(1, 1, Math.sqrt(2))]
            },
            {
                surface: new Ellipsoid(1, 2, 3),
                points: [V3(1, 0, 0), V3(0, 2, 0), V3(0, 0, 3)]
            }
        ];

        for (const { surface, points } of surfaces) {
            for (const p of points) {
                const { tangent0, tangent1, normal } = surface.getFrame(p);

                // Orthonormal.
                expect(length(tangent0)).toBeCloseTo(1, 12);
                expect(length(tangent1)).toBeCloseTo(1, 12);
                expect(length(normal)).toBeCloseTo(1, 12);
                expect(dot(tangent0, tangent1)).toBeCloseTo(0, 12);
                expect(dot(tangent0, normal)).toBeCloseTo(0, 12);
                expect(dot(tangent1, normal)).toBeCloseTo(0, 12);

                // Right-handed: tangent0 x tangent1 = normal.
                const c = cross(tangent0, tangent1);
                expect(length(sub(c, normal))).toBeCloseTo(0, 12);

                // The normal is the normalized gradient.
                const g = surface.getGradient(p);
                const gLength = length(g);
                for (let i = 0; i < 3; ++i) {
                    expect(normal.values[i]).toBeCloseTo(g.values[i] / gLength, 12);
                }
            }
        }
    });

    it('computes principal curvatures 1/r everywhere on a sphere of radius r', () => {
        for (const r of [0.5, 1, 2, 7.25]) {
            const sphere = new Sphere(r);
            // Points on the sphere from spherical coordinates.
            for (let i = 1; i <= 4; ++i) {
                for (let j = 0; j < 5; ++j) {
                    const phi = Math.PI * i / 5;
                    const theta = 2 * Math.PI * j / 5;
                    const p = V3(
                        r * Math.sin(phi) * Math.cos(theta),
                        r * Math.sin(phi) * Math.sin(theta),
                        r * Math.cos(phi));
                    expect(sphere.isOnSurface(p, 1e-10)).toBe(true);

                    const info = sphere.getPrincipalInformation(p);
                    expect(info.valid).toBe(true);
                    expect(info.curvature0).toBeCloseTo(1 / r, 10);
                    expect(info.curvature1).toBeCloseTo(1 / r, 10);

                    // Every direction is principal; the returned pair must
                    // still be an orthonormal set in the tangent plane.
                    expect(length(info.direction0)).toBeCloseTo(1, 10);
                    expect(length(info.direction1)).toBeCloseTo(1, 10);
                    expect(dot(info.direction0, info.direction1)).toBeCloseTo(0, 10);
                    const n = sphere.getFrame(p).normal;
                    expect(dot(info.direction0, n)).toBeCloseTo(0, 10);
                    expect(dot(info.direction1, n)).toBeCloseTo(0, 10);
                }
            }
        }
    });

    it('matches the analytic ellipsoid curvatures at the axis points', () => {
        const a = 1, b = 2, c = 3;
        const ellipsoid = new Ellipsoid(a, b, c);

        // At (a,0,0) the principal curvatures are a/b^2 and a/c^2.
        let info = ellipsoid.getPrincipalInformation(V3(a, 0, 0));
        expect(info.valid).toBe(true);
        const kx = [a / (c * c), a / (b * b)].sort((p, q) => p - q);
        expect(info.curvature0).toBeCloseTo(kx[0], 12);
        expect(info.curvature1).toBeCloseTo(kx[1], 12);

        // At (0,b,0) they are b/a^2 and b/c^2.
        info = ellipsoid.getPrincipalInformation(V3(0, b, 0));
        const ky = [b / (a * a), b / (c * c)].sort((p, q) => p - q);
        expect(info.curvature0).toBeCloseTo(ky[0], 12);
        expect(info.curvature1).toBeCloseTo(ky[1], 12);

        // At (0,0,c) they are c/a^2 and c/b^2.
        info = ellipsoid.getPrincipalInformation(V3(0, 0, c));
        const kz = [c / (a * a), c / (b * b)].sort((p, q) => p - q);
        expect(info.curvature0).toBeCloseTo(kz[0], 12);
        expect(info.curvature1).toBeCloseTo(kz[1], 12);
    });

    it('returns curvatures sorted increasingly with orthonormal directions', () => {
        const ellipsoid = new Ellipsoid(1, 2, 3);
        const points = [
            V3(1, 0, 0), V3(0, 2, 0), V3(0, 0, 3),
            V3(Math.SQRT1_2, Math.SQRT2, 0),
            V3(0.5, Math.sqrt(3), 0)
        ];
        for (const p of points) {
            const info = ellipsoid.getPrincipalInformation(p);
            expect(info.valid).toBe(true);
            expect(info.curvature0).toBeLessThanOrEqual(info.curvature1);
            expect(length(info.direction0)).toBeCloseTo(1, 10);
            expect(length(info.direction1)).toBeCloseTo(1, 10);
            expect(dot(info.direction0, info.direction1)).toBeCloseTo(0, 10);
            const n = ellipsoid.getFrame(p).normal;
            expect(dot(info.direction0, n)).toBeCloseTo(0, 10);
            expect(dot(info.direction1, n)).toBeCloseTo(0, 10);
        }
    });

    it('detects the special point where the gradient is zero', () => {
        const sphere = new Sphere(2);
        const info = sphere.getPrincipalInformation(V3(0, 0, 0));
        expect(info.valid).toBe(false);
        expect(info.curvature0).toBe(0);
        expect(info.curvature1).toBe(0);
        expect(info.direction0.values).toEqual([0, 0, 0]);
        expect(info.direction1.values).toEqual([0, 0, 0]);
    });

    it('has mean and Gaussian curvature consistent with the shape operator', () => {
        // For an implicit surface, k0 + k1 is the trace of the shape
        // operator and k0*k1 is its determinant. Cross-check against the
        // classical formulas
        //   H = (grad F . Hess F . grad F - |grad F|^2 trace(Hess F))
        //       / (2 |grad F|^3)
        //   K = (grad F . adj(Hess F) . grad F) / |grad F|^4
        // written for the sign convention of the port (positive curvature
        // for a sphere with an outward gradient).
        const ellipsoid = new Ellipsoid(1, 2, 3);
        const points = [
            V3(1, 0, 0), V3(0, 2, 0), V3(0, 0, 3),
            V3(Math.SQRT1_2, Math.SQRT2, 0)
        ];
        for (const p of points) {
            const info = ellipsoid.getPrincipalInformation(p);
            const g = ellipsoid.getGradient(p);
            const H = ellipsoid.getHessian(p);
            const gLen = length(g);

            // gT H g
            let gHg = 0;
            for (let r = 0; r < 3; ++r) {
                for (let c = 0; c < 3; ++c) {
                    gHg += g.values[r] * H.get(r, c) * g.values[c];
                }
            }
            const traceH = H.get(0, 0) + H.get(1, 1) + H.get(2, 2);
            const twoH = (gLen * gLen * traceH - gHg) / (gLen * gLen * gLen);
            expect(info.curvature0 + info.curvature1).toBeCloseTo(twoH, 10);

            // Gaussian curvature via the adjugate of the Hessian.
            const m = (r: number, c: number) => H.get(r, c);
            const adj = [
                [m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1),
                    m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2),
                    m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)],
                [m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2),
                    m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0),
                    m(0, 2) * m(1, 0) - m(0, 0) * m(1, 2)],
                [m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0),
                    m(0, 1) * m(2, 0) - m(0, 0) * m(2, 1),
                    m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0)]
            ];
            let gAg = 0;
            for (let r = 0; r < 3; ++r) {
                for (let c = 0; c < 3; ++c) {
                    gAg += g.values[r] * adj[r][c] * g.values[c];
                }
            }
            const K = gAg / Math.pow(gLen, 4);
            expect(info.curvature0 * info.curvature1).toBeCloseTo(K, 10);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification block (V16): properties that would catch a translation error
// in the port of GetGradient/GetHessian/GetFrame/GetPrincipalInformation.
// ---------------------------------------------------------------------------

// A general quadric F(p) = p^T M p + b . p + c with M symmetric, so the
// gradient (2 M p + b) and the Hessian (2 M) exercise every off-diagonal
// entry of the port's Hessian assembly.
class Quadric extends ImplicitSurface3 {
    constructor(readonly M: number[][], readonly b: number[],
        readonly c: number) {
        super();
    }

    f(p: Vector): number {
        const v = p.values;
        let s = this.c;
        for (let i = 0; i < 3; ++i) {
            s += this.b[i] * v[i];
            for (let j = 0; j < 3; ++j) {
                s += v[i] * this.M[i][j] * v[j];
            }
        }
        return s;
    }

    private g(p: Vector, i: number): number {
        const v = p.values;
        let s = this.b[i];
        for (let j = 0; j < 3; ++j) {
            s += 2 * this.M[i][j] * v[j];
        }
        return s;
    }

    fx(p: Vector): number { return this.g(p, 0); }
    fy(p: Vector): number { return this.g(p, 1); }
    fz(p: Vector): number { return this.g(p, 2); }
    fxx(_p: Vector): number { return 2 * this.M[0][0]; }
    fxy(_p: Vector): number { return 2 * this.M[0][1]; }
    fxz(_p: Vector): number { return 2 * this.M[0][2]; }
    fyy(_p: Vector): number { return 2 * this.M[1][1]; }
    fyz(_p: Vector): number { return 2 * this.M[1][2]; }
    fzz(_p: Vector): number { return 2 * this.M[2][2]; }
}

// A symmetric 3x3 matrix, a linear term and a constant, all at a moderate
// scale. wellScaled snaps tiny magnitudes to zero, so no entry is subnormal;
// the algorithm normalizes the gradient, so this is enough conditioning.
const quadricArb = fc.tuple(
    fc.array(wellScaled(-3, 3), { minLength: 6, maxLength: 6 }),
    fc.array(wellScaled(-3, 3), { minLength: 3, maxLength: 3 }),
    wellScaled(-3, 3))
    .map(([m, b, c]) => {
        const M = [
            [m[0], m[1], m[2]],
            [m[1], m[3], m[4]],
            [m[2], m[4], m[5]]
        ];
        return new Quadric(M, b, c);
    });

// A quadric together with a point whose gradient is comfortably nonzero.
const quadricAndPoint = fc.tuple(quadricArb, wellScaledVector(3, -3, 3))
    .filter(([q, p]) => length(q.getGradient(p)) > 1e-2);

function matVec3(M: Matrix, v: Vector): Vector {
    const r = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        let s = 0;
        for (let j = 0; j < 3; ++j) { s += M.get(i, j) * v.values[j]; }
        r.values[i] = s;
    }
    return r;
}

function scaleVector(v: Vector, s: number): Vector {
    return Vector.fromArray(v.values.map(x => x * s));
}

describe('ImplicitSurface3 verification', () => {
    it('assembles the gradient and the symmetric Hessian from the partials', () => {
        check(quadricAndPoint, ([q, p]) => {
            const g = q.getGradient(p);
            expect(g.size).toBe(3);
            expect(g.values[0]).toBe(q.fx(p));
            expect(g.values[1]).toBe(q.fy(p));
            expect(g.values[2]).toBe(q.fz(p));

            const H = q.getHessian(p);
            expect(H.numRows).toBe(3);
            expect(H.numCols).toBe(3);
            expect(H.get(0, 0)).toBe(q.fxx(p));
            expect(H.get(1, 1)).toBe(q.fyy(p));
            expect(H.get(2, 2)).toBe(q.fzz(p));
            // The upstream Matrix3x3 initializer list writes fxy, fxz and fyz
            // twice, once in each triangle.
            expect(H.get(0, 1)).toBe(q.fxy(p));
            expect(H.get(1, 0)).toBe(q.fxy(p));
            expect(H.get(0, 2)).toBe(q.fxz(p));
            expect(H.get(2, 0)).toBe(q.fxz(p));
            expect(H.get(1, 2)).toBe(q.fyz(p));
            expect(H.get(2, 1)).toBe(q.fyz(p));
        });
    });

    it('returns a right-handed orthonormal frame whose normal is the unit gradient', () => {
        check(quadricAndPoint, ([q, p]) => {
            const { tangent0, tangent1, normal } = q.getFrame(p);
            expectClose(length(tangent0), 1, 1e-12, 1e-12);
            expectClose(length(tangent1), 1, 1e-12, 1e-12);
            expectClose(length(normal), 1, 1e-12, 1e-12);
            expectClose(dot(tangent0, tangent1), 0, 1e-12, 1e-12);
            expectClose(dot(tangent0, normal), 0, 1e-12, 1e-12);
            expectClose(dot(tangent1, normal), 0, 1e-12, 1e-12);
            expectVectorClose(cross(tangent0, tangent1), normal, 1e-12, 1e-12);

            const g = q.getGradient(p);
            expectVectorClose(normal, scaleVector(g, 1 / length(g)),
                1e-12, 1e-12);
        });
    });

    it('satisfies the shape-operator eigen-equation of the referenced PDF', () => {
        // Upstream solves the eigensystem of barA = J^T A J with A = H/|grad|
        // and J = [U | V] the tangent frame, then returns direction_i = J*w_i.
        // Because J J^T = I - N N^T, that is equivalent to: the tangential
        // part of A*direction_i equals curvature_i * direction_i. The identity
        // pins the matrix products, the eigenvalue ordering and the J*w
        // back-substitution at once.
        check(quadricAndPoint, ([q, p]) => {
            const info = q.getPrincipalInformation(p);
            expect(info.valid).toBe(true);
            expect(info.curvature0).toBeLessThanOrEqual(info.curvature1);

            const g = q.getGradient(p);
            const gLen = length(g);
            const N = scaleVector(g, 1 / gLen);
            const H = q.getHessian(p);
            const A = Matrix.fromArray(3, 3, H.values.map(x => x / gLen));

            const pairs: Array<[number, Vector]> = [
                [info.curvature0, info.direction0],
                [info.curvature1, info.direction1]
            ];
            for (const [k, d] of pairs) {
                expectClose(length(d), 1, 1e-9, 1e-9);
                expectClose(dot(d, N), 0, 1e-9, 1e-9);

                const Ad = matVec3(A, d);
                const tangential = sub(Ad, scaleVector(N, dot(N, Ad)));
                // The residual is relative to the scale of A = H/|grad|, which
                // the curvature bounds from below; both are moderate here.
                const scale = Math.max(1, Math.abs(k));
                expectVectorClose(tangential, scaleVector(d, k),
                    1e-8 * scale, 1e-8);
            }
            expectClose(dot(info.direction0, info.direction1), 0, 1e-8, 1e-8);
        });
    });

    it('matches the trace and determinant of the tangential Hessian', () => {
        // Sum and product of the eigenvalues of barA = J^T A J. Because
        // J J^T = I - N N^T, trace(barA) = trace(A) - N^T A N.
        check(quadricAndPoint, ([q, p]) => {
            const info = q.getPrincipalInformation(p);
            const g = q.getGradient(p);
            const gLen = length(g);
            const N = scaleVector(g, 1 / gLen);
            const H = q.getHessian(p);
            const A = Matrix.fromArray(3, 3, H.values.map(x => x / gLen));

            const traceA = A.get(0, 0) + A.get(1, 1) + A.get(2, 2);
            const nAn = dot(N, matVec3(A, N));
            expectClose(info.curvature0 + info.curvature1, traceA - nAn,
                1e-8, 1e-8);

            // The determinant, computed from the same frame the port uses.
            const { tangent0: U, tangent1: V } = q.getFrame(p);
            const s00 = dot(U, matVec3(A, U));
            const s01 = dot(U, matVec3(A, V));
            const s10 = dot(V, matVec3(A, U));
            const s11 = dot(V, matVec3(A, V));
            const avr = 0.5 * (s01 + s10);
            expectClose(info.curvature0 * info.curvature1,
                s00 * s11 - avr * avr, 1e-8, 1e-8);
        });
    });

    it('is invariant under a rigid motion of the surface', () => {
        // F'(p) = F(R^T (p - t)) describes the rigidly moved surface, so the
        // principal curvatures of F' at R p + t equal those of F at p.
        check(fc.tuple(quadricAndPoint, rotationFrame(3),
            wellScaledVector(3, -2, 2)), ([[q, p], R, t]) => {
            // R has the frame vectors as its columns.
            const rot = (x: Vector): Vector => {
                const r = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    r.values[i] = R[0].values[i] * x.values[0]
                        + R[1].values[i] * x.values[1]
                        + R[2].values[i] * x.values[2];
                }
                return r;
            };
            // M' = R M R^T, b' = R b - (M' + M'^T) t,
            // c' = c - (R b) . t + t^T M' t.
            const Mp: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
            for (let i = 0; i < 3; ++i) {
                for (let j = 0; j < 3; ++j) {
                    let s = 0;
                    for (let a = 0; a < 3; ++a) {
                        for (let b2 = 0; b2 < 3; ++b2) {
                            s += R[a].values[i] * q.M[a][b2] * R[b2].values[j];
                        }
                    }
                    Mp[i][j] = s;
                }
            }
            const Rb = rot(Vector.fromArray(q.b));
            const bp: number[] = [0, 0, 0];
            for (let i = 0; i < 3; ++i) {
                let s = Rb.values[i];
                for (let j = 0; j < 3; ++j) {
                    s -= (Mp[i][j] + Mp[j][i]) * t.values[j];
                }
                bp[i] = s;
            }
            let cp = q.c;
            for (let i = 0; i < 3; ++i) {
                cp -= Rb.values[i] * t.values[i];
                for (let j = 0; j < 3; ++j) {
                    cp += t.values[i] * Mp[i][j] * t.values[j];
                }
            }
            const moved = new Quadric(Mp, bp, cp);
            const pp = add(rot(p), t);
            // Sanity: the moved surface takes the same value at the moved
            // point. The rotation and translation lose a few digits.
            expectClose(moved.f(pp), q.f(p), 1e-6, 1e-8);

            const a = q.getPrincipalInformation(p);
            const b = moved.getPrincipalInformation(pp);
            expect(b.valid).toBe(true);
            expectClose(a.curvature0, b.curvature0, 1e-7, 1e-7);
            expectClose(a.curvature1, b.curvature1, 1e-7, 1e-7);
        });
    });

    it('reports invalid with zero outputs wherever the gradient vanishes', () => {
        // The gradient 2 M p + b vanishes at p when b = -2 M p, so build that
        // quadric rather than inverting M.
        check(fc.tuple(quadricArb, wellScaledVector(3, -3, 3)), ([q, p]) => {
            const b: number[] = [0, 0, 0];
            for (let i = 0; i < 3; ++i) {
                let s = 0;
                for (let j = 0; j < 3; ++j) { s += 2 * q.M[i][j] * p.values[j]; }
                b[i] = -s;
            }
            const critical = new Quadric(q.M, b, q.c);
            expect(length(critical.getGradient(p))).toBe(0);
            const info = critical.getPrincipalInformation(p);
            expect(info.valid).toBe(false);
            expect(info.curvature0).toBe(0);
            expect(info.curvature1).toBe(0);
            expect(info.direction0.values).toEqual([0, 0, 0]);
            expect(info.direction1.values).toEqual([0, 0, 0]);
        });
    });

    it('agrees with finite differences for the gradient and the Hessian', () => {
        check(quadricAndPoint, ([q, p]) => {
            const h = 1e-4;
            expectVectorClose(q.getGradient(p), fdGradient(q, p, h), 1e-6, 1e-6);
            const H = q.getHessian(p);
            const fd = fdHessian(q, p, h);
            for (let i = 0; i < 3; ++i) {
                for (let j = 0; j < 3; ++j) {
                    expectClose(H.get(i, j), fd[i][j], 1e-4, 1e-4);
                }
            }
        });
    });

    it('is unchanged when the implicit function is scaled by a nonzero factor', () => {
        // s*F = 0 defines the same surface and the port divides the Hessian by
        // |grad|, so the curvatures pick up sign(s) only (and swap when the
        // sign flips, because they are returned in increasing order).
        check(fc.tuple(quadricAndPoint, nonzero(-5, 5, 0.1)), ([[q, p], s]) => {
            const scaled = new Quadric(q.M.map(row => row.map(x => x * s)),
                q.b.map(x => x * s), q.c * s);
            const a = q.getPrincipalInformation(p);
            const b = scaled.getPrincipalInformation(p);
            expect(b.valid).toBe(true);
            if (s > 0) {
                expectClose(a.curvature0, b.curvature0, 1e-9, 1e-9);
                expectClose(a.curvature1, b.curvature1, 1e-9, 1e-9);
            } else {
                expectClose(a.curvature0, -b.curvature1, 1e-9, 1e-9);
                expectClose(a.curvature1, -b.curvature0, 1e-9, 1e-9);
            }
        });
    });
});
