import { describe, it, expect } from 'vitest';
import { ImplicitSurface3 } from '../src/ImplicitSurface3';
import { Vector, dot, length, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

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
