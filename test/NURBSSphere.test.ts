import { describe, it, expect } from 'vitest';
import {
    NURBSEighthSphereDegree4, NURBSHalfSphereDegree3, NURBSFullSphereDegree3,
    createEighthSphereValues
} from '../src/NURBSSphere.js';
import { NURBSSurface } from '../src/NURBSSurface.js';
import { ParametricSurface } from '../src/ParametricSurface.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function positionOf(surface: NURBSSurface, u: number, v: number): Vector {
    const jet = surface.createJet();
    surface.evaluate(u, v, 0, jet);
    return jet[0];
}

function norm(v: Vector): number {
    return Math.hypot(v.values[0], v.values[1], v.values[2]);
}

describe('NURBSEighthSphereDegree4', () => {
    const sphere = new NURBSEighthSphereDegree4();

    it('allocates six 3D output values', () => {
        const values = createEighthSphereValues();
        expect(values.length).toBe(6);
        for (const v of values) {
            expect(v.size).toBe(3);
            expect(v.values).toEqual([0, 0, 0]);
        }
    });

    it('interpolates the three corners of the triangular domain', () => {
        const values = createEighthSphereValues();
        sphere.evaluate(0, 0, 0, values);
        expect(values[0].values[0]).toBeCloseTo(0, 12);
        expect(values[0].values[1]).toBeCloseTo(0, 12);
        expect(values[0].values[2]).toBeCloseTo(1, 12);

        sphere.evaluate(1, 0, 0, values);
        expect(values[0].values[0]).toBeCloseTo(1, 12);
        expect(values[0].values[1]).toBeCloseTo(0, 12);
        expect(values[0].values[2]).toBeCloseTo(0, 12);

        sphere.evaluate(0, 1, 0, values);
        expect(values[0].values[0]).toBeCloseTo(0, 12);
        expect(values[0].values[1]).toBeCloseTo(1, 12);
        expect(values[0].values[2]).toBeCloseTo(0, 12);
    });

    it('produces points very near the unit sphere in the first octant', () => {
        const values = createEighthSphereValues();
        let maxError = 0;
        for (let i = 0; i <= 10; ++i) {
            for (let j = 0; i + j <= 10; ++j) {
                const u = i / 10, v = j / 10;
                sphere.evaluate(u, v, 0, values);
                const p = values[0];
                maxError = Math.max(maxError, Math.abs(norm(p) - 1));
                expect(p.values[0]).toBeGreaterThanOrEqual(-1e-12);
                expect(p.values[1]).toBeGreaterThanOrEqual(-1e-12);
                expect(p.values[2]).toBeGreaterThanOrEqual(-1e-12);
            }
        }
        // The degree-4 rational triangular patch represents the sphere
        // octant EXACTLY (section 3.1.2 of the upstream document); the
        // residual is pure floating-point rounding. The original bound of
        // 1e-3 here was far weaker than the construction warrants.
        expect(maxError).toBeLessThan(1e-14);
        expect(maxError).toBeGreaterThan(0);
    });

    it('is symmetric under swapping u and v (x <-> y)', () => {
        const a = createEighthSphereValues();
        const b = createEighthSphereValues();
        for (const [u, v] of [[0.1, 0.3], [0.25, 0.5], [0.4, 0.15]]) {
            sphere.evaluate(u, v, 0, a);
            sphere.evaluate(v, u, 0, b);
            expect(a[0].values[0]).toBeCloseTo(b[0].values[1], 10);
            expect(a[0].values[1]).toBeCloseTo(b[0].values[0], 10);
            expect(a[0].values[2]).toBeCloseTo(b[0].values[2], 10);
        }
    });

    it('computes first derivatives that match finite differences', () => {
        const values = createEighthSphereValues();
        const p = createEighthSphereValues();
        const h = 1e-6;
        for (const [u, v] of [[0.2, 0.3], [0.4, 0.25], [0.15, 0.55]]) {
            sphere.evaluate(u, v, 1, values);
            const du: number[] = [];
            const dv: number[] = [];
            sphere.evaluate(u + h, v, 0, p);
            const xu1 = p[0].values.slice();
            sphere.evaluate(u - h, v, 0, p);
            const xu0 = p[0].values.slice();
            sphere.evaluate(u, v + h, 0, p);
            const xv1 = p[0].values.slice();
            sphere.evaluate(u, v - h, 0, p);
            const xv0 = p[0].values.slice();
            for (let k = 0; k < 3; ++k) {
                du.push((xu1[k] - xu0[k]) / (2 * h));
                dv.push((xv1[k] - xv0[k]) / (2 * h));
            }
            for (let k = 0; k < 3; ++k) {
                expect(values[1].values[k]).toBeCloseTo(du[k], 5);
                expect(values[2].values[k]).toBeCloseTo(dv[k], 5);
            }
        }
    });

    it('computes second derivatives that match finite differences', () => {
        const values = createEighthSphereValues();
        const p = createEighthSphereValues();
        const h = 1e-4;
        const u = 0.3, v = 0.25;
        sphere.evaluate(u, v, 2, values);
        const at = (a: number, b: number): number[] => {
            sphere.evaluate(a, b, 0, p);
            return p[0].values.slice();
        };
        const c = at(u, v);
        const up = at(u + h, v), um = at(u - h, v);
        const vp = at(u, v + h), vm = at(u, v - h);
        const pp = at(u + h, v + h), pm = at(u + h, v - h);
        const mp = at(u - h, v + h), mm = at(u - h, v - h);
        for (let k = 0; k < 3; ++k) {
            const duu = (up[k] - 2 * c[k] + um[k]) / (h * h);
            const dvv = (vp[k] - 2 * c[k] + vm[k]) / (h * h);
            const duv = (pp[k] - pm[k] - mp[k] + mm[k]) / (4 * h * h);
            expect(values[3].values[k]).toBeCloseTo(duu, 2);
            expect(values[4].values[k]).toBeCloseTo(duv, 2);
            expect(values[5].values[k]).toBeCloseTo(dvv, 2);
        }
    });

    it('exposes its controls and weights', () => {
        const controls = sphere.getControls();
        const weights = sphere.getWeights();
        expect(controls.length).toBe(5);
        expect(weights.length).toBe(5);
        // P004 = (0,0,1) with weight 4*sqrt(3)*(sqrt(3)-1).
        expect(controls[0][0].values).toEqual([0, 0, 1]);
        expect(weights[0][0]).toBeCloseTo(
            4 * Math.sqrt(3) * (Math.sqrt(3) - 1), 12);
        // P400 = (1,0,0) with the same weight.
        expect(controls[4][0].values).toEqual([1, 0, 0]);
        expect(weights[4][0]).toBeCloseTo(weights[0][0], 12);
        // w022 = 4.
        expect(weights[0][2]).toBe(4);
    });
});

describe('NURBSHalfSphereDegree3', () => {
    const sphere = new NURBSHalfSphereDegree3();

    it('is a 3D NURBS surface on the unit square', () => {
        expect(sphere).toBeInstanceOf(NURBSSurface);
        expect(sphere.isConstructed()).toBe(true);
        expect(sphere.getDimension()).toBe(3);
        expect(sphere.getNumControls(0)).toBe(4);
        expect(sphere.getNumControls(1)).toBe(4);
        expect(sphere.getUMin()).toBe(0);
        expect(sphere.getUMax()).toBe(1);
        expect(sphere.getVMin()).toBe(0);
        expect(sphere.getVMax()).toBe(1);
    });

    it('places every point exactly on the unit sphere', () => {
        for (let i = 0; i <= 12; ++i) {
            for (let j = 0; j <= 12; ++j) {
                const p = positionOf(sphere, i / 12, j / 12);
                expect(norm(p)).toBeCloseTo(1, 10);
            }
        }
    });

    it('collapses the v = 0 and v = 1 seams to the poles', () => {
        for (let i = 0; i <= 8; ++i) {
            const north = positionOf(sphere, i / 8, 0);
            expect(north.values[0]).toBeCloseTo(0, 12);
            expect(north.values[1]).toBeCloseTo(0, 12);
            expect(north.values[2]).toBeCloseTo(1, 12);
            const south = positionOf(sphere, i / 8, 1);
            expect(south.values[0]).toBeCloseTo(0, 12);
            expect(south.values[1]).toBeCloseTo(0, 12);
            expect(south.values[2]).toBeCloseTo(-1, 12);
        }
    });

    it('covers the half sphere y >= 0', () => {
        let maxY = -1;
        for (let i = 0; i <= 10; ++i) {
            for (let j = 0; j <= 10; ++j) {
                const p = positionOf(sphere, i / 10, j / 10);
                expect(p.values[1]).toBeGreaterThan(-1e-10);
                maxY = Math.max(maxY, p.values[1]);
            }
        }
        expect(maxY).toBeCloseTo(1, 8);
    });

    it('passes through the equator point (1,0,0) at u=0, v=1/2', () => {
        const p = positionOf(sphere, 0, 0.5);
        expect(p.values[0]).toBeCloseTo(1, 12);
        expect(p.values[1]).toBeCloseTo(0, 12);
        expect(p.values[2]).toBeCloseTo(0, 12);
        const q = positionOf(sphere, 1, 0.5);
        expect(q.values[0]).toBeCloseTo(-1, 12);
        expect(q.values[1]).toBeCloseTo(0, 12);
        expect(q.values[2]).toBeCloseTo(0, 12);
    });

    it('has a surface normal parallel to the position (a sphere)', () => {
        const jet = sphere.createJet();
        for (const [u, v] of [[0.3, 0.4], [0.7, 0.25], [0.5, 0.8]]) {
            sphere.evaluate(u, v, 1, jet);
            const X = jet[0], Xu = jet[1], Xv = jet[2];
            // The tangent vectors are perpendicular to the radius vector.
            let dotU = 0, dotV = 0;
            for (let k = 0; k < 3; ++k) {
                dotU += X.values[k] * Xu.values[k];
                dotV += X.values[k] * Xv.values[k];
            }
            expect(dotU).toBeCloseTo(0, 9);
            expect(dotV).toBeCloseTo(0, 9);
        }
    });
});

describe('NURBSFullSphereDegree3', () => {
    const sphere = new NURBSFullSphereDegree3();

    it('is a 3D NURBS surface with 4-by-7 control points', () => {
        expect(sphere).toBeInstanceOf(ParametricSurface);
        expect(sphere.isConstructed()).toBe(true);
        expect(sphere.getNumControls(0)).toBe(4);
        expect(sphere.getNumControls(1)).toBe(7);
    });

    it('places every point exactly on the unit sphere', () => {
        for (let i = 0; i <= 12; ++i) {
            for (let j = 0; j <= 12; ++j) {
                const p = positionOf(sphere, i / 12, j / 12);
                expect(norm(p)).toBeCloseTo(1, 10);
            }
        }
    });

    it('collapses the u = 0 and u = 1 seams to the poles', () => {
        for (let j = 0; j <= 8; ++j) {
            const north = positionOf(sphere, 0, j / 8);
            expect(north.values[0]).toBeCloseTo(0, 12);
            expect(north.values[1]).toBeCloseTo(0, 12);
            expect(north.values[2]).toBeCloseTo(1, 12);
            const south = positionOf(sphere, 1, j / 8);
            expect(south.values[0]).toBeCloseTo(0, 12);
            expect(south.values[1]).toBeCloseTo(0, 12);
            expect(south.values[2]).toBeCloseTo(-1, 12);
        }
    });

    it('closes up in v: X(u,0) equals X(u,1)', () => {
        for (let i = 1; i < 8; ++i) {
            const u = i / 8;
            const a = positionOf(sphere, u, 0);
            const b = positionOf(sphere, u, 1);
            for (let k = 0; k < 3; ++k) {
                expect(a.values[k]).toBeCloseTo(b.values[k], 10);
            }
        }
    });

    it('reaches both signs of y (a full sphere, unlike the half sphere)', () => {
        let minY = 1, maxY = -1;
        for (let i = 0; i <= 12; ++i) {
            for (let j = 0; j <= 12; ++j) {
                const p = positionOf(sphere, i / 12, j / 12);
                minY = Math.min(minY, p.values[1]);
                maxY = Math.max(maxY, p.values[1]);
            }
        }
        expect(maxY).toBeGreaterThan(0.9);
        expect(minY).toBeLessThan(-0.9);
    });

    it('agrees with the half sphere on the v in [0,1/2] half', () => {
        // The two surfaces transpose the roles of the parameters: for the
        // full sphere u is latitude and v is longitude, while for the half
        // sphere u is longitude and v is latitude. The full sphere's
        // v in [0,1/2] is the same Bezier patch (knots 0,0,0,0,1/2,1/2,1/2)
        // as the half sphere's whole longitude range.
        const half = new NURBSHalfSphereDegree3();
        for (let i = 0; i <= 8; ++i) {
            for (let j = 0; j <= 8; ++j) {
                const latitude = i / 8, longitude = j / 8;
                const a = positionOf(sphere, latitude, longitude / 2);
                const b = positionOf(half, longitude, latitude);
                for (let k = 0; k < 3; ++k) {
                    expect(a.values[k]).toBeCloseTo(b.values[k], 9);
                }
            }
        }
    });
});


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). NURBSSphere.h was read line by
// line against src/NURBSSphere.ts: every control point, every weight and every
// Bernstein/derivative table entry of NURBSEighthSphereDegree4 was compared
// term by term. The properties below check the constructions against the
// geometry they claim rather than against the formulas: the patches lie on the
// unit sphere, the triangular patch degenerates to the correct coordinate
// planes on its three edges, and the derivative tables agree with numerical
// differentiation of the position.
import {
    check, fc, expectClose
} from './helpers/arbitraries.js';

// Barycentric-like parameters of the triangular domain u >= 0, v >= 0,
// u + v <= 1.
const triangleParam = fc.tuple(
    fc.double({ min: 0, max: 1, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0, max: 1, noNaN: true, noDefaultInfinity: true }))
    .map(([a, b]) => (a + b <= 1 ? [a, b] : [1 - a, 1 - b]));

const square = fc.tuple(
    fc.double({ min: 0, max: 1, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0, max: 1, noNaN: true, noDefaultInfinity: true }));

describe('NURBSSphere verification', () => {
    const eighth = new NURBSEighthSphereDegree4();

    it('the degree-4 triangular patch is exactly the unit sphere octant', () => {
        check(triangleParam, ([u, v]) => {
            const values = createEighthSphereValues();
            eighth.evaluate(u, v, 0, values);
            const p = values[0];
            expectClose(norm(p), 1, 1e-13, 1e-13);
            for (let k = 0; k < 3; ++k) {
                expect(p.values[k]).toBeGreaterThanOrEqual(-1e-13);
                expect(p.values[k]).toBeLessThanOrEqual(1 + 1e-13);
            }
        });
    });

    it('the three edges of the triangular patch lie in the coordinate planes', () => {
        // On v = 0 only the j0 = 0 controls contribute and all of them have
        // y = 0; on u = 0 only the j1 = 0 controls contribute and all have
        // x = 0; on u + v = 1 only the controls with j0 + j1 = 4 contribute
        // and all have z = 0. Each edge is therefore an exact quarter circle.
        check(fc.double({ min: 0, max: 1, noNaN: true,
            noDefaultInfinity: true }), t => {
            const values = createEighthSphereValues();
            eighth.evaluate(t, 0, 0, values);
            expect(values[0].values[1]).toBe(0);
            expectClose(norm(values[0]), 1, 1e-13, 1e-13);

            eighth.evaluate(0, t, 0, values);
            expect(values[0].values[0]).toBe(0);
            expectClose(norm(values[0]), 1, 1e-13, 1e-13);

            eighth.evaluate(t, 1 - t, 0, values);
            expect(values[0].values[2]).toBe(0);
            expectClose(norm(values[0]), 1, 1e-13, 1e-13);
        });
    });

    it('the triangular patch is symmetric in u and v', () => {
        // The control net is symmetric under the reflection that swaps x and
        // y, so X(v,u) is X(u,v) with the first two components swapped.
        check(triangleParam, ([u, v]) => {
            const a = createEighthSphereValues();
            const b = createEighthSphereValues();
            eighth.evaluate(u, v, 0, a);
            eighth.evaluate(v, u, 0, b);
            expectClose(a[0].values[0], b[0].values[1], 1e-14, 1e-14);
            expectClose(a[0].values[1], b[0].values[0], 1e-14, 1e-14);
            expectClose(a[0].values[2], b[0].values[2], 1e-14, 1e-14);
        });
    });

    it('the triangular patch derivative tables agree with differencing', () => {
        // Five-point stencils at two step sizes; the tolerance is calibrated
        // from the observed convergence, which keeps the property meaningful
        // near the corners where the rational patch varies quickly.
        check(triangleParam, ([u, v]) => {
            // Keep the whole stencil strictly inside the triangle.
            const h = 5e-3;
            if (u < 3 * h || v < 3 * h || u + v > 1 - 3 * h) { return; }
            const values = createEighthSphereValues();
            eighth.evaluate(u, v, 2, values);
            const at = (a: number, b: number, order: number): number[] => {
                const p = createEighthSphereValues();
                eighth.evaluate(a, b, order, p);
                return p[order].values.slice();
            };
            const stencil = (step: number, k: number, du: boolean,
                order: number): number => {
                const p1 = du ? at(u + step, v, order) : at(u, v + step, order);
                const m1 = du ? at(u - step, v, order) : at(u, v - step, order);
                const p2 = du ? at(u + 2 * step, v, order)
                    : at(u, v + 2 * step, order);
                const m2 = du ? at(u - 2 * step, v, order)
                    : at(u, v - 2 * step, order);
                return (8 * (p1[k] - m1[k]) - (p2[k] - m2[k])) / (12 * step);
            };
            for (let k = 0; k < 3; ++k) {
                for (const [du, index] of [[true, 1], [false, 2]] as
                    [boolean, number][]) {
                    const coarse = stencil(h, k, du, 0);
                    const fine = stencil(0.5 * h, k, du, 0);
                    expectClose(fine, values[index].values[k],
                        1e-9 + Math.abs(coarse - fine), 1e-9);
                }
                // d2X/du2 and d2X/dudv from differencing dX/du; d2X/dv2 from
                // differencing dX/dv.
                for (const [du, index, order] of [[true, 3, 1], [false, 4, 1],
                    [false, 5, 2]] as [boolean, number, number][]) {
                    const coarse = stencil(h, k, du, order);
                    const fine = stencil(0.5 * h, k, du, order);
                    expectClose(fine, values[index].values[k],
                        1e-8 + Math.abs(coarse - fine), 1e-8);
                }
            }
        }, 40);
    });

    it('the degree-3 half and full spheres lie on the unit sphere', () => {
        const half = new NURBSHalfSphereDegree3();
        const full = new NURBSFullSphereDegree3();
        check(square, ([u, v]) => {
            for (const surface of [half, full]) {
                const p = positionOf(surface, u, v);
                expectClose(norm(p), 1, 1e-13, 1e-13);
            }
            // The half sphere covers y >= 0 only.
            expect(positionOf(half, u, v).values[1])
                .toBeGreaterThanOrEqual(-1e-13);
        });
    });

    it('the degree-3 sphere normals are parallel to the position', () => {
        // A sphere centered at the origin has an outward normal along the
        // position vector, so Cross(dX/du, dX/dv) must be parallel to X.
        const full = new NURBSFullSphereDegree3();
        check(square, ([u, v]) => {
            const jet = full.createJet();
            full.evaluate(u, v, 1, jet);
            const n = cross(jet[1], jet[2]);
            const len = norm(n);
            if (len < 1e-6) { return; }   // the poles are seam points
            const p = jet[0];
            for (let k = 0; k < 3; ++k) {
                expectClose(n.values[k] / len,
                    Math.sign(dot(n, p)) * p.values[k], 1e-8, 1e-8);
            }
        });
    });

    it('the full sphere reproduces the half sphere on its first half', () => {
        // The v knot vector of the full sphere has the interior knot 0.5 with
        // multiplicity 3, so v in [0,1/2] is the half sphere reparameterized
        // by v -> v/2.
        const half = new NURBSHalfSphereDegree3();
        const full = new NURBSFullSphereDegree3();
        // The two surfaces transpose the roles of the parameters: for the
        // full sphere u is latitude and v is longitude, while for the half
        // sphere u is longitude and v is latitude.
        check(square, ([latitude, longitude]) => {
            const a = positionOf(full, latitude, 0.5 * longitude);
            const b = positionOf(half, longitude, latitude);
            for (let k = 0; k < 3; ++k) {
                expectClose(a.values[k], b.values[k], 1e-12, 1e-12);
            }
        });
    });

    it('the degree-3 seams collapse to the poles', () => {
        const half = new NURBSHalfSphereDegree3();
        const full = new NURBSFullSphereDegree3();
        check(fc.double({ min: 0, max: 1, noNaN: true,
            noDefaultInfinity: true }), t => {
            // Half sphere: v = 0 and v = 1 are the two poles.
            expectClose(norm(sub(positionOf(half, t, 0),
                Vector.fromArray([0, 0, 1]))), 0, 1e-13, 1e-13);
            expectClose(norm(sub(positionOf(half, t, 1),
                Vector.fromArray([0, 0, -1]))), 0, 1e-13, 1e-13);
            // Full sphere: u = 0 and u = 1 are the poles, and v wraps.
            expectClose(norm(sub(positionOf(full, 0, t),
                Vector.fromArray([0, 0, 1]))), 0, 1e-13, 1e-13);
            expectClose(norm(sub(positionOf(full, 1, t),
                Vector.fromArray([0, 0, -1]))), 0, 1e-13, 1e-13);
            expectClose(norm(sub(positionOf(full, t, 0),
                positionOf(full, t, 1))), 0, 1e-12, 1e-12);
        });
    });
});
