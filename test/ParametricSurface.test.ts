import { describe, it, expect } from 'vitest';
import { ParametricSurface } from '../src/ParametricSurface';
import { Vector, dot, length, sub, div } from '../src/Vector';

// ---------------------------------------------------------------------------
// Concrete subclasses used to exercise the abstract base. The real surface
// classes (BSplineSurface, NURBSSurface, ...) arrive in later batches; these
// stand in for them and follow the same Evaluate contract: only the jet
// entries required by 'order' are written, and an unconstructed surface (or
// an out-of-range order) produces an all-zero jet.
// ---------------------------------------------------------------------------

// X(u,v) = P + u * U + v * V, an affine patch with exact derivatives.
class PlanePatch extends ParametricSurface {
    private readonly p: Vector;
    private readonly u: Vector;
    private readonly v: Vector;

    constructor(p: Vector, u: Vector, v: Vector, umin: number, umax: number,
        vmin: number, vmax: number, rectangular: boolean = true) {
        super(p.size, umin, umax, vmin, vmax, rectangular);
        this.p = p;
        this.u = u;
        this.v = v;
        this.mConstructed = true;
    }

    evaluate(u: number, v: number, order: number, jet: Vector[]): void {
        if (!this.mConstructed || order >= ParametricSurface.SUP_ORDER) {
            for (let i = 0; i < ParametricSurface.SUP_ORDER; ++i) {
                jet[i].makeZero();
            }
            return;
        }

        for (let i = 0; i < this.mDimension; ++i) {
            jet[0].set(i, this.p.get(i) + u * this.u.get(i) + v * this.v.get(i));
        }
        if (order >= 1) {
            // Deliberately assigns the stored direction vectors into the jet
            // (upstream C++ would copy them). The base class must not damage
            // them when it normalizes the tangents.
            jet[1] = this.u;
            jet[2] = this.v;
            if (order >= 2) {
                jet[3].makeZero();
                jet[4].makeZero();
                jet[5].makeZero();
            }
        }
    }
}

// The unit sphere patch X(u,v) = (cos(u) sin(v), sin(u) sin(v), cos(v)),
// u in [0, 2*pi] (longitude), v in [0, pi] (colatitude).
class SpherePatch extends ParametricSurface {
    constructor() {
        super(3, 0, 2 * Math.PI, 0, Math.PI, true);
        this.mConstructed = true;
    }

    evaluate(u: number, v: number, order: number, jet: Vector[]): void {
        if (!this.mConstructed || order >= ParametricSurface.SUP_ORDER) {
            for (let i = 0; i < ParametricSurface.SUP_ORDER; ++i) {
                jet[i].makeZero();
            }
            return;
        }

        const cu = Math.cos(u), su = Math.sin(u);
        const cv = Math.cos(v), sv = Math.sin(v);

        jet[0] = Vector.fromArray([cu * sv, su * sv, cv]);
        if (order >= 1) {
            jet[1] = Vector.fromArray([-su * sv, cu * sv, 0]);
            jet[2] = Vector.fromArray([cu * cv, su * cv, -sv]);
            if (order >= 2) {
                jet[3] = Vector.fromArray([-cu * sv, -su * sv, 0]);
                jet[4] = Vector.fromArray([-su * cv, cu * cv, 0]);
                jet[5] = Vector.fromArray([-cu * sv, -su * sv, -cv]);
            }
        }
    }
}

// A surface that is never marked as constructed.
class BrokenPatch extends ParametricSurface {
    constructor() {
        super(3, 0, 1, 0, 1, false);
    }

    evaluate(_u: number, _v: number, _order: number, jet: Vector[]): void {
        if (!this.mConstructed) {
            for (let i = 0; i < ParametricSurface.SUP_ORDER; ++i) {
                jet[i].makeZero();
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers.
// ---------------------------------------------------------------------------

function cross3(a: Vector, b: Vector): Vector {
    return Vector.fromArray([
        a.get(1) * b.get(2) - a.get(2) * b.get(1),
        a.get(2) * b.get(0) - a.get(0) * b.get(2),
        a.get(0) * b.get(1) - a.get(1) * b.get(0),
    ]);
}

function expectVectorClose(actual: Vector, expected: readonly number[],
    tolerance = 1e-12): void {
    expect(actual.size).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(actual.get(i)).toBeCloseTo(expected[i], -Math.log10(tolerance));
    }
}

const P = Vector.fromArray([1, 2, 3]);
const U = Vector.fromArray([2, 0, 0]);
const V = Vector.fromArray([0, 3, 0]);

function makePlane(rectangular = true): PlanePatch {
    return new PlanePatch(P, U, V, -1, 4, 0, 5, rectangular);
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

describe('ParametricSurface', () => {
    describe('domain accessors', () => {
        it('returns the umin/umax/vmin/vmax passed to the constructor', () => {
            const surface = makePlane();
            expect(surface.getUMin()).toBe(-1);
            expect(surface.getUMax()).toBe(4);
            expect(surface.getVMin()).toBe(0);
            expect(surface.getVMax()).toBe(5);
        });

        it('reports a rectangular domain', () => {
            expect(makePlane(true).isRectangular()).toBe(true);
        });

        it('reports a triangular domain', () => {
            expect(makePlane(false).isRectangular()).toBe(false);
        });

        it('exposes the runtime dimension that replaces the N template', () => {
            expect(makePlane().getDimension()).toBe(3);
            const planar2D = new PlanePatch(Vector.fromArray([0, 0]),
                Vector.fromArray([1, 0]), Vector.fromArray([0, 1]), 0, 1, 0, 1);
            expect(planar2D.getDimension()).toBe(2);
            expect(planar2D.getPosition(0.25, 0.5).values).toEqual([0.25, 0.5]);
        });
    });

    describe('construction validity (the port of operator bool)', () => {
        it('is true once the derived constructor succeeds', () => {
            expect(makePlane().isConstructed()).toBe(true);
        });

        it('is false when the derived constructor did not complete', () => {
            const broken = new BrokenPatch();
            expect(broken.isConstructed()).toBe(false);
            expect(broken.getPosition(0.5, 0.5).values).toEqual([0, 0, 0]);
        });
    });

    describe('createJet', () => {
        it('has SUP_ORDER equal to 6 (position, 2 first, 3 second order)', () => {
            expect(ParametricSurface.SUP_ORDER).toBe(6);
        });

        it('allocates SUP_ORDER zero vectors of the surface dimension', () => {
            const jet = makePlane().createJet();
            expect(jet.length).toBe(6);
            for (const entry of jet) {
                expect(entry.size).toBe(3);
                expect(entry.values).toEqual([0, 0, 0]);
            }
        });

        it('allocates independent vectors', () => {
            const jet = makePlane().createJet();
            jet[0].set(0, 1);
            expect(jet[1].values).toEqual([0, 0, 0]);
        });
    });

    describe('evaluate jet ordering', () => {
        it('order 0 fills only the position and leaves the rest untouched', () => {
            const surface = makePlane();
            const jet = surface.createJet();
            surface.evaluate(2, 1, 0, jet);
            expectVectorClose(jet[0], [1 + 2 * 2, 2 + 1 * 3, 3]);
            for (let i = 1; i < 6; ++i) {
                expect(jet[i].values).toEqual([0, 0, 0]);
            }
        });

        it('order 1 fills dX/du in slot 1 and dX/dv in slot 2', () => {
            const surface = makePlane();
            const jet = surface.createJet();
            surface.evaluate(2, 1, 1, jet);
            expectVectorClose(jet[0], [5, 5, 3]);
            expectVectorClose(jet[1], [2, 0, 0]);
            expectVectorClose(jet[2], [0, 3, 0]);
            for (let i = 3; i < 6; ++i) {
                expect(jet[i].values).toEqual([0, 0, 0]);
            }
        });

        it('order 2 fills d2X/du2, d2X/dudv, d2X/dv2 in slots 3, 4, 5', () => {
            const surface = new SpherePatch();
            const jet = surface.createJet();
            const u = 0.7, v = 1.1;
            surface.evaluate(u, v, 2, jet);
            const cu = Math.cos(u), su = Math.sin(u);
            const cv = Math.cos(v), sv = Math.sin(v);
            expectVectorClose(jet[0], [cu * sv, su * sv, cv]);
            expectVectorClose(jet[1], [-su * sv, cu * sv, 0]);
            expectVectorClose(jet[2], [cu * cv, su * cv, -sv]);
            expectVectorClose(jet[3], [-cu * sv, -su * sv, 0]);
            expectVectorClose(jet[4], [-su * cv, cu * cv, 0]);
            expectVectorClose(jet[5], [-cu * sv, -su * sv, -cv]);
        });

        it('matches central finite differences of the position function', () => {
            const surface = new SpherePatch();
            const jet = surface.createJet();
            const u = 1.3, v = 0.8, h = 1e-5;
            surface.evaluate(u, v, 2, jet);

            const dU = div(sub(surface.getPosition(u + h, v),
                surface.getPosition(u - h, v)), 2 * h);
            const dV = div(sub(surface.getPosition(u, v + h),
                surface.getPosition(u, v - h)), 2 * h);
            expectVectorClose(dU, jet[1].values, 1e-8);
            expectVectorClose(dV, jet[2].values, 1e-8);

            // d2X/dudv from a central difference of dX/du in v.
            const jetPlus = surface.createJet();
            const jetMinus = surface.createJet();
            surface.evaluate(u, v + h, 1, jetPlus);
            surface.evaluate(u, v - h, 1, jetMinus);
            const dUV = div(sub(jetPlus[1], jetMinus[1]), 2 * h);
            expectVectorClose(dUV, jet[4].values, 1e-8);
        });

        it('produces a zero jet when order is at or beyond SUP_ORDER', () => {
            const surface = new SpherePatch();
            const jet = surface.createJet();
            surface.evaluate(0.5, 0.5, ParametricSurface.SUP_ORDER, jet);
            for (const entry of jet) {
                expect(entry.values).toEqual([0, 0, 0]);
            }
        });
    });

    describe('getPosition', () => {
        it('returns jet[0] for the affine patch', () => {
            const surface = makePlane();
            expectVectorClose(surface.getPosition(0, 0), [1, 2, 3]);
            expectVectorClose(surface.getPosition(1.5, -2), [4, -4, 3]);
        });

        it('returns unit-length points on the sphere patch', () => {
            const surface = new SpherePatch();
            expectVectorClose(surface.getPosition(0, Math.PI / 2), [1, 0, 0]);
            expectVectorClose(surface.getPosition(Math.PI / 2, Math.PI / 2), [0, 1, 0]);
            expectVectorClose(surface.getPosition(0, 0), [0, 0, 1]);
            expect(length(surface.getPosition(2.4, 1.9))).toBeCloseTo(1, 12);
        });

        it('returns a copy the caller may modify freely', () => {
            const surface = makePlane();
            const position = surface.getPosition(1, 1);
            position.set(0, 100);
            expectVectorClose(surface.getPosition(1, 1), [3, 5, 3]);
        });
    });

    describe('getUTangent / getVTangent', () => {
        it('returns the normalized first-order derivatives', () => {
            const surface = makePlane();
            expectVectorClose(surface.getUTangent(2, 3), [1, 0, 0]);
            expectVectorClose(surface.getVTangent(2, 3), [0, 1, 0]);
        });

        it('does not alias or damage vectors stored by the derived class', () => {
            const surface = makePlane();
            surface.getUTangent(0, 0);
            surface.getVTangent(0, 0);
            // U and V would have been shortened to unit length by an
            // aliasing implementation.
            expect(U.values).toEqual([2, 0, 0]);
            expect(V.values).toEqual([0, 3, 0]);
            expectVectorClose(surface.getUTangent(0, 0), [1, 0, 0]);
        });

        it('produces unit-length tangents on the sphere patch', () => {
            const surface = new SpherePatch();
            const u = 2.1, v = 1.4;
            const tu = surface.getUTangent(u, v);
            const tv = surface.getVTangent(u, v);
            expect(length(tu)).toBeCloseTo(1, 12);
            expect(length(tv)).toBeCloseTo(1, 12);
            // The sphere's parameterization is orthogonal, and both tangents
            // are perpendicular to the (radial) position.
            expect(dot(tu, tv)).toBeCloseTo(0, 12);
            expect(dot(tu, surface.getPosition(u, v))).toBeCloseTo(0, 12);
            expect(dot(tv, surface.getPosition(u, v))).toBeCloseTo(0, 12);
        });
    });

    describe('normal orientation', () => {
        it('uses cross(uTangent, vTangent) consistent with the jet ordering', () => {
            const surface = makePlane();
            const normal = cross3(surface.getUTangent(1, 1), surface.getVTangent(1, 1));
            expectVectorClose(normal, [0, 0, 1]);
        });

        it('flips when the two parameter directions are swapped', () => {
            const swapped = new PlanePatch(P, V, U, -1, 4, 0, 5);
            const normal = cross3(swapped.getUTangent(1, 1), swapped.getVTangent(1, 1));
            expectVectorClose(normal, [0, 0, -1]);
        });

        it('points inward for the (longitude, colatitude) sphere patch', () => {
            const surface = new SpherePatch();
            const u = 0.9, v = 1.2;
            const normal = cross3(surface.getUTangent(u, v), surface.getVTangent(u, v));
            const position = surface.getPosition(u, v);
            // cross(Xu, Xv) = -sin(v) * X, so the unit normal is -X.
            expectVectorClose(normal, position.values.map((x) => -x));
            expect(length(normal)).toBeCloseTo(1, 12);
        });
    });
});
