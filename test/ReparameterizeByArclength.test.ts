import { describe, it, expect } from 'vitest';
import { ReparameterizeByArclength } from '../src/ReparameterizeByArclength.js';
import { BezierCurve } from '../src/BezierCurve.js';
import { ParametricCurve } from '../src/ParametricCurve.js';
import { Vector } from '../src/Vector.js';

function vec(...values: number[]): Vector {
    const v = new Vector(values.length);
    for (let i = 0; i < values.length; ++i) {
        v.values[i] = values[i];
    }
    return v;
}

// A circle of radius r traversed on [0, 2*pi]; its speed is the constant r,
// so arclength s corresponds exactly to t = s / r.
class Circle2 extends ParametricCurve {
    constructor(private r: number) {
        super(2, 0, 2 * Math.PI);
        this.mConstructed = true;
    }

    override evaluate(t: number, order: number, jet: Vector[]): void {
        const c = Math.cos(t), s = Math.sin(t), r = this.r;
        jet[0] = vec(r * c, r * s);
        if (order >= 1) { jet[1] = vec(-r * s, r * c); }
        if (order >= 2) { jet[2] = vec(-r * c, -r * s); }
        if (order >= 3) { jet[3] = vec(r * s, -r * c); }
    }
}

describe('ReparameterizeByArclength', () => {
    it('exposes the curve, its domain and its total arclength', () => {
        const curve = new Circle2(3);
        const repar = new ReparameterizeByArclength(curve);
        expect(repar.getCurve()).toBe(curve);
        expect(repar.getTMin()).toBe(0);
        expect(repar.getTMax()).toBeCloseTo(2 * Math.PI, 14);
        expect(repar.getTotalArclength()).toBeCloseTo(6 * Math.PI, 8);
    });

    it('clamps s at or below zero to tMin', () => {
        const repar = new ReparameterizeByArclength(new Circle2(2));
        for (const s of [0, -1, -1e6]) {
            const out = repar.getT(s, true);
            expect(out.t).toBe(0);
            expect(out.f).toBe(0);
            expect(out.numIterations).toBe(0);
        }
    });

    it('clamps s at or above the total arclength to tMax', () => {
        const repar = new ReparameterizeByArclength(new Circle2(2));
        const L = repar.getTotalArclength();
        for (const s of [L, L + 1, 1e6]) {
            const out = repar.getT(s, true);
            expect(out.t).toBeCloseTo(2 * Math.PI, 14);
            expect(out.f).toBe(0);
            expect(out.numIterations).toBe(0);
        }
    });

    it('gives uniform angles for uniform arclength samples on a circle', () => {
        for (const r of [1, 2.5, 7]) {
            const repar = new ReparameterizeByArclength(new Circle2(r));
            const L = repar.getTotalArclength();
            const n = 16;
            for (const useBisection of [true, false]) {
                let previous = -1;
                for (let i = 1; i < n; ++i) {
                    const s = L * i / n;
                    const out = repar.getT(s, useBisection);
                    // Constant speed r means t = s / r.
                    expect(out.t).toBeCloseTo(s / r, 8);
                    expect(out.t).toBeGreaterThan(previous);
                    previous = out.t;
                    expect(Math.abs(out.f)).toBeLessThan(1e-8);
                }
            }
        }
    });

    it('gives equal angular increments between consecutive samples', () => {
        const repar = new ReparameterizeByArclength(new Circle2(4));
        const L = repar.getTotalArclength();
        const n = 12;
        const ts: number[] = [];
        for (let i = 0; i <= n; ++i) {
            ts.push(repar.getT(L * i / n, false).t);
        }
        const expectedStep = 2 * Math.PI / n;
        for (let i = 1; i <= n; ++i) {
            expect(ts[i] - ts[i - 1]).toBeCloseTo(expectedStep, 7);
        }
    });

    it('inverts the arclength of a non-uniform-speed Bezier curve', () => {
        const curve = new BezierCurve(2, 3,
            [vec(0, 0), vec(0.2, 3), vec(2.8, -2), vec(4, 1)]);
        const repar = new ReparameterizeByArclength(curve);
        const L = repar.getTotalArclength();
        expect(L).toBeGreaterThan(0);
        for (const useBisection of [true, false]) {
            for (let i = 1; i < 10; ++i) {
                const s = L * i / 10;
                const out = repar.getT(s, useBisection);
                expect(out.t).toBeGreaterThan(0);
                expect(out.t).toBeLessThan(1);
                // The residual F(t,s) = Arclength(tMin,t) - s is tiny.
                expect(Math.abs(curve.getLength(0, out.t) - s))
                    .toBeLessThan(1e-9);
            }
        }
    });

    it('agrees between bisection and the Newton hybrid', () => {
        const curve = new BezierCurve(3, 4,
            [vec(0, 0, 0), vec(1, 2, -1), vec(3, -1, 2), vec(4, 2, 1),
                vec(6, 0, 3)]);
        const repar = new ReparameterizeByArclength(curve);
        const L = repar.getTotalArclength();
        for (let i = 1; i < 8; ++i) {
            const s = L * i / 8;
            const a = repar.getT(s, true);
            const b = repar.getT(s, false);
            expect(a.t).toBeCloseTo(b.t, 9);
        }
    });

    it('uses far fewer iterations for the Newton hybrid than for bisection',
        () => {
            const curve = new BezierCurve(2, 3,
                [vec(0, 0), vec(1, 3), vec(3, -1), vec(4, 2)]);
            const repar = new ReparameterizeByArclength(curve);
            const s = 0.5 * repar.getTotalArclength();
            const bisect = repar.getT(s, true);
            const newton = repar.getT(s, false);
            expect(bisect.numIterations).toBeGreaterThan(10);
            expect(newton.numIterations).toBeLessThan(bisect.numIterations);
            expect(newton.numIterations).toBeGreaterThan(0);
        });

    it('is monotone in s', () => {
        const curve = new BezierCurve(2, 3,
            [vec(0, 0), vec(1, 3), vec(3, -1), vec(4, 2)]);
        const repar = new ReparameterizeByArclength(curve);
        const L = repar.getTotalArclength();
        let previous = -1;
        for (let i = 0; i <= 20; ++i) {
            const out = repar.getT(L * i / 20, false);
            expect(out.t).toBeGreaterThanOrEqual(previous);
            previous = out.t;
        }
        expect(previous).toBeCloseTo(1, 12);
    });

    it('is consistent with ParametricCurve.getTime', () => {
        const curve = new BezierCurve(2, 3,
            [vec(0, 0), vec(1, 3), vec(3, -1), vec(4, 2)]);
        const repar = new ReparameterizeByArclength(curve);
        const L = repar.getTotalArclength();
        for (let i = 1; i < 6; ++i) {
            const s = L * i / 6;
            expect(repar.getT(s, false).t).toBeCloseTo(curve.getTime(s), 6);
        }
    });
});


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). ReparameterizeByArclength.h
// was read line by line against src/ReparameterizeByArclength.ts. The class is
// a root finder for F(t,s) = Arclength(tMin,t) - s, so the properties below
// check the returned root against an independently known arclength function,
// check that the two solver modes agree, and check the unit-speed property
// that motivates the class.
import {
    check, fc, expectClose, wellScaledVector
} from './helpers/arbitraries.js';

// X(t) = (t^3, 0) on [-1,1]: the speed 3t^2 vanishes at t = 0, which is the
// F'(t) == 0 case the hybrid solver has to handle, and the arclength from -1
// to t is exactly t^3 + 1 (a quadratic speed that Romberg integrates
// exactly).
class CuspCurve extends ParametricCurve {
    constructor() {
        super(2, -1, 1);
        this.mConstructed = true;
    }

    override evaluate(t: number, order: number, jet: Vector[]): void {
        jet[0] = vec(t * t * t, 0);
        if (order >= 1) { jet[1] = vec(3 * t * t, 0); }
        if (order >= 2) { jet[2] = vec(6 * t, 0); }
        if (order >= 3) { jet[3] = vec(6, 0); }
    }
}

const bezierCurve = fc.integer({ min: 2, max: 4 }).chain(degree =>
    fc.array(wellScaledVector(2, -4, 4),
        { minLength: degree + 1, maxLength: degree + 1 }))
    .map(controls => new BezierCurve(2, controls.length - 1, controls))
    .filter(curve => {
        // The speed |X'(t)| must stay bounded away from zero: near a
        // stationary point the speed has a corner, Romberg integration of it
        // converges slowly, and the arclength values the properties below
        // compare are then dominated by quadrature error rather than by the
        // root finder.
        const L = curve.getTotalLength();
        if (!(L > 0.5)) { return false; }
        for (let i = 0; i <= 20; ++i) {
            if (curve.getSpeed(i / 20) < 0.25 * L) { return false; }
        }
        return true;
    });

const fraction = fc.double({ min: 0.02, max: 0.98, noNaN: true,
    noDefaultInfinity: true });

describe('ReparameterizeByArclength verification', () => {
    it('returns a root of F(t,s) for both solver modes', () => {
        check(fc.tuple(bezierCurve, fraction), ([curve, r]) => {
            const repar = new ReparameterizeByArclength(curve);
            const L = repar.getTotalArclength();
            const s = r * L;
            const bisect = repar.getT(s, true);
            const newton = repar.getT(s, false);
            for (const out of [bisect, newton]) {
                expect(out.t).toBeGreaterThanOrEqual(repar.getTMin());
                expect(out.t).toBeLessThanOrEqual(repar.getTMax());
                // The residual is limited by the Romberg quadrature inside
                // GetLength, not by the root finder.
                expect(Math.abs(curve.getLength(repar.getTMin(), out.t) - s))
                    .toBeLessThan(1e-9 * Math.max(1, L));
                expect(Math.abs(out.f)).toBeLessThan(1e-9 * Math.max(1, L));
                expect(out.numIterations).toBeGreaterThan(0);
            }
            expectClose(bisect.t, newton.t, 1e-8, 1e-8);
        }, 25);
    });

    it('inverts an arclength known in closed form, stationary point included',
        () => {
            // Arclength(-1, t) = t^3 + 1, so the root of F(t,s) is
            // t = (s-1)^(1/3). The solver has to cross t = 0, where the speed
            // (and hence F') is exactly zero.
            const curve = new CuspCurve();
            const repar = new ReparameterizeByArclength(curve);
            expectClose(repar.getTotalArclength(), 2, 1e-12, 1e-12);
            check(fc.tuple(fc.double({ min: 0.01, max: 1.99, noNaN: true,
                noDefaultInfinity: true }), fc.boolean()),
            ([s, useBisection]) => {
                const out = repar.getT(s, useBisection);
                const expected = Math.cbrt(s - 1);
                // Near t = 0 the map t -> t^3 flattens, so a small residual in
                // s corresponds to a much larger error in t; compare the
                // arclength instead, which is what the solver actually
                // controls.
                expectClose(out.t ** 3 + 1, s, 1e-9, 1e-9);
                expect(Number.isFinite(out.t)).toBe(true);
                if (Math.abs(expected) > 1e-3) {
                    expect(Math.sign(out.t)).toBe(Math.sign(expected));
                    expectClose(out.t, expected, 1e-6, 1e-6);
                }
            }, 40);
        });

    it('is monotone in the arclength', () => {
        check(fc.tuple(bezierCurve, fraction, fraction, fc.boolean()),
            ([curve, r0, r1, useBisection]) => {
                const repar = new ReparameterizeByArclength(curve);
                const L = repar.getTotalArclength();
                const lo = Math.min(r0, r1) * L;
                const hi = Math.max(r0, r1) * L;
                const a = repar.getT(lo, useBisection);
                const b = repar.getT(hi, useBisection);
                expect(a.t).toBeLessThanOrEqual(b.t + 1e-12);
            }, 25);
    });

    it('produces a unit-speed reparameterization', () => {
        // Sampling s uniformly and mapping through getT gives points that are
        // equally spaced in arclength along the curve; that is the whole
        // point of the class.
        check(fc.tuple(bezierCurve, fc.integer({ min: 3, max: 6 })),
            ([curve, n]) => {
                const repar = new ReparameterizeByArclength(curve);
                const L = repar.getTotalArclength();
                const ts: number[] = [];
                for (let i = 0; i <= n; ++i) {
                    ts.push(repar.getT((L * i) / n, false).t);
                }
                for (let i = 0; i <= n; ++i) {
                    // The cumulative arclength is exactly what the solver
                    // controls, so this identity holds to the residual of the
                    // root finder.
                    expectClose(curve.getLength(repar.getTMin(), ts[i]),
                        (L * i) / n, 1e-9 * Math.max(1, L), 1e-9);
                }
                for (let i = 1; i <= n; ++i) {
                    // The length of a single subinterval is a different
                    // Romberg integration from the difference of the two
                    // cumulative ones, so the agreement here is limited by the
                    // quadrature rather than by the solver.
                    expectClose(curve.getLength(ts[i - 1], ts[i]), L / n,
                        1e-4 * Math.max(1, L), 1e-4);
                }
            }, 20);
    });

    it('clamps the arclength to the ends of the domain', () => {
        check(fc.tuple(bezierCurve, fc.double({ min: 0, max: 1e6,
            noNaN: true, noDefaultInfinity: true }), fc.boolean()),
        ([curve, d, useBisection]) => {
            const repar = new ReparameterizeByArclength(curve);
            const L = repar.getTotalArclength();
            const lo = repar.getT(-d, useBisection);
            expect(lo.t).toBe(repar.getTMin());
            expect(lo.f).toBe(0);
            expect(lo.numIterations).toBe(0);
            const hi = repar.getT(L + d, useBisection);
            expect(hi.t).toBe(repar.getTMax());
            expect(hi.f).toBe(0);
            expect(hi.numIterations).toBe(0);
        }, 25);
    });
});
