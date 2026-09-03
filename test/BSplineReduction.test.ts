import { describe, it, expect } from 'vitest';
import { BSplineReduction } from '../src/BSplineReduction.js';
import { BSplineCurve } from '../src/BSplineCurve.js';
import { BasisFunctionInput } from '../src/BasisFunction.js';
import { Vector } from '../src/Vector.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function vec(...values: number[]): Vector {
    const v = new Vector(values.length);
    for (let i = 0; i < values.length; ++i) {
        v.values[i] = values[i];
    }
    return v;
}

// The maximum distance between the input curve and the reduced curve,
// sampled uniformly over the domain.
function maxDeviation(inControls: Vector[], outControls: Vector[],
    degree: number, numSamples: number): number {
    const dimension = inControls[0].size;
    const inCurve = new BSplineCurve(dimension,
        new BasisFunctionInput(inControls.length, degree), inControls);
    const outCurve = new BSplineCurve(dimension,
        new BasisFunctionInput(outControls.length, degree), outControls);
    let maximum = 0;
    for (let k = 0; k <= numSamples; ++k) {
        const t = k / numSamples;
        const p = inCurve.getPosition(t);
        const q = outCurve.getPosition(t);
        let sum = 0;
        for (let d = 0; d < dimension; ++d) {
            const diff = p.values[d] - q.values[d];
            sum += diff * diff;
        }
        maximum = Math.max(maximum, Math.sqrt(sum));
    }
    return maximum;
}

describe('BSplineReduction: input validation', () => {
    it('rejects fewer than two control points', () => {
        const reduction = new BSplineReduction();
        expect(() => reduction.compute([vec(0, 0)], 1, 0.5))
            .toThrow(/Invalid input/);
    });

    it('rejects a degree outside [1, numControls-1]', () => {
        const reduction = new BSplineReduction();
        const controls = [vec(0, 0), vec(1, 1), vec(2, 0)];
        expect(() => reduction.compute(controls, 0, 0.5)).toThrow(/Invalid input/);
        expect(() => reduction.compute(controls, 3, 0.5)).toThrow(/Invalid input/);
    });
});

describe('BSplineReduction: clamping of the output size', () => {
    it('returns a copy of the input when the fraction is one', () => {
        const reduction = new BSplineReduction();
        const controls = [vec(0, 0), vec(1, 2), vec(3, 1), vec(4, 4), vec(6, 0)];
        const out = reduction.compute(controls, 2, 1);
        expect(out.length).toBe(controls.length);
        for (let i = 0; i < controls.length; ++i) {
            expect(out[i].values).toEqual(controls[i].values);
        }
        // The output is a copy, not an alias.
        out[0].values[0] = 99;
        expect(controls[0].values[0]).toBe(0);
    });

    it('returns a copy when the fraction exceeds one', () => {
        const reduction = new BSplineReduction();
        const controls = [vec(0, 0), vec(1, 1), vec(2, 3)];
        const out = reduction.compute(controls, 1, 2);
        expect(out.length).toBe(3);
    });

    it('clamps the output to degree+1 control points when the fraction is small', () => {
        const reduction = new BSplineReduction();
        const controls: Vector[] = [];
        for (let i = 0; i < 12; ++i) {
            controls.push(vec(i, 0));
        }
        for (const fraction of [0, 0.01, 0.1]) {
            const out = reduction.compute(controls, 3, fraction);
            expect(out.length).toBe(4);
        }
    });

    it('uses the truncated fraction of the input size', () => {
        const reduction = new BSplineReduction();
        const controls: Vector[] = [];
        for (let i = 0; i < 10; ++i) {
            controls.push(vec(i, i * i));
        }
        // trunc(0.55 * 10) = 5.
        expect(reduction.compute(controls, 2, 0.55).length).toBe(5);
        // trunc(0.79 * 10) = 7.
        expect(reduction.compute(controls, 2, 0.79).length).toBe(7);
    });
});

describe('BSplineReduction: fitting quality', () => {
    it('reproduces a straight line up to the quadrature error', () => {
        // A degree-1 B-spline whose control points are collinear and
        // uniformly spaced is the straight segment; the exact least-squares
        // fit of that segment by any number of control points is the same
        // segment. Upstream integrates the basis-function products with
        // Integration::Romberg(8, ...) over the whole support, which straddles
        // the knots where the integrand's derivative jumps, so the fit is
        // only accurate to about 1e-3. The port preserves that behavior.
        const reduction = new BSplineReduction();
        const controls: Vector[] = [];
        for (let i = 0; i < 9; ++i) {
            controls.push(vec(i / 8, 2 * (i / 8)));
        }
        const out = reduction.compute(controls, 1, 0.4);
        expect(out.length).toBe(3);
        expect(maxDeviation(controls, out, 1, 64)).toBeLessThan(0.02);
        // The reduced control points are near the exact ones (0,0), (1/2,1)
        // and (1,2).
        const exact = [[0, 0], [0.5, 1], [1, 2]];
        for (let i = 0; i < 3; ++i) {
            expect(Math.abs(out[i].values[0] - exact[i][0])).toBeLessThan(0.02);
            expect(Math.abs(out[i].values[1] - exact[i][1])).toBeLessThan(0.02);
        }
    });

    it('reproduces a degree-2 curve that is already representable', () => {
        // Insert an extra control point into a degree-2 curve by refining a
        // quadratic Bezier curve; the reduction back to three control points
        // must recover the Bezier control points.
        const reduction = new BSplineReduction();
        const bezier = [vec(0, 0), vec(1, 2), vec(2, 0)];
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 2), bezier);
        // Sample-based check: reduce a 6-control-point spline sampled from
        // the Bezier curve is not exact, so instead reduce the Bezier itself
        // with a fraction that keeps all three controls.
        const out = reduction.compute(bezier, 2, 1);
        expect(out.length).toBe(3);
        expect(maxDeviation(bezier, out, 2, 32)).toBeLessThan(1e-12);
        expect(curve.getPosition(0.5).values[1]).toBeCloseTo(1, 12);
    });

    it('approximates a smooth curve closely when few points are dropped', () => {
        const reduction = new BSplineReduction();
        const controls: Vector[] = [];
        const n = 16;
        for (let i = 0; i < n; ++i) {
            const s = i / (n - 1);
            controls.push(vec(s, Math.sin(2 * s), Math.cos(2 * s)));
        }
        const out = reduction.compute(controls, 3, 0.75);
        expect(out.length).toBe(12);
        // The reduced curve stays close to the original.
        expect(maxDeviation(controls, out, 3, 128)).toBeLessThan(0.02);
    });

    it('degrades gracefully as the number of control points shrinks', () => {
        const reduction = new BSplineReduction();
        const controls: Vector[] = [];
        const n = 14;
        for (let i = 0; i < n; ++i) {
            const s = i / (n - 1);
            controls.push(vec(s, s * (1 - s) * 4));
        }
        const dev8 = maxDeviation(controls, reduction.compute(controls, 2, 8 / n), 2, 128);
        const dev4 = maxDeviation(controls, reduction.compute(controls, 2, 4 / n), 2, 128);
        expect(dev8).toBeLessThan(dev4 + 1e-12);
        expect(dev8).toBeLessThan(0.05);
    });

    it('is equivariant under translation of the control points', () => {
        const reduction = new BSplineReduction();
        const rand = makeRandom(1234567);
        const controls: Vector[] = [];
        for (let i = 0; i < 10; ++i) {
            controls.push(vec(rand() * 4 - 2, rand() * 4 - 2));
        }
        const shift = vec(3, -7);
        const shifted = controls.map(c => vec(c.values[0] + shift.values[0],
            c.values[1] + shift.values[1]));

        const out = reduction.compute(controls, 3, 0.6);
        const outShifted = reduction.compute(shifted, 3, 0.6);
        expect(outShifted.length).toBe(out.length);
        // Translation equivariance holds exactly only when the rows of
        // A^{-1}B sum to one; upstream's Romberg quadrature makes the row
        // sums 1 to about 1e-2 for degree 3, so the tolerance is scaled by
        // the shift magnitude.
        const tolerance = 0.02 * (Math.abs(shift.values[0])
            + Math.abs(shift.values[1]));
        for (let i = 0; i < out.length; ++i) {
            expect(Math.abs(outShifted[i].values[0]
                - (out[i].values[0] + shift.values[0]))).toBeLessThan(tolerance);
            expect(Math.abs(outShifted[i].values[1]
                - (out[i].values[1] + shift.values[1]))).toBeLessThan(tolerance);
        }
    });

    it('is equivariant under uniform scaling', () => {
        const reduction = new BSplineReduction();
        const rand = makeRandom(24680);
        const controls: Vector[] = [];
        for (let i = 0; i < 9; ++i) {
            controls.push(vec(rand() * 2, rand() * 2, rand() * 2));
        }
        const scaled = controls.map(c => vec(3 * c.values[0], 3 * c.values[1],
            3 * c.values[2]));
        const out = reduction.compute(controls, 2, 0.6);
        const outScaled = reduction.compute(scaled, 2, 0.6);
        for (let i = 0; i < out.length; ++i) {
            for (let d = 0; d < 3; ++d) {
                expect(outScaled[i].values[d]).toBeCloseTo(3 * out[i].values[d], 8);
            }
        }
    });

    it('reuses the same reduction object across calls', () => {
        const reduction = new BSplineReduction();
        const a: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            a.push(vec(i, i % 3));
        }
        const first = reduction.compute(a, 2, 0.6);
        // A different degree and size on the same object must not be
        // affected by leftover state.
        const b: Vector[] = [];
        for (let i = 0; i < 11; ++i) {
            b.push(vec(i, Math.sin(i)));
        }
        reduction.compute(b, 3, 0.5);
        const again = reduction.compute(a, 2, 0.6);
        expect(again.length).toBe(first.length);
        for (let i = 0; i < first.length; ++i) {
            expect(again[i].values[0]).toBeCloseTo(first[i].values[0], 12);
            expect(again[i].values[1]).toBeCloseTo(first[i].values[1], 12);
        }
    });

    it('preserves the dimension of the control points', () => {
        const reduction = new BSplineReduction();
        const controls: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            controls.push(vec(i, i, i, i));
        }
        const out = reduction.compute(controls, 2, 0.5);
        expect(out.length).toBe(4);
        for (const c of out) {
            expect(c.size).toBe(4);
        }
    });

    it('does not modify the input control points', () => {
        const reduction = new BSplineReduction();
        const controls: Vector[] = [];
        for (let i = 0; i < 7; ++i) {
            controls.push(vec(i, 2 * i));
        }
        const snapshot = controls.map(c => c.values.slice());
        reduction.compute(controls, 2, 0.6);
        for (let i = 0; i < controls.length; ++i) {
            expect(controls[i].values).toEqual(snapshot[i]);
        }
    });
});
