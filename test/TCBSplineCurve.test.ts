import { describe, it, expect } from 'vitest';
import { TCBSplineCurve } from '../src/TCBSplineCurve.js';
import { Vector, length as vectorLength } from '../src/Vector.js';

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

function zeros(n: number): number[] {
    return new Array<number>(n).fill(0);
}

function positionOf(curve: TCBSplineCurve, t: number): Vector {
    const jet = curve.createJet();
    curve.evaluate(t, 0, jet);
    return jet[0];
}

// The classical uniform Catmull-Rom spline: on segment k the tangents are
// m_k = (P_{k+1} - P_{k-1})/2 with the boundaries using duplicated points.
function catmullRom(points: readonly Vector[], u: number, key: number): Vector {
    const n = points.length;
    const tangent = (k: number): Vector => {
        const prev = points[Math.max(k - 1, 0)];
        const next = points[Math.min(k + 1, n - 1)];
        const m = new Vector(prev.size);
        for (let i = 0; i < m.size; ++i) {
            m.values[i] = 0.5 * (next.values[i] - prev.values[i]);
        }
        return m;
    };
    const P0 = points[key], P1 = points[key + 1];
    const m0 = tangent(key), m1 = tangent(key + 1);
    const h00 = 2 * u * u * u - 3 * u * u + 1;
    const h10 = u * u * u - 2 * u * u + u;
    const h01 = -2 * u * u * u + 3 * u * u;
    const h11 = u * u * u - u * u;
    const result = new Vector(P0.size);
    for (let i = 0; i < result.size; ++i) {
        result.values[i] = h00 * P0.values[i] + h10 * m0.values[i]
            + h01 * P1.values[i] + h11 * m1.values[i];
    }
    return result;
}

describe('TCBSplineCurve', () => {
    const points = [vec(0, 0), vec(1, 2), vec(3, 1), vec(4, -1), vec(6, 2)];
    const times = [0, 1, 2, 3, 4];

    function zeroTCB(): TCBSplineCurve {
        const n = points.length;
        return new TCBSplineCurve(2, points, times, zeros(n), zeros(n),
            zeros(n), []);
    }

    it('reports its key frames and domain', () => {
        const curve = zeroTCB();
        expect(curve.isConstructed()).toBe(true);
        expect(curve.getNumKeyFrames()).toBe(5);
        expect(curve.getNumSegments()).toBe(4);
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBe(4);
        expect(curve.getPoints().length).toBe(5);
        expect(curve.getTensions()).toEqual(zeros(5));
        expect(curve.getContinuities()).toEqual(zeros(5));
        expect(curve.getBiases()).toEqual(zeros(5));
        expect(curve.getLambdas()).toEqual([]);
        expect(curve.getInTangents().length).toBe(5);
        expect(curve.getOutTangents().length).toBe(5);
    });

    it('rejects mismatched input sizes', () => {
        expect(() => new TCBSplineCurve(2, [vec(0, 0)], [0], [0], [0], [0], []))
            .toThrow(/Invalid size/);
        expect(() => new TCBSplineCurve(2, points, [0, 1, 2, 3], zeros(5),
            zeros(5), zeros(5), [])).toThrow(/Invalid size/);
        expect(() => new TCBSplineCurve(2, points, times, zeros(4), zeros(5),
            zeros(5), [])).toThrow(/Invalid size/);
        expect(() => new TCBSplineCurve(2, points, times, zeros(5), zeros(5),
            zeros(5), [1, 2])).toThrow(/Invalid size/);
    });

    it('copies its inputs (C++ value semantics)', () => {
        const p = points.map(q => q.clone());
        const curve = new TCBSplineCurve(2, p, times, zeros(5), zeros(5),
            zeros(5), []);
        p[1].values[0] = 100;
        expect(curve.getPoints()[1].values[0]).toBe(1);
    });

    it('passes through the key frame points at the key times', () => {
        const curve = zeroTCB();
        for (let i = 0; i < points.length; ++i) {
            const p = positionOf(curve, times[i]);
            expect(p.values[0]).toBeCloseTo(points[i].values[0], 12);
            expect(p.values[1]).toBeCloseTo(points[i].values[1], 12);
        }
    });

    it('reproduces Catmull-Rom for zero tension, continuity and bias', () => {
        const curve = zeroTCB();
        for (let key = 0; key < 4; ++key) {
            for (let s = 0; s <= 8; ++s) {
                const u = s / 8;
                const t = times[key] + u;
                const actual = positionOf(curve, t);
                const expected = catmullRom(points, u, key);
                expect(actual.values[0]).toBeCloseTo(expected.values[0], 11);
                expect(actual.values[1]).toBeCloseTo(expected.values[1], 11);
            }
        }
    });

    it('has the Catmull-Rom tangents at the key frames', () => {
        const curve = zeroTCB();
        const inTan = curve.getInTangents();
        const outTan = curve.getOutTangents();
        // Interior: (P_{k+1} - P_{k-1}) / 2, with unit time spacing.
        for (let k = 1; k <= 3; ++k) {
            for (let i = 0; i < 2; ++i) {
                const expected = 0.5 * (points[k + 1].values[i]
                    - points[k - 1].values[i]);
                expect(inTan[k].values[i]).toBeCloseTo(expected, 12);
                expect(outTan[k].values[i]).toBeCloseTo(expected, 12);
            }
        }
        // Boundaries: half the one-sided difference.
        for (let i = 0; i < 2; ++i) {
            expect(outTan[0].values[i]).toBeCloseTo(
                0.5 * (points[1].values[i] - points[0].values[i]), 12);
            expect(inTan[4].values[i]).toBeCloseTo(
                0.5 * (points[4].values[i] - points[3].values[i]), 12);
        }
    });

    it('clamps evaluation outside the time domain', () => {
        const curve = zeroTCB();
        const before = positionOf(curve, -5);
        const after = positionOf(curve, 100);
        expect(before.values).toEqual(points[0].values);
        expect(after.values[0]).toBeCloseTo(points[4].values[0], 11);
        expect(after.values[1]).toBeCloseTo(points[4].values[1], 11);
    });

    it('honors user-specified boundary tangents', () => {
        const first = vec(5, -3);
        const last = vec(-2, 7);
        const curve = new TCBSplineCurve(2, points, times, zeros(5), zeros(5),
            zeros(5), [], first, last);
        expect(curve.getOutTangents()[0].values).toEqual([5, -3]);
        expect(curve.getInTangents()[0].values).toEqual([5, -3]);
        expect(curve.getInTangents()[4].values).toEqual([-2, 7]);
        expect(curve.getOutTangents()[4].values).toEqual([-2, 7]);

        // The first derivative at t=0 is the specified outgoing tangent.
        const jet = curve.createJet();
        curve.evaluate(0, 1, jet);
        expect(jet[1].values[0]).toBeCloseTo(5, 11);
        expect(jet[1].values[1]).toBeCloseTo(-3, 11);
    });

    it('has derivatives that match finite differences', () => {
        const curve = new TCBSplineCurve(2, points, times, [0.2, 0.2, 0.2, 0.2, 0.2],
            [-0.3, -0.3, -0.3, -0.3, -0.3], [0.1, 0.1, 0.1, 0.1, 0.1], []);
        const jet = curve.createJet();
        const h = 1e-5;
        for (const t of [0.3, 1.4, 2.5, 3.6]) {
            curve.evaluate(t, 3, jet);
            const p = positionOf(curve, t);
            const pp = positionOf(curve, t + h);
            const pm = positionOf(curve, t - h);
            for (let k = 0; k < 2; ++k) {
                const d1 = (pp.values[k] - pm.values[k]) / (2 * h);
                const d2 = (pp.values[k] - 2 * p.values[k] + pm.values[k])
                    / (h * h);
                expect(jet[1].values[k]).toBeCloseTo(d1, 6);
                expect(jet[2].values[k]).toBeCloseTo(d2, 3);
            }
            // The third derivative is constant on a segment: 6*D/delta^3.
            expect(Number.isFinite(jet[3].values[0])).toBe(true);
        }
    });

    it('makes the curve C1 across interior key frames when lambda is unset',
        () => {
            const curve = zeroTCB();
            const jet = curve.createJet();
            const h = 1e-7;
            for (const t of [1, 2, 3]) {
                curve.evaluate(t - h, 1, jet);
                const left = [jet[1].values[0], jet[1].values[1]];
                curve.evaluate(t + h, 1, jet);
                expect(jet[1].values[0]).toBeCloseTo(left[0], 5);
                expect(jet[1].values[1]).toBeCloseTo(left[1], 5);
            }
        });

    it('scales the tangents by lambda for speed continuity', () => {
        const lambda = [1, 2, 3, 4, 5];
        const curve = new TCBSplineCurve(2, points, times, zeros(5), zeros(5),
            zeros(5), lambda);
        // At an interior key frame the in- and out-tangent lengths become
        // equal (both are 2*lambda*inLen*outLen/(inLen+outLen)).
        for (let k = 1; k <= 3; ++k) {
            const inLen = vectorLength(curve.getInTangents()[k]);
            const outLen = vectorLength(curve.getOutTangents()[k]);
            expect(inLen).toBeCloseTo(outLen, 11);
        }
        // Boundary tangents are simply scaled by lambda.
        expect(vectorLength(curve.getOutTangents()[0])).toBeCloseTo(
            lambda[0] * 0.5 * Math.hypot(1, 2), 11);
        expect(vectorLength(curve.getInTangents()[4])).toBeCloseTo(
            lambda[4] * 0.5 * Math.hypot(2, 3), 11);
    });

    it('handles the two-point case (a cubic Hermite segment)', () => {
        const curve = new TCBSplineCurve(3, [vec(0, 0, 0), vec(1, 1, 1)],
            [0, 2], zeros(2), zeros(2), zeros(2), []);
        expect(curve.getNumSegments()).toBe(1);
        const p0 = positionOf(curve, 0);
        const p1 = positionOf(curve, 2);
        expect(p0.values).toEqual([0, 0, 0]);
        for (let k = 0; k < 3; ++k) {
            expect(p1.values[k]).toBeCloseTo(1, 11);
        }
    });

    it('handles nonuniform time spacing', () => {
        const nonuniform = [0, 0.5, 2.5, 3, 7];
        const curve = new TCBSplineCurve(2, points, nonuniform, zeros(5),
            zeros(5), zeros(5), []);
        for (let i = 0; i < 5; ++i) {
            const p = positionOf(curve, nonuniform[i]);
            expect(p.values[0]).toBeCloseTo(points[i].values[0], 11);
            expect(p.values[1]).toBeCloseTo(points[i].values[1], 11);
        }
    });

    it('reduces to a line for collinear points with zero TCB', () => {
        const rand = makeRandom(0x7ea51);
        const line: Vector[] = [];
        for (let i = 0; i < 6; ++i) {
            line.push(vec(i, 2 * i, -i));
        }
        const curve = new TCBSplineCurve(3, line, [0, 1, 2, 3, 4, 5],
            zeros(6), zeros(6), zeros(6), []);
        for (let s = 0; s < 20; ++s) {
            const t = 5 * rand();
            const p = positionOf(curve, t);
            expect(p.values[1]).toBeCloseTo(2 * p.values[0], 10);
            expect(p.values[2]).toBeCloseTo(-p.values[0], 10);
        }
    });
});
