import { describe, it, expect } from 'vitest';
import {
    perp, unitPerp, dotPerp, computeOrthogonalComplement2,
    computeBarycentrics2, IntrinsicsVector2
} from '../src/Vector2';
import { Vector, add, dot, length, mul } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

describe('perp, unitPerp, dotPerp', () => {
    it('perp((x0,x1)) = (x1,-x0)', () => {
        expect(perp(v2(3, 7)).values).toEqual([7, -3]);
        expect(perp(v2(1, 0)).values).toEqual([0, -1]);
    });

    it('perp is orthogonal to its input and perp(perp(v)) = -v', () => {
        const rng = makeRng(1234);
        for (let trial = 0; trial < 10; ++trial) {
            const v = v2((rng() - 0.5) * 10, (rng() - 0.5) * 10);
            expect(dot(v, perp(v))).toBe(0);
            const pp = perp(perp(v));
            expect(pp.values[0]).toBe(-v.values[0]);
            expect(pp.values[1]).toBe(-v.values[1]);
            // perp is a clockwise rotation, so it preserves length.
            expect(length(perp(v))).toBeCloseTo(length(v), 14);
        }
    });

    it('perp requires size 2', () => {
        expect(() => perp(Vector.fromArray([1, 2, 3]))).toThrow('size 2');
    });

    it('unitPerp is the normalized perpendicular', () => {
        const u = unitPerp(v2(3, 4));
        expect(u.values[0]).toBeCloseTo(0.8, 15);
        expect(u.values[1]).toBeCloseTo(-0.6, 15);
        expect(length(u)).toBeCloseTo(1, 15);

        const ur = unitPerp(v2(3e200, 4e200), true);
        expect(ur.values[0]).toBeCloseTo(0.8, 15);
        expect(ur.values[1]).toBeCloseTo(-0.6, 15);
    });

    it('unitPerp of the zero vector is the zero vector', () => {
        expect(unitPerp(v2(0, 0)).values).toEqual([0, 0]);
        expect(unitPerp(v2(0, 0), true).values).toEqual([0, 0]);
    });

    it('dotPerp((x0,x1),(y0,y1)) = x0*y1 - x1*y0', () => {
        expect(dotPerp(v2(1, 2), v2(3, 4))).toBe(1 * 4 - 2 * 3);
        expect(dotPerp(v2(1, 0), v2(0, 1))).toBe(1);
        expect(dotPerp(v2(0, 1), v2(1, 0))).toBe(-1);
        // dotPerp(v, v) = 0 (the determinant with repeated column).
        expect(dotPerp(v2(5, -7), v2(5, -7))).toBe(0);
    });
});

describe('computeOrthogonalComplement2', () => {
    it('builds a right-handed orthonormal basis from one vector', () => {
        const v = [v2(2, 1), new Vector(2)];
        const minLength = computeOrthogonalComplement2(1, v);
        expect(minLength).toBeGreaterThan(0);
        expect(length(v[0])).toBeCloseTo(1, 14);
        expect(length(v[1])).toBeCloseTo(1, 14);
        expect(dot(v[0], v[1])).toBeCloseTo(0, 14);
        // Right-handed: det[v0 v1] = dotPerp(v0, v1) = +1.
        expect(dotPerp(v[0], v[1])).toBeCloseTo(1, 14);
        // v[0] keeps the input direction.
        expect(v[0].values[0]).toBeCloseTo(2 / Math.sqrt(5), 14);
        expect(v[0].values[1]).toBeCloseTo(1 / Math.sqrt(5), 14);
    });

    it('returns 0 for numInputs other than 1', () => {
        const v = [v2(1, 0), v2(0, 1)];
        expect(computeOrthogonalComplement2(2, v)).toBe(0);
        expect(computeOrthogonalComplement2(0, v)).toBe(0);
    });

    it('returns a small length for a nearly zero input', () => {
        const v = [v2(1e-200, 0), new Vector(2)];
        const minLength = computeOrthogonalComplement2(1, v, true);
        expect(minLength).toBeLessThanOrEqual(1e-200);
    });
});

describe('computeBarycentrics2', () => {
    const t0 = v2(0, 0);
    const t1 = v2(1, 0);
    const t2 = v2(0, 1);

    it('computes known coordinates for the standard triangle', () => {
        const { valid, bary } = computeBarycentrics2(v2(0.25, 0.25), t0, t1, t2);
        expect(valid).toBe(true);
        expect(bary[0]).toBeCloseTo(0.5, 14);
        expect(bary[1]).toBeCloseTo(0.25, 14);
        expect(bary[2]).toBeCloseTo(0.25, 14);
    });

    it('vertices map to the canonical coordinates', () => {
        expect(computeBarycentrics2(t0, t0, t1, t2).bary[0]).toBeCloseTo(1, 14);
        expect(computeBarycentrics2(t1, t0, t1, t2).bary[1]).toBeCloseTo(1, 14);
        expect(computeBarycentrics2(t2, t0, t1, t2).bary[2]).toBeCloseTo(1, 14);
    });

    it('satisfies the partition of unity and reconstructs P', () => {
        const rng = makeRng(5678);
        const w0 = v2(-1, 2);
        const w1 = v2(3, 0.5);
        const w2 = v2(0.5, -4);
        for (let trial = 0; trial < 10; ++trial) {
            const p = v2((rng() - 0.5) * 8, (rng() - 0.5) * 8);
            const { valid, bary } = computeBarycentrics2(p, w0, w1, w2);
            expect(valid).toBe(true);
            expect(bary[0] + bary[1] + bary[2]).toBeCloseTo(1, 13);
            const q = add(add(mul(w0, bary[0]), mul(w1, bary[1])),
                mul(w2, bary[2]));
            expect(q.values[0]).toBeCloseTo(p.values[0], 12);
            expect(q.values[1]).toBeCloseTo(p.values[1], 12);
        }
    });

    it('returns invalid with zeroed coordinates for a degenerate triangle', () => {
        const { valid, bary } = computeBarycentrics2(v2(1, 1),
            v2(0, 0), v2(1, 1), v2(2, 2));
        expect(valid).toBe(false);
        expect(bary).toEqual([0, 0, 0]);
    });

    it('epsilon widens the degeneracy test', () => {
        const p = v2(0.1, 0.1);
        // A very thin triangle: det = 1e-6.
        const r0 = v2(0, 0);
        const r1 = v2(1, 0);
        const r2 = v2(1, 1e-6);
        expect(computeBarycentrics2(p, r0, r1, r2).valid).toBe(true);
        expect(computeBarycentrics2(p, r0, r1, r2, 1e-3).valid).toBe(false);
    });
});

describe('IntrinsicsVector2', () => {
    it('classifies a single point as dimension 0', () => {
        const info = new IntrinsicsVector2([v2(1, 2)], 0);
        expect(info.dimension).toBe(0);
        expect(info.origin.values).toEqual([1, 2]);
        expect(info.maxRange).toBe(0);
        expect(info.extreme).toEqual([0, 0, 0]);
    });

    it('classifies coincident points within epsilon as dimension 0', () => {
        const pts = [v2(0, 0), v2(1e-9, -1e-9), v2(-1e-9, 1e-9)];
        const info = new IntrinsicsVector2(pts, 1e-6);
        expect(info.dimension).toBe(0);
    });

    it('classifies collinear points as dimension 1 with the line direction', () => {
        const pts = [v2(0, 0), v2(1, 0), v2(3, 0), v2(2, 0)];
        const info = new IntrinsicsVector2(pts, 0);
        expect(info.dimension).toBe(1);
        expect(info.maxRange).toBe(3);
        expect(info.min).toEqual([0, 0]);
        expect(info.max).toEqual([3, 0]);
        expect(info.origin.values).toEqual([0, 0]);
        expect(info.extreme[0]).toBe(0);
        expect(info.extreme[1]).toBe(2);
        expect(info.extreme[2]).toBe(2);
        expect(info.direction[0].values[0]).toBeCloseTo(1, 15);
        expect(info.direction[0].values[1]).toBeCloseTo(0, 15);
        // direction[1] spans the orthogonal complement.
        expect(dot(info.direction[0], info.direction[1])).toBeCloseTo(0, 15);
    });

    it('classifies a genuine 2D set as dimension 2 with CCW extremes', () => {
        const pts = [v2(0, 0), v2(1, 0), v2(0, 1)];
        const info = new IntrinsicsVector2(pts, 0);
        expect(info.dimension).toBe(2);
        expect(info.extreme[0]).toBe(0);
        expect(info.extreme[1]).toBe(1);
        expect(info.extreme[2]).toBe(2);
        // (0,0) -> (1,0) -> (0,1) is counterclockwise.
        expect(info.extremeCCW).toBe(true);
        expect(length(info.direction[0])).toBeCloseTo(1, 14);
        expect(length(info.direction[1])).toBeCloseTo(1, 14);
    });

    it('detects clockwise extreme order', () => {
        const pts = [v2(0, 0), v2(1, 0), v2(0.5, -1)];
        const info = new IntrinsicsVector2(pts, 0);
        expect(info.dimension).toBe(2);
        expect(info.extremeCCW).toBe(false);
    });

    it('keeps zero defaults for invalid input', () => {
        const empty = new IntrinsicsVector2([], 0);
        expect(empty.dimension).toBe(0);
        expect(empty.maxRange).toBe(0);
        expect(empty.origin.values).toEqual([0, 0]);

        const badEps = new IntrinsicsVector2([v2(1, 2), v2(3, 4)], -1);
        expect(badEps.dimension).toBe(0);
        expect(badEps.origin.values).toEqual([0, 0]);
    });
});
