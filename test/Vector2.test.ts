import { describe, it, expect } from 'vitest';
import {
    check, vector, expectClose, expectVectorClose, fc
} from './helpers/arbitraries.js';
import {
    perp, unitPerp, dotPerp, computeOrthogonalComplement2,
    computeBarycentrics2, IntrinsicsVector2
} from '../src/Vector2.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';

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

// ---------------------------------------------------------------------------
// Verification wave: property-based checks against upstream Vector2.h.
// ---------------------------------------------------------------------------

describe('Vector2 verification', () => {
    it('perp is the exact 90-degree clockwise rotation', () => {
        check(vector(2), v => {
            const p = perp(v);
            // perp = (x1, -x0): exact component moves, so these are exact.
            expect(p.values).toEqual([v.get(1), -v.get(0)]);
            expect(dot(v, p)).toBe(0);
            // perp(perp(v)) = -v.
            expect(perp(p).values).toEqual([-v.get(0), -v.get(1)]);
            expect(length(p)).toBe(length(v));
        });
    });

    it('dotPerp is the 2x2 determinant and is antisymmetric', () => {
        check(fc.tuple(vector(2), vector(2)), ([a, b]) => {
            // Dot(v0, Perp(v1)) = x0*y1 - x1*y0 with v1 = (y0, y1).
            expect(dotPerp(a, b)).toBe(a.get(0) * b.get(1) - a.get(1) * b.get(0));
            // IEEE subtraction is sign-symmetric, so the swap is exact.
            expect(dotPerp(b, a) + dotPerp(a, b)).toBe(0);
            expect(dotPerp(a, a)).toBe(0);
        });
    });

    it('unitPerp is a unit vector orthogonal to the input', () => {
        check(fc.tuple(vector(2).filter(v => length(v) > 1e-2), fc.boolean()),
            ([v, robust]) => {
                const u = unitPerp(v, robust);
                expectClose(length(u), 1, 1e-12, 1e-12);
                // The two components are scaled by 1/L independently, so the
                // cancellation is only up to rounding of |v|.
                expectClose(dot(u, v), 0, 1e-12 * length(v), 0);
                // The perpendicular points to the "right" of v.
                expect(dotPerp(v, u)).toBeLessThan(0);
            });
        // Upstream Normalize sets the zero vector to zero.
        expect(unitPerp(Vector.fromArray([0, 0])).values).toEqual([0, 0]);
    });

    it('computeOrthogonalComplement2 builds a right-handed frame', () => {
        check(fc.tuple(vector(2).filter(v => length(v) > 1e-2), fc.boolean()),
            ([v0, robust]) => {
                const v = [v0.clone()];
                const minLength = computeOrthogonalComplement2(1, v, robust);
                expectClose(minLength, length(v0), 1e-9, 1e-9);
                expect(v.length).toBe(2);
                expectClose(length(v[0]), 1, 1e-9, 1e-9);
                expectClose(length(v[1]), 1, 1e-9, 1e-9);
                expectClose(dot(v[0], v[1]), 0, 1e-9, 1e-9);
                // v[1] = -Perp(v[0]) = (-y, x), so det[v0 v1] = +1.
                expectClose(dotPerp(v[0], v[1]), 1, 1e-9, 1e-9);
            });
        // Any numInputs other than 1 is rejected with a 0 return.
        expect(computeOrthogonalComplement2(0, [Vector.fromArray([1, 0])]))
            .toBe(0);
        expect(computeOrthogonalComplement2(2, [Vector.fromArray([1, 0])]))
            .toBe(0);
    });

    it('computeBarycentrics2 reconstructs the point and sums to one', () => {
        check(fc.tuple(vector(2), vector(2), vector(2), vector(2)),
            ([p, a, b, c]) => {
                const r = computeBarycentrics2(p, a, b, c);
                if (!r.valid) {
                    expect(r.bary).toEqual([0, 0, 0]);
                    return;
                }
                const det = Math.abs(dotPerp(sub(a, c), sub(b, c)));
                if (det < 1e-3) {
                    return;   // ill-conditioned: the division amplifies error
                }
                expectClose(r.bary[0] + r.bary[1] + r.bary[2], 1, 1e-9, 1e-9);
                const q = add(add(mul(a, r.bary[0]), mul(b, r.bary[1])),
                    mul(c, r.bary[2]));
                expectVectorClose(q, p, 1e-6, 1e-6);
            });
    });

    it('computeBarycentrics2 gives the unit coordinates at the vertices', () => {
        const a = Vector.fromArray([0, 0]);
        const b = Vector.fromArray([4, 0]);
        const c = Vector.fromArray([0, 2]);
        // '+ 0' normalizes the signed zeros that the divisions produce.
        const bary = (p: Vector) =>
            computeBarycentrics2(p, a, b, c).bary.map(x => x + 0);
        expect(bary(a)).toEqual([1, 0, 0]);
        expect(bary(b)).toEqual([0, 1, 0]);
        expect(bary(c)).toEqual([0, 0, 1]);
        // Degenerate (collinear) triangle: invalid with zeroed coordinates.
        const d = Vector.fromArray([8, 0]);
        const r = computeBarycentrics2(a, a, b, d);
        expect(r.valid).toBe(false);
        expect(r.bary).toEqual([0, 0, 0]);
        // The epsilon test is |det| <= epsilon, i.e. strict rejection.
        expect(computeBarycentrics2(a, a, b, c, 8).valid).toBe(false);
        expect(computeBarycentrics2(a, a, b, c, 7.9).valid).toBe(true);
    });

    it('IntrinsicsVector2 classifies the intrinsic dimension', () => {
        // d = 0: all points identical.
        const same = [Vector.fromArray([3, -1]), Vector.fromArray([3, -1])];
        const i0 = new IntrinsicsVector2(same, 0);
        expect(i0.dimension).toBe(0);
        expect(i0.extreme).toEqual([0, 0, 0]);
        expect(i0.origin.values).toEqual([3, -1]);

        // d = 1: collinear points.
        const line = [Vector.fromArray([0, 0]), Vector.fromArray([1, 1]),
            Vector.fromArray([2, 2]), Vector.fromArray([-1, -1])];
        const i1 = new IntrinsicsVector2(line, 0);
        expect(i1.dimension).toBe(1);
        expect(i1.extreme[2]).toBe(i1.extreme[1]);

        // d = 2: a genuine triangle; extremeCCW records the orientation of
        // (V[e0], V[e1], V[e2]) relative to direction[1] = -Perp(direction[0]).
        const tri = [Vector.fromArray([0, 0]), Vector.fromArray([4, 0]),
            Vector.fromArray([0, 3])];
        const i2 = new IntrinsicsVector2(tri, 0);
        expect(i2.dimension).toBe(2);
        expect(i2.extremeCCW).toBe(
            dot(i2.direction[1], sub(tri[i2.extreme[2]], i2.origin)) > 0);
        expect(i2.maxRange).toBe(4);
    });

    it('IntrinsicsVector2 members agree with a brute-force computation', () => {
        check(fc.array(vector(2), { minLength: 1, maxLength: 10 }), vs => {
            const ins = new IntrinsicsVector2(vs, 0);
            // The bounding box is exact.
            for (let j = 0; j < 2; ++j) {
                let lo = vs[0].get(j), hi = vs[0].get(j);
                for (const v of vs) {
                    lo = Math.min(lo, v.get(j));
                    hi = Math.max(hi, v.get(j));
                }
                expect(ins.min[j]).toBe(lo);
                expect(ins.max[j]).toBe(hi);
            }
            expect(ins.maxRange).toBe(Math.max(ins.max[0] - ins.min[0],
                ins.max[1] - ins.min[1]));
            expect(ins.origin.values).toEqual(vs[ins.extreme[0]].values);
            if (ins.dimension === 2) {
                // direction[0] and direction[1] are an orthonormal frame.
                expectClose(length(ins.direction[0]), 1, 1e-9, 1e-9);
                expectClose(length(ins.direction[1]), 1, 1e-9, 1e-9);
                expect(dot(ins.direction[0], ins.direction[1])).toBe(0);
            }
        });
    });

    it('IntrinsicsVector2 leaves the defaults for invalid inputs', () => {
        // Upstream returns without touching the members when numVectors == 0
        // or epsilon < 0.
        for (const ins of [new IntrinsicsVector2([], 0),
            new IntrinsicsVector2([Vector.fromArray([1, 2])], -1)]) {
            expect(ins.dimension).toBe(0);
            expect(ins.maxRange).toBe(0);
            expect(ins.min).toEqual([0, 0]);
            expect(ins.max).toEqual([0, 0]);
            expect(ins.origin.values).toEqual([0, 0]);
            expect(ins.extreme).toEqual([0, 0, 0]);
            expect(ins.extremeCCW).toBe(false);
        }
    });

    it('the size-2 requirement is enforced', () => {
        expect(() => perp(Vector.fromArray([1, 2, 3])))
            .toThrow('must have size 2');
        expect(() => new IntrinsicsVector2([Vector.fromArray([1, 2, 3])], 0))
            .toThrow('must have size 2');
    });
});
