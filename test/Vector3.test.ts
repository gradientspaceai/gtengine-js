import { describe, it, expect } from 'vitest';
import {
    cross, unitCross, dotCross, computeOrthogonalComplement3,
    fastComputeOrthogonalComplement, computeBarycentrics3, IntrinsicsVector3
} from '../src/Vector3.js';
import { Vector, add, dot, length, mul } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

function randomV3(rng: () => number): Vector {
    return v3((rng() - 0.5) * 10, (rng() - 0.5) * 10, (rng() - 0.5) * 10);
}

describe('cross, unitCross, dotCross', () => {
    it('reproduces the right-handed basis products', () => {
        const e0 = v3(1, 0, 0);
        const e1 = v3(0, 1, 0);
        const e2 = v3(0, 0, 1);
        expect(cross(e0, e1).values).toEqual([0, 0, 1]);
        expect(cross(e1, e2).values).toEqual([1, 0, 0]);
        expect(cross(e2, e0).values).toEqual([0, 1, 0]);
    });

    it('computes the formal determinant for known values', () => {
        const c = cross(v3(1, 2, 3), v3(4, 5, 6));
        expect(c.values).toEqual([2 * 6 - 3 * 5, 3 * 4 - 1 * 6, 1 * 5 - 2 * 4]);
    });

    it('is anticommutative and orthogonal to both inputs', () => {
        const rng = makeRng(42);
        for (let trial = 0; trial < 10; ++trial) {
            const a = randomV3(rng);
            const b = randomV3(rng);
            const c = cross(a, b);
            expect(dot(c, a)).toBeCloseTo(0, 12);
            expect(dot(c, b)).toBeCloseTo(0, 12);
            const cba = cross(b, a);
            for (let i = 0; i < 3; ++i) {
                expect(cba.values[i]).toBeCloseTo(-c.values[i], 12);
            }
            // Lagrange identity: |a x b|^2 = |a|^2 |b|^2 - (a.b)^2.
            expect(dot(c, c)).toBeCloseTo(
                dot(a, a) * dot(b, b) - dot(a, b) * dot(a, b), 9);
        }
    });

    it('supports affine 4-tuples with w = 0', () => {
        const a = Vector.fromArray([1, 0, 0, 0]);
        const b = Vector.fromArray([0, 1, 0, 0]);
        const c = cross(a, b);
        expect(c.size).toBe(4);
        expect(c.values).toEqual([0, 0, 1, 0]);
        // The w components of the inputs are ignored; the result has w = 0.
        const d = cross(Vector.fromArray([1, 0, 0, 5]),
            Vector.fromArray([0, 1, 0, 7]));
        expect(d.values).toEqual([0, 0, 1, 0]);
    });

    it('rejects invalid dimensions', () => {
        expect(() => cross(Vector.fromArray([1, 2]), Vector.fromArray([3, 4])))
            .toThrow('Dimension must be 3 or 4.');
        expect(() => cross(v3(1, 2, 3), Vector.fromArray([1, 2, 3, 0])))
            .toThrow('mismatched sizes');
    });

    it('unitCross is the normalized cross product', () => {
        const u = unitCross(v3(2, 0, 0), v3(0, 3, 0));
        expect(u.values).toEqual([0, 0, 1]);
        // The cross product magnitude 6e300 overflows the non-robust
        // squared-length computation; robust normalization handles it.
        const ur = unitCross(v3(2e150, 0, 0), v3(0, 3e150, 0), true);
        expect(ur.values[2]).toBeCloseTo(1, 15);
        // Parallel inputs produce the zero vector.
        expect(unitCross(v3(1, 2, 3), v3(2, 4, 6)).values).toEqual([0, 0, 0]);
    });

    it('dotCross is the triple scalar product (a determinant)', () => {
        expect(dotCross(v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1))).toBe(1);
        expect(dotCross(v3(0, 1, 0), v3(1, 0, 0), v3(0, 0, 1))).toBe(-1);
        // Cyclic invariance: [a,b,c] = [b,c,a] = [c,a,b].
        const rng = makeRng(99);
        const a = randomV3(rng);
        const b = randomV3(rng);
        const c = randomV3(rng);
        const abc = dotCross(a, b, c);
        expect(dotCross(b, c, a)).toBeCloseTo(abc, 10);
        expect(dotCross(c, a, b)).toBeCloseTo(abc, 10);
        // Repeated argument gives zero.
        expect(dotCross(a, a, b)).toBeCloseTo(0, 10);
    });
});

function expectRightHandedOrthonormal(v: Vector[]): void {
    for (let i = 0; i < 3; ++i) {
        expect(length(v[i])).toBeCloseTo(1, 13);
        for (let j = i + 1; j < 3; ++j) {
            expect(dot(v[i], v[j])).toBeCloseTo(0, 13);
        }
    }
    expect(dotCross(v[0], v[1], v[2])).toBeCloseTo(1, 13);
}

describe('computeOrthogonalComplement3', () => {
    it('builds a right-handed orthonormal basis from one vector', () => {
        const v = [v3(1, 2, 3), new Vector(3), new Vector(3)];
        const minLength = computeOrthogonalComplement3(1, v);
        expect(minLength).toBeGreaterThan(0);
        expectRightHandedOrthonormal(v);
        // v[0] keeps the input direction.
        const invLen = 1 / Math.sqrt(14);
        expect(v[0].values[0]).toBeCloseTo(1 * invLen, 14);
        expect(v[0].values[1]).toBeCloseTo(2 * invLen, 14);
        expect(v[0].values[2]).toBeCloseTo(3 * invLen, 14);
    });

    it('covers both branches of the initial perpendicular choice', () => {
        // |x| > |y| branch.
        const va = [v3(5, 1, 2), new Vector(3), new Vector(3)];
        computeOrthogonalComplement3(1, va);
        expectRightHandedOrthonormal(va);
        // |x| <= |y| branch.
        const vb = [v3(1, 5, 2), new Vector(3), new Vector(3)];
        computeOrthogonalComplement3(1, vb);
        expectRightHandedOrthonormal(vb);
    });

    it('builds the basis from two vectors', () => {
        const v = [v3(1, 0, 0), v3(1, 1, 0), new Vector(3)];
        const minLength = computeOrthogonalComplement3(2, v);
        expect(minLength).toBeGreaterThan(0);
        expectRightHandedOrthonormal(v);
        // The third vector completes {(1,0,0),(0,1,0)} right-handedly.
        expect(v[2].values[2]).toBeCloseTo(1, 14);
    });

    it('returns 0 for numInputs other than 1 or 2', () => {
        const v = [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        expect(computeOrthogonalComplement3(0, v)).toBe(0);
        expect(computeOrthogonalComplement3(3, v)).toBe(0);
    });
});

describe('fastComputeOrthogonalComplement', () => {
    it('produces a right-handed orthonormal basis (v2[2] >= 0 branch)', () => {
        const rng = makeRng(7);
        for (let trial = 0; trial < 10; ++trial) {
            const w = randomV3(rng);
            w.values[2] = Math.abs(w.values[2]);
            const len = length(w);
            const unit = mul(w, 1 / len);
            const { v0, v1 } = fastComputeOrthogonalComplement(unit);
            expectRightHandedOrthonormal([v0, v1, unit]);
            const c = cross(v0, v1);
            for (let i = 0; i < 3; ++i) {
                expect(c.values[i]).toBeCloseTo(unit.values[i], 13);
            }
        }
    });

    it('produces a right-handed orthonormal basis (v2[2] < 0 branch)', () => {
        const rng = makeRng(11);
        for (let trial = 0; trial < 10; ++trial) {
            const w = randomV3(rng);
            w.values[2] = -Math.abs(w.values[2]) - 0.1;
            const len = length(w);
            const unit = mul(w, 1 / len);
            const { v0, v1 } = fastComputeOrthogonalComplement(unit);
            expectRightHandedOrthonormal([v0, v1, unit]);
        }
    });

    it('handles the poles exactly', () => {
        const north = fastComputeOrthogonalComplement(v3(0, 0, 1));
        // Some components are -0; compare numerically.
        expect(north.v0.values[0]).toBe(1);
        expect(Math.abs(north.v0.values[1])).toBe(0);
        expect(Math.abs(north.v0.values[2])).toBe(0);
        expect(Math.abs(north.v1.values[0])).toBe(0);
        expect(north.v1.values[1]).toBe(1);
        expect(Math.abs(north.v1.values[2])).toBe(0);
        const south = fastComputeOrthogonalComplement(v3(0, 0, -1));
        const c = cross(south.v0, south.v1);
        expect(c.values[2]).toBeCloseTo(-1, 15);
    });
});

describe('computeBarycentrics3', () => {
    const t0 = v3(0, 0, 0);
    const t1 = v3(1, 0, 0);
    const t2 = v3(0, 1, 0);
    const t3 = v3(0, 0, 1);

    it('computes known coordinates for the standard tetrahedron', () => {
        const { valid, bary } =
            computeBarycentrics3(v3(0.25, 0.25, 0.25), t0, t1, t2, t3);
        expect(valid).toBe(true);
        expect(bary[0]).toBeCloseTo(0.25, 14);
        expect(bary[1]).toBeCloseTo(0.25, 14);
        expect(bary[2]).toBeCloseTo(0.25, 14);
        expect(bary[3]).toBeCloseTo(0.25, 14);
    });

    it('vertices map to the canonical coordinates', () => {
        expect(computeBarycentrics3(t0, t0, t1, t2, t3).bary[0])
            .toBeCloseTo(1, 14);
        expect(computeBarycentrics3(t1, t0, t1, t2, t3).bary[1])
            .toBeCloseTo(1, 14);
        expect(computeBarycentrics3(t2, t0, t1, t2, t3).bary[2])
            .toBeCloseTo(1, 14);
        expect(computeBarycentrics3(t3, t0, t1, t2, t3).bary[3])
            .toBeCloseTo(1, 14);
    });

    it('satisfies the partition of unity and reconstructs P', () => {
        const rng = makeRng(31415);
        const w0 = v3(-1, 2, 0.5);
        const w1 = v3(3, 0.5, -1);
        const w2 = v3(0.5, -4, 2);
        const w3 = v3(1, 1, 5);
        for (let trial = 0; trial < 10; ++trial) {
            const p = randomV3(rng);
            const { valid, bary } = computeBarycentrics3(p, w0, w1, w2, w3);
            expect(valid).toBe(true);
            expect(bary[0] + bary[1] + bary[2] + bary[3]).toBeCloseTo(1, 12);
            const q = add(add(mul(w0, bary[0]), mul(w1, bary[1])),
                add(mul(w2, bary[2]), mul(w3, bary[3])));
            for (let i = 0; i < 3; ++i) {
                expect(q.values[i]).toBeCloseTo(p.values[i], 11);
            }
        }
    });

    it('returns invalid with zeroed coordinates for coplanar vertices', () => {
        const { valid, bary } = computeBarycentrics3(v3(0.1, 0.2, 0),
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0));
        expect(valid).toBe(false);
        expect(bary).toEqual([0, 0, 0, 0]);
    });

    it('epsilon widens the degeneracy test', () => {
        // A very flat tetrahedron: det = 1e-6.
        const r3 = v3(0, 0, 1e-6);
        const p = v3(0.1, 0.1, 0);
        expect(computeBarycentrics3(p, t0, t1, t2, r3).valid).toBe(true);
        expect(computeBarycentrics3(p, t0, t1, t2, r3, 1e-3).valid).toBe(false);
    });
});

describe('IntrinsicsVector3', () => {
    it('classifies a single point as dimension 0', () => {
        const info = new IntrinsicsVector3([v3(1, 2, 3)], 0);
        expect(info.dimension).toBe(0);
        expect(info.origin.values).toEqual([1, 2, 3]);
        expect(info.maxRange).toBe(0);
        expect(info.extreme).toEqual([0, 0, 0, 0]);
    });

    it('classifies collinear points as dimension 1', () => {
        const pts = [v3(0, 0, 0), v3(0, 0, 1), v3(0, 0, 5), v3(0, 0, 3)];
        const info = new IntrinsicsVector3(pts, 0);
        expect(info.dimension).toBe(1);
        expect(info.maxRange).toBe(5);
        expect(info.origin.values).toEqual([0, 0, 0]);
        expect(info.extreme[0]).toBe(0);
        expect(info.extreme[1]).toBe(2);
        expect(info.extreme[2]).toBe(2);
        expect(info.extreme[3]).toBe(2);
        expect(Math.abs(info.direction[0].values[2])).toBeCloseTo(1, 15);
        // The remaining directions span the orthogonal complement.
        expect(dot(info.direction[0], info.direction[1])).toBeCloseTo(0, 14);
        expect(dot(info.direction[0], info.direction[2])).toBeCloseTo(0, 14);
    });

    it('classifies coplanar points as dimension 2', () => {
        const pts = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0),
            v3(0.3, 0.7, 0)];
        const info = new IntrinsicsVector3(pts, 0);
        expect(info.dimension).toBe(2);
        // direction[2] is normal to the plane z = 0.
        expect(Math.abs(info.direction[2].values[2])).toBeCloseTo(1, 14);
        expect(Math.abs(info.direction[2].values[0])).toBeCloseTo(0, 14);
        expect(Math.abs(info.direction[2].values[1])).toBeCloseTo(0, 14);
        expect(info.extreme[3]).toBe(info.extreme[2]);
    });

    it('classifies a full-dimensional set as dimension 3', () => {
        const pts = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0),
            v3(0, 0, 1)];
        const info = new IntrinsicsVector3(pts, 0);
        expect(info.dimension).toBe(3);
        expect(info.extreme[0]).toBe(0);
        expect(info.extreme[1]).toBe(1);
        expect(info.extreme[2]).toBe(2);
        expect(info.extreme[3]).toBe(4);
        // (0,0,0),(1,0,0),(0,1,0),(0,0,1) is counterclockwise.
        expect(info.extremeCCW).toBe(true);
        expectRightHandedOrthonormalDirections(info.direction);
    });

    it('detects clockwise extreme order', () => {
        const pts = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0),
            v3(0, 0, -1)];
        const info = new IntrinsicsVector3(pts, 0);
        expect(info.dimension).toBe(3);
        expect(info.extremeCCW).toBe(false);
    });

    it('uses epsilon*maxRange as the relative flatness threshold', () => {
        // A plane-like cloud with tiny z-jitter is planar for a loose epsilon
        // but 3D for epsilon = 0.
        const pts = [v3(0, 0, 0), v3(10, 0, 0), v3(0, 10, 0), v3(10, 10, 1e-8)];
        expect(new IntrinsicsVector3(pts, 0).dimension).toBe(3);
        expect(new IntrinsicsVector3(pts, 1e-6).dimension).toBe(2);
    });

    it('keeps zero defaults for invalid input', () => {
        const empty = new IntrinsicsVector3([], 0);
        expect(empty.dimension).toBe(0);
        expect(empty.origin.values).toEqual([0, 0, 0]);

        const badEps = new IntrinsicsVector3([v3(1, 2, 3), v3(4, 5, 6)], -1);
        expect(badEps.dimension).toBe(0);
        expect(badEps.origin.values).toEqual([0, 0, 0]);
    });
});

function expectRightHandedOrthonormalDirections(
    d: [Vector, Vector, Vector]): void {
    for (let i = 0; i < 3; ++i) {
        expect(length(d[i])).toBeCloseTo(1, 13);
        for (let j = i + 1; j < 3; ++j) {
            expect(dot(d[i], d[j])).toBeCloseTo(0, 13);
        }
    }
    expect(dotCross(d[0], d[1], d[2])).toBeCloseTo(1, 13);
}
