import { describe, it, expect } from 'vitest';
import { DistOrientedBox3OrientedBox3 } from '../src/DistOrientedBox3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function frame(w: Vector): [Vector, Vector, Vector] {
    const n = w.clone();
    normalize(n);
    let u = Math.abs(n.values[0]) > 0.5
        ? v(-n.values[1], n.values[0], 0)
        : v(0, -n.values[2], n.values[1]);
    normalize(u);
    return [u, cross(n, u), n];
}

function box(center: Vector, axis: Vector[], extent: Vector): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(center, axis, extent);
}

function axisAlignedBox(center: Vector, extent: Vector): OrientedBox {
    return box(center, [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], extent);
}

// The distance from a point to a solid oriented box, computed by clamping in
// the box frame.
function distPointBox(p: Vector, b: OrientedBox): number {
    const delta = sub(p, b.center);
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const y = dot(delta, b.axis[i]);
        const e = b.extent.values[i];
        const d = y < -e ? -e - y : (y > e ? y - e : 0);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

// Brute-force minimum of the distance from the surface of box0 to solid box1.
// The closest pair of two solid convex bodies always has a point on each
// boundary (or the bodies overlap and the distance is zero), so sampling the
// six faces of box0 converges to the true distance.
function bruteForce(box0: OrientedBox, box1: OrientedBox, n: number): number {
    let best = Infinity;
    for (let face = 0; face < 3; ++face) {
        const i0 = face, i1 = (face + 1) % 3, i2 = (face + 2) % 3;
        for (const sign of [-1, 1]) {
            const base = add(box0.center,
                mul(sign * box0.extent.values[i0], box0.axis[i0]));
            for (let i = 0; i <= n; ++i) {
                const e1 = box0.extent.values[i1];
                const p0 = add(base, mul(-e1 + 2 * e1 * i / n, box0.axis[i1]));
                for (let j = 0; j <= n; ++j) {
                    const e2 = box0.extent.values[i2];
                    const p = add(p0,
                        mul(-e2 + 2 * e2 * j / n, box0.axis[i2]));
                    const d = distPointBox(p, box1);
                    if (d < best) {
                        best = d;
                    }
                }
            }
        }
    }
    return best;
}

describe('DistOrientedBox3OrientedBox3', () => {
    it('reports the gap between two axis-aligned boxes', () => {
        const b0 = axisAlignedBox(v(0, 0, 0), v(1, 1, 1));
        const b1 = axisAlignedBox(v(5, 0, 0), v(1, 1, 1));
        const result = new DistOrientedBox3OrientedBox3().compute(b0, b1);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.sqrDistance).toBeCloseTo(9, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(4, 12);
    });

    it('reports zero for touching boxes', () => {
        const b0 = axisAlignedBox(v(0, 0, 0), v(1, 1, 1));
        const b1 = axisAlignedBox(v(2, 0, 0), v(1, 1, 1));
        const result = new DistOrientedBox3OrientedBox3().compute(b0, b1);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(length(sub(result.closest[0], result.closest[1])))
            .toBeCloseTo(0, 12);
    });

    it('reports zero for overlapping boxes', () => {
        const b0 = axisAlignedBox(v(0, 0, 0), v(1, 1, 1));
        const b1 = axisAlignedBox(v(1, 1, 1), v(1, 1, 1));
        const result = new DistOrientedBox3OrientedBox3().compute(b0, b1);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
    });

    it('measures a corner-to-corner separation', () => {
        // The boxes are diagonally offset, so the closest points are the
        // corners (1,1,1) and (3,3,3).
        const b0 = axisAlignedBox(v(0, 0, 0), v(1, 1, 1));
        const b1 = axisAlignedBox(v(4, 4, 4), v(1, 1, 1));
        const result = new DistOrientedBox3OrientedBox3().compute(b0, b1);
        expect(result.distance).toBeCloseTo(2 * Math.sqrt(3), 12);
        expect(result.closest[0].values).toEqual([1, 1, 1]);
        expect(result.closest[1].values).toEqual([3, 3, 3]);
    });

    it('accounts for a rotated box', () => {
        // Box1 is a cube rotated 45 degrees about the z-axis, so its extreme
        // point toward -x is at x = 5 - sqrt(2).
        const s = Math.SQRT1_2;
        const b0 = axisAlignedBox(v(0, 0, 0), v(1, 1, 1));
        const b1 = box(v(5, 0, 0), [v(s, s, 0), v(-s, s, 0), v(0, 0, 1)],
            v(1, 1, 1));
        const result = new DistOrientedBox3OrientedBox3().compute(b0, b1);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(5 - Math.SQRT2, 12);
    });

    it('is symmetric in its arguments', () => {
        const [u0, u1, u2] = frame(v(1, -2, 0.5));
        const b0 = box(v(0.25, -0.5, 1), [u0, u1, u2], v(1, 0.5, 2));
        const [w0, w1, w2] = frame(v(-3, 1, 2));
        const b1 = box(v(5, 2, -1), [w0, w1, w2], v(0.75, 1.5, 0.5));
        const query = new DistOrientedBox3OrientedBox3();
        const forward = query.compute(b0, b1);
        const backward = query.compute(b1, b0);
        expect(backward.distance).toBeCloseTo(forward.distance, 10);
        expect(length(sub(backward.closest[1], forward.closest[0])))
            .toBeCloseTo(0, 8);
        expect(length(sub(backward.closest[0], forward.closest[1])))
            .toBeCloseTo(0, 8);
    });

    it('agrees with brute-force sampling and reports consistent points', () => {
        let seed = 555444333;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistOrientedBox3OrientedBox3();
        for (let trial = 0; trial < 25; ++trial) {
            const mk = () => {
                const w = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
                if (length(w) < 1e-3) {
                    return null;
                }
                const [a0, a1, a2] = frame(w);
                return box(v(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3),
                    [a0, a1, a2],
                    v(0.3 + rnd(), 0.3 + rnd(), 0.3 + rnd()));
            };
            const b0 = mk(), b1 = mk();
            if (!b0 || !b1) {
                continue;
            }

            const result = query.compute(b0, b1);
            const brute = Math.min(bruteForce(b0, b1, 24),
                bruteForce(b1, b0, 24));
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.25);
            expect(result.sqrDistance)
                .toBeCloseTo(result.distance * result.distance, 9);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            expect(distPointBox(result.closest[0], b0)).toBeCloseTo(0, 8);
            expect(distPointBox(result.closest[1], b1)).toBeCloseTo(0, 8);
        }
    });
});
