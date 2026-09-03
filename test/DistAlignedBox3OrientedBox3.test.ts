import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistAlignedBox3OrientedBox3 } from '../src/DistAlignedBox3OrientedBox3.js';
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

function distPointOrientedBox(p: Vector, b: OrientedBox): number {
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

function distPointAlignedBox(p: Vector, b: AlignedBox): number {
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const x = p.values[i];
        const lo = b.min.values[i], hi = b.max.values[i];
        const d = x < lo ? lo - x : (x > hi ? x - hi : 0);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

// Brute-force minimum from the surface of the aligned box to the solid
// oriented box.
function bruteForce(box0: AlignedBox, box1: OrientedBox, n: number): number {
    let best = Infinity;
    for (let face = 0; face < 3; ++face) {
        const i0 = face, i1 = (face + 1) % 3, i2 = (face + 2) % 3;
        for (const which of [box0.min, box0.max]) {
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    const p = new Vector(3);
                    p.values[i0] = which.values[i0];
                    p.values[i1] = box0.min.values[i1] + i / n *
                        (box0.max.values[i1] - box0.min.values[i1]);
                    p.values[i2] = box0.min.values[i2] + j / n *
                        (box0.max.values[i2] - box0.min.values[i2]);
                    const d = distPointOrientedBox(p, box1);
                    if (d < best) {
                        best = d;
                    }
                }
            }
        }
    }
    return best;
}

describe('DistAlignedBox3OrientedBox3', () => {
    it('reports the gap between an aligned box and an aligned oriented box',
        () => {
            const box0 = AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));
            const box1 = OrientedBox.fromCenterAxisExtent(v(5, 0, 0),
                [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(1, 1, 1));
            const result =
                new DistAlignedBox3OrientedBox3().compute(box0, box1);
            expect(result.distance).toBeCloseTo(3, 12);
            expect(result.closest[0].values[0]).toBeCloseTo(1, 12);
            expect(result.closest[1].values[0]).toBeCloseTo(4, 12);
        });

    it('handles an aligned box that is not centered at the origin', () => {
        const box0 = AlignedBox.fromMinMax(v(10, 20, 30), v(12, 24, 36));
        const box1 = OrientedBox.fromCenterAxisExtent(v(11, 22, 40),
            [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(0.5, 0.5, 1));
        const result = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(36, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(39, 12);
        expect(distPointAlignedBox(result.closest[0], box0)).toBeCloseTo(0, 8);
    });

    it('accounts for the rotation of the oriented box', () => {
        const s = Math.SQRT1_2;
        const box0 = AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));
        const box1 = OrientedBox.fromCenterAxisExtent(v(5, 0, 0),
            [v(s, s, 0), v(-s, s, 0), v(0, 0, 1)], v(1, 1, 1));
        const result = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 12);
    });

    it('reports zero when the boxes overlap', () => {
        const s = Math.SQRT1_2;
        const box0 = AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));
        const box1 = OrientedBox.fromCenterAxisExtent(v(1, 1, 0),
            [v(s, s, 0), v(-s, s, 0), v(0, 0, 1)], v(1, 1, 1));
        const result = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
    });

    it('matches the box-box query on the equivalent oriented box', () => {
        let seed = 246813579;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistAlignedBox3OrientedBox3();
        const reference = new DistOrientedBox3OrientedBox3();
        for (let trial = 0; trial < 25; ++trial) {
            const lo = v(rnd() * 4 - 3, rnd() * 4 - 3, rnd() * 4 - 3);
            const box0 = AlignedBox.fromMinMax(lo,
                add(lo, v(0.5 + rnd() * 2, 0.5 + rnd() * 2, 0.5 + rnd() * 2)));

            const w = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (length(w) < 1e-3) {
                continue;
            }
            const [a0, a1, a2] = frame(w);
            const box1 = OrientedBox.fromCenterAxisExtent(
                v(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3), [a0, a1, a2],
                v(0.3 + rnd(), 0.3 + rnd(), 0.3 + rnd()));

            const result = query.compute(box0, box1);

            // The same query, expressed with box0 converted by hand.
            const centered = box0.getCenteredForm();
            const obox0 = OrientedBox.fromCenterAxisExtent(centered.center,
                [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], centered.extent);
            const expected = reference.compute(obox0, box1);
            expect(result.distance).toBe(expected.distance);

            // Cross-check against brute-force sampling of the aligned box.
            const brute = bruteForce(box0, box1, 24);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.25);

            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            expect(distPointAlignedBox(result.closest[0], box0))
                .toBeCloseTo(0, 8);
            expect(distPointOrientedBox(result.closest[1], box1))
                .toBeCloseTo(0, 8);
        }
    });

    it('handles a degenerate aligned box that is a single point', () => {
        // A degenerate aligned box (a point) still produces a valid query.
        const box0 = AlignedBox.fromMinMax(v(2, 2, 2), v(2, 2, 2));
        const box1 = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(1, 1, 1));
        const result = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        expect(result.distance).toBeCloseTo(Math.sqrt(3), 12);
        expect(result.closest[0].values).toEqual([2, 2, 2]);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });
});
