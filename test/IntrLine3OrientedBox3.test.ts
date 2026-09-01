import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { OrientedBox } from '../src/OrientedBox';
import { Line } from '../src/Line';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';
import { computeOrthogonalComplement3 } from '../src/Vector3';
import {
    IntrLine3AlignedBox3TI,
    IntrLine3AlignedBox3FI
} from '../src/IntrLine3AlignedBox3';
import {
    IntrLine3OrientedBox3TI,
    IntrLine3OrientedBox3FI
} from '../src/IntrLine3OrientedBox3';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function line(p: number[], d: number[]): Line {
    const dir = vec(d);
    normalize(dir);
    return Line.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Build an oriented box whose first axis is the given (unnormalized) vector.
function orientedBox(center: number[], axis0: number[], extent: number[]):
    OrientedBox {
    const a0 = vec(axis0);
    normalize(a0);
    const basis = [a0, new Vector(3), new Vector(3)];
    computeOrthogonalComplement3(1, basis);
    return OrientedBox.fromCenterAxisExtent(vec(center), basis, vec(extent));
}

function insideBox(box: OrientedBox, x: Vector): boolean {
    const diff = sub(x, box.center);
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(dot(diff, box.axis[i])) > box.extent.values[i] + 1e-12) {
            return false;
        }
    }
    return true;
}

describe('IntrLine3OrientedBox3', () => {
    const ti = new IntrLine3OrientedBox3TI();
    const fi = new IntrLine3OrientedBox3FI();

    it('matches the aligned-box query when the axes are the standard basis', () => {
        const abox = AlignedBox.fromMinMax(vec([-1, -2, -3]), vec([1, 2, 3]));
        const obox = OrientedBox.fromCenterAxisExtent(
            vec([0, 0, 0]),
            [vec([1, 0, 0]), vec([0, 1, 0]), vec([0, 0, 1])],
            vec([1, 2, 3]));
        const aTI = new IntrLine3AlignedBox3TI();
        const aFI = new IntrLine3AlignedBox3FI();
        const rand = makeRandom(12345);
        for (let trial = 0; trial < 200; ++trial) {
            const l = line(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const expectedTI = aTI.test(l, abox);
            const expectedFI = aFI.find(l, abox);
            expect(ti.test(l, obox).intersect).toBe(expectedTI.intersect);
            const got = fi.find(l, obox);
            expect(got.intersect).toBe(expectedFI.intersect);
            expect(got.numIntersections).toBe(expectedFI.numIntersections);
            if (got.intersect) {
                expect(got.parameter[0]).toBeCloseTo(expectedFI.parameter[0], 12);
                expect(got.parameter[1]).toBeCloseTo(expectedFI.parameter[1], 12);
            }
        }
    });

    it('finds the crossing of a rotated box along its own first axis', () => {
        // The box axis[0] is (1,1,0)/sqrt(2) with extent 2, so the line
        // through the center along that axis enters at t = -2 and exits at
        // t = 2.
        const box = orientedBox([1, 2, 3], [1, 1, 0], [2, 1, 1]);
        const l = Line.fromOriginDirection(box.center.clone(),
            box.axis[0].clone());
        const result = fi.find(l, box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
        expect(ti.test(l, box).intersect).toBe(true);
    });

    it('reports no intersection for a line that misses a rotated box', () => {
        const box = orientedBox([0, 0, 0], [1, 1, 1], [1, 1, 1]);
        const l = line([10, 0, 0], [0, 0, 1]);
        expect(ti.test(l, box).intersect).toBe(false);
        expect(fi.find(l, box).intersect).toBe(false);
    });

    it('reports a single point for a line tangent to a box face', () => {
        // The box is axis aligned with extent 1; the line lies in the plane
        // x = 1 and is parallel to the z-axis, so it touches the face.
        const box = OrientedBox.fromCenterAxisExtent(
            vec([0, 0, 0]),
            [vec([1, 0, 0]), vec([0, 1, 0]), vec([0, 0, 1])],
            vec([1, 1, 1]));
        const l = line([1, 1, 0], [0, 0, 1]);
        const result = fi.find(l, box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
    });

    it('agrees with dense sampling along the line', () => {
        const rand = makeRandom(777);
        const box = orientedBox([0.5, -0.25, 1], [2, -1, 0.5], [1.5, 1, 0.75]);
        for (let trial = 0; trial < 100; ++trial) {
            const l = line(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(l, box);
            expect(ti.test(l, box).intersect).toBe(result.intersect);

            // Sample the line and find the first/last inside parameters.
            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = -10 + (20 * k) / n;
                const x = add(l.origin, mul(t, l.direction));
                if (insideBox(box, x)) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                // The sampled bracket is inside the exact interval and within
                // one sample spacing of it.
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(2e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(2e-3);
            }
        }
    });
});
