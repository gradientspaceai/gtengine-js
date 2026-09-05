import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import {
    IntrLine3AlignedBox3TI,
    IntrLine3AlignedBox3FI
} from '../src/IntrLine3AlignedBox3.js';
import {
    IntrLine3OrientedBox3TI,
    IntrLine3OrientedBox3FI
} from '../src/IntrLine3OrientedBox3.js';

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

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, expectClose,
    expectVectorClose, seededRandom
} from './helpers/arbitraries.js';
import { length as vlength } from '../src/Vector.js';

const lineBox3 = fc.tuple(wellScaledVector(3, -6, 6), unitVector(3),
    wellScaledVector(3, -5, 5), rotationFrame(3),
    fc.double({ min: 0.2, max: 4, noNaN: true }),
    fc.double({ min: 0.2, max: 4, noNaN: true }),
    fc.double({ min: 0.2, max: 4, noNaN: true }))
    .map(([o, d, c, R, e0, e1, e2]) => ({
        line: Line.fromOriginDirection(o, d),
        box: OrientedBox.fromCenterAxisExtent(c, R,
            Vector.fromArray([e0, e1, e2]))
    }));

function boxCoords3(box: OrientedBox, p: Vector): number[] {
    const diff = sub(p, box.center);
    return [dot(diff, box.axis[0]), dot(diff, box.axis[1]),
        dot(diff, box.axis[2])];
}

describe('IntrLine3OrientedBox3 verification', () => {
    const tiq = new IntrLine3OrientedBox3TI();
    const fiq = new IntrLine3OrientedBox3FI();

    it('TI and FI agree on intersect', () => {
        check(lineBox3, ({ line: l, box: b }) => {
            expect(tiq.test(l, b).intersect).toBe(fiq.find(l, b).intersect);
        });
    });

    it('matches the aligned-box query on the box-frame line', () => {
        const atiq = new IntrLine3AlignedBox3TI();
        const afiq = new IntrLine3AlignedBox3FI();
        check(lineBox3, ({ line: l, box: b }) => {
            const diff = sub(l.origin, b.center);
            const local = Line.fromOriginDirection(
                Vector.fromArray(b.axis.map(a => dot(diff, a))),
                Vector.fromArray(b.axis.map(a => dot(l.direction, a))));
            const e = b.extent.values;
            const alignedBox3 = AlignedBox.fromMinMax(
                Vector.fromArray([-e[0], -e[1], -e[2]]),
                Vector.fromArray([e[0], e[1], e[2]]));
            expect(tiq.test(l, b).intersect)
                .toBe(atiq.test(local, alignedBox3).intersect);
            const af = afiq.find(local, alignedBox3);
            const f = fiq.find(l, b);
            expect(f.intersect).toBe(af.intersect);
            expect(f.numIntersections).toBe(af.numIntersections);
            if (f.intersect) {
                // Bit-identical: both go through the same DoQuery.
                expect(f.parameter[0]).toBe(af.parameter[0]);
                expect(f.parameter[1]).toBe(af.parameter[1]);
            }
        });
    });

    it('the reported points are on the line and inside the box', () => {
        check(lineBox3, ({ line: l, box: b }) => {
            const f = fiq.find(l, b);
            if (!f.intersect) {
                return;
            }
            expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
            // Upstream fills both entries whenever intersect is true.
            for (let i = 0; i < 2; ++i) {
                expectVectorClose(f.point[i],
                    add(l.origin, mul(f.parameter[i], l.direction)), 0, 0);
                const q = boxCoords3(b, f.point[i]);
                for (let k = 0; k < 3; ++k) {
                    expect(Math.abs(q[k]))
                        .toBeLessThanOrEqual(b.extent.values[k] + 1e-9);
                }
            }
        });
    });

    it('a line through a sampled box point always hits', () => {
        check(fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
            fc.double({ min: 0.2, max: 4, noNaN: true }),
            fc.double({ min: 0.2, max: 4, noNaN: true }),
            fc.double({ min: 0.2, max: 4, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }), unitVector(3)),
            ([c, R, e0, e1, e2, s0, s1, s2, d]) => {
                const b = OrientedBox.fromCenterAxisExtent(c, R,
                    Vector.fromArray([e0, e1, e2]));
                const target = add(c, add(mul(s0 * e0, R[0]),
                    add(mul(s1 * e1, R[1]), mul(s2 * e2, R[2]))));
                const l = Line.fromOriginDirection(target, d);
                expect(tiq.test(l, b).intersect).toBe(true);
                const f = fiq.find(l, b);
                expect(f.intersect).toBe(true);
                expect(f.parameter[0]).toBeLessThanOrEqual(1e-9);
                expect(f.parameter[1]).toBeGreaterThanOrEqual(-1e-9);
            });
    });

    it('a fine sweep of the line agrees with the reported interval', () => {
        const rnd = seededRandom(0x3c1f77a1);
        for (let trial = 0; trial < 120; ++trial) {
            const ang = [rnd() * 6.28, rnd() * 6.28, rnd() * 6.28];
            const R = [
                Vector.fromArray([Math.cos(ang[0]), Math.sin(ang[0]), 0]),
                Vector.fromArray([-Math.sin(ang[0]), Math.cos(ang[0]), 0]),
                Vector.fromArray([0, 0, 1])];
            const b = OrientedBox.fromCenterAxisExtent(
                Vector.fromArray([rnd() * 4 - 2, rnd() * 4 - 2,
                    rnd() * 4 - 2]),
                R, Vector.fromArray([0.5 + rnd(), 0.5 + rnd(), 0.5 + rnd()]));
            const dir = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (vlength(dir) < 0.3) {
                continue;
            }
            normalize(dir);
            const l = Line.fromOriginDirection(
                Vector.fromArray([rnd() * 8 - 4, rnd() * 8 - 4,
                    rnd() * 8 - 4]), dir);
            const f = fiq.find(l, b);
            for (let k = 0; k <= 500; ++k) {
                const t = -10 + (20 * k) / 500;
                const q = boxCoords3(b, add(l.origin, mul(t, dir)));
                const inside = q.every((qi, i) =>
                    Math.abs(qi) < b.extent.values[i] - 1e-6);
                if (inside) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-9);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-9);
                }
            }
        }
    }, 30000);

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineBox3, rotationFrame(3), wellScaledVector(3, -5, 5)),
            ([{ line: l, box: b }, R, tr]) => {
                const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                    add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
                const xf = (p: Vector): Vector => add(tr, rot(p));
                const l2 = Line.fromOriginDirection(xf(l.origin),
                    rot(l.direction));
                const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                    b.axis.map(a => rot(a)), b.extent);
                const f1 = fiq.find(l, b);
                const f2 = fiq.find(l2, b2);
                // A grazing hit legitimately flips under perturbation.
                const graze = (f: typeof f1): boolean => f.intersect
                    && f.parameter[1] - f.parameter[0] < 1e-6;
                if (graze(f1) || graze(f2)) {
                    return;
                }
                expect(f1.intersect).toBe(f2.intersect);
                if (!f1.intersect) {
                    return;
                }
                expectClose(f1.parameter[0], f2.parameter[0], 1e-8, 1e-9);
                expectClose(f1.parameter[1], f2.parameter[1], 1e-8, 1e-9);
            });
    });
});
