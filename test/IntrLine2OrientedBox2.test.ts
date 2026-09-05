import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Line } from '../src/Line.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, mul, sub, dot, normalize } from '../src/Vector.js';
import {
    IntrLine2AlignedBox2TI,
    IntrLine2AlignedBox2FI
} from '../src/IntrLine2AlignedBox2.js';
import {
    IntrLine2OrientedBox2TI,
    IntrLine2OrientedBox2FI
} from '../src/IntrLine2OrientedBox2.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(px: number, py: number, dx: number, dy: number): Line {
    return Line.fromOriginDirection(v2(px, py), v2(dx, dy));
}

// An oriented box with the given center, rotation angle and extents.
function obox(cx: number, cy: number, angle: number, e0: number,
    e1: number): OrientedBox {
    const c = Math.cos(angle), s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(v2(cx, cy),
        [v2(c, s), v2(-s, c)], v2(e0, e1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Brute-force containment in the solid oriented box.
function inBox(p: Vector, box: OrientedBox): boolean {
    const diff = sub(p, box.center);
    for (let i = 0; i < 2; ++i) {
        if (Math.abs(dot(diff, box.axis[i])) > box.extent.values[i] + 1e-12) {
            return false;
        }
    }
    return true;
}

describe('IntrLine2OrientedBox2', () => {
    const ti = new IntrLine2OrientedBox2TI();
    const fi = new IntrLine2OrientedBox2FI();
    const abTI = new IntrLine2AlignedBox2TI();
    const abFI = new IntrLine2AlignedBox2FI();

    it('finds the entry and exit parameters for an axis-aligned box', () => {
        // A zero-angle oriented box is the aligned box [-1,1]^2.
        const result = fi.find(line(0, 0, 1, 0), obox(0, 0, 0, 1, 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a rotated box', () => {
        // A unit square rotated 45 degrees has corners at distance sqrt(2)
        // along the coordinate axes.
        const box = obox(0, 0, Math.PI / 4, 1, 1);
        const result = fi.find(line(0, 0, 1, 0), box);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-Math.SQRT2, 12);
        expect(result.parameter[1]).toBeCloseTo(Math.SQRT2, 12);

        // A horizontal line above the top corner misses.
        expect(ti.test(line(0, 1.5, 1, 0), box).intersect).toBe(false);
        expect(ti.test(line(0, 1.4, 1, 0), box).intersect).toBe(true);
    });

    it('reports a corner contact for a tangent line', () => {
        // A rotated box touched by a horizontal line at its topmost corner.
        // The topmost corner of the box with axes (0.6,0.8) and (-0.8,0.6)
        // and unit extents is (-0.2, 1.4).
        const box = OrientedBox.fromCenterAxisExtent(v2(0, 0),
            [v2(0.6, 0.8), v2(-0.8, 0.6)], v2(1, 1));
        const result = fi.find(line(0, 1.4, 1, 0), box);
        expect(result.intersect).toBe(true);
        for (let i = 0; i < result.numIntersections; ++i) {
            expect(result.point[i].values[0]).toBeCloseTo(-0.2, 9);
            expect(result.point[i].values[1]).toBeCloseTo(1.4, 9);
        }

        // Just above the corner there is no intersection.
        expect(ti.test(line(0, 1.4 + 1e-6, 1, 0), box).intersect).toBe(false);
    });

    it('matches the aligned-box query when the axes are the standard basis',
        () => {
            const rnd = makeRandom(5150);
            const ab = AlignedBox.fromMinMax(v2(-1, -2), v2(3, 1));
            const { center, extent } = ab.getCenteredForm();
            const ob = OrientedBox.fromCenterAxisExtent(center,
                [v2(1, 0), v2(0, 1)], extent);
            for (let k = 0; k < 300; ++k) {
                const l = line(rnd() * 8 - 4, rnd() * 8 - 4,
                    rnd() * 2 - 1, rnd() * 2 - 1);
                if (dot(l.direction, l.direction) < 1e-8) {
                    continue;
                }
                expect(ti.test(l, ob).intersect)
                    .toBe(abTI.test(l, ab).intersect);
                const r0 = fi.find(l, ob);
                const r1 = abFI.find(l, ab);
                expect(r0.numIntersections).toBe(r1.numIntersections);
                for (let i = 0; i < r0.numIntersections; ++i) {
                    expect(r0.parameter[i]).toBeCloseTo(r1.parameter[i], 9);
                }
            }
        });

    it('agrees with the TI query and a brute-force sampling', () => {
        const rnd = makeRandom(606);
        let numHit = 0, numMiss = 0;
        for (let k = 0; k < 300; ++k) {
            const box = obox(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * Math.PI,
                0.2 + rnd(), 0.2 + rnd());
            const d = v2(rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d, d) < 1e-8) {
                continue;
            }
            normalize(d);
            const l = Line.fromOriginDirection(
                v2(rnd() * 6 - 3, rnd() * 6 - 3), d);

            const tiResult = ti.test(l, box).intersect;
            const fiResult = fi.find(l, box);
            expect(tiResult).toBe(fiResult.intersect);

            // Every reported intersection point must be in the box.
            for (let i = 0; i < fiResult.numIntersections; ++i) {
                expect(inBox(fiResult.point[i], box)).toBe(true);
            }

            // Dense sampling of the line: a sample inside the box forces the
            // query to report an intersection.
            let sampleHit = false;
            for (let s = -600; s <= 600 && !sampleHit; ++s) {
                if (inBox(add(l.origin, mul(s * 0.01, l.direction)), box)) {
                    sampleHit = true;
                }
            }
            if (sampleHit) {
                expect(tiResult).toBe(true);
                ++numHit;
            } else if (!tiResult) {
                ++numMiss;
            }
        }
        expect(numHit).toBeGreaterThan(0);
        expect(numMiss).toBeGreaterThan(0);
    });

    it('handles a degenerate box with zero extents', () => {
        // The box degenerates to the point (1,1).
        const box = obox(1, 1, 0.3, 0, 0);
        expect(ti.test(line(0, 0, 1, 1), box).intersect).toBe(true);
        expect(ti.test(line(0, 0.5, 1, 1), box).intersect).toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, expectClose,
    expectVectorClose
} from './helpers/arbitraries.js';
import { length } from '../src/Vector.js';

// An oriented box built from a rotation frame, so no axis component is
// subnormal, together with a line whose direction is unit length.
const lineBox2 = fc.tuple(wellScaledVector(2, -5, 5), unitVector(2),
    wellScaledVector(2, -5, 5), rotationFrame(2),
    fc.double({ min: 0.1, max: 4, noNaN: true }),
    fc.double({ min: 0.1, max: 4, noNaN: true }))
    .map(([o, d, c, R, e0, e1]) => ({
        line: Line.fromOriginDirection(o, d),
        box: OrientedBox.fromCenterAxisExtent(c, R, Vector.fromArray([e0, e1]))
    }));

// The box coordinates of a point.
function boxCoords(box: OrientedBox, p: Vector): number[] {
    const diff = sub(p, box.center);
    return [dot(diff, box.axis[0]), dot(diff, box.axis[1])];
}

describe('IntrLine2OrientedBox2 verification', () => {
    const tiq = new IntrLine2OrientedBox2TI();
    const fiq = new IntrLine2OrientedBox2FI();

    it('TI and FI agree on intersect', () => {
        check(lineBox2, ({ line: l, box: b }) => {
            expect(tiq.test(l, b).intersect).toBe(fiq.find(l, b).intersect);
        });
    });

    it('matches the aligned-box query on the box-frame line', () => {
        const atiq = new IntrLine2AlignedBox2TI();
        const afiq = new IntrLine2AlignedBox2FI();
        check(lineBox2, ({ line: l, box: b }) => {
            const diff = sub(l.origin, b.center);
            const local = Line.fromOriginDirection(
                Vector.fromArray([dot(diff, b.axis[0]), dot(diff, b.axis[1])]),
                Vector.fromArray([dot(l.direction, b.axis[0]),
                    dot(l.direction, b.axis[1])]));
            const alignedBox = AlignedBox.fromMinMax(
                Vector.fromArray([-b.extent.values[0], -b.extent.values[1]]),
                Vector.fromArray([b.extent.values[0], b.extent.values[1]]));
            const at = atiq.test(local, alignedBox);
            const af = afiq.find(local, alignedBox);
            const t = tiq.test(l, b);
            const f = fiq.find(l, b);
            expect(t.intersect).toBe(at.intersect);
            expect(f.intersect).toBe(af.intersect);
            expect(f.numIntersections).toBe(af.numIntersections);
            for (let i = 0; i < f.numIntersections; ++i) {
                // Bit-identical: both go through the same DoQuery.
                expect(f.parameter[i]).toBe(af.parameter[i]);
            }
        });
    });

    it('the reported points are on the line and inside the box', () => {
        check(lineBox2, ({ line: l, box: b }) => {
            const f = fiq.find(l, b);
            if (!f.intersect) {
                return;
            }
            expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
            for (let i = 0; i < f.numIntersections; ++i) {
                const p = add(l.origin, mul(f.parameter[i], l.direction));
                expectVectorClose(f.point[i], p, 0, 0);
                const q = boxCoords(b, f.point[i]);
                expect(Math.abs(q[0]))
                    .toBeLessThanOrEqual(b.extent.values[0] + 1e-9);
                expect(Math.abs(q[1]))
                    .toBeLessThanOrEqual(b.extent.values[1] + 1e-9);
            }
        });
    });

    it('a line through a sampled box point always hits', () => {
        check(fc.tuple(wellScaledVector(2, -5, 5), rotationFrame(2),
            fc.double({ min: 0.1, max: 4, noNaN: true }),
            fc.double({ min: 0.1, max: 4, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }), unitVector(2)),
            ([c, R, e0, e1, s0, s1, d]) => {
                const b = OrientedBox.fromCenterAxisExtent(c, R,
                    Vector.fromArray([e0, e1]));
                const target = add(c, add(mul(s0 * e0, R[0]),
                    mul(s1 * e1, R[1])));
                const l = Line.fromOriginDirection(target, d);
                const t = tiq.test(l, b);
                const f = fiq.find(l, b);
                expect(t.intersect).toBe(true);
                expect(f.intersect).toBe(true);
                expect(f.parameter[0]).toBeLessThanOrEqual(0 + 1e-9);
                expect(f.parameter[1]).toBeGreaterThanOrEqual(0 - 1e-9);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineBox2, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([{ line: l, box: b }, R, tr]) => {
            const rot = (p: Vector): Vector =>
                add(mul(p.values[0], R[0]), mul(p.values[1], R[1]));
            const xf = (p: Vector): Vector => add(tr, rot(p));
            const l2 = Line.fromOriginDirection(xf(l.origin), rot(l.direction));
            const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                [rot(b.axis[0]), rot(b.axis[1])], b.extent);
            const f1 = fiq.find(l, b);
            const f2 = fiq.find(l2, b2);
            // A grazing hit legitimately flips under perturbation; require a
            // robustly transverse interval.
            if (f1.intersect && f1.parameter[1] - f1.parameter[0] < 1e-6) {
                return;
            }
            if (f2.intersect && f2.parameter[1] - f2.parameter[0] < 1e-6) {
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

    it('a zero-extent box behaves as a point on the line', () => {
        check(fc.tuple(wellScaledVector(2, -5, 5), unitVector(2),
            rotationFrame(2), fc.double({ min: -3, max: 3, noNaN: true })),
            ([o, d, R, t]) => {
                const center = add(o, mul(t, d));
                const b = OrientedBox.fromCenterAxisExtent(center, R,
                    Vector.fromArray([0, 0]));
                const l = Line.fromOriginDirection(o, d);
                const f = fiq.find(l, b);
                for (let i = 0; i < f.numIntersections; ++i) {
                    expect(length(sub(f.point[i], center)))
                        .toBeLessThan(1e-6);
                }
            });
    });
});
