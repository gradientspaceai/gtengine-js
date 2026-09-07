import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import { Vector, mul } from '../src/Vector.js';
import { expectClose, seededRandom } from './helpers/arbitraries.js';
import {
    AreaEllipse2Ellipse2,
    AreaEllipse2Ellipse2Configuration as Cfg,
    defaultAreaEllipse2Ellipse2Result
} from '../src/IntrAreaEllipse2Ellipse2.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function ellipse(cx: number, cy: number, angle: number, a: number,
    b: number): Hyperellipsoid {
    const c = Math.cos(angle), s = Math.sin(angle);
    return Hyperellipsoid.fromCenterAxisExtent(v2(cx, cy),
        [v2(c, s), v2(-s, c)], v2(a, b));
}

function circle(cx: number, cy: number, r: number): Hyperellipsoid {
    return ellipse(cx, cy, 0, r, r);
}

// The symmetric matrix M of the algebraic form (X-C)^T*M*(X-C) = 1, returned
// as [M00, M01, M11].
function standardForm(e: Hyperellipsoid): [number, number, number] {
    const [u, v] = e.axis;
    const a2 = e.extent.values[0] * e.extent.values[0];
    const b2 = e.extent.values[1] * e.extent.values[1];
    const m00 = u.values[0] * u.values[0] / a2 + v.values[0] * v.values[0] / b2;
    const m01 = u.values[0] * u.values[1] / a2 + v.values[0] * v.values[1] / b2;
    const m11 = u.values[1] * u.values[1] / a2 + v.values[1] * v.values[1] / b2;
    return [m00, m01, m11];
}

// The y-interval of the ellipse for the given x, or null when the vertical
// line misses the ellipse. The ellipse points satisfy
//   m00*dx^2 + 2*m01*dx*dy + m11*dy^2 = 1,
// a quadratic in dy = y - cy.
function yInterval(e: Hyperellipsoid, x: number): [number, number] | null {
    const [m00, m01, m11] = standardForm(e);
    const dx = x - e.center.values[0];
    const disc = m01 * m01 * dx * dx - m11 * (m00 * dx * dx - 1);
    if (disc <= 0) {
        return null;
    }
    const root = Math.sqrt(disc);
    const cy = e.center.values[1];
    return [cy + (-m01 * dx - root) / m11, cy + (-m01 * dx + root) / m11];
}

// The x-extent of the ellipse.
function xInterval(e: Hyperellipsoid): [number, number] {
    const [m00, m01, m11] = standardForm(e);
    const det = m00 * m11 - m01 * m01;
    const halfWidth = Math.sqrt(m11 / det);
    return [e.center.values[0] - halfWidth, e.center.values[0] + halfWidth];
}

// A high-accuracy numeric evaluation of the area of intersection: the length
// of the vertical chord common to both ellipses, integrated in x by the
// trapezoid rule.
function numericArea(e0: Hyperellipsoid, e1: Hyperellipsoid,
    numSamples: number = 200000): number {
    const [lo0, hi0] = xInterval(e0);
    const [lo1, hi1] = xInterval(e1);
    const lo = Math.max(lo0, lo1);
    const hi = Math.min(hi0, hi1);
    if (hi <= lo) {
        return 0;
    }
    const h = (hi - lo) / numSamples;
    let sum = 0;
    for (let i = 0; i <= numSamples; ++i) {
        const x = lo + i * h;
        const s0 = yInterval(e0, x);
        const s1 = yInterval(e1, x);
        let chord = 0;
        if (s0 !== null && s1 !== null) {
            chord = Math.max(0, Math.min(s0[1], s1[1]) - Math.max(s0[0], s1[0]));
        }
        sum += (i === 0 || i === numSamples ? 0.5 : 1) * chord;
    }
    return sum * h;
}

// The exact area of the lens common to two circles of radii r0 and r1 whose
// centers are a distance d apart.
function lensArea(r0: number, r1: number, d: number): number {
    if (d >= r0 + r1) {
        return 0;
    }
    if (d <= Math.abs(r0 - r1)) {
        return Math.PI * Math.min(r0, r1) * Math.min(r0, r1);
    }
    const a0 = Math.acos((d * d + r0 * r0 - r1 * r1) / (2 * d * r0));
    const a1 = Math.acos((d * d + r1 * r1 - r0 * r0) / (2 * d * r1));
    const tri = 0.5 * Math.sqrt(
        (-d + r0 + r1) * (d + r0 - r1) * (d - r0 + r1) * (d + r0 + r1));
    return r0 * r0 * a0 + r1 * r1 * a1 - tri;
}

describe('IntrAreaEllipse2Ellipse2', () => {
    it('default result is invalid with zero area', () => {
        const result = defaultAreaEllipse2Ellipse2Result();
        expect(result.configuration).toBe(Cfg.INVALID);
        expect(result.area).toBe(0);
        expect(result.findResult.intersect).toBe(false);
    });

    it('reports equal ellipses with the full area', () => {
        const query = new AreaEllipse2Ellipse2();
        const e = ellipse(1, -2, 0.4, 3, 2);
        const result = query.compute(e, e);
        expect(result.configuration).toBe(Cfg.ELLIPSES_ARE_EQUAL);
        expect(result.area).toBeCloseTo(Math.PI * 6, 12);
    });

    it('reports separated ellipses with zero area', () => {
        const query = new AreaEllipse2Ellipse2();
        const result = query.compute(circle(0, 0, 1), circle(5, 0, 1));
        expect(result.configuration).toBe(Cfg.ELLIPSES_ARE_SEPARATED);
        expect(result.area).toBe(0);
    });

    it('reports containment with the area of the contained ellipse', () => {
        const query = new AreaEllipse2Ellipse2();

        const r0 = query.compute(circle(0, 0, 5), ellipse(0.5, 0, 0.3, 2, 1));
        expect(r0.configuration).toBe(Cfg.E0_CONTAINS_E1);
        expect(r0.area).toBeCloseTo(Math.PI * 2, 10);

        const r1 = query.compute(ellipse(0.5, 0, 0.3, 2, 1), circle(0, 0, 5));
        expect(r1.configuration).toBe(Cfg.E1_CONTAINS_E0);
        expect(r1.area).toBeCloseTo(Math.PI * 2, 10);
    });

    it('matches the closed-form lens area for two circles', () => {
        const query = new AreaEllipse2Ellipse2();
        for (const d of [0.25, 0.5, 1, 1.5, 1.9, 2.5]) {
            const result = query.compute(circle(0, 0, 1), circle(d, 0, 1.5));
            expect(result.area).toBeCloseTo(lensArea(1, 1.5, d), 9);
        }
    });

    it('is invariant under a rigid motion of both ellipses', () => {
        const query = new AreaEllipse2Ellipse2();
        const base = query.compute(ellipse(0, 0, 0, 3, 1),
            ellipse(1, 0.5, 1.1, 2, 1.5));
        // Rotate both ellipses by t about the origin and translate.
        const t = 0.7, c = Math.cos(t), s = Math.sin(t);
        const rot = (e: Hyperellipsoid): Hyperellipsoid => {
            const rc = (x: number, y: number): Vector =>
                v2(c * x - s * y, s * x + c * y);
            const center = rc(e.center.values[0], e.center.values[1]);
            center.values[0] += 4;
            center.values[1] -= 3;
            return Hyperellipsoid.fromCenterAxisExtent(center, [
                rc(e.axis[0].values[0], e.axis[0].values[1]),
                rc(e.axis[1].values[0], e.axis[1].values[1])
            ], e.extent.clone());
        };
        const moved = query.compute(rot(ellipse(0, 0, 0, 3, 1)),
            rot(ellipse(1, 0.5, 1.1, 2, 1.5)));
        expect(moved.configuration).toBe(base.configuration);
        expect(moved.area).toBeCloseTo(base.area, 8);
    });

    it('matches numeric integration for transverse configurations', () => {
        const query = new AreaEllipse2Ellipse2();
        const cases: [Hyperellipsoid, Hyperellipsoid][] = [
            // Two-point (one chord region) configurations.
            [circle(0, 0, 1), circle(1.2, 0, 1)],
            [ellipse(0, 0, 0, 3, 1), ellipse(2, 0.5, 0.4, 2, 1.5)],
            [ellipse(0, 0, 0.2, 2, 1), ellipse(0, 2.2, 1.3, 2.5, 0.8)],
            // Four-point (four chord region) configurations.
            [ellipse(0, 0, 0, 3, 1), ellipse(0, 0, 1.4, 2.5, 1.1)],
            [ellipse(0, 0, 0, 2, 1), ellipse(0.1, 0.05, 1.2, 2, 1)],
            [ellipse(0, 0, 0, 4, 1), ellipse(0.3, 0, 1.4, 3.5, 1.2)]
        ];
        for (const [e0, e1] of cases) {
            const result = query.compute(e0, e1);
            const expected = numericArea(e0, e1);
            expect(result.area).toBeGreaterThan(0);
            expect(Math.abs(result.area - expected)).toBeLessThan(
                1e-4 * Math.max(1, expected));
            // The configuration must agree with the number of intersections.
            if (result.findResult.numPoints === 4) {
                expect(result.configuration).toBe(Cfg.FOUR_CHORD_REGION);
            }
            else {
                expect(result.configuration).toBe(Cfg.ONE_CHORD_REGION);
            }
            // The area of the intersection never exceeds the smaller ellipse.
            const area0 = Math.PI * e0.extent.values[0] * e0.extent.values[1];
            const area1 = Math.PI * e1.extent.values[0] * e1.extent.values[1];
            expect(result.area).toBeLessThanOrEqual(
                Math.min(area0, area1) * (1 + 1e-12));
        }
    });

    it('matches numeric integration for randomized ellipse pairs', () => {
        // A deterministic linear congruential generator keeps the test
        // reproducible.
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const query = new AreaEllipse2Ellipse2();
        let numTransverse = 0;
        for (let trial = 0; trial < 40; ++trial) {
            const e0 = ellipse(2 * rand() - 1, 2 * rand() - 1,
                Math.PI * rand(), 0.5 + 2 * rand(), 0.5 + 2 * rand());
            const e1 = ellipse(2 * rand() - 1, 2 * rand() - 1,
                Math.PI * rand(), 0.5 + 2 * rand(), 0.5 + 2 * rand());
            const result = query.compute(e0, e1);
            const expected = numericArea(e0, e1, 20000);
            expect(Math.abs(result.area - expected)).toBeLessThan(
                2e-3 * Math.max(1, expected));
            if (result.configuration === Cfg.ONE_CHORD_REGION ||
                result.configuration === Cfg.FOUR_CHORD_REGION) {
                ++numTransverse;
            }
        }
        // The random sampling must exercise the chord-region code paths.
        expect(numTransverse).toBeGreaterThan(10);
    });

    it('is symmetric in its arguments', () => {
        const query = new AreaEllipse2Ellipse2();
        const pairs: [Hyperellipsoid, Hyperellipsoid][] = [
            [circle(0, 0, 1), circle(1.2, 0.3, 1.5)],
            [ellipse(0, 0, 0, 3, 1), ellipse(0.2, 0, 1.3, 2.5, 1.1)],
            [circle(0, 0, 4), circle(0.5, 0.5, 1)]
        ];
        for (const [e0, e1] of pairs) {
            const a = query.compute(e0, e1).area;
            const b = query.compute(e1, e0).area;
            expect(b).toBeCloseTo(a, 8);
        }
    });

    it('rejects ellipses of the wrong dimension', () => {
        const query = new AreaEllipse2Ellipse2();
        const e3 = new Hyperellipsoid(3);
        expect(() => query.compute(e3, e3)).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrAreaEllipse2Ellipse2.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('AreaEllipse2Ellipse2 verification', () => {
    const q = new AreaEllipse2Ellipse2();

    it('uses initialized constants (upstream leaves five members unset)', () => {
        // Upstream declares 'T mZero, mOne, mTwo, mPi, mTwoPi;' with no
        // constructor and never assigns them, so every area, every pi scaling
        // and the 'dtheta <= mPi' branch read indeterminate values
        // (gtengine-js #301 item 1). The port uses literals and the shared
        // pi constants; these exact values pin them.
        const unit = circle(0, 0, 1);
        const same = q.compute(unit, circle(0, 0, 1));
        expect(same.configuration).toBe(Cfg.ELLIPSES_ARE_EQUAL);
        expect(same.area).toBe(Math.PI);

        const inner = circle(0, 0, 0.5);
        const contained = q.compute(unit, inner);
        expect(contained.configuration).toBe(Cfg.E0_CONTAINS_E1);
        expect(contained.area).toBe(Math.PI * 0.25);

        const away = q.compute(unit, circle(10, 0, 1));
        expect(away.configuration).toBe(Cfg.ELLIPSES_ARE_SEPARATED);
        expect(away.area).toBe(0);

        // A half-overlap of two unit circles whose centers are 1 apart has the
        // closed-form lens area; the constant 2*pi enters through the angle
        // wrapping of ComputeAreaChordRegion.
        const lens = q.compute(unit, circle(1, 0, 1));
        expect(lens.configuration).toBe(Cfg.ONE_CHORD_REGION);
        expectClose(lens.area, lensArea(1, 1, 1), 1e-12, 1e-12);
    });

    it('matches numeric integration over many random pairs', () => {
        const rnd = seededRandom(0xa2ea11);
        const counts = new Map<number, number>();
        let compared = 0;
        for (let iter = 0; iter < 220; ++iter) {
            const e0 = ellipse(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * Math.PI,
                0.5 + rnd() * 1.5, 0.5 + rnd() * 1.5);
            const e1 = ellipse(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * Math.PI,
                0.5 + rnd() * 1.5, 0.5 + rnd() * 1.5);
            let result;
            try {
                result = q.compute(e0, e1);
            } catch {
                // IntrEllipse2Ellipse2 throws "Unexpected condition." for
                // nearly concentric ellipses (gtengine-js #458 item 6);
                // preserved from upstream.
                continue;
            }
            counts.set(result.configuration,
                (counts.get(result.configuration) ?? 0) + 1);
            const reference = numericArea(e0, e1, 20000);
            // The chord integrand has a square-root singularity at the ends
            // of the x-range, so the trapezoid rule converges like h^1.5; at
            // 20000 samples over a range of at most 6 that is a few 1e-5.
            expectClose(result.area, reference, 3e-4, 3e-4);
            ++compared;
        }
        expect(compared).toBeGreaterThan(180);
        // Every branch of AreaDispatch is exercised.
        expect(counts.get(Cfg.ONE_CHORD_REGION) ?? 0).toBeGreaterThan(20);
        expect(counts.get(Cfg.FOUR_CHORD_REGION) ?? 0).toBeGreaterThan(5);
        expect((counts.get(Cfg.E0_CONTAINS_E1) ?? 0)
            + (counts.get(Cfg.E1_CONTAINS_E0) ?? 0)
            + (counts.get(Cfg.ELLIPSES_ARE_SEPARATED) ?? 0))
            .toBeGreaterThan(12);
    }, 60000);

    it('the reported configuration matches the geometry', () => {
        const rnd = seededRandom(0xc0f16);
        for (let iter = 0; iter < 250; ++iter) {
            const e0 = ellipse(rnd() * 3 - 1.5, rnd() * 3 - 1.5,
                rnd() * Math.PI, 0.4 + rnd() * 1.6, 0.4 + rnd() * 1.6);
            const e1 = ellipse(rnd() * 3 - 1.5, rnd() * 3 - 1.5,
                rnd() * Math.PI, 0.4 + rnd() * 1.6, 0.4 + rnd() * 1.6);
            let r;
            try {
                r = q.compute(e0, e1);
            } catch {
                continue;
            }
            const area0 = Math.PI * e0.extent.values[0] * e0.extent.values[1];
            const area1 = Math.PI * e1.extent.values[0] * e1.extent.values[1];
            switch (r.configuration) {
                case Cfg.ELLIPSES_ARE_SEPARATED:
                    expect(r.area).toBe(0);
                    expectClose(numericArea(e0, e1, 4000), 0, 1e-6, 1e-6);
                    break;
                case Cfg.E0_CONTAINS_E1:
                    expectClose(r.area, area1, 1e-12, 1e-12);
                    expect(area1).toBeLessThanOrEqual(area0 * (1 + 1e-12));
                    break;
                case Cfg.E1_CONTAINS_E0:
                    expectClose(r.area, area0, 1e-12, 1e-12);
                    expect(area0).toBeLessThanOrEqual(area1 * (1 + 1e-12));
                    break;
                case Cfg.ELLIPSES_ARE_EQUAL:
                    expectClose(r.area, area0, 1e-12, 1e-12);
                    break;
                default:
                    // A chord region: the area is strictly between 0 and the
                    // area of the smaller ellipse.
                    expect(r.area).toBeGreaterThan(-1e-12);
                    expect(r.area).toBeLessThanOrEqual(
                        Math.min(area0, area1) * (1 + 1e-9));
                    break;
            }
            expect(Number.isNaN(r.area)).toBe(false);
        }
    }, 30000);

    it('is symmetric under swapping the ellipses', () => {
        const rnd = seededRandom(0x5a11);
        for (let iter = 0; iter < 200; ++iter) {
            const e0 = ellipse(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * Math.PI,
                0.5 + rnd(), 0.5 + rnd());
            const e1 = ellipse(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * Math.PI,
                0.5 + rnd(), 0.5 + rnd());
            let a;
            let b;
            try {
                a = q.compute(e0, e1);
                b = q.compute(e1, e0);
            } catch {
                continue;
            }
            expectClose(a.area, b.area, 1e-9, 1e-9);
        }
    }, 30000);

    it('is invariant when the ellipse axes are not unit length', () => {
        // Upstream's operator() comment says "The ellipse axes are not
        // required to be normalized", but the polar angles used by Area4 and
        // ComputeAreaChordRegion are not length invariant. The port
        // normalizes private copies of the axes, which makes the documented
        // contract true (gtengine-js #301 item 2).
        //
        // The invariance holds for a COMMON scale of the two axes. The
        // FIQuery this class consumes builds its standard form as
        // (U*U^T/a^2 + V*V^T/b^2) / |U|^2, which divides both terms by the
        // squared length of the FIRST axis, so independently scaled axes
        // are not supported by the query pipeline as a whole.
        const rnd = seededRandom(0x5ca1e);
        for (let iter = 0; iter < 120; ++iter) {
            const e0 = ellipse(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * Math.PI,
                0.5 + rnd(), 0.5 + rnd());
            const e1 = ellipse(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * Math.PI,
                0.5 + rnd(), 0.5 + rnd());
            let base;
            try {
                base = q.compute(e0, e1);
            } catch {
                continue;
            }
            const scale = 0.25 + rnd() * 4;
            const scaled = Hyperellipsoid.fromCenterAxisExtent(e0.center,
                [mul(scale, e0.axis[0]), mul(scale, e0.axis[1])],
                e0.extent);
            const scaledResult = q.compute(scaled, e1);
            expect(scaledResult.configuration).toBe(base.configuration);
            expectClose(scaledResult.area, base.area, 1e-9, 1e-9);
        }
    }, 30000);

    it('is equivariant under a rigid motion of both ellipses', () => {
        const rnd = seededRandom(0x8181d);
        for (let iter = 0; iter < 150; ++iter) {
            const e0 = ellipse(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * Math.PI,
                0.5 + rnd(), 0.5 + rnd());
            const e1 = ellipse(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * Math.PI,
                0.5 + rnd(), 0.5 + rnd());
            let base;
            try {
                base = q.compute(e0, e1);
            } catch {
                continue;
            }
            const a = rnd() * 2 * Math.PI;
            const ca = Math.cos(a);
            const sa = Math.sin(a);
            const tx = rnd() * 4 - 2;
            const ty = rnd() * 4 - 2;
            const map = (p: Vector) => v2(
                ca * p.values[0] - sa * p.values[1] + tx,
                sa * p.values[0] + ca * p.values[1] + ty);
            const rot = (p: Vector) => v2(
                ca * p.values[0] - sa * p.values[1],
                sa * p.values[0] + ca * p.values[1]);
            const move = (e: Hyperellipsoid) =>
                Hyperellipsoid.fromCenterAxisExtent(map(e.center),
                    [rot(e.axis[0]), rot(e.axis[1])], e.extent);
            let moved;
            try {
                moved = q.compute(move(e0), move(e1));
            } catch {
                continue;
            }
            // The chord-region areas are sums of elliptic sector integrals
            // whose arguments drift under the rotation.
            expectClose(moved.area, base.area, 1e-7, 1e-7);
        }
    }, 30000);

    it('reports containment for concentric ellipses of different sizes', () => {
        const outer = ellipse(0.25, -0.5, 0.4, 2, 1);
        const inner = ellipse(0.25, -0.5, 0.4, 1, 0.5);
        const r = q.compute(outer, inner);
        expect(r.configuration).toBe(Cfg.E0_CONTAINS_E1);
        expect(r.area).toBe(Math.PI * 0.5);
        const swapped = q.compute(inner, outer);
        expect(swapped.configuration).toBe(Cfg.E1_CONTAINS_E0);
        expect(swapped.area).toBe(Math.PI * 0.5);
    });

    it('consumes the four-point FI result in counterclockwise order', () => {
        // Two congruent ellipses rotated 90 degrees about a common center
        // cross in four points; the FOUR_CHORD_REGION branch orders them by
        // the polar angle in the frame of E0 (a std::multimap upstream, a
        // stable sort here).
        const e0 = ellipse(0, 0, 0, 2, 1);
        const e1 = ellipse(0, 0, Math.PI / 2, 2, 1);
        const r = q.compute(e0, e1);
        expect(r.configuration).toBe(Cfg.FOUR_CHORD_REGION);
        expect(r.findResult.numPoints).toBe(4);
        expectClose(r.area, numericArea(e0, e1, 400000), 1e-7, 1e-7);
        // The four crossings are at (+-t, +-t) with t = 2/sqrt(5).
        const t = 2 / Math.sqrt(5);
        for (let i = 0; i < 4; ++i) {
            // The quartic solver of IntrEllipse2Ellipse2 places the four
            // crossings to about 1e-8 for this symmetric configuration.
            expectClose(Math.abs(r.findResult.points[i].values[0]), t, 1e-7,
                1e-7);
            expectClose(Math.abs(r.findResult.points[i].values[1]), t, 1e-7,
                1e-7);
        }
    });
});
