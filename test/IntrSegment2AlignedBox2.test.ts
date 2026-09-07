import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { inContainerAlignedBox } from '../src/ContAlignedBox.js';
import {
    IntrSegment2AlignedBox2TI,
    IntrSegment2AlignedBox2FI
} from '../src/IntrSegment2AlignedBox2.js';
import {
    check, expectClose, expectVectorClose, fc, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrSegment2AlignedBox2', () => {
    const ti = new IntrSegment2AlignedBox2TI();
    const fi = new IntrSegment2AlignedBox2FI();
    const unit = box([-1, -1], [1, 1]);

    it('clips a segment crossing the box', () => {
        const s = segment([-3, 0], [3, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The endpoint-form parameters for the crossing at x = -1 and x = 1.
        expect(result.parameter[0]).toBeCloseTo(1 / 3, 12);
        expect(result.parameter[1]).toBeCloseTo(2 / 3, 12);
        // The centered form of the segment has extent 3 and direction (1,0).
        expect(result.cdeParameter[0]).toBeCloseTo(-1, 12);
        expect(result.cdeParameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports the whole segment when it is inside the box', () => {
        const s = segment([-0.5, -0.25], [0.5, 0.25]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0.5, 12);
    });

    it('clips a segment with one endpoint inside the box', () => {
        const s = segment([0, 0], [4, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports 2 intersections when the segment touches only a corner', () => {
        // The segment touches the box only at the corner (1,1).
        const s = segment([0, 2], [2, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        // Upstream promotes a single-point overlap to 2 intersections so the
        // caller computes both (identical) points.
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('misses when the segment stops short of the box', () => {
        const s = segment([-5, 0], [-2, 0]);
        expect(ti.test(s, unit).intersect).toBe(false);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('misses when the supporting line misses the box', () => {
        const s = segment([-5, 3], [5, 3]);
        expect(ti.test(s, unit).intersect).toBe(false);
        expect(fi.find(s, unit).intersect).toBe(false);
    });

    it('handles a degenerate segment inside the box', () => {
        const s = segment([0.25, -0.5], [0.25, -0.5]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.cdeParameter).toEqual([0, 0]);
        expect(result.point[0].values[0]).toBeCloseTo(0.25, 12);
        expect(result.point[1].values[1]).toBeCloseTo(-0.5, 12);
    });

    it('handles a degenerate segment outside the box', () => {
        const s = segment([5, 5], [5, 5]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('handles a translated box', () => {
        const shifted = box([9, 9], [11, 11]);
        const s = segment([7, 10], [13, 10]);
        expect(ti.test(s, shifted).intersect).toBe(true);
        const result = fi.find(s, shifted);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(9, 12);
        expect(result.point[1].values[0]).toBeCloseTo(11, 12);
        expect(result.point[0].values[1]).toBeCloseTo(10, 12);
    });

    it('agrees with brute-force sampling on random configurations', () => {
        const rnd = makeRandom(31415);
        let tiFiMismatch = 0;
        let sampleMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;
        const samples = 3000;

        for (let trial = 0; trial < 250; ++trial) {
            const cx = 2 * rnd() - 1, cy = 2 * rnd() - 1;
            const ex = 0.2 + rnd(), ey = 0.2 + rnd();
            const b = box([cx - ex, cy - ey], [cx + ex, cy + ey]);
            const p0 = vec(6 * rnd() - 3, 6 * rnd() - 3);
            const p1 = add(p0, add(mul(2, sub(vec(cx, cy), p0)),
                vec(3 * rnd() - 1.5, 3 * rnd() - 1.5)));
            const d = sub(p1, p0);
            if (dot(d, d) < 1e-6) {
                continue;
            }
            const s = Segment.fromEndpoints(p0, p1);

            const tiResult = ti.test(s, b);
            const fiResult = fi.find(s, b);
            if (tiResult.intersect !== fiResult.intersect) {
                ++tiFiMismatch;
            }

            let sampled = false;
            for (let k = 0; k <= samples; ++k) {
                if (inContainerAlignedBox(add(p0, mul(k / samples, d)), b)) {
                    sampled = true;
                    break;
                }
            }
            if (sampled && !fiResult.intersect) {
                ++sampleMismatch;
            }

            if (fiResult.intersect) {
                ++hits;
                for (let i = 0; i < fiResult.numIntersections; ++i) {
                    const t = fiResult.parameter[i];
                    if (t < -1e-9 || t > 1 + 1e-9) {
                        ++pointMismatch;
                    }
                    // The reported point must match the endpoint-form
                    // evaluation of the reported parameter.
                    const expected = add(p0, mul(t, d));
                    const diff = sub(fiResult.point[i], expected);
                    if (Math.sqrt(dot(diff, diff)) > 1e-9) {
                        ++pointMismatch;
                    }
                    // ... and it must be in the box (up to rounding).
                    const grown = AlignedBox.fromMinMax(
                        sub(b.min, vec(1e-9, 1e-9)),
                        add(b.max, vec(1e-9, 1e-9)));
                    if (!inContainerAlignedBox(expected, grown)) {
                        ++pointMismatch;
                    }
                }
            }
        }

        expect(hits).toBeGreaterThan(20);
        expect([tiFiMismatch, sampleMismatch, pointMismatch]).toEqual([0, 0, 0]);
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrSegment2AlignedBox2.h.
// ---------------------------------------------------------------------------

describe('IntrSegment2AlignedBox2 verification', () => {
    const ti = new IntrSegment2AlignedBox2TI();
    const fi = new IntrSegment2AlignedBox2FI();

    // Moderately scaled segments and boxes: the shared segment()/alignedBox()
    // generators draw from finite(), which emits subnormal coordinates that
    // make the containment comparisons below meaningless.
    const arbSegment = fc.tuple(wellScaledVector(2, -8, 8),
        wellScaledVector(2, -8, 8))
        .filter(([a, b]) => {
            const d = sub(b, a);
            return dot(d, d) > 1e-2;
        })
        .map(([a, b]) => Segment.fromEndpoints(a, b));
    const arbBox = fc.tuple(wellScaledVector(2, -6, 6),
        wellScaledVector(2, -6, 6))
        .map(([a, b]) => {
            const lo = new Vector(2), hi = new Vector(2);
            for (let i = 0; i < 2; ++i) {
                lo.set(i, Math.min(a.get(i), b.get(i)));
                hi.set(i, Math.max(a.get(i), b.get(i)));
            }
            return AlignedBox.fromMinMax(lo, hi);
        });
    const arbPair = fc.tuple(arbSegment, arbBox);

    // A sound separation certificate: both endpoints are strictly beyond the
    // same box face, so no point of the segment is in the box.
    function separatedByAxis(s: Segment, b: AlignedBox): boolean {
        for (let d = 0; d < 2; ++d) {
            const tol = 1e-6 * (1 + Math.abs(b.max.get(d))
                + Math.abs(b.min.get(d)));
            if (s.p[0].get(d) > b.max.get(d) + tol
                && s.p[1].get(d) > b.max.get(d) + tol) {
                return true;
            }
            if (s.p[0].get(d) < b.min.get(d) - tol
                && s.p[1].get(d) < b.min.get(d) - tol) {
                return true;
            }
        }
        return false;
    }

    it('both queries report no intersection when an axis separates the'
        + ' segment from the box', () => {
        check(arbPair, ([s, b]) => {
            if (separatedByAxis(s, b)) {
                expect(ti.test(s, b).intersect).toBe(false);
                expect(fi.find(s, b).intersect).toBe(false);
            }
        });
    });

    it('TI and FI agree on intersect away from tangential configurations',
        () => {
            // The TI query uses the method of separating axes while the FI
            // query uses Liang-Barsky clipping followed by an interval
            // intersection. For a segment that only touches the box the two
            // computations differ by an ulp, so they can disagree; see the
            // deterministic case below. The property therefore requires
            // agreement only when the configuration is robustly one-sided.
            check(arbPair, ([s, b]) => {
                let insideMargin = false;
                for (let k = 0; k <= 256 && !insideMargin; ++k) {
                    const t = k / 256;
                    const p = add(mul(s.p[0], 1 - t), mul(s.p[1], t));
                    let all = true;
                    for (let d = 0; d < 2 && all; ++d) {
                        const tol = 1e-6 * (1 + Math.abs(p.get(d)));
                        all = p.get(d) > b.min.get(d) + tol
                            && p.get(d) < b.max.get(d) - tol;
                    }
                    insideMargin = all;
                }
                if (insideMargin) {
                    expect(ti.test(s, b).intersect).toBe(true);
                    expect(fi.find(s, b).intersect).toBe(true);
                }
                else if (separatedByAxis(s, b)) {
                    expect(ti.test(s, b).intersect).toBe(false);
                    expect(fi.find(s, b).intersect).toBe(false);
                }
            });
        });

    it('TI and FI can disagree for a segment that only touches the box'
        + ' (structural, matches upstream)', () => {
        // Degenerate box: the vertical segment from (0,-6) to (0,0). The query
        // segment ends exactly on the box corner (0,0). The separating-axis
        // test of the TI query finds |segOrigin[1]| == boxExtent[1] +
        // segExtent, so it reports contact, while the clipped t-interval of
        // the FI query starts one ulp beyond +segExtent, so it reports none.
        const s = segment([0, 1.9974098846622947], [0, 0]);
        const b = box([0, -5.9999999999990825], [0, 0]);
        expect(ti.test(s, b).intersect).toBe(true);
        expect(fi.find(s, b).intersect).toBe(false);
    });

    it('FI points are on the segment, inside the box and consistent with the'
        + ' two parameter conventions', () => {
        check(arbPair, ([s, b]) => {
            const res = fi.find(s, b);
            if (!res.intersect) {
                expect(res.numIntersections).toBe(0);
                return;
            }
            const { center, direction, extent } = s.getCenteredForm();
            for (let i = 0; i < res.numIntersections; ++i) {
                const t = res.parameter[i];
                const sParam = res.cdeParameter[i];
                expect(Number.isFinite(t)).toBe(true);
                expect(Number.isFinite(sParam)).toBe(true);
                // The endpoint-form parameter is in [0,1] and the centered
                // form parameter is in [-e,e]; the two are related by
                // t = (s/e + 1)/2.
                expect(t).toBeGreaterThanOrEqual(-1e-12);
                expect(t).toBeLessThanOrEqual(1 + 1e-12);
                expect(Math.abs(sParam)).toBeLessThanOrEqual(
                    extent * (1 + 1e-12) + 1e-12);
                expectClose(t, (sParam / extent + 1) * 0.5, 1e-12, 1e-12);
                // point[i] = C + s*D and, to rounding, (1-t)*P0 + t*P1.
                expectVectorClose(res.point[i],
                    add(center, mul(sParam, direction)), 1e-9, 1e-9);
                expectVectorClose(res.point[i],
                    add(mul(s.p[0], 1 - t), mul(s.p[1], t)), 1e-8, 1e-8);
                // The point is in the (slightly grown) box.
                for (let d = 0; d < 2; ++d) {
                    const tol = 1e-8 * (1 + Math.abs(res.point[i].get(d)));
                    expect(res.point[i].get(d))
                        .toBeGreaterThanOrEqual(b.min.get(d) - tol);
                    expect(res.point[i].get(d))
                        .toBeLessThanOrEqual(b.max.get(d) + tol);
                }
            }
        });
    });

    it('intersect is true whenever a finely sampled segment point is strictly'
        + ' inside the box', () => {
        check(arbPair, ([s, b]) => {
            let inside = false;
            for (let k = 0; k <= 512 && !inside; ++k) {
                const t = k / 512;
                const p = add(mul(s.p[0], 1 - t), mul(s.p[1], t));
                let all = true;
                for (let d = 0; d < 2 && all; ++d) {
                    const tol = 1e-9 * (1 + Math.abs(p.get(d)));
                    all = p.get(d) > b.min.get(d) + tol
                        && p.get(d) < b.max.get(d) - tol;
                }
                inside = all;
            }
            if (inside) {
                expect(ti.test(s, b).intersect).toBe(true);
                expect(fi.find(s, b).intersect).toBe(true);
            }
        });
    });

    it('a degenerate segment intersects exactly when its point is in the box',
        () => {
            check(fc.tuple(wellScaledVector(2, -6, 6), arbBox), ([p, b]) => {
                const s = Segment.fromEndpoints(p, p.clone());
                const res = fi.find(s, b);
                const contained = inContainerAlignedBox(p, b);
                expect(res.intersect).toBe(contained);
                if (contained) {
                    expect(res.numIntersections).toBe(1);
                    expect(res.parameter).toEqual([0, 0]);
                    expect(res.cdeParameter).toEqual([0, 0]);
                    expect(res.point[0].equals(p)).toBe(true);
                    expect(res.point[1].equals(p)).toBe(true);
                }
            });
        });

    it('a segment touching the box only at a corner reports two coincident'
        + ' points, as upstream forces numIntersections to 2', () => {
        // The upstream DoQuery replaces numIntersections == 1 by 2 so that the
        // caller always fills both point[] slots.
        const s = segment([-1, 1], [1, -1]);
        const b = box([0, 0], [2, 2]);
        const res = fi.find(s, b);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(2);
        expect(res.point[0].equals(vec(0, 0))).toBe(true);
        expect(res.point[1].equals(vec(0, 0))).toBe(true);
        expect(res.cdeParameter[0]).toBe(res.cdeParameter[1]);
    });

    it('a contained segment reports its own endpoints', () => {
        const s = segment([-1, -1], [1, 1]);
        const b = box([-4, -4], [4, 4]);
        const res = fi.find(s, b);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(2);
        expectClose(res.parameter[0], 0);
        expectClose(res.parameter[1], 1);
        expectVectorClose(res.point[0], vec(-1, -1));
        expectVectorClose(res.point[1], vec(1, 1));
    });

    it('the result does not alias the input segment or box', () => {
        const p0 = vec(-1, -1), p1 = vec(1, 1);
        const s = Segment.fromEndpoints(p0, p1);
        const b = box([-4, -4], [4, 4]);
        const res = fi.find(s, b);
        res.point[0].set(0, 99);
        expect(s.p[0].get(0)).toBe(-1);
        expect(p0.get(0)).toBe(-1);
        expect(res.point[1].get(0)).toBe(1);
    });
});
