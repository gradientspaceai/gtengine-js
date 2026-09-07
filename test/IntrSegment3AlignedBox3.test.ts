import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { inContainerAlignedBox } from '../src/ContAlignedBox.js';
import {
    IntrSegment3AlignedBox3TI,
    IntrSegment3AlignedBox3FI,
    defaultIntrSegment3AlignedBox3FIResult,
    intrSegment3AlignedBox3FIDoQuery
} from '../src/IntrSegment3AlignedBox3.js';
import { Line } from '../src/Line.js';
import {
    IntrLine3AlignedBox3TI, IntrLine3AlignedBox3FI
} from '../src/IntrLine3AlignedBox3.js';
import {
    check, expectClose, expectVectorClose, fc, positive, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
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

describe('IntrSegment3AlignedBox3', () => {
    const ti = new IntrSegment3AlignedBox3TI();
    const fi = new IntrSegment3AlignedBox3FI();
    const unit = box([-1, -1, -1], [1, 1, 1]);

    it('clips a segment crossing the box', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The parameters are relative to the centered form, extent 3.
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports the whole segment when it lies inside the box', () => {
        const s = segment([-0.5, -0.25, 0.1], [0.5, 0.25, -0.1]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(-0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0.5, 12);
    });

    it('clips a diagonal segment through the box', () => {
        const s = segment([-3, -3, -3], [3, 3, 3]);
        const result = fi.find(s, unit);
        expect(result.numIntersections).toBe(2);
        for (let j = 0; j < 3; ++j) {
            expect(result.point[0].values[j]).toBeCloseTo(-1, 12);
            expect(result.point[1].values[j]).toBeCloseTo(1, 12);
        }
    });

    it('reports a single point when the segment touches a corner', () => {
        // The segment lies in the plane z = 1 and touches the box only at
        // the corner (1,1,1).
        const s = segment([0, 2, 1], [2, 0, 1]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[2]).toBeCloseTo(1, 12);
    });

    it('misses when the segment stops short of the box', () => {
        const s = segment([-5, 0, 0], [-2, 0, 0]);
        expect(ti.test(s, unit).intersect).toBe(false);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('misses when the supporting line misses the box', () => {
        const s = segment([-5, 3, 0], [5, 3, 0]);
        expect(ti.test(s, unit).intersect).toBe(false);
        expect(fi.find(s, unit).intersect).toBe(false);
    });

    it('handles a degenerate segment', () => {
        // A zero-length segment has zero direction after normalize; the
        // separating-axis test then compares |origin[i]| with boxExtent[i].
        const inside = segment([0.25, -0.5, 0], [0.25, -0.5, 0]);
        expect(ti.test(inside, unit).intersect).toBe(true);
        const outside = segment([5, 5, 5], [5, 5, 5]);
        expect(ti.test(outside, unit).intersect).toBe(false);
        expect(fi.find(outside, unit).intersect).toBe(false);
    });

    it('handles a translated box', () => {
        const shifted = box([9, 9, 9], [11, 11, 11]);
        const s = segment([7, 10, 10], [13, 10, 10]);
        expect(ti.test(s, shifted).intersect).toBe(true);
        const result = fi.find(s, shifted);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(9, 12);
        expect(result.point[1].values[0]).toBeCloseTo(11, 12);
        expect(result.point[0].values[1]).toBeCloseTo(10, 12);
    });

    it('the exported FI DoQuery matches the class query', () => {
        const s = segment([-4, 0.3, -0.2], [4, -0.1, 0.5]);
        const { center: boxCenter, extent: boxExtent } = unit.getCenteredForm();
        const shifted = Segment.fromEndpoints(sub(s.p[0], boxCenter),
            sub(s.p[1], boxCenter));
        const cf = shifted.getCenteredForm();
        const direct = defaultIntrSegment3AlignedBox3FIResult();
        intrSegment3AlignedBox3FIDoQuery(cf.center, cf.direction, cf.extent,
            boxExtent, direct);
        const viaClass = fi.find(s, unit);
        expect(direct.intersect).toBe(viaClass.intersect);
        expect(direct.numIntersections).toBe(viaClass.numIntersections);
        expect(direct.parameter[0]).toBeCloseTo(viaClass.parameter[0], 12);
        expect(direct.parameter[1]).toBeCloseTo(viaClass.parameter[1], 12);
    });

    it('agrees with brute-force sampling on random configurations', () => {
        const rnd = makeRandom(271828);
        let tiFiMismatch = 0;
        let sampleMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;
        const samples = 2000;

        for (let trial = 0; trial < 200; ++trial) {
            const c = vec(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            const e = vec(0.2 + rnd(), 0.2 + rnd(), 0.2 + rnd());
            const b = AlignedBox.fromMinMax(sub(c, e), add(c, e));
            const p0 = vec(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const p1 = add(p0, add(mul(2, sub(c, p0)),
                vec(3 * rnd() - 1.5, 3 * rnd() - 1.5, 3 * rnd() - 1.5)));
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
                const cf = s.getCenteredForm();
                const eps = vec(1e-9, 1e-9, 1e-9);
                const grown = AlignedBox.fromMinMax(sub(b.min, eps),
                    add(b.max, eps));
                for (let i = 0; i < fiResult.numIntersections; ++i) {
                    const t = fiResult.parameter[i];
                    if (Math.abs(t) > cf.extent + 1e-9) {
                        ++pointMismatch;
                    }
                    const expected = add(cf.center, mul(t, cf.direction));
                    const diff = sub(fiResult.point[i], expected);
                    if (Math.sqrt(dot(diff, diff)) > 1e-9) {
                        ++pointMismatch;
                    }
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
// IntrSegment3AlignedBox3.h.
// ---------------------------------------------------------------------------

describe('IntrSegment3AlignedBox3 verification', () => {
    const ti = new IntrSegment3AlignedBox3TI();
    const fi = new IntrSegment3AlignedBox3FI();
    const lineTi = new IntrLine3AlignedBox3TI();
    const lineFi = new IntrLine3AlignedBox3FI();

    const arbSegment = fc.tuple(wellScaledVector(3, -8, 8),
        wellScaledVector(3, -8, 8))
        .filter(([a, b]) => {
            const d = sub(b, a);
            return dot(d, d) > 1e-2;
        })
        .map(([a, b]) => Segment.fromEndpoints(a, b));
    const arbBox = fc.tuple(wellScaledVector(3, -6, 6),
        wellScaledVector(3, -6, 6))
        .map(([a, b]) => {
            const lo = new Vector(3), hi = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                lo.set(i, Math.min(a.get(i), b.get(i)));
                hi.set(i, Math.max(a.get(i), b.get(i)));
            }
            return AlignedBox.fromMinMax(lo, hi);
        });
    const arbPair = fc.tuple(arbSegment, arbBox);

    // A sound separation certificate: both endpoints lie strictly beyond the
    // same box face.
    function separatedByAxis(s: Segment, b: AlignedBox): boolean {
        for (let d = 0; d < 3; ++d) {
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

    it('both queries agree with a robust brute-force classification', () => {
        // The TI query uses separating axes while the FI query clips the
        // t-interval; the two can differ by an ulp for a segment that only
        // touches the box, so the property asserts agreement only when the
        // configuration is robustly one-sided.
        check(arbPair, ([s, b]) => {
            let insideMargin = false;
            for (let k = 0; k <= 256 && !insideMargin; ++k) {
                const t = k / 256;
                const p = add(mul(s.p[0], 1 - t), mul(s.p[1], t));
                let all = true;
                for (let d = 0; d < 3 && all; ++d) {
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

    it('FI reports centered-form parameters in [-e,e] and points in the box',
        () => {
            check(arbPair, ([s, b]) => {
                const res = fi.find(s, b);
                if (!res.intersect) {
                    expect(res.numIntersections).toBe(0);
                    return;
                }
                const { center, direction, extent } = s.getCenteredForm();
                for (let i = 0; i < res.numIntersections; ++i) {
                    const t = res.parameter[i];
                    expect(Number.isFinite(t)).toBe(true);
                    expect(Math.abs(t)).toBeLessThanOrEqual(
                        extent * (1 + 1e-12) + 1e-12);
                    // The reported points use the centered form C + t*D.
                    expectVectorClose(res.point[i],
                        add(center, mul(t, direction)), 1e-9, 1e-9);
                    for (let d = 0; d < 3; ++d) {
                        const tol = 1e-8 * (1 + Math.abs(res.point[i].get(d)));
                        expect(res.point[i].get(d))
                            .toBeGreaterThanOrEqual(b.min.get(d) - tol);
                        expect(res.point[i].get(d))
                            .toBeLessThanOrEqual(b.max.get(d) + tol);
                    }
                }
                expect(res.parameter[0]).toBeLessThanOrEqual(res.parameter[1]);
            });
        });

    // Boxes centered at the origin: the FI query translates the segment by the
    // box center before taking its centered form, so for an off-center box the
    // world-frame centered form used by this property would differ from the
    // query's by rounding.
    const arbCenteredBox = fc.array(positive(6, 1e-2),
        { minLength: 3, maxLength: 3 })
        .map(e => AlignedBox.fromMinMax(
            Vector.fromArray(e.map(x => -x)), Vector.fromArray(e)));

    it('the segment result is the line result clipped to [-e,e]', () => {
        check(fc.tuple(arbSegment, arbCenteredBox), ([s, b]) => {
            const { center, direction, extent } = s.getCenteredForm();
            const lineRes = lineFi.find(
                Line.fromOriginDirection(center, direction), b);
            const segRes = fi.find(s, b);
            if (!lineRes.intersect) {
                expect(segRes.intersect).toBe(false);
                return;
            }
            const c0 = Math.max(lineRes.parameter[0], -extent);
            const c1 = Math.min(lineRes.parameter[1], extent);
            if (c0 > c1) {
                expect(segRes.intersect).toBe(false);
                return;
            }
            expect(segRes.intersect).toBe(true);
            expect(segRes.parameter[0] + 0).toBe(c0 + 0);
            expect(segRes.parameter[1] + 0).toBe(c1 + 0);
            expect(segRes.numIntersections).toBe(c0 < c1 ? 2 : 1);
        });
    });

    it('the TI query never reports more than the line TI query', () => {
        // The separating-axis prefilter can only remove intersections the line
        // query would have reported. The comparison uses the same box-relative
        // centered form the query itself builds; the world-frame form differs
        // by rounding, which flips the answer for segments that merely touch a
        // degenerate box.
        check(fc.tuple(arbSegment, arbCenteredBox), ([s, b]) => {
            const { center, direction } = s.getCenteredForm();
            if (ti.test(s, b).intersect) {
                expect(lineTi.test(
                    Line.fromOriginDirection(center, direction), b).intersect)
                    .toBe(true);
            }
        });
    });

    // The query classifies with exact comparisons, so a segment that only
    // touches the box can change its answer under an arbitrarily small
    // perturbation. Equivariance is asserted only for configurations that
    // robustly cross the box interior.
    function crossesInterior(s: Segment, b: AlignedBox): boolean {
        for (let k = 0; k <= 256; ++k) {
            const t = k / 256;
            const p = add(mul(s.p[0], 1 - t), mul(s.p[1], t));
            let all = true;
            for (let d = 0; d < 3 && all; ++d) {
                const tol = 1e-4 * (1 + Math.abs(p.get(d)));
                all = p.get(d) > b.min.get(d) + tol
                    && p.get(d) < b.max.get(d) - tol;
            }
            if (all) {
                return true;
            }
        }
        return false;
    }

    it('is equivariant under translation', () => {
        check(fc.tuple(arbPair, wellScaledVector(3, -5, 5)),
            ([[s, b], shift]) => {
                if (!crossesInterior(s, b)) {
                    return;
                }
                const s2 = Segment.fromEndpoints(add(s.p[0], shift),
                    add(s.p[1], shift));
                const b2 = AlignedBox.fromMinMax(add(b.min, shift),
                    add(b.max, shift));
                const a = fi.find(s, b), c = fi.find(s2, b2);
                expect(c.intersect).toBe(a.intersect);
                if (a.intersect) {
                    expect(c.numIntersections).toBe(a.numIntersections);
                    for (let i = 0; i < a.numIntersections; ++i) {
                        expectClose(c.parameter[i], a.parameter[i],
                            1e-8, 1e-8);
                        expectVectorClose(c.point[i], add(a.point[i], shift),
                            1e-7, 1e-7);
                    }
                }
            });
    });

    it('reports the whole segment for a contained segment', () => {
        const s = segment([-1, -1, -1], [1, 1, 1]);
        const b = box([-4, -4, -4], [4, 4, 4]);
        const res = fi.find(s, b);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(2);
        const e = s.getCenteredForm().extent;
        expectClose(res.parameter[0], -e);
        expectClose(res.parameter[1], e);
        expectVectorClose(res.point[0], vec(-1, -1, -1));
        expectVectorClose(res.point[1], vec(1, 1, 1));
    });

    it('the FI DoQuery resets a result whose interval misses the segment',
        () => {
            const res = defaultIntrSegment3AlignedBox3FIResult();
            // Centered segment [-1,1] along +x at x = 5, box extent 1 at the
            // origin: the line hits the box for t in [-6,-4].
            intrSegment3AlignedBox3FIDoQuery(vec(5, 0, 0), vec(1, 0, 0), 1,
                vec(1, 1, 1), res);
            expect(res.intersect).toBe(false);
            expect(res.numIntersections).toBe(0);
            expect(res.parameter).toEqual([0, 0]);
            expect(res.point[0].equals(Vector.zero(3))).toBe(true);
        });
});
