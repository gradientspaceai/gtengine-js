import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrSegment3AlignedBox3TI,
    IntrSegment3AlignedBox3FI
} from '../src/IntrSegment3AlignedBox3.js';
import {
    IntrSegment3OrientedBox3TI,
    IntrSegment3OrientedBox3FI,
    defaultIntrSegment3OrientedBox3TIResult,
    defaultIntrSegment3OrientedBox3FIResult
} from '../src/IntrSegment3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

// An orthonormal frame: a rotation about z by 'a' followed by a rotation
// about the rotated x axis by 'b'.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        vec(ca, sa, 0),
        vec(-sa * cb, ca * cb, sb),
        vec(sa * sb, -ca * sb, cb)
    ];
}

function obb(center: number[], axes: Vector[], extent: number[]):
    OrientedBox {
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center), axes,
        Vector.fromArray(extent));
}

// The signed "outside" amount of a point relative to the box: negative when
// strictly inside, zero on the boundary, positive when outside.
function boxSignedDepth(p: Vector, box: OrientedBox): number {
    const delta = sub(p, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 3; ++i) {
        const value = Math.abs(dot(box.axis[i], delta)) -
            box.extent.values[i];
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

// The most deeply contained sampled point of the segment, expressed as the
// signed depth above; negative means some sample is strictly inside.
function bestSampledDepth(s: Segment, box: OrientedBox): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= 4000; ++i) {
        const p = add(s.p[0], mul(i / 4000, sub(s.p[1], s.p[0])));
        const value = boxSignedDepth(p, box);
        if (value < best) {
            best = value;
        }
    }
    return best;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrSegment3OrientedBox3', () => {
    const ti = new IntrSegment3OrientedBox3TI();
    const fi = new IntrSegment3OrientedBox3FI();
    const identityAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];
    const unit = obb([0, 0, 0], identityAxes, [1, 1, 1]);

    it('default results report no intersection', () => {
        const tiResult = defaultIntrSegment3OrientedBox3TIResult();
        expect(tiResult.intersect).toBe(false);
        const fiResult = defaultIntrSegment3OrientedBox3FIResult();
        expect(fiResult.intersect).toBe(false);
        expect(fiResult.numIntersections).toBe(0);
        expect(fiResult.parameter).toEqual([0, 0]);
    });

    it('clips a segment crossing an axis-aligned box', () => {
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

    it('clips a segment with only one endpoint inside the box', () => {
        // The segment straddles the face x = 1.
        const s = segment([0, 0, 0], [4, 0, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports a tangent (single-point) contact on a face', () => {
        // The segment lies in the plane x = 1 and grazes the box face.
        const s = segment([1, -3, 0], [1, 3, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        for (const p of result.point) {
            expect(p.values[0]).toBeCloseTo(1, 12);
        }
    });

    it('reports a single-point contact at a corner', () => {
        // The segment touches the corner (1,1,1) and immediately leaves.
        const s = segment([1, 1, 1], [3, 4, 5]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        for (let i = 0; i < 3; ++i) {
            expect(result.point[0].values[i]).toBeCloseTo(1, 10);
        }
    });

    it('reports no intersection for a segment that stops short of the box',
        () => {
            const s = segment([2, 0, 0], [5, 0, 0]);
            expect(ti.test(s, unit).intersect).toBe(false);
            expect(fi.find(s, unit).intersect).toBe(false);
        });

    it('reports no intersection when the line hits but the segment does not',
        () => {
            // The line through these points passes through the box, but the
            // overlap with the segment t-interval is empty.
            const s = segment([3, 0, 0], [6, 0, 0]);
            expect(ti.test(s, unit).intersect).toBe(false);
            expect(fi.find(s, unit).intersect).toBe(false);
        });

    it('is sensitive to the box orientation', () => {
        // A box rotated 45 degrees about z reaches out to x = sqrt(2), so a
        // point at x = 1.2 on the x axis is inside it but outside the
        // axis-aligned unit box.
        const c = Math.SQRT1_2;
        const box45 = obb([0, 0, 0], [vec(c, c, 0), vec(-c, c, 0),
            vec(0, 0, 1)], [1, 1, 1]);
        const s = segment([1.2, 0, 0], [3, 0, 0]);
        expect(ti.test(s, box45).intersect).toBe(true);
        expect(ti.test(s, unit).intersect).toBe(false);
        const result = fi.find(s, box45);
        expect(result.intersect).toBe(true);
        expect(result.point[0].values[0]).toBeCloseTo(1.2, 10);
        expect(result.point[1].values[0]).toBeCloseTo(Math.SQRT2, 10);
    });

    it('handles a degenerate zero-length segment', () => {
        const inside = segment([0.25, -0.5, 0.75], [0.25, -0.5, 0.75]);
        expect(ti.test(inside, unit).intersect).toBe(true);
        const insideResult = fi.find(inside, unit);
        expect(insideResult.intersect).toBe(true);
        for (let i = 0; i < 2; ++i) {
            expect(insideResult.point[i].values[0]).toBeCloseTo(0.25, 12);
            expect(insideResult.point[i].values[1]).toBeCloseTo(-0.5, 12);
            expect(insideResult.point[i].values[2]).toBeCloseTo(0.75, 12);
        }

        const outside = segment([3, 4, 5], [3, 4, 5]);
        expect(ti.test(outside, unit).intersect).toBe(false);
        expect(fi.find(outside, unit).intersect).toBe(false);
    });

    it('matches the aligned-box query when the axes are the identity', () => {
        const abTI = new IntrSegment3AlignedBox3TI();
        const abFI = new IntrSegment3AlignedBox3FI();
        const random = makeRandom(90134);
        for (let trial = 0; trial < 300; ++trial) {
            const min = [3 * random() - 3, 3 * random() - 3,
                3 * random() - 3];
            const max = [min[0] + 0.2 + 3 * random(),
                min[1] + 0.2 + 3 * random(), min[2] + 0.2 + 3 * random()];
            const aligned = AlignedBox.fromMinMax(Vector.fromArray(min),
                Vector.fromArray(max));
            const { center, extent } = aligned.getCenteredForm();
            const oriented = OrientedBox.fromCenterAxisExtent(center,
                identityAxes, extent);
            const s = segment(
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3],
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3]);

            expect(ti.test(s, oriented).intersect)
                .toBe(abTI.test(s, aligned).intersect);
            const expected = abFI.find(s, aligned);
            const actual = fi.find(s, oriented);
            expect(actual.intersect).toBe(expected.intersect);
            if (actual.intersect) {
                expect(actual.numIntersections)
                    .toBe(expected.numIntersections);
                for (let i = 0; i < 2; ++i) {
                    expect(actual.parameter[i])
                        .toBeCloseTo(expected.parameter[i], 9);
                    for (let j = 0; j < 3; ++j) {
                        expect(actual.point[i].values[j])
                            .toBeCloseTo(expected.point[i].values[j], 9);
                    }
                }
            }
        }
    });

    it('agrees with sampling, and TI agrees with FI, on random configurations',
        () => {
            const random = makeRandom(90135);
            let intersections = 0;
            for (let trial = 0; trial < 400; ++trial) {
                const axes = frame(2 * Math.PI * random(), Math.PI * random());
                const box = obb(
                    [3 * random() - 1.5, 3 * random() - 1.5,
                        3 * random() - 1.5],
                    axes,
                    [0.2 + 1.5 * random(), 0.2 + 1.5 * random(),
                        0.2 + 1.5 * random()]);
                const s = segment(
                    [6 * random() - 3, 6 * random() - 3, 6 * random() - 3],
                    [6 * random() - 3, 6 * random() - 3, 6 * random() - 3]);

                const tiResult = ti.test(s, box);
                const fiResult = fi.find(s, box);
                expect(tiResult.intersect).toBe(fiResult.intersect);

                const depth = bestSampledDepth(s, box);
                if (depth < -1e-6) {
                    // A sampled segment point is strictly inside the box.
                    expect(tiResult.intersect).toBe(true);
                }
                else if (depth > 1e-2) {
                    // Every segment point is well outside the box; the dense
                    // sampling cannot have missed a thin crossing.
                    expect(tiResult.intersect).toBe(false);
                }

                if (fiResult.intersect) {
                    ++intersections;
                    // The reported points are on the segment and on or in
                    // the box, and the parameters are ordered.
                    expect(fiResult.parameter[0])
                        .toBeLessThanOrEqual(fiResult.parameter[1] + 1e-12);
                    const { center, direction, extent } =
                        s.getCenteredForm();
                    for (let i = 0; i < 2; ++i) {
                        expect(Math.abs(fiResult.parameter[i]))
                            .toBeLessThanOrEqual(extent + 1e-9);
                        const onSegment = add(center,
                            mul(fiResult.parameter[i], direction));
                        for (let j = 0; j < 3; ++j) {
                            expect(fiResult.point[i].values[j])
                                .toBeCloseTo(onSegment.values[j], 9);
                        }
                        expect(boxSignedDepth(fiResult.point[i], box))
                            .toBeLessThan(1e-9);
                    }
                }
            }
            // Guard against a degenerate test that never intersects.
            expect(intersections).toBeGreaterThan(20);
        });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrSegment3OrientedBox3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { IntrLine3OrientedBox3FI } from '../src/IntrLine3OrientedBox3.js';
import { Line } from '../src/Line.js';
import { length, normalize } from '../src/Vector.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no frame component is a
// subnormal that would underflow when squared.
const angle3 = () => wellScaled(-Math.PI, Math.PI);
const extent3 = () => fc.double(
    { min: 0.25, max: 3, noNaN: true, noDefaultInfinity: true });

// R = Rz(a)*Ry(b)*Rx(c); the columns are the box axes.
function rotFrame3(a: number, b: number, c: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const cc = Math.cos(c), sc = Math.sin(c);
    return [
        vec(ca * cb, sa * cb, -sb),
        vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc, cb * sc),
        vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc, cb * cc)];
}

const segBox3 = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    angle3(), angle3(), angle3(),
    extent3(), extent3(), extent3()
).filter(([x0, y0, z0, x1, y1, z1]) =>
    (x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2 > 1e-4)
    .map(([x0, y0, z0, x1, y1, z1, cx, cy, cz, a, b, c, e0, e1, e2]) => ({
        segment: Segment.fromEndpoints(vec(x0, y0, z0), vec(x1, y1, z1)),
        box: OrientedBox.fromCenterAxisExtent(vec(cx, cy, cz),
            rotFrame3(a, b, c), vec(e0, e1, e2))
    }));

// boxSignedDepth is a maximum of convex functions, hence convex along the
// segment; a ternary search finds its minimum reliably. A value near zero
// means the segment grazes the box, where the TI separating-axis test and the
// FI interval clip round differently.
function minSegBoxDepth(box: OrientedBox, s: Segment): number {
    const e = sub(s.p[1], s.p[0]);
    let lo = 0, hi = 1;
    for (let k = 0; k < 200; ++k) {
        const p = lo + (hi - lo) / 3, q = hi - (hi - lo) / 3;
        const fp = boxSignedDepth(add(s.p[0], mul(p, e)), box);
        const fq = boxSignedDepth(add(s.p[0], mul(q, e)), box);
        if (fp < fq) { hi = q; } else { lo = p; }
    }
    return boxSignedDepth(add(s.p[0], mul(0.5 * (lo + hi), e)), box);
}

describe('IntrSegment3OrientedBox3 verification', () => {
    const tiq = new IntrSegment3OrientedBox3TI();
    const fiq = new IntrSegment3OrientedBox3FI();
    const lbq = new IntrLine3OrientedBox3FI();

    it('TI and FI agree away from grazing configurations', () => {
        check(segBox3, ({ segment: s, box: b }) => {
            if (Math.abs(minSegBoxDepth(b, s)) < 1e-6) {
                return;
            }
            expect(tiq.test(s, b).intersect).toBe(fiq.find(s, b).intersect);
        });
    });

    it('the segment hit is the line hit clipped to |t| <= extent', () => {
        check(segBox3, ({ segment: s, box: b }) => {
            const cf = s.getCenteredForm();
            const line = Line.fromOriginDirection(cf.center, cf.direction);
            const lb = lbq.find(line, b);
            const f = fiq.find(s, b);
            if (!lb.intersect) {
                expect(f.intersect).toBe(false);
                return;
            }
            const t0 = Math.max(lb.parameter[0], -cf.extent);
            const t1 = Math.min(lb.parameter[1], cf.extent);
            if (t0 > t1) {
                expect(f.intersect).toBe(false);
                return;
            }
            expect(f.intersect).toBe(true);
            expect(f.numIntersections).toBe(t0 < t1 ? 2 : 1);
            // Normalize the -0/+0 tie (toBe uses Object.is).
            expect(f.parameter[0] + 0).toBe(t0 + 0);
        });
    });

    it('reported points lie on the segment and in the closed box', () => {
        check(segBox3, ({ segment: s, box: b }) => {
            const cf = s.getCenteredForm();
            const f = fiq.find(s, b);
            if (!f.intersect) {
                expect(f.point[0].values).toEqual([0, 0, 0]);
                expect(f.point[1].values).toEqual([0, 0, 0]);
                return;
            }
            // Upstream fills both entries whenever intersect is true.
            for (let i = 0; i < 2; ++i) {
                // The parameter is in the centered convention, C + t*D with
                // |t| <= e.
                expect(Math.abs(f.parameter[i]))
                    .toBeLessThanOrEqual(cf.extent + 1e-12);
                // The point is the world-space evaluation of the centered
                // form; upstream reconstructs it through the box frame, which
                // must agree.
                expectVectorClose(f.point[i],
                    add(cf.center, mul(f.parameter[i], cf.direction)),
                    1e-9, 1e-9);
                expect(boxSignedDepth(f.point[i], b)).toBeLessThanOrEqual(1e-8);
            }
        });
    });

    it('a segment strictly inside the box is clipped to its own endpoints',
        () => {
            check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3),
                wellScaled(-3, 3), angle3(), angle3(), angle3(),
                extent3(), extent3(), extent3(),
                fc.array(fc.double({ min: -0.8, max: 0.8, noNaN: true,
                    noDefaultInfinity: true }), { minLength: 6, maxLength: 6 })),
            ([cx, cy, cz, a, b, c, e0, e1, e2, u]) => {
                const box = OrientedBox.fromCenterAxisExtent(vec(cx, cy, cz),
                    rotFrame3(a, b, c), vec(e0, e1, e2));
                const at = (i: number): Vector => add(box.center,
                    add(mul(u[i] * e0, box.axis[0]),
                        add(mul(u[i + 1] * e1, box.axis[1]),
                            mul(u[i + 2] * e2, box.axis[2]))));
                const p0 = at(0), p1 = at(3);
                if (length(sub(p1, p0)) < 1e-3) { return; }
                const s = Segment.fromEndpoints(p0, p1);
                const cf = s.getCenteredForm();
                expect(tiq.test(s, box).intersect).toBe(true);
                const f = fiq.find(s, box);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(2);
                expectClose(f.parameter[0], -cf.extent, 1e-12, 1e-12);
                expectClose(f.parameter[1], cf.extent, 1e-12, 1e-12);
                expectVectorClose(f.point[0], p0, 1e-8, 1e-9);
                expectVectorClose(f.point[1], p1, 1e-8, 1e-9);
            });
        });

    it('a fine sweep of the segment agrees with the reported interval', () => {
        const rnd = seededRandom(0x6b1de44);
        for (let trial = 0; trial < 150; ++trial) {
            const box = OrientedBox.fromCenterAxisExtent(
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rotFrame3(rnd() * 6, rnd() * 6, rnd() * 6),
                vec(0.3 + rnd() * 2, 0.3 + rnd() * 2, 0.3 + rnd() * 2));
            const s = Segment.fromEndpoints(
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5),
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5));
            if (length(sub(s.p[1], s.p[0])) < 1e-2) { continue; }
            const cf = s.getCenteredForm();
            const f = fiq.find(s, box);
            for (let k = 0; k <= 500; ++k) {
                const t = -cf.extent + (2 * cf.extent * k) / 500;
                const p = add(cf.center, mul(t, cf.direction));
                if (boxSignedDepth(p, box) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(segBox3, angle3(), angle3(), angle3(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2)),
        ([{ segment: s, box: b }, a1, a2, a3, tx, ty, tz]) => {
            if (Math.abs(minSegBoxDepth(b, s)) < 1e-6) {
                return;   // grazing face, edge or corner
            }
            const fr = rotFrame3(a1, a2, a3);
            const rot = (v: Vector): Vector => add(mul(v.get(0), fr[0]),
                add(mul(v.get(1), fr[1]), mul(v.get(2), fr[2])));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty, tz));
            const s2 = Segment.fromEndpoints(xf(s.p[0]), xf(s.p[1]));
            const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                [rot(b.axis[0]), rot(b.axis[1]), rot(b.axis[2])], b.extent);
            const f0 = fiq.find(s, b);
            const f1 = fiq.find(s2, b2);
            expect(f1.intersect).toBe(f0.intersect);
            if (!f0.intersect) { return; }
            expect(f1.numIntersections).toBe(f0.numIntersections);
            for (let i = 0; i < 2; ++i) {
                expectVectorClose(f1.point[i], xf(f0.point[i]), 1e-7, 1e-8);
            }
        });
    });

    it('reconstructs world-space points for every rotation of the box', () => {
        // The query evaluates the intersection in the box frame and rotates
        // it back through the axes; a shot through the box center must come
        // back on the box boundary in world space.
        check(fc.tuple(angle3(), angle3(), angle3(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2),
            extent3(), extent3(), extent3(), angle3(), angle3()),
        ([a, b, c, cx, cy, cz, e0, e1, e2, th, ph]) => {
            const box = OrientedBox.fromCenterAxisExtent(vec(cx, cy, cz),
                rotFrame3(a, b, c), vec(e0, e1, e2));
            const d = vec(Math.cos(th) * Math.cos(ph),
                Math.sin(th) * Math.cos(ph), Math.sin(ph));
            normalize(d);
            const s = Segment.fromEndpoints(sub(box.center, mul(30, d)),
                add(box.center, mul(30, d)));
            const f = fiq.find(s, box);
            expect(f.intersect).toBe(true);
            for (let i = 0; i < 2; ++i) {
                expectClose(boxSignedDepth(f.point[i], box), 0, 1e-8, 1e-9);
            }
        });
    });
});
