import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import { IntrLine2AlignedBox2FI } from '../src/IntrLine2AlignedBox2.js';
import {
    IntrRay2AlignedBox2TI,
    IntrRay2AlignedBox2FI
} from '../src/IntrRay2AlignedBox2.js';
import {
    intrRay2AlignedBox2TIDoQuery,
    intrRay2AlignedBox2FIDoQuery
} from '../src/IntrRay2AlignedBox2.js';
import {
    defaultIntrLine2AlignedBox2TIResult,
    defaultIntrLine2AlignedBox2FIResult
} from '../src/IntrLine2AlignedBox2.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(vec(min), vec(max));
}

function ray(p: number[], d: number[]): Ray {
    const dir = vec(d);
    normalize(dir);
    return Ray.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function insideBox(b: AlignedBox, x: Vector): boolean {
    for (let i = 0; i < 2; ++i) {
        if (x.values[i] < b.min.values[i] - 1e-12 ||
            x.values[i] > b.max.values[i] + 1e-12) {
            return false;
        }
    }
    return true;
}

describe('IntrRay2AlignedBox2', () => {
    const ti = new IntrRay2AlignedBox2TI();
    const fi = new IntrRay2AlignedBox2FI();
    const unit = box([-1, -1], [1, 1]);

    it('finds both crossings for a ray aimed at the box', () => {
        const r = ray([-5, 0], [1, 0]);
        expect(ti.test(r, unit).intersect).toBe(true);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips the near end when the ray origin is inside the box', () => {
        const r = ray([0.5, 0], [1, 0]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
    });

    it('reports no intersection when the ray points away', () => {
        const r = ray([-5, 0], [-1, 0]);
        expect(ti.test(r, unit).intersect).toBe(false);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports a single touching point on the exit edge', () => {
        const r = ray([1, 0], [1, 0]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a ray that grazes a corner', () => {
        // The ray runs along x = 1, touching the right edge of the box.
        const r = ray([1, -5], [0, 1]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
    });

    it('is the line query clipped to t >= 0', () => {
        const rand = makeRandom(20260901);
        const lineFI = new IntrLine2AlignedBox2FI();
        const b = box([-1, 0.5], [2, 3]);
        for (let trial = 0; trial < 400; ++trial) {
            const r = ray([8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1]);
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lineResult = lineFI.find(l, b);
            const rayResult = fi.find(r, b);
            expect(ti.test(r, b).intersect).toBe(rayResult.intersect);

            if (!lineResult.intersect || lineResult.parameter[1] < 0) {
                expect(rayResult.intersect).toBe(false);
            }
            else {
                expect(rayResult.intersect).toBe(true);
                expect(rayResult.parameter[0]).toBeCloseTo(
                    Math.max(lineResult.parameter[0], 0), 12);
                expect(rayResult.parameter[1]).toBeCloseTo(
                    lineResult.parameter[1], 12);
            }
        }
    });

    it('agrees with dense sampling along the ray', () => {
        const rand = makeRandom(31337);
        const b = box([-1, 0.5], [2, 3]);
        for (let trial = 0; trial < 80; ++trial) {
            const r = ray([6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(r, b);

            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = (12 * k) / n;
                const x = add(r.origin, mul(t, r.direction));
                if (insideBox(b, x)) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(2e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(2e-3);
            }
        }
    });
});

describe('intrRay2AlignedBox2 DoQuery helpers', () => {
    // The helpers take the ray in the box-centered coordinate system. The box
    // below is already centered at the origin, so no translation is needed.
    const b = box([-2, -1], [2, 1]);
    const extent = vec([2, 1]);

    it('the TI helper agrees with the class query', () => {
        const cases: Array<[number, number, number, number]> = [
            [-5, 0, 1, 0],
            [5, 0, 1, 0],
            [0, 0, 0, 1],
            [-5, 3, 1, 0],
            [-4, -3, 1, 1]
        ];
        for (const [px, py, dx, dy] of cases) {
            const d = vec([dx, dy]);
            normalize(d);
            const result = defaultIntrLine2AlignedBox2TIResult();
            intrRay2AlignedBox2TIDoQuery(vec([px, py]), d, extent, result);
            const expected = new IntrRay2AlignedBox2TI().test(
                Ray.fromOriginDirection(vec([px, py]), d), b);
            expect(result.intersect).toBe(expected.intersect);
        }
    });

    it('the FI helper clips to the ray t-interval and leaves points alone', () => {
        const result = defaultIntrLine2AlignedBox2FIResult();
        intrRay2AlignedBox2FIDoQuery(vec([0, 0]), vec([1, 0]), extent, result);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The ray starts inside the box, so t0 is clamped to 0.
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('the FI helper rejects a ray pointing away from the box', () => {
        const result = defaultIntrLine2AlignedBox2FIResult();
        intrRay2AlignedBox2FIDoQuery(vec([-5, 0]), vec([-1, 0]), extent,
            result);
        expect(result.intersect).toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, unitVector, wellScaledVector, expectClose,
    expectVectorClose, seededRandom
} from './helpers/arbitraries.js';
import { sub } from '../src/Vector.js';

// A well-scaled aligned box: the shared 'alignedBox' generator draws its
// corners from fc.double, which emits subnormal coordinates (1e-323) whose
// centered form has no significant digits left, so the boundary tests are
// meaningless there.
const box2Arb = fc.tuple(wellScaledVector(2, -4, 4),
    fc.double({ min: 0.1, max: 4, noNaN: true }),
    fc.double({ min: 0.1, max: 4, noNaN: true }))
    .map(([lo, w, h]) => AlignedBox.fromMinMax(lo,
        Vector.fromArray([lo.values[0] + w, lo.values[1] + h])));

const rayBox2 = fc.tuple(wellScaledVector(2, -6, 6), unitVector(2), box2Arb)
    .map(([o, d, b]) => ({ ray: Ray.fromOriginDirection(o, d), box: b }));

function insideBox2(box: AlignedBox, p: Vector, tol: number): boolean {
    for (let i = 0; i < 2; ++i) {
        if (p.values[i] < box.min.values[i] - tol
            || p.values[i] > box.max.values[i] + tol) {
            return false;
        }
    }
    return true;
}

describe('IntrRay2AlignedBox2 verification', () => {
    const tiq = new IntrRay2AlignedBox2TI();
    const fiq = new IntrRay2AlignedBox2FI();

    it('TI and FI agree on intersect', () => {
        check(rayBox2, ({ ray: r, box: b }) => {
            expect(tiq.test(r, b).intersect).toBe(fiq.find(r, b).intersect);
        });
    });

    it('a ray hit is the line hit clipped to nonnegative parameters', () => {
        const lfi = new IntrLine2AlignedBox2FI();
        check(rayBox2, ({ ray: r, box: b }) => {
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lf = lfi.find(l, b);
            const f = fiq.find(r, b);
            if (!lf.intersect) {
                expect(f.intersect).toBe(false);
                return;
            }
            const lo = Math.max(lf.parameter[0], 0);
            const hi = lf.parameter[1];
            expect(f.intersect).toBe(hi >= 0);
            if (!f.intersect) {
                return;
            }
            expectClose(f.parameter[0], lo, 0, 0);
            expectClose(f.parameter[1], hi, 0, 0);
            expect(f.numIntersections).toBe(lo < hi ? 2 : 1);
        });
    });

    it('the reported points are on the ray and inside the box', () => {
        check(rayBox2, ({ ray: r, box: b }) => {
            const f = fiq.find(r, b);
            if (!f.intersect) {
                return;
            }
            expect(f.parameter[0]).toBeGreaterThanOrEqual(0);
            expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
            for (let i = 0; i < f.numIntersections; ++i) {
                expectVectorClose(f.point[i],
                    add(r.origin, mul(f.parameter[i], r.direction)), 0, 0);
                expect(insideBox2(b, f.point[i], 1e-9)).toBe(true);
            }
        });
    });

    it('a fine sweep of the ray agrees with the reported interval', () => {
        const rnd = seededRandom(0x1a2b3c4d);
        for (let trial = 0; trial < 200; ++trial) {
            const lo = Vector.fromArray([rnd() * 4 - 4, rnd() * 4 - 4]);
            const hi = Vector.fromArray([lo.values[0] + 0.5 + rnd() * 3,
                lo.values[1] + 0.5 + rnd() * 3]);
            const b = AlignedBox.fromMinMax(lo, hi);
            const a = rnd() * 2 * Math.PI;
            const d = Vector.fromArray([Math.cos(a), Math.sin(a)]);
            const r = Ray.fromOriginDirection(
                Vector.fromArray([rnd() * 8 - 4, rnd() * 8 - 4]), d);
            const f = fiq.find(r, b);
            let anyInside = false;
            for (let k = 0; k <= 600; ++k) {
                const t = (12 * k) / 600;
                const p = add(r.origin, mul(t, d));
                if (insideBox2(b, p, -1e-6)) {
                    anyInside = true;
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-9);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-9);
                }
            }
            if (anyInside) {
                expect(f.numIntersections).toBe(2);
            }
        }
    }, 30000);

    it('a ray whose origin is inside the box hits at parameter 0', () => {
        check(fc.tuple(box2Arb, unitVector(2),
            fc.double({ min: 0.05, max: 0.95, noNaN: true }),
            fc.double({ min: 0.05, max: 0.95, noNaN: true })),
            ([b, d, s0, s1]) => {
                const o = Vector.fromArray([
                    b.min.values[0]
                        + s0 * (b.max.values[0] - b.min.values[0]),
                    b.min.values[1]
                        + s1 * (b.max.values[1] - b.min.values[1])]);
                const r = Ray.fromOriginDirection(o, d);
                const f = fiq.find(r, b);
                expect(tiq.test(r, b).intersect).toBe(true);
                expect(f.intersect).toBe(true);
                expectClose(f.parameter[0], 0, 1e-12, 0);
                expect(f.parameter[1]).toBeGreaterThanOrEqual(0);
            });
    });

    it('the exported DoQuery helpers reproduce the class results', () => {
        check(rayBox2, ({ ray: r, box: b }) => {
            const cf = b.getCenteredForm();
            const localOrigin = sub(r.origin, cf.center);
            const tres = defaultIntrLine2AlignedBox2TIResult();
            intrRay2AlignedBox2TIDoQuery(localOrigin, r.direction, cf.extent,
                tres);
            expect(tres.intersect).toBe(tiq.test(r, b).intersect);

            const fres = defaultIntrLine2AlignedBox2FIResult();
            intrRay2AlignedBox2FIDoQuery(localOrigin, r.direction, cf.extent,
                fres);
            const f = fiq.find(r, b);
            expect(fres.intersect).toBe(f.intersect);
            expect(fres.numIntersections).toBe(f.numIntersections);
            for (let i = 0; i < fres.numIntersections; ++i) {
                expect(fres.parameter[i]).toBe(f.parameter[i]);
            }
            // The helper leaves point[] untouched; the class fills it.
            if (fres.numIntersections === 0) {
                expect(fres.point[0].values).toEqual([0, 0]);
            }
        });
    });
});
