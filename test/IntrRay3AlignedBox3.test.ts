import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import {
    IntrLine3AlignedBox3FI
} from '../src/IntrLine3AlignedBox3.js';
import { Line } from '../src/Line.js';
import {
    IntrRay3AlignedBox3TI,
    IntrRay3AlignedBox3FI
} from '../src/IntrRay3AlignedBox3.js';
import {
    intrRay3AlignedBox3TIDoQuery,
    intrRay3AlignedBox3FIDoQuery
} from '../src/IntrRay3AlignedBox3.js';
import {
    defaultIntrLine3AlignedBox3TIResult,
    defaultIntrLine3AlignedBox3FIResult
} from '../src/IntrLine3AlignedBox3.js';

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
    for (let i = 0; i < 3; ++i) {
        if (x.values[i] < b.min.values[i] - 1e-12 ||
            x.values[i] > b.max.values[i] + 1e-12) {
            return false;
        }
    }
    return true;
}

describe('IntrRay3AlignedBox3', () => {
    const ti = new IntrRay3AlignedBox3TI();
    const fi = new IntrRay3AlignedBox3FI();
    const unit = box([-1, -1, -1], [1, 1, 1]);

    it('clips the near end when the ray origin is inside the box', () => {
        const r = ray([0, 0, 0], [0, 0, 1]);
        expect(ti.test(r, unit).intersect).toBe(true);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values[2]).toBeCloseTo(1, 12);
    });

    it('finds both crossings for a ray that starts outside', () => {
        const r = ray([0, 0, -5], [0, 0, 1]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[2]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports no intersection when the ray points away from the box', () => {
        const r = ray([0, 0, -5], [0, 0, -1]);
        expect(ti.test(r, unit).intersect).toBe(false);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('reports the touching point when the ray origin is on the far face', () => {
        // The origin is on the +z face and the ray points outward, so the
        // overlap with [0,+infinity) is the single parameter 0.
        const r = ray([0, 0, 1], [0, 0, 1]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(0, 12);
        expect(ti.test(r, unit).intersect).toBe(true);
    });

    it('handles a ray grazing an edge of the box', () => {
        const r = ray([1, 1, -5], [0, 0, 1]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
    });

    it('is the line query clipped to t >= 0', () => {
        const rand = makeRandom(4242);
        const lineFI = new IntrLine3AlignedBox3FI();
        const b = box([-1, -2, 0.5], [2, 1, 3]);
        for (let trial = 0; trial < 300; ++trial) {
            const r = ray(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
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
        const rand = makeRandom(99);
        const b = box([-1, -2, 0.5], [2, 1, 3]);
        for (let trial = 0; trial < 60; ++trial) {
            const r = ray(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
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

describe('intrRay3AlignedBox3 DoQuery helpers', () => {
    // The helpers take the ray in the box-centered coordinate system. The box
    // below is centered at the origin, so no translation is needed.
    const b = box([-2, -1, -3], [2, 1, 3]);
    const extent = vec([2, 1, 3]);

    it('the TI helper agrees with the class query', () => {
        const cases: Array<[number[], number[]]> = [
            [[-5, 0, 0], [1, 0, 0]],
            [[5, 0, 0], [1, 0, 0]],
            [[0, 0, 0], [0, 1, 0]],
            [[-5, 3, 0], [1, 0, 0]],
            [[-4, -3, 0], [1, 1, 0]]
        ];
        for (const [p, d] of cases) {
            const r = ray(p, d);
            const result = defaultIntrLine3AlignedBox3TIResult();
            intrRay3AlignedBox3TIDoQuery(r.origin, r.direction, extent,
                result);
            const expected = new IntrRay3AlignedBox3TI().test(r, b);
            expect(result.intersect).toBe(expected.intersect);
        }
    });

    it('the FI helper clips to the ray t-interval and leaves points alone', () => {
        const result = defaultIntrLine3AlignedBox3FIResult();
        intrRay3AlignedBox3FIDoQuery(vec([0, 0, 0]), vec([1, 0, 0]), extent,
            result);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The ray starts inside the box, so t0 is clamped to 0.
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([0, 0, 0]);
    });

    it('the FI helper rejects a ray pointing away from the box', () => {
        const result = defaultIntrLine3AlignedBox3FIResult();
        intrRay3AlignedBox3FIDoQuery(vec([-5, 0, 0]), vec([-1, 0, 0]), extent,
            result);
        expect(result.intersect).toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, unitVector, wellScaledVector, expectClose, expectVectorClose,
    seededRandom
} from './helpers/arbitraries.js';
import { sub } from '../src/Vector.js';

// A well-scaled aligned box: the shared 'alignedBox' generator draws its
// corners from fc.double, which emits subnormal coordinates whose centered
// form has no significant digits left.
const box3Arb = fc.tuple(wellScaledVector(3, -4, 4),
    fc.double({ min: 0.2, max: 4, noNaN: true }),
    fc.double({ min: 0.2, max: 4, noNaN: true }),
    fc.double({ min: 0.2, max: 4, noNaN: true }))
    .map(([lo, w, h, d]) => AlignedBox.fromMinMax(lo,
        Vector.fromArray([lo.values[0] + w, lo.values[1] + h,
            lo.values[2] + d])));

const rayBox3 = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3), box3Arb)
    .map(([o, d, b]) => ({ ray: Ray.fromOriginDirection(o, d), box: b }));

function insideBox3(box: AlignedBox, p: Vector, tol: number): boolean {
    for (let i = 0; i < 3; ++i) {
        if (p.values[i] < box.min.values[i] - tol
            || p.values[i] > box.max.values[i] + tol) {
            return false;
        }
    }
    return true;
}

describe('IntrRay3AlignedBox3 verification', () => {
    const tiq = new IntrRay3AlignedBox3TI();
    const fiq = new IntrRay3AlignedBox3FI();

    it('TI and FI agree on intersect', () => {
        check(rayBox3, ({ ray: r, box: b }) => {
            expect(tiq.test(r, b).intersect).toBe(fiq.find(r, b).intersect);
        });
    });

    it('a ray hit is the line hit clipped to nonnegative parameters', () => {
        const lfi = new IntrLine3AlignedBox3FI();
        check(rayBox3, ({ ray: r, box: b }) => {
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lf = lfi.find(l, b);
            const f = fiq.find(r, b);
            if (!lf.intersect) {
                expect(f.intersect).toBe(false);
                return;
            }
            // Upstream clips against the semi-infinite interval [0,+inf).
            expect(f.intersect).toBe(lf.parameter[1] >= 0);
            if (!f.intersect) {
                // The 'result = Result{}' reset restores every field.
                expect(f.numIntersections).toBe(0);
                expect(f.parameter).toEqual([0, 0]);
                expect(f.point[0].values).toEqual([0, 0, 0]);
                expect(f.point[1].values).toEqual([0, 0, 0]);
                return;
            }
            expectClose(f.parameter[0], Math.max(lf.parameter[0], 0), 0, 0);
            expectClose(f.parameter[1], lf.parameter[1], 0, 0);
        });
    });

    it('the reported points are on the ray and inside the box', () => {
        check(rayBox3, ({ ray: r, box: b }) => {
            const f = fiq.find(r, b);
            if (!f.intersect) {
                return;
            }
            expect(f.parameter[0]).toBeGreaterThanOrEqual(0);
            expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
            // Upstream fills both point entries whenever intersect is true.
            for (let i = 0; i < 2; ++i) {
                expectVectorClose(f.point[i],
                    add(r.origin, mul(f.parameter[i], r.direction)), 0, 0);
                expect(insideBox3(b, f.point[i], 1e-9)).toBe(true);
            }
        });
    });

    it('a fine sweep of the ray agrees with the reported interval', () => {
        const rnd = seededRandom(0x2b7d54e9);
        for (let trial = 0; trial < 120; ++trial) {
            const lo = Vector.fromArray([rnd() * 4 - 4, rnd() * 4 - 4,
                rnd() * 4 - 4]);
            const b = AlignedBox.fromMinMax(lo, Vector.fromArray([
                lo.values[0] + 0.5 + rnd() * 3,
                lo.values[1] + 0.5 + rnd() * 3,
                lo.values[2] + 0.5 + rnd() * 3]));
            const dir = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            const n = Math.hypot(dir.values[0], dir.values[1], dir.values[2]);
            if (n < 0.3) {
                continue;
            }
            normalize(dir);
            const r = Ray.fromOriginDirection(
                Vector.fromArray([rnd() * 10 - 5, rnd() * 10 - 5,
                    rnd() * 10 - 5]), dir);
            const f = fiq.find(r, b);
            for (let k = 0; k <= 500; ++k) {
                const t = (16 * k) / 500;
                if (insideBox3(b, add(r.origin, mul(t, dir)), -1e-6)) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-9);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-9);
                }
            }
        }
    }, 30000);

    it('a ray whose origin is inside the box hits at parameter 0', () => {
        check(fc.tuple(box3Arb, unitVector(3),
            fc.double({ min: 0.05, max: 0.95, noNaN: true }),
            fc.double({ min: 0.05, max: 0.95, noNaN: true }),
            fc.double({ min: 0.05, max: 0.95, noNaN: true })),
            ([b, d, s0, s1, s2]) => {
                const s = [s0, s1, s2];
                const o = Vector.fromArray([0, 1, 2].map(i =>
                    b.min.values[i]
                        + s[i] * (b.max.values[i] - b.min.values[i])));
                const r = Ray.fromOriginDirection(o, d);
                expect(tiq.test(r, b).intersect).toBe(true);
                const f = fiq.find(r, b);
                expect(f.intersect).toBe(true);
                expectClose(f.parameter[0], 0, 1e-12, 0);
            });
    });

    it('the exported DoQuery helpers reproduce the class results', () => {
        check(rayBox3, ({ ray: r, box: b }) => {
            const cf = b.getCenteredForm();
            const localOrigin = sub(r.origin, cf.center);
            const tres = defaultIntrLine3AlignedBox3TIResult();
            intrRay3AlignedBox3TIDoQuery(localOrigin, r.direction, cf.extent,
                tres);
            expect(tres.intersect).toBe(tiq.test(r, b).intersect);

            const fres = defaultIntrLine3AlignedBox3FIResult();
            intrRay3AlignedBox3FIDoQuery(localOrigin, r.direction, cf.extent,
                fres);
            const f = fiq.find(r, b);
            expect(fres.intersect).toBe(f.intersect);
            expect(fres.numIntersections).toBe(f.numIntersections);
            if (fres.intersect) {
                expect(fres.parameter[0]).toBe(f.parameter[0]);
                expect(fres.parameter[1]).toBe(f.parameter[1]);
            }
            // The helper leaves point[] untouched; the class fills it.
            expect(fres.point[0].values).toEqual([0, 0, 0]);
        });
    });
});
