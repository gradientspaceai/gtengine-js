import { describe, expect, it } from 'vitest';
import { DistPoint3Frustum3 } from '../src/DistPoint3Frustum3.js';
import { Frustum3 } from '../src/Frustum3.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The canonical frustum: origin at (0,0,0), D = +z, U = +y, R = +x.
function canonical(dMin: number, dMax: number, uBound: number,
    rBound: number): Frustum3 {
    return Frustum3.fromParameters(v(0, 0, 0), v(0, 0, 1), v(0, 1, 0),
        v(1, 0, 0), dMin, dMax, uBound, rBound);
}

describe('DistPoint3Frustum3', () => {
    const query = new DistPoint3Frustum3();

    it('reports zero distance for a point inside the frustum', () => {
        const f = canonical(1, 3, 1, 1);
        const result = query.compute(v(0, 0, 2), f);
        expect(result.distance).toBe(0);
        expect(result.closest[1].values[2]).toBeCloseTo(2, 12);
    });

    it('measures a point beyond the far plane', () => {
        const f = canonical(1, 3, 1, 1);
        const result = query.compute(v(0, 0, 7), f);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(3, 12);
    });

    it('measures a point behind the near plane', () => {
        const f = canonical(1, 3, 1, 1);
        const result = query.compute(v(0, 0, -2), f);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('clamps to the far-plane rectangle corner', () => {
        const f = canonical(1, 3, 1, 1);
        // rmax = umax = 3 at the far plane z = 3.
        const result = query.compute(v(6, 7, 3), f);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(3, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(3, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(3, 12);
    });

    it('is symmetric under reflection in the R and U directions', () => {
        const f = canonical(1, 4, 0.75, 1.5);
        const rnd = makeRandom(555);
        for (let trial = 0; trial < 40; ++trial) {
            const p = v(10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5);
            const r0 = query.compute(p, f);
            const r1 = query.compute(
                v(-p.values[0], p.values[1], p.values[2]), f);
            const r2 = query.compute(
                v(p.values[0], -p.values[1], p.values[2]), f);
            expect(r1.distance).toBeCloseTo(r0.distance, 10);
            expect(r2.distance).toBeCloseTo(r0.distance, 10);
            expect(r1.closest[1].values[0]).toBeCloseTo(
                -r0.closest[1].values[0], 10);
            expect(r2.closest[1].values[1]).toBeCloseTo(
                -r0.closest[1].values[1], 10);
        }
    });

    it('is invariant to a rigid motion of point and frustum', () => {
        const rnd = makeRandom(8080);
        const a = 0.9;
        const ca = Math.cos(a), sa = Math.sin(a);
        const rot = (p: Vector): Vector => v(
            ca * p.values[0] - sa * p.values[2],
            p.values[1],
            sa * p.values[0] + ca * p.values[2]);
        const shift = v(2, -3, 1);

        const f0 = canonical(1, 3, 1, 2);
        const f1 = Frustum3.fromParameters(add(shift, v(0, 0, 0)),
            rot(v(0, 0, 1)), rot(v(0, 1, 0)), rot(v(1, 0, 0)), 1, 3, 1, 2);

        for (let trial = 0; trial < 40; ++trial) {
            const p = v(10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5);
            const r0 = query.compute(p, f0);
            const r1 = query.compute(add(rot(p), shift), f1);
            expect(r1.distance).toBeCloseTo(r0.distance, 10);
            const expected = add(rot(r0.closest[1]), shift);
            for (let i = 0; i < 3; ++i) {
                expect(r1.closest[1].values[i]).toBeCloseTo(
                    expected.values[i], 9);
            }
        }
    });

    it('agrees with a dense sampling of the solid frustum', () => {
        const rnd = makeRandom(13579);
        const origin = v(0.5, -1, 0.25);
        const dVec = v(1, 1, 2);
        normalize(dVec);
        // A right-handed orthonormal frame containing dVec.
        const uVec = v(-1, 1, 0);
        normalize(uVec);
        const rVec = v(
            uVec.values[1] * dVec.values[2] - uVec.values[2] * dVec.values[1],
            uVec.values[2] * dVec.values[0] - uVec.values[0] * dVec.values[2],
            uVec.values[0] * dVec.values[1] - uVec.values[1] * dVec.values[0]);
        normalize(rVec);

        const dMin = 1, dMax = 2.5, uBound = 0.8, rBound = 1.2;
        const f = Frustum3.fromParameters(origin, dVec, uVec, rVec,
            dMin, dMax, uBound, rBound);

        for (let trial = 0; trial < 25; ++trial) {
            const p = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const result = query.compute(p, f);

            // The reported closest point is in the solid frustum.
            const delta = sub(result.closest[1], origin);
            const z = dot(delta, dVec);
            const x = dot(delta, rVec);
            const y = dot(delta, uVec);
            expect(z).toBeGreaterThanOrEqual(dMin - 1e-9);
            expect(z).toBeLessThanOrEqual(dMax + 1e-9);
            expect(Math.abs(x)).toBeLessThanOrEqual(
                rBound * z / dMin + 1e-9);
            expect(Math.abs(y)).toBeLessThanOrEqual(
                uBound * z / dMin + 1e-9);

            // The reported closest point realizes the reported distance.
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 9);

            // No sampled solid-frustum point is closer.
            const n = 24;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                const zz = dMin + (dMax - dMin) * i / n;
                const scale = zz / dMin;
                for (let j = 0; j <= n; ++j) {
                    const xx = rBound * scale * (2 * j / n - 1);
                    for (let k = 0; k <= n; ++k) {
                        const yy = uBound * scale * (2 * k / n - 1);
                        const q = add(origin, add(mul(zz, dVec),
                            add(mul(xx, rVec), mul(yy, uVec))));
                        const g = sub(p, q);
                        best = Math.min(best, dot(g, g));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});
