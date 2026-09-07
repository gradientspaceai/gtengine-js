import { describe, it, expect } from 'vitest';
import { Ray } from '../src/Ray.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Line } from '../src/Line.js';
import { IntrLine3Sphere3FI } from '../src/IntrLine3Sphere3.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';
import { Vector, add, mul, sub, dot, normalize } from '../src/Vector.js';
import {
    IntrRay3Sphere3TI,
    IntrRay3Sphere3FI,
    defaultIntrRay3Sphere3FIResult,
    intrRay3Sphere3DoQuery
} from '../src/IntrRay3Sphere3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(o: number[], d: number[]): Ray {
    const dir = Vector.fromArray(d);
    normalize(dir);
    return Ray.fromOriginDirection(Vector.fromArray(o), dir);
}

function sphere(c: number[], r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(c), r);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrRay3Sphere3', () => {
    const ti = new IntrRay3Sphere3TI();
    const fi = new IntrRay3Sphere3FI();
    const unit = sphere([0, 0, 0], 1);

    it('finds the two crossings of a ray through the sphere center', () => {
        const r = ray([-5, 0, 0], [1, 0, 0]);
        expect(ti.test(r, unit).intersect).toBe(true);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips the near root when the origin is inside the sphere', () => {
        const r = ray([0, 0, 0], [0, 1, 0]);
        expect(ti.test(r, unit).intersect).toBe(true);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        // The line interval is [-1,1]; clipping to [0,+inf) gives [0,1].
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[1].values[1]).toBeCloseTo(1, 12);
    });

    it('reports a single point for a tangent ray', () => {
        const r = ray([-5, 1, 0], [1, 0, 0]);
        expect(ti.test(r, unit).intersect).toBe(true);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(5, 12);
        expect(result.parameter[1]).toBeCloseTo(5, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('reports a single point when the origin is on the sphere and the ray leaves', () => {
        const r = ray([1, 0, 0], [1, 0, 0]);
        // The origin is on the sphere, so a0 == 0 and the TI query reports an
        // intersection because a0 <= 0.
        expect(ti.test(r, unit).intersect).toBe(true);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
    });

    it('misses when the sphere is behind the ray origin', () => {
        const r = ray([5, 0, 0], [1, 0, 0]);
        expect(ti.test(r, unit).intersect).toBe(false);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('misses when the line misses the sphere entirely', () => {
        const r = ray([-5, 3, 0], [1, 0, 0]);
        expect(ti.test(r, unit).intersect).toBe(false);
        expect(fi.find(r, unit).intersect).toBe(false);
    });

    it('handles a zero-radius sphere on the ray', () => {
        const degenerate = sphere([2, 0, 0], 0);
        const r = ray([0, 0, 0], [1, 0, 0]);
        expect(ti.test(r, degenerate).intersect).toBe(true);
        const result = fi.find(r, degenerate);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(2, 12);
    });

    it('handles a zero-radius sphere off the ray', () => {
        const degenerate = sphere([2, 1, 0], 0);
        const r = ray([0, 0, 0], [1, 0, 0]);
        expect(ti.test(r, degenerate).intersect).toBe(false);
        expect(fi.find(r, degenerate).intersect).toBe(false);
    });

    it('the exported DoQuery matches the class query', () => {
        const r = ray([-3, 0.5, 0.25], [1, 0.1, -0.2]);
        const direct = defaultIntrRay3Sphere3FIResult();
        intrRay3Sphere3DoQuery(r.origin, r.direction, unit, direct);
        const viaClass = fi.find(r, unit);
        expect(direct.intersect).toBe(viaClass.intersect);
        expect(direct.numIntersections).toBe(viaClass.numIntersections);
        expect(direct.parameter[0]).toBeCloseTo(viaClass.parameter[0], 12);
        expect(direct.parameter[1]).toBeCloseTo(viaClass.parameter[1], 12);
    });

    it('agrees with a dense sampling of the ray on random configurations', () => {
        const rnd = makeRandom(20250901);
        let tiFiMismatch = 0;
        let sampleMismatch = 0;
        let parameterMismatch = 0;
        let hits = 0;
        const tMax = 12;
        const samples = 800;

        for (let trial = 0; trial < 250; ++trial) {
            const o = vec(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const s = sphere([2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1],
                0.25 + 1.5 * rnd());
            // Aim roughly at the sphere so that a useful fraction of the
            // trials produce hits, then jitter to also produce misses.
            const d = add(sub(s.center, o),
                vec(3 * rnd() - 1.5, 3 * rnd() - 1.5, 3 * rnd() - 1.5));
            if (dot(d, d) < 1e-6) {
                continue;
            }
            normalize(d);
            const r = Ray.fromOriginDirection(o, d);

            const tiResult = ti.test(r, s);
            const fiResult = fi.find(r, s);
            if (tiResult.intersect !== fiResult.intersect) {
                ++tiFiMismatch;
            }

            // Brute-force sampling: does any sampled ray point lie in the
            // sphere?
            let sampled = false;
            for (let k = 0; k <= samples; ++k) {
                const t = (tMax * k) / samples;
                const p = add(o, mul(t, d));
                const diff = sub(p, s.center);
                if (dot(diff, diff) <= s.radius * s.radius) {
                    sampled = true;
                    break;
                }
            }
            if (sampled && !fiResult.intersect) {
                ++sampleMismatch;
            }

            if (fiResult.intersect) {
                ++hits;
                // Both reported parameters must be nonnegative and their
                // points must be on the sphere (or the endpoint clipped to
                // t = 0 must be inside).
                for (let i = 0; i < fiResult.numIntersections; ++i) {
                    const t = fiResult.parameter[i];
                    if (t < -1e-12) {
                        ++parameterMismatch;
                        continue;
                    }
                    const p = add(o, mul(t, d));
                    const diff = sub(p, s.center);
                    const onSphere =
                        Math.abs(Math.sqrt(dot(diff, diff)) - s.radius) < 1e-9;
                    const clippedOrigin = Math.abs(t) < 1e-12 &&
                        dot(diff, diff) <= s.radius * s.radius + 1e-9;
                    if (!onSphere && !clippedOrigin) {
                        ++parameterMismatch;
                    }
                    if (!fiResult.point[i].equals(p) &&
                        Math.abs(fiResult.point[i].values[0] - p.values[0]) > 1e-9) {
                        ++parameterMismatch;
                    }
                }
            }
        }

        expect(hits).toBeGreaterThan(20);
        expect([tiFiMismatch, sampleMismatch, parameterMismatch])
            .toEqual([0, 0, 0]);
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrRay3Sphere3.h.
// ---------------------------------------------------------------------------

describe('IntrRay3Sphere3 verification', () => {
    const ti = new IntrRay3Sphere3TI();
    const fi = new IntrRay3Sphere3FI();
    const lineFi = new IntrLine3Sphere3FI();

    // Rays and spheres with moderate, non-subnormal scales. The default
    // vector() generator emits 1e-160-scale components, which make the
    // relative comparisons below meaningless.
    const arbRay = fc.tuple(wellScaledVector(3, -6, 6), unitVector(3))
        .map(([o, d]) => Ray.fromOriginDirection(o, d));
    const arbSphere = fc.tuple(wellScaledVector(3, -6, 6), positive(4, 0.25))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));
    const arbPair = fc.tuple(arbRay, arbSphere);

    it('TI and FI agree on intersect (both treat the sphere as a solid)', () => {
        check(arbPair, ([r, s]) => {
            expect(fi.find(r, s).intersect).toBe(ti.test(r, s).intersect);
        });
    });

    it('FI reports parameters in the ray domain and points on the ray', () => {
        check(arbPair, ([r, s]) => {
            const res = fi.find(r, s);
            if (!res.intersect) {
                expect(res.numIntersections).toBe(0);
                return;
            }
            expect(res.numIntersections === 1 || res.numIntersections === 2)
                .toBe(true);
            for (let i = 0; i < res.numIntersections; ++i) {
                const t = res.parameter[i];
                expect(Number.isFinite(t)).toBe(true);
                expect(t).toBeGreaterThanOrEqual(0);
                // point[i] = origin + t * direction.
                const p = add(r.origin, mul(t, r.direction));
                expectVectorClose(res.point[i], p, 1e-12, 1e-12);
                // The point is on the sphere unless the ray origin is inside
                // the solid sphere and the interval was clipped to t = 0.
                const diff = sub(res.point[i], s.center);
                const onSurface = Math.abs(dot(diff, diff)
                    - s.radius * s.radius) <= 1e-8 * (1 + s.radius * s.radius);
                const clipped = t === 0
                    && dot(diff, diff) <= s.radius * s.radius + 1e-9;
                expect(onSurface || clipped).toBe(true);
            }
            expect(res.parameter[0]).toBeLessThanOrEqual(res.parameter[1]);
        });
    });

    it('the ray result is the line result clipped to t >= 0', () => {
        check(arbPair, ([r, s]) => {
            const lineRes = lineFi.find(
                Line.fromOriginDirection(r.origin, r.direction), s);
            const rayRes = fi.find(r, s);
            if (!lineRes.intersect) {
                expect(rayRes.intersect).toBe(false);
                return;
            }
            const t0 = lineRes.parameter[0], t1 = lineRes.parameter[1];
            const c0 = Math.max(t0, 0), c1 = t1;
            if (c0 > c1) {
                expect(rayRes.intersect).toBe(false);
                return;
            }
            expect(rayRes.intersect).toBe(true);
            expect(rayRes.parameter[0] + 0).toBe(c0 + 0);
            expect(rayRes.parameter[1] + 0).toBe(c1 + 0);
            expect(rayRes.numIntersections).toBe(c0 < c1 ? 2 : 1);
        });
    });

    it('intersect is true whenever a sampled ray point is inside the sphere',
        () => {
            check(arbPair, ([r, s]) => {
                const rsqr = s.radius * s.radius;
                let inside = false;
                for (let k = 0; k <= 400 && !inside; ++k) {
                    const t = (k * 20) / 400;
                    const diff = sub(add(r.origin, mul(t, r.direction)),
                        s.center);
                    // Shrink the radius slightly so that points within
                    // rounding distance of the surface do not force a
                    // conclusion the exact predicate need not share.
                    if (dot(diff, diff) < rsqr - 1e-9) {
                        inside = true;
                    }
                }
                if (inside) {
                    expect(ti.test(r, s).intersect).toBe(true);
                }
            });
        });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(arbPair, rotationFrame(3), wellScaledVector(3, -5, 5)),
            ([[r, s], frame, shift]) => {
                const rot = (v: Vector): Vector => {
                    const out = new Vector(3);
                    for (let i = 0; i < 3; ++i) {
                        out.set(i, dot(frame[i], v));
                    }
                    return out;
                };
                const map = (v: Vector): Vector => add(rot(v), shift);
                const r2 = Ray.fromOriginDirection(map(r.origin),
                    rot(r.direction));
                const s2 = Hypersphere.fromCenterRadius(map(s.center),
                    s.radius);
                const a = fi.find(r, s), b = fi.find(r2, s2);
                expect(b.intersect).toBe(a.intersect);
                if (a.intersect) {
                    expect(b.numIntersections).toBe(a.numIntersections);
                    for (let i = 0; i < a.numIntersections; ++i) {
                        expectClose(b.parameter[i], a.parameter[i], 1e-8, 1e-8);
                        expectVectorClose(b.point[i], map(a.point[i]),
                            1e-7, 1e-7);
                    }
                }
            });
    });

    it('is tangent-exact for a ray grazing the unit sphere', () => {
        const r = Ray.fromOriginDirection(vec(-5, 1, 0), vec(1, 0, 0));
        const s = sphere([0, 0, 0], 1);
        const res = fi.find(r, s);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(1);
        expect(res.parameter[0]).toBe(5);
        expect(res.parameter[1]).toBe(5);
        expect(ti.test(r, s).intersect).toBe(true);
    });

    it('reports no intersection for a ray pointing away from the sphere', () => {
        const r = Ray.fromOriginDirection(vec(5, 0, 0), vec(1, 0, 0));
        const s = sphere([0, 0, 0], 1);
        expect(ti.test(r, s).intersect).toBe(false);
        const res = fi.find(r, s);
        expect(res.intersect).toBe(false);
        expect(res.numIntersections).toBe(0);
        expect(res.parameter).toEqual([0, 0]);
    });

    it('the DoQuery helper leaves a non-intersecting result fully reset', () => {
        const res = defaultIntrRay3Sphere3FIResult();
        // A line that meets the sphere only for t < 0.
        intrRay3Sphere3DoQuery(vec(5, 0, 0), vec(1, 0, 0),
            sphere([0, 0, 0], 1), res);
        expect(res.intersect).toBe(false);
        expect(res.numIntersections).toBe(0);
        expect(res.parameter).toEqual([0, 0]);
        expect(res.point[0].equals(Vector.zero(3))).toBe(true);
        expect(res.point[1].equals(Vector.zero(3))).toBe(true);
    });
});
