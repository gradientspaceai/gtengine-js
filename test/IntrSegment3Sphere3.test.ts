import { describe, it, expect } from 'vitest';
import { Segment } from '../src/Segment.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, add, mul, sub, dot } from '../src/Vector.js';
import {
    IntrSegment3Sphere3TI,
    IntrSegment3Sphere3FI,
    defaultIntrSegment3Sphere3FIResult,
    intrSegment3Sphere3DoQuery
} from '../src/IntrSegment3Sphere3.js';
import { Line } from '../src/Line.js';
import { IntrLine3Sphere3FI } from '../src/IntrLine3Sphere3.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
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

describe('IntrSegment3Sphere3', () => {
    const ti = new IntrSegment3Sphere3TI();
    const fi = new IntrSegment3Sphere3FI();
    const unit = sphere([0, 0, 0], 1);

    it('finds the two crossings of a segment through the sphere center', () => {
        // Center (0,0,0), unit direction (1,0,0), extent 3.
        const s = segment([-3, 0, 0], [3, 0, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // Parameters are relative to the centered form.
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips to the segment when one endpoint is inside the sphere', () => {
        // Center (0.5,0,0), unit direction (1,0,0), extent 2.5.
        const s = segment([-2, 0, 0], [3, 0, 0]);
        const result = fi.find(s, unit);
        expect(ti.test(s, unit).intersect).toBe(true);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports the whole segment when it lies inside the sphere', () => {
        const s = segment([-0.25, 0, 0], [0.25, 0, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-0.25, 12);
        expect(result.parameter[1]).toBeCloseTo(0.25, 12);

        // Upstream bug (FIXED; see upstream-bug issue (B71)): when the whole
        // segment is strictly inside the solid sphere, Q(-e) < 0 and Q(e) < 0
        // and upstream's final test 'qm > 0 && |a1| < e' returned false,
        // disagreeing with the FI query above. The port tests the endpoints
        // for containment first, as IntrRay3Sphere3TI does.
        expect(ti.test(s, unit).intersect).toBe(true);
    });

    it('reports a single point for a tangent segment', () => {
        const s = segment([-2, 1, 0], [2, 1, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
    });

    it('reports the endpoint when the segment just touches the sphere', () => {
        // The endpoint (1,0,0) is on the sphere; the segment goes outward.
        const s = segment([1, 0, 0], [4, 0, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
    });

    it('misses when the segment stops short of the sphere', () => {
        const s = segment([-5, 0, 0], [-2, 0, 0]);
        expect(ti.test(s, unit).intersect).toBe(false);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('misses when the supporting line misses the sphere', () => {
        const s = segment([-5, 3, 0], [5, 3, 0]);
        expect(ti.test(s, unit).intersect).toBe(false);
        expect(fi.find(s, unit).intersect).toBe(false);
    });

    it('handles a degenerate (zero-length) segment', () => {
        const inside = segment([0.25, 0, 0], [0.25, 0, 0]);
        // The extent is 0, so Q(-e) = Q(e) = a0 < 0: the single point is
        // inside the solid sphere and both queries report an intersection.
        expect(ti.test(inside, unit).intersect).toBe(true);
        expect(fi.find(inside, unit).intersect).toBe(true);
        const outside = segment([5, 0, 0], [5, 0, 0]);
        expect(ti.test(outside, unit).intersect).toBe(false);
        expect(fi.find(outside, unit).intersect).toBe(false);
    });

    it('handles a zero-radius sphere', () => {
        const degenerate = sphere([1, 0, 0], 0);
        const hit = segment([0, 0, 0], [2, 0, 0]);
        expect(ti.test(hit, degenerate).intersect).toBe(true);
        const result = fi.find(hit, degenerate);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);

        const miss = segment([0, 1, 0], [2, 1, 0]);
        expect(ti.test(miss, degenerate).intersect).toBe(false);
        expect(fi.find(miss, degenerate).intersect).toBe(false);
    });

    it('the exported DoQuery matches the class query', () => {
        const s = segment([-2, 0.3, -0.4], [2.5, 0.1, 0.6]);
        const cf = s.getCenteredForm();
        const direct = defaultIntrSegment3Sphere3FIResult();
        intrSegment3Sphere3DoQuery(cf.center, cf.direction, cf.extent, unit,
            direct);
        const viaClass = fi.find(s, unit);
        expect(direct.intersect).toBe(viaClass.intersect);
        expect(direct.numIntersections).toBe(viaClass.numIntersections);
        expect(direct.parameter[0]).toBeCloseTo(viaClass.parameter[0], 12);
        expect(direct.parameter[1]).toBeCloseTo(viaClass.parameter[1], 12);
    });

    it('TI and FI agree on random segments versus a sphere, including contained ones', () => {
        const rnd = makeRandom(20260901);
        let mismatch = 0;
        let containedCases = 0;
        let crossingCases = 0;
        let missCases = 0;

        for (let trial = 0; trial < 600; ++trial) {
            const s = sphere([2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1],
                0.5 + 1.5 * rnd());
            let p0: Vector;
            let p1: Vector;
            const mode = trial % 3;
            if (mode === 0) {
                // Force both endpoints strictly inside the sphere: this is
                // the case the upstream TI query got wrong.
                const r = 0.9 * s.radius;
                p0 = add(s.center, vec(2 * rnd() - 1, 2 * rnd() - 1,
                    2 * rnd() - 1));
                p1 = add(s.center, vec(2 * rnd() - 1, 2 * rnd() - 1,
                    2 * rnd() - 1));
                const scale = (p: Vector): Vector => {
                    const d = sub(p, s.center);
                    const len = Math.sqrt(dot(d, d));
                    return len > r
                        ? add(s.center, mul(r / len, d))
                        : p;
                };
                p0 = scale(p0);
                p1 = scale(p1);
                ++containedCases;
            }
            else if (mode === 1) {
                // A segment aimed through the sphere.
                p0 = vec(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
                p1 = add(p0, mul(2, sub(s.center, p0)));
                ++crossingCases;
            }
            else {
                // A mostly-random segment; many of these miss.
                p0 = vec(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
                p1 = vec(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
                ++missCases;
            }

            const d = sub(p1, p0);
            if (dot(d, d) < 1e-8) {
                continue;
            }
            const seg = Segment.fromEndpoints(p0, p1);
            if (ti.test(seg, s).intersect !== fi.find(seg, s).intersect) {
                ++mismatch;
            }
        }

        expect(containedCases).toBeGreaterThan(50);
        expect(crossingCases).toBeGreaterThan(50);
        expect(missCases).toBeGreaterThan(50);
        expect(mismatch).toBe(0);
    });

    it('agrees with a dense sampling of the segment on random configurations', () => {
        const rnd = makeRandom(987654321);
        let tiFiMismatch = 0;
        let sampleMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;
        const samples = 2000;

        for (let trial = 0; trial < 200; ++trial) {
            const p0 = vec(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const s = sphere([2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1],
                0.25 + 1.25 * rnd());
            // Aim the segment roughly at the sphere, then jitter.
            const p1 = add(p0, add(mul(2, sub(s.center, p0)),
                vec(3 * rnd() - 1.5, 3 * rnd() - 1.5, 3 * rnd() - 1.5)));
            const seg = Segment.fromEndpoints(p0, p1);
            const d = sub(p1, p0);
            if (dot(d, d) < 1e-8) {
                continue;
            }

            const tiResult = ti.test(seg, s);
            const fiResult = fi.find(seg, s);
            if (tiResult.intersect !== fiResult.intersect) {
                ++tiFiMismatch;
            }

            // Brute force: sample the segment and test sphere containment.
            let sampled = false;
            for (let k = 0; k <= samples; ++k) {
                const p = add(p0, mul(k / samples, d));
                const diff = sub(p, s.center);
                if (dot(diff, diff) <= s.radius * s.radius) {
                    sampled = true;
                    break;
                }
            }
            if (sampled !== fiResult.intersect) {
                // A sampled hit must be reported; a reported hit whose
                // overlap is shorter than the sample spacing may be missed by
                // the sampling, so only flag the first direction.
                if (sampled) {
                    ++sampleMismatch;
                }
            }

            if (fiResult.intersect) {
                ++hits;
                const cf = seg.getCenteredForm();
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
                    // The reported point is on the sphere unless the
                    // parameter was clipped to a segment endpoint.
                    const rel = sub(expected, s.center);
                    const onSphere =
                        Math.abs(Math.sqrt(dot(rel, rel)) - s.radius) < 1e-8;
                    const clipped = Math.abs(Math.abs(t) - cf.extent) < 1e-9;
                    if (!onSphere && !clipped) {
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
// IntrSegment3Sphere3.h.
// ---------------------------------------------------------------------------

describe('IntrSegment3Sphere3 verification', () => {
    const ti = new IntrSegment3Sphere3TI();
    const fi = new IntrSegment3Sphere3FI();
    const lineFi = new IntrLine3Sphere3FI();

    const arbSegment = fc.tuple(wellScaledVector(3, -6, 6),
        wellScaledVector(3, -6, 6))
        .filter(([a, b]) => {
            const d = sub(b, a);
            return dot(d, d) > 1e-2;
        })
        .map(([a, b]) => Segment.fromEndpoints(a, b));
    const arbSphere = fc.tuple(wellScaledVector(3, -5, 5), positive(4, 0.25))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));
    const arbPair = fc.tuple(arbSegment, arbSphere);

    it('TI and FI agree on intersect (both treat the sphere as a solid)', () => {
        // Upstream's TI query misses segments wholly inside the sphere; the
        // port fixes that (upstream-bug issue #203), so the two queries agree
        // for every configuration.
        check(arbPair, ([s, sph]) => {
            expect(ti.test(s, sph).intersect).toBe(fi.find(s, sph).intersect);
        });
    });

    it('TI matches a finely sampled solid-sphere containment test', () => {
        check(arbPair, ([s, sph]) => {
            const rsqr = sph.radius * sph.radius;
            let inside = false;
            for (let k = 0; k <= 512 && !inside; ++k) {
                const t = k / 512;
                const p = add(mul(s.p[0], 1 - t), mul(s.p[1], t));
                const diff = sub(p, sph.center);
                inside = dot(diff, diff) < rsqr - 1e-9;
            }
            if (inside) {
                expect(ti.test(s, sph).intersect).toBe(true);
            }
        });
    });

    it('a segment with both endpoints strictly inside the sphere intersects it',
        () => {
            // This is the case upstream's 'qm * qp <= 0' test misses: Q(-e)
            // and Q(e) are both negative, so their product is positive.
            check(arbPair, ([s, sph]) => {
                const inside = (p: Vector): boolean => {
                    const d = sub(p, sph.center);
                    return dot(d, d) < sph.radius * sph.radius;
                };
                if (inside(s.p[0]) && inside(s.p[1])) {
                    expect(ti.test(s, sph).intersect).toBe(true);
                    const res = fi.find(s, sph);
                    expect(res.intersect).toBe(true);
                    expect(res.numIntersections).toBe(2);
                    // The whole segment is the intersection set.
                    const e = s.getCenteredForm().extent;
                    expectClose(res.parameter[0], -e, 1e-9, 1e-9);
                    expectClose(res.parameter[1], e, 1e-9, 1e-9);
                }
            });
        });

    it('FI reports centered-form parameters in [-e,e] and points on the'
        + ' segment', () => {
        check(arbPair, ([s, sph]) => {
            const res = fi.find(s, sph);
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
                expectVectorClose(res.point[i], add(center, mul(t, direction)),
                    1e-12, 1e-12);
                // The point is on the sphere unless the interval was clipped
                // to a segment endpoint inside the solid sphere.
                const diff = sub(res.point[i], sph.center);
                const rsqr = sph.radius * sph.radius;
                const onSurface = Math.abs(dot(diff, diff) - rsqr)
                    <= 1e-8 * (1 + rsqr);
                const clipped = Math.abs(Math.abs(t) - extent) <= 1e-12 * extent
                    && dot(diff, diff) <= rsqr + 1e-9;
                expect(onSurface || clipped).toBe(true);
            }
            expect(res.parameter[0]).toBeLessThanOrEqual(res.parameter[1]);
        });
    });

    it('the segment result is the line result clipped to [-e,e]', () => {
        check(arbPair, ([s, sph]) => {
            const { center, direction, extent } = s.getCenteredForm();
            const lineRes = lineFi.find(
                Line.fromOriginDirection(center, direction), sph);
            const segRes = fi.find(s, sph);
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

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(arbPair, rotationFrame(3), wellScaledVector(3, -4, 4)),
            ([[s, sph], frame, shift]) => {
                const rot = (v: Vector): Vector => Vector.fromArray([
                    dot(frame[0], v), dot(frame[1], v), dot(frame[2], v)]);
                const map = (v: Vector): Vector => add(rot(v), shift);
                const s2 = Segment.fromEndpoints(map(s.p[0]), map(s.p[1]));
                const sph2 = Hypersphere.fromCenterRadius(map(sph.center),
                    sph.radius);
                const a = fi.find(s, sph), b = fi.find(s2, sph2);
                expect(b.intersect).toBe(a.intersect);
                if (a.intersect) {
                    expect(b.numIntersections).toBe(a.numIntersections);
                    for (let i = 0; i < a.numIntersections; ++i) {
                        expectClose(b.parameter[i], a.parameter[i], 1e-7, 1e-7);
                        expectVectorClose(b.point[i], map(a.point[i]),
                            1e-6, 1e-6);
                    }
                }
            });
    });

    it('reports a contained segment (upstream-bug regression)', () => {
        // Upstream's TI query concludes "no intersection" for a segment whose
        // endpoints are both strictly inside the solid sphere, because
        // Q(-e) * Q(e) > 0 there and the follow-up test requires Q(-e) > 0.
        const s = segment([-0.5, 0, 0], [0.5, 0, 0]);
        const unitSphere = sphere([0, 0, 0], 1);
        expect(ti.test(s, unitSphere).intersect).toBe(true);
        const res = fi.find(s, unitSphere);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(2);
        expect(res.parameter[0]).toBe(-0.5);
        expect(res.parameter[1]).toBe(0.5);
        expectVectorClose(res.point[0], vec(-0.5, 0, 0));
        expectVectorClose(res.point[1], vec(0.5, 0, 0));
    });

    it('reports a tangent segment as a single point', () => {
        const s = segment([-2, 1, 0], [2, 1, 0]);
        const unitSphere = sphere([0, 0, 0], 1);
        expect(ti.test(s, unitSphere).intersect).toBe(true);
        const res = fi.find(s, unitSphere);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(1);
        expectVectorClose(res.point[0], vec(0, 1, 0));
    });

    it('the DoQuery helper resets a result outside the segment domain', () => {
        const res = defaultIntrSegment3Sphere3FIResult();
        // Centered segment of extent 1 at x = 5 along +x; the unit sphere at
        // the origin is met by the line for t in [-6,-4].
        intrSegment3Sphere3DoQuery(vec(5, 0, 0), vec(1, 0, 0), 1,
            sphere([0, 0, 0], 1), res);
        expect(res.intersect).toBe(false);
        expect(res.numIntersections).toBe(0);
        expect(res.parameter).toEqual([0, 0]);
    });
});
