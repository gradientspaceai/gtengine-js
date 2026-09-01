import { describe, it, expect } from 'vitest';
import { Segment } from '../src/Segment';
import { Hypersphere } from '../src/Hypersphere';
import { Vector, add, mul, sub, dot } from '../src/Vector';
import {
    IntrSegment3Sphere3TI,
    IntrSegment3Sphere3FI,
    defaultIntrSegment3Sphere3FIResult,
    intrSegment3Sphere3DoQuery
} from '../src/IntrSegment3Sphere3';

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

    it('reports the whole segment when it lies inside the sphere (FI), with the upstream TI quirk', () => {
        const s = segment([-0.25, 0, 0], [0.25, 0, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-0.25, 12);
        expect(result.parameter[1]).toBeCloseTo(0.25, 12);

        // Upstream bug (preserved): when the whole segment is strictly inside
        // the solid sphere, Q(-e) < 0 and Q(e) < 0, and the TI query's final
        // test 'qm > 0 && |a1| < e' returns false. The comment above that
        // test asserts the segment does not intersect the sphere, but the
        // sphere is documented as a solid, so it does. Compare
        // IntrRay3Sphere3TI, which explicitly returns true when the ray
        // origin is inside the sphere.
        expect(ti.test(s, unit).intersect).toBe(false);
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
        // The extent is 0, so Q(-e) = Q(e) = a0 < 0 and the TI query hits
        // the same "entirely inside the solid sphere" quirk described above.
        expect(ti.test(inside, unit).intersect).toBe(false);
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
            // The TI and FI queries agree except for the preserved upstream
            // quirk where the whole segment is strictly inside the sphere.
            const d0 = sub(p0, s.center);
            const d1 = sub(p1, s.center);
            const rsqr = s.radius * s.radius;
            const contained = dot(d0, d0) < rsqr && dot(d1, d1) < rsqr;
            if (tiResult.intersect !== fiResult.intersect && !contained) {
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
