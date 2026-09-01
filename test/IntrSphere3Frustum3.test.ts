import { describe, it, expect } from 'vitest';
import { DistPoint3Frustum3 } from '../src/DistPoint3Frustum3';
import { Frustum3 } from '../src/Frustum3';
import { Hypersphere } from '../src/Hypersphere';
import {
    IntrSphere3Frustum3TI,
    defaultIntrSphere3Frustum3TIResult
} from '../src/IntrSphere3Frustum3';
import { Vector, add, mul, normalize } from '../src/Vector';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function sphere(center: Vector, radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(center, radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrSphere3Frustum3', () => {
    const query = new IntrSphere3Frustum3TI();

    // The frustum with apex at the origin looking down +z: the near face is
    // the square |x| <= 1, |y| <= 1 at z = 1 and the far face is the square
    // |x| <= 3, |y| <= 3 at z = 3. The four side planes are x = +-z and
    // y = +-z.
    const frustum = Frustum3.fromParameters(vec(0, 0, 0), vec(0, 0, 1),
        vec(0, 1, 0), vec(1, 0, 0), 1, 3, 1, 1);

    // Assert that a sphere at 'center' whose closest frustum point is at
    // distance 'd' touches the frustum at radius d: it misses for a slightly
    // smaller radius and intersects for a slightly larger one. When the
    // geometry is exact in binary floating point (axis-aligned faces),
    // 'exact' asks for the verdict at radius exactly d as well.
    function expectTouchesAt(center: Vector, d: number,
        exact: boolean = false): void {
        const eps = exact ? 1e-9 : 1e-6;
        expect(query.test(sphere(center, d * (1 - eps)), frustum).intersect)
            .toBe(false);
        expect(query.test(sphere(center, d * (1 + eps)), frustum).intersect)
            .toBe(true);
        if (exact) {
            // The query uses 'distance <= radius', so an exact touch counts
            // as an intersection.
            expect(query.test(sphere(center, d), frustum).intersect)
                .toBe(true);
        }
        // Cross-check the assumed distance against the point-frustum query.
        expect(new DistPoint3Frustum3().compute(center, frustum).distance)
            .toBeCloseTo(d, 9);
    }

    it('default results report no intersection', () => {
        expect(defaultIntrSphere3Frustum3TIResult().intersect).toBe(false);
    });

    it('reports an intersection for a sphere inside the frustum', () => {
        expect(query.test(sphere(vec(0, 0, 2), 0.25), frustum).intersect)
            .toBe(true);
        // A degenerate zero-radius sphere is the frustum-contains-point test.
        expect(query.test(sphere(vec(0, 0, 2), 0), frustum).intersect)
            .toBe(true);
        expect(query.test(sphere(vec(0, 0, 0.5), 0), frustum).intersect)
            .toBe(false);
    });

    it('reports an intersection for a sphere containing the frustum', () => {
        expect(query.test(sphere(vec(0, 0, 2), 100), frustum).intersect)
            .toBe(true);
    });

    it('reports no intersection for a sphere well outside the frustum', () => {
        expect(query.test(sphere(vec(0, 0, 10), 1), frustum).intersect)
            .toBe(false);
        expect(query.test(sphere(vec(20, 0, 2), 1), frustum).intersect)
            .toBe(false);
        expect(query.test(sphere(vec(0, 0, -5), 1), frustum).intersect)
            .toBe(false);
    });

    it('touches the near face', () => {
        // The nearest frustum point to (0,0,0.5) is (0,0,1).
        expectTouchesAt(vec(0, 0, 0.5), 0.5, true);
    });

    it('touches the far face', () => {
        // The nearest frustum point to (0,0,3.5) is (0,0,3).
        expectTouchesAt(vec(0, 0, 3.5), 0.5, true);
    });

    it('touches each of the four side faces', () => {
        const h = 0.5;
        const s = Math.SQRT1_2;
        // Interior points of the four side faces at z = 2, with the outward
        // unit normals of the planes x = z, x = -z, y = z, y = -z.
        const cases: Array<[Vector, Vector]> = [
            [vec(2, 0, 2), vec(s, 0, -s)],
            [vec(-2, 0, 2), vec(-s, 0, -s)],
            [vec(0, 2, 2), vec(0, s, -s)],
            [vec(0, -2, 2), vec(0, -s, -s)]
        ];
        for (const [facePoint, outward] of cases) {
            expectTouchesAt(add(facePoint, mul(h, outward)), h);
        }
    });

    it('touches a side edge of the frustum', () => {
        // The edge shared by the planes x = z and y = z is the ray x = y = z;
        // (2,2,2) is an interior point of that edge. The outward bisector of
        // the two face normals is (1,1,-2)/sqrt(6).
        const outward = vec(1, 1, -2);
        normalize(outward);
        const h = 0.75;
        expectTouchesAt(add(vec(2, 2, 2), mul(h, outward)), h);
    });

    it('touches a near-face corner', () => {
        // (1,1,1) is a corner of the near face; approaching it along -D puts
        // the closest frustum point exactly at that corner.
        expectTouchesAt(vec(1, 1, 1 - 0.4), 0.4, true);
    });

    it('touches a far-face corner', () => {
        // (3,3,3) is a corner of the far face. The normal cone there is
        // spanned by the far normal (0,0,1) and the two side normals; their
        // sum is an interior direction of that cone.
        const s = Math.SQRT1_2;
        const outward = vec(s, s, 1 - Math.SQRT2);
        normalize(outward);
        const h = 0.6;
        expectTouchesAt(add(vec(3, 3, 3), mul(h, outward)), h);
    });

    it('touches the apex region in front of the near face', () => {
        // The frustum apex is at the origin but the solid starts at z = 1,
        // so the distance from the origin to the frustum is 1.
        expectTouchesAt(vec(0, 0, 0), 1, true);
    });

    it('is invariant under a rigid motion of sphere and frustum', () => {
        // The same frustum expressed in a rotated, translated frame, with
        // D, U, R rotated 90 degrees about the y axis: D -> (1,0,0),
        // U -> (0,1,0), R -> (0,0,-1), and the origin moved to (5,-2,7).
        const moved = Frustum3.fromParameters(vec(5, -2, 7), vec(1, 0, 0),
            vec(0, 1, 0), vec(0, 0, -1), 1, 3, 1, 1);
        const mapPoint = (p: Vector) => vec(
            5 + p.values[2], -2 + p.values[1], 7 - p.values[0]);

        const random = makeRandom(90043);
        let hits = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const p = vec(8 * random() - 4, 8 * random() - 4,
                8 * random() - 1);
            const radius = 2 * random();
            const expected = query.test(sphere(p, radius), frustum).intersect;
            const actual = query.test(sphere(mapPoint(p), radius), moved)
                .intersect;
            expect(actual).toBe(expected);
            if (expected) {
                ++hits;
            }
        }
        expect(hits).toBeGreaterThan(20);
        expect(hits).toBeLessThan(280);
    });

    it('agrees with a brute-force sampling of the frustum', () => {
        // Sample the frustum solid densely and compare the minimum distance
        // from the sphere center to the samples against the query verdict.
        const samples: Vector[] = [];
        const n = 24;
        for (let i = 0; i <= n; ++i) {
            const z = 1 + 2 * (i / n);
            for (let j = 0; j <= n; ++j) {
                for (let k = 0; k <= n; ++k) {
                    samples.push(vec(z * (2 * (j / n) - 1),
                        z * (2 * (k / n) - 1), z));
                }
            }
        }

        const random = makeRandom(90044);
        for (let trial = 0; trial < 200; ++trial) {
            const center = vec(10 * random() - 5, 10 * random() - 5,
                8 * random() - 1);
            let best = Number.MAX_VALUE;
            for (const s of samples) {
                const dx = center.values[0] - s.values[0];
                const dy = center.values[1] - s.values[1];
                const dz = center.values[2] - s.values[2];
                const value = dx * dx + dy * dy + dz * dz;
                if (value < best) {
                    best = value;
                }
            }
            // The sampling is an upper bound on the true distance; it is
            // never smaller, so a radius above it must intersect and a
            // radius far below the sampled distance must not.
            const sampled = Math.sqrt(best);
            expect(query.test(sphere(center, sampled), frustum).intersect)
                .toBe(true);
            if (sampled > 0.5) {
                expect(query.test(sphere(center, sampled - 0.5), frustum)
                    .intersect).toBe(false);
            }
        }
    });
});
