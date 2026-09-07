import { describe, it, expect } from 'vitest';
import { DistPoint3Frustum3 } from '../src/DistPoint3Frustum3.js';
import { Frustum3 } from '../src/Frustum3.js';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    IntrSphere3Frustum3TI,
    defaultIntrSphere3Frustum3TIResult
} from '../src/IntrSphere3Frustum3.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';

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

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrSphere3Frustum3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { dot, length, sub } from '../src/Vector.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no frame component is a
// subnormal that would underflow when squared.
const angleF = () => wellScaled(-Math.PI, Math.PI);

// R = Rz(a)*Ry(b)*Rx(c); the columns are an orthonormal frame.
function rotFrameF(a: number, b: number, c: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const cc = Math.cos(c), sc = Math.sin(c);
    return [
        vec(ca * cb, sa * cb, -sb),
        vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc, cb * sc),
        vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc, cb * cc)];
}

// A point of the solid frustum at depth d in [dMin, dMax] and normalized
// bounds (u, r) in [-1, 1]. The frustum widens linearly with depth.
function frustumPoint(f: Frustum3, d: number, u: number, r: number): Vector {
    const s = d / f.dMin;
    return add(f.origin, add(mul(d, f.dVector),
        add(mul(u * s * f.uBound, f.uVector),
            mul(r * s * f.rBound, f.rVector))));
}

// Membership in the solid frustum, from the definition in Frustum3.h. The
// tolerance matters because both callers below feed it points that are on the
// boundary by construction (a sampled face point, or the closest frustum
// point of the distance query), where an exact comparison is a coin flip.
function inFrustum(f: Frustum3, p: Vector, eps = 0): boolean {
    const diff = sub(p, f.origin);
    const d = dot(diff, f.dVector);
    if (d < f.dMin - eps || d > f.dMax + eps) { return false; }
    const s = Math.max(d, f.dMin) / f.dMin;
    return Math.abs(dot(diff, f.uVector)) <= s * f.uBound + eps
        && Math.abs(dot(diff, f.rVector)) <= s * f.rBound + eps;
}

const frustumArb = fc.tuple(
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    angleF(), angleF(), angleF(),
    fc.double({ min: 0.3, max: 2, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 1.2, max: 5, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.2, max: 2, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.2, max: 2, noNaN: true, noDefaultInfinity: true })
).map(([ox, oy, oz, a, b, c, dMin, dSpan, u, r]) => {
    const fr = rotFrameF(a, b, c);
    return Frustum3.fromParameters(vec(ox, oy, oz), fr[0], fr[1], fr[2],
        dMin, dMin + dSpan, u, r);
});

const sphereFrustum = fc.tuple(frustumArb,
    wellScaled(-6, 6), wellScaled(-6, 6), wellScaled(-6, 6),
    fc.double({ min: 0, max: 3, noNaN: true, noDefaultInfinity: true }))
    .map(([f, cx, cy, cz, rad]) => ({
        frustum: f,
        sphere: Hypersphere.fromCenterRadius(vec(cx, cy, cz), rad)
    }));

describe('IntrSphere3Frustum3 verification', () => {
    const q = new IntrSphere3Frustum3TI();
    const dq = new DistPoint3Frustum3();

    it('agrees with the point-frustum distance it is defined by', () => {
        check(sphereFrustum, ({ frustum: f, sphere: s }) => {
            const d = dq.compute(s.center, f);
            expect(q.test(s, f).intersect).toBe(d.distance <= s.radius);
            // The closest point really is in the frustum and at the reported
            // distance from the sphere center.
            expect(inFrustum(f, d.closest[1], 1e-9)
                || Math.abs(d.distance) < 1e-9).toBe(true);
            expectClose(length(sub(d.closest[1], s.center)), d.distance,
                1e-8, 1e-9);
        });
    });

    it('a sphere centred on a sampled frustum point always intersects', () => {
        // The sample is kept strictly interior (a 5% margin on each of the
        // three frustum coordinates): a point exactly on a face or corner
        // with radius 0 is a knife edge, where the distance query can return
        // a few-ulp positive distance and the answer is genuinely ambiguous.
        check(fc.tuple(frustumArb,
            fc.double({ min: 0.05, max: 0.95, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -0.95, max: 0.95, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -0.95, max: 0.95, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: 0, max: 2, noNaN: true,
                noDefaultInfinity: true })),
        ([f, t, u, r, rad]) => {
            const d = f.dMin + t * (f.dMax - f.dMin);
            const p = frustumPoint(f, d, u, r);
            expect(inFrustum(f, p, 1e-9)).toBe(true);
            expect(q.test(sphere(p, rad), f).intersect).toBe(true);
        });
    });

    it('any sampled frustum point inside the sphere forces intersect', () => {
        const rnd = seededRandom(0x5f30ab1);
        let witnessed = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const fr = rotFrameF(rnd() * 6, rnd() * 6, rnd() * 6);
            const dMin = 0.3 + rnd() * 1.5;
            const f = Frustum3.fromParameters(
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                fr[0], fr[1], fr[2], dMin, dMin + 0.5 + rnd() * 3,
                0.2 + rnd() * 1.5, 0.2 + rnd() * 1.5);
            const s = Hypersphere.fromCenterRadius(
                vec(rnd() * 8 - 4, rnd() * 8 - 4, rnd() * 8 - 4),
                0.2 + rnd() * 2);
            let common = false;
            for (let k = 0; k < 300 && !common; ++k) {
                const p = frustumPoint(f, f.dMin + rnd() * (f.dMax - f.dMin),
                    2 * rnd() - 1, 2 * rnd() - 1);
                if (length(sub(p, s.center)) <= s.radius) { common = true; }
            }
            if (common) {
                expect(q.test(s, f).intersect).toBe(true);
                ++witnessed;
            }
        }
        expect(witnessed).toBeGreaterThan(20);
    }, 60000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(sphereFrustum, angleF(), angleF(), angleF(),
            wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3)),
        ([{ frustum: f, sphere: s }, a1, a2, a3, tx, ty, tz]) => {
            const fr = rotFrameF(a1, a2, a3);
            const rot = (v: Vector): Vector => add(mul(v.get(0), fr[0]),
                add(mul(v.get(1), fr[1]), mul(v.get(2), fr[2])));
            const t = vec(tx, ty, tz);
            const xf = (v: Vector): Vector => add(rot(v), t);
            const f2 = Frustum3.fromParameters(xf(f.origin), rot(f.dVector),
                rot(f.uVector), rot(f.rVector), f.dMin, f.dMax, f.uBound,
                f.rBound);
            const s2 = Hypersphere.fromCenterRadius(xf(s.center), s.radius);
            const a = q.test(s, f).intersect;
            const b = q.test(s2, f2).intersect;
            if (a === b) { return; }
            // Only a grazing sphere may flip: shrinking and growing the
            // radius by 1e-6 must change the original answer too.
            const d = dq.compute(s.center, f).distance;
            expect(Math.abs(d - s.radius)).toBeLessThan(1e-6);
        });
    });

    it('never intersects a sphere entirely behind the near plane', () => {
        check(fc.tuple(frustumArb,
            fc.double({ min: 0.05, max: 3, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: 0, max: 0.5, noNaN: true,
                noDefaultInfinity: true })),
        ([f, gap, rad]) => {
            // A center at depth dMin - gap - rad along -D from the near
            // plane, on the frustum axis: its whole ball is at depth below
            // dMin, so it cannot meet the solid frustum.
            const depth = f.dMin - gap - rad;
            const c = add(f.origin, mul(depth, f.dVector));
            const s = Hypersphere.fromCenterRadius(c, rad);
            expect(q.test(s, f).intersect).toBe(false);
        });
    });

    it('a sphere large enough to swallow the frustum always intersects',
        () => {
            check(sphereFrustum, ({ frustum: f, sphere: s }) => {
                let far = 0;
                for (const v of f.computeVertices()) {
                    const d = length(sub(v, s.center));
                    if (d > far) { far = d; }
                }
                const big = Hypersphere.fromCenterRadius(s.center, far + 1);
                expect(q.test(big, f).intersect).toBe(true);
            });
        });

    it('a zero-radius sphere is the frustum containment test', () => {
        check(fc.tuple(frustumArb, wellScaled(-6, 6), wellScaled(-6, 6),
            wellScaled(-6, 6)),
        ([f, cx, cy, cz]) => {
            const p = vec(cx, cy, cz);
            const inside = inFrustum(f, p);
            const d = dq.compute(p, f).distance;
            if (Math.abs(d) < 1e-9) { return; }   // on the boundary
            expect(q.test(sphere(p, 0), f).intersect).toBe(inside);
        });
    });
});
