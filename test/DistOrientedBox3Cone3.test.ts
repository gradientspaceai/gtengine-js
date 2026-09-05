import { describe, expect, it } from 'vitest';
import { Cone } from '../src/Cone.js';
import { DistOrientedBox3Cone3 } from '../src/DistOrientedBox3Cone3.js';
import type { DistOrientedBox3Cone3Result }
    from '../src/DistOrientedBox3Cone3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { seededRandom } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function obox(center: number[], axis: Vector[], extent: number[]):
    OrientedBox {
    return OrientedBox.fromCenterAxisExtent(v(...center), axis, v(...extent));
}

// An orthonormal frame parameterized by two angles.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        v(ca, sa, 0),
        v(-sa * cb, ca * cb, sb),
        v(sa * sb, -ca * sb, cb)
    ];
}

function frustum(vertex: number[], direction: number[], angle: number,
    hMin: number, hMax: number): Cone {
    const d = v(...direction);
    const ray = Ray.fromOriginDirection(v(...vertex), d);
    return Cone.fromRayAngleMinMaxHeight(ray, angle, hMin, hMax);
}

// The exact distance from a point to a solid oriented box.
function pointBoxDistance(p: Vector, box: OrientedBox): number {
    const delta = sub(p, box.center);
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const s = dot(box.axis[i], delta);
        const d = Math.abs(s) - box.extent.values[i];
        if (d > 0) {
            sqr += d * d;
        }
    }
    return Math.sqrt(sqr);
}

// A sampling of the boundary of the solid cone frustum. Because the distance
// from each sample to the box is exact, the result is an upper bound for the
// box-frustum distance.
function sampledDistance(box: OrientedBox, cone: Cone, nH: number,
    nR: number, nA: number): number {
    // A frame with the cone axis as its third vector.
    const d = cone.ray.direction;
    let w0 = Math.abs(d.values[0]) > Math.abs(d.values[1])
        ? v(-d.values[2], 0, d.values[0])
        : v(0, d.values[2], -d.values[1]);
    const len = Math.sqrt(dot(w0, w0));
    w0 = mul(1 / len, w0);
    const w1 = v(
        d.values[1] * w0.values[2] - d.values[2] * w0.values[1],
        d.values[2] * w0.values[0] - d.values[0] * w0.values[2],
        d.values[0] * w0.values[1] - d.values[1] * w0.values[0]);

    const hMin = cone.getMinHeight();
    const hMax = cone.getMaxHeight();
    let best = Number.MAX_VALUE;
    for (let ih = 0; ih <= nH; ++ih) {
        const h = hMin + ((hMax - hMin) * ih) / nH;
        const radius = h * cone.tanAngle;
        const axisPoint = add(cone.ray.origin, mul(h, d));
        for (let ir = 0; ir <= nR; ++ir) {
            const r = (radius * ir) / nR;
            for (let ia = 0; ia < nA; ++ia) {
                const t = (2 * Math.PI * ia) / nA;
                const p = add(axisPoint,
                    add(mul(r * Math.cos(t), w0), mul(r * Math.sin(t), w1)));
                best = Math.min(best, pointBoxDistance(p, box));
            }
        }
    }
    return best;
}

function expectConsistent(box: OrientedBox, cone: Cone,
    result: DistOrientedBox3Cone3Result): void {
    // The reported distance is the length of the segment joining the closest
    // points.
    const delta = sub(result.boxClosestPoint, result.coneClosestPoint);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 6);

    // The box closest point is in the box.
    const boxDelta = sub(result.boxClosestPoint, box.center);
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(dot(box.axis[i], boxDelta)))
            .toBeLessThanOrEqual(box.extent.values[i] + 1e-6);
    }

    // The cone closest point is in the frustum: its height is in
    // [hMin,hMax] and its radial distance is at most h * tan(angle).
    const coneDelta = sub(result.coneClosestPoint, cone.ray.origin);
    const h = dot(cone.ray.direction, coneDelta);
    expect(h).toBeGreaterThanOrEqual(cone.getMinHeight() - 1e-6);
    expect(h).toBeLessThanOrEqual(cone.getMaxHeight() + 1e-6);
    const radial = sub(coneDelta, mul(h, cone.ray.direction));
    expect(Math.sqrt(dot(radial, radial)))
        .toBeLessThanOrEqual(h * cone.tanAngle + 1e-6);
}

describe('DistOrientedBox3Cone3', () => {
    const query = new DistOrientedBox3Cone3();
    const angle = Math.PI / 6;
    const tan = Math.tan(angle);

    it('computes the distance to a box below the frustum cap', () => {
        // The frustum spans z in [1,3] with radius h*tan(30 degrees). A unit
        // box centered at (0,0,-5) has its top face at z = -4, and the
        // frustum cap at z = 1 lies over that face, so the gap is 5.
        const cone = frustum([0, 0, 0], [0, 0, 1], angle, 1, 3);
        const box = obox([0, 0, -5], frame(0, 0), [1, 1, 1]);
        const result = query.compute(box, cone);
        expect(result.distance).toBeCloseTo(5, 6);
        expect(result.boxClosestPoint.values[2]).toBeCloseTo(-4, 6);
        expect(result.coneClosestPoint.values[2]).toBeCloseTo(1, 6);
        expectConsistent(box, cone, result);
    });

    it('computes the distance to a box beside the frustum', () => {
        // The widest frustum point at the box height range z in [1,3] is
        // (3*tan(30), 0, 3), and the box face nearest it is x = 9.
        const cone = frustum([0, 0, 0], [0, 0, 1], angle, 1, 3);
        const box = obox([10, 0, 2], frame(0, 0), [1, 1, 1]);
        const result = query.compute(box, cone);
        expect(result.distance).toBeCloseTo(9 - 3 * tan, 5);
        expectConsistent(box, cone, result);
    });

    it('reports zero distance when the box contains the frustum', () => {
        const cone = frustum([0, 0, 0], [0, 0, 1], angle, 1, 3);
        const box = obox([0, 0, 2], frame(0, 0), [5, 5, 5]);
        const result = query.compute(box, cone);
        expect(result.distance).toBeLessThan(1e-6);
        expectConsistent(box, cone, result);
    });

    it('is invariant under rigid motion', () => {
        const c = Math.cos(0.8), s = Math.sin(0.8);
        const rot = (p: Vector) => v(
            c * p.values[0] - s * p.values[2],
            p.values[1],
            s * p.values[0] + c * p.values[2]);
        const shift = v(2, -3, 1);

        const cone = frustum([0, 0, 0], [0, 0, 1], angle, 1, 3);
        const box = obox([4, 1, 2], frame(0.3, 0.6), [0.7, 1.2, 0.5]);
        const r0 = query.compute(box, cone);

        const movedCone = Cone.fromRayAngleMinMaxHeight(
            Ray.fromOriginDirection(add(rot(cone.ray.origin), shift),
                rot(cone.ray.direction)),
            angle, 1, 3);
        const movedBox = OrientedBox.fromCenterAxisExtent(
            add(rot(box.center), shift),
            [rot(box.axis[0]), rot(box.axis[1]), rot(box.axis[2])],
            box.extent);
        const r1 = query.compute(movedBox, movedCone);
        expect(r1.distance).toBeCloseTo(r0.distance, 5);
    });

    it('accepts custom minimizer controls', () => {
        const cone = frustum([0, 0, 0], [0, 0, 1], angle, 1, 3);
        const box = obox([6, 0, 2], frame(0, 0), [1, 1, 1]);
        const coarse = query.compute(box, cone,
            { maxSubdivisions: 4, maxBisections: 32, epsilon: 1e-6,
                tolerance: 1e-3 });
        const fine = query.compute(box, cone);
        expect(coarse.distance).toBeCloseTo(fine.distance, 4);
    });

    it('rejects an infinite cone', () => {
        const cone = Cone.fromRayAngle(
            Ray.fromOriginDirection(v(0, 0, 0), v(0, 0, 1)), angle);
        const box = obox([6, 0, 2], frame(0, 0), [1, 1, 1]);
        expect(() => query.compute(box, cone)).toThrow();
    });

    it('matches a sampling of the frustum on assorted configurations', () => {
        let seed = 20261107;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const cone = frustum([0, 0, 0], [0, 0, 1], angle, 1, 4);
        let misses = 0;
        for (let trial = 0; trial < 10; ++trial) {
            const box = obox(
                [rand() * 10 - 5, rand() * 10 - 5, rand() * 10 - 3],
                frame(rand() * Math.PI, rand() * Math.PI),
                [rand() + 0.2, rand() + 0.2, rand() + 0.2]);
            const result = query.compute(box, cone);
            const sampled = sampledDistance(box, cone, 24, 6, 48);
            // The reported pair is always a valid pair of points, so its
            // distance is never below the true minimum by more than the
            // sampling error.
            expectConsistent(box, cone, result);
            if (Math.abs(result.distance - sampled) > 0.2) {
                ++misses;
            }
        }
        // The upstream minimizer misses the global minimum for one of these
        // configurations; see the limitation test below.
        expect(misses).toBeLessThanOrEqual(1);
    });

    it('records the known upstream minimizer limitation', () => {
        // This configuration defeats Minimize1: F(-pi/2) and F(+pi/2) are
        // equal to within round-off (both angles describe the same
        // quadrilateral), so the interpolating parabola has its vertex a few
        // ulps away from the midpoint of the bracket. Minimize1 then takes
        // its asymmetric branch, collapses the bracket and stops at the
        // midpoint value instead of the true minimum near angle 0.22. The
        // port preserves upstream behavior; update this test if Minimize1 is
        // ever made robust to a near-midpoint parabola vertex.
        const cone = frustum([0, 0, 0], [0, 0, 1], angle, 1, 4);
        const box = obox(
            [-1.3862013816833496, 2.672417163848877, 3.3069558143615723],
            frame(0.21659464178479224, 0.5082338445396278),
            [0.9043597102165222, 0.38772029876708985, 0.6854440987110137]);
        const result = query.compute(box, cone);
        const sampled = sampledDistance(box, cone, 48, 10, 120);

        // The returned pair is a genuine box point and frustum point.
        expectConsistent(box, cone, result);

        // But it is not the closest such pair.
        expect(sampled).toBeLessThan(0.1);
        expect(result.distance).toBeCloseTo(0.362774519113053, 6);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): cross-checks of the port against
// upstream DistOrientedBox3Cone3.h.
//
// This query is the least robust in the group and the properties below are
// chosen accordingly. Two upstream defects are involved:
//
//  * Minimize1 collapses its bracket when F(-pi/2) and F(+pi/2) agree to
//    within round-off, which they always do here (both angles describe the
//    same quadrilateral). The reported angle, and therefore the reported
//    distance, is then whatever the truncated search happened to reach. This
//    is recorded by the 'records the known upstream minimizer limitation'
//    test above and in the port notes.
//  * The 5-variable box-quadrilateral subproblem has a rank-deficient
//    Hessian (the Gram matrix of five vectors in R^3), so its 10-dimensional
//    LCP is degenerate. For a small fraction of configurations Lemke's
//    algorithm reports HAS_NONTRIVIAL_SOLUTION with a (w,z) pair that does
//    not satisfy w = q + M*z, and the query then returns points that lie
//    outside the box and outside the frustum. See the regression test at the
//    end of this block.
//
// Both defects make the query sensitive to tiny input perturbations, so
// rigid-motion and uniform-scaling equivariance are NOT properties of this
// query (they fail for roughly 10-20% of random configurations) and are not
// asserted here. What does hold unconditionally is the internal consistency
// of the returned triple, which is what these properties check, over a fixed
// seeded sample so that the outcome is deterministic.
// ---------------------------------------------------------------------------

describe('DistOrientedBox3Cone3 verification', () => {
    const query = new DistOrientedBox3Cone3();

    // A deterministic sample of (box, frustum) configurations.
    function sampleConfigurations(seed: number, count: number):
        Array<{ box: OrientedBox, cone: Cone }> {
        const rng = seededRandom(seed);
        const out: Array<{ box: OrientedBox, cone: Cone }> = [];
        for (let k = 0; k < count; ++k) {
            const dir = v(rng() - 0.5, rng() - 0.5, rng() - 0.5);
            normalize(dir);
            const h0 = 0.3 + 1.5 * rng();
            const h1 = h0 + 0.5 + 3 * rng();
            const cone = Cone.fromRayAngleMinMaxHeight(
                Ray.fromOriginDirection(
                    v(4 * rng() - 2, 4 * rng() - 2, 4 * rng() - 2), dir),
                0.15 + rng(), h0, h1);
            const box = OrientedBox.fromCenterAxisExtent(
                v(10 * rng() - 5, 10 * rng() - 5, 10 * rng() - 5),
                frame(2 * Math.PI * rng(), 2 * Math.PI * rng()),
                v(0.1 + 1.8 * rng(), 0.1 + 1.8 * rng(), 0.1 + 1.8 * rng()));
            out.push({ box, cone });
        }
        return out;
    }

    const sample = sampleConfigurations(0x51ce0001, 120);

    it('always reports a finite distance equal to the pair separation', () => {
        for (const { box, cone } of sample) {
            const r = query.compute(box, cone);
            expect(Number.isFinite(r.distance)).toBe(true);
            expect(r.distance).toBeGreaterThanOrEqual(0);
            const d = sub(r.boxClosestPoint, r.coneClosestPoint);
            expect(Math.sqrt(dot(d, d))).toBeCloseTo(r.distance, 8);
        }
    }, 30000);

    it('reports points on the box and on the frustum for the sample', () => {
        // Over wider random sampling roughly one configuration in 400 breaks
        // this because of the degenerate-LCP defect described above; the
        // fixed seed keeps this test deterministic and the known bad
        // configuration is pinned by the regression test below.
        for (const { box, cone } of sample) {
            const r = query.compute(box, cone);
            const bd = sub(r.boxClosestPoint, box.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(box.axis[i], bd)))
                    .toBeLessThanOrEqual(box.extent.values[i] + 1e-6);
            }
            const cd = sub(r.coneClosestPoint, cone.ray.origin);
            const h = dot(cone.ray.direction, cd);
            expect(h).toBeGreaterThanOrEqual(cone.getMinHeight() - 1e-6);
            expect(h).toBeLessThanOrEqual(cone.getMaxHeight() + 1e-6);
            const radial = sub(cd, mul(h, cone.ray.direction));
            expect(Math.sqrt(dot(radial, radial)))
                .toBeLessThanOrEqual(h * cone.tanAngle + 1e-6);
        }
    }, 30000);

    // The frustum is contained in the ball centered on its mid-height axis
    // point whose radius reaches both cap rims, so the box-to-ball distance
    // is a lower bound for the true minimum and therefore for any distance
    // the query can legitimately report.
    it('is never below the box-to-bounding-ball distance', () => {
        for (const { box, cone } of sample) {
            const h0 = cone.getMinHeight(), h1 = cone.getMaxHeight();
            const hc = 0.5 * (h0 + h1);
            const radiusSqr = Math.max(
                (h0 - hc) ** 2 + (h0 * cone.tanAngle) ** 2,
                (h1 - hc) ** 2 + (h1 * cone.tanAngle) ** 2);
            const center = add(cone.ray.origin, mul(hc, cone.ray.direction));
            const lower = pointBoxDistance(center, box) - Math.sqrt(radiusSqr);
            const r = query.compute(box, cone);
            expect(r.distance)
                .toBeGreaterThanOrEqual(Math.max(lower, 0) - 1e-6);
        }
    }, 30000);

    it('reports zero when the box swallows the whole frustum', () => {
        for (const { cone } of sample) {
            const h0 = cone.getMinHeight(), h1 = cone.getMaxHeight();
            const p0 = add(cone.ray.origin, mul(h0, cone.ray.direction));
            const p1 = add(cone.ray.origin, mul(h1, cone.ray.direction));
            const e = 2 * (h1 - h0) + 2 * h1 * cone.tanAngle;
            const box = OrientedBox.fromCenterAxisExtent(
                mul(0.5, add(p0, p1)),
                [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(e, e, e));
            expect(query.compute(box, cone).distance).toBeLessThan(1e-6);
        }
    }, 30000);

    // The LCP solver is a member reused across calls, so a result that
    // aliased its scratch state would change when another query ran.
    it('does not leak solver state between calls', () => {
        const cone = frustum([0, 0, 0], [0, 0, 1], Math.PI / 6, 1, 4);
        const box0 = obox([4, 1, 2], frame(0.3, 0.6), [0.7, 1.2, 0.5]);
        const box1 = obox([-3, 2, 5], frame(1.1, 0.2), [0.4, 0.9, 1.3]);
        const a = query.compute(box0, cone);
        const savedDistance = a.distance;
        const savedBox = a.boxClosestPoint.clone();
        const savedCone = a.coneClosestPoint.clone();
        query.compute(box1, cone);
        expect(a.distance).toBe(savedDistance);
        expect(a.boxClosestPoint.equals(savedBox)).toBe(true);
        expect(a.coneClosestPoint.equals(savedCone)).toBe(true);
        expect(query.compute(box0, cone).distance).toBe(savedDistance);
    });

    it('records the known upstream degenerate-LCP failure', () => {
        // For this configuration the LCP built for quadrilateral angle -pi/4
        // converges to w = 0 for the constraint row z[1] <= 2*extent[1] even
        // though q + M*z gives -0.71 there: the returned (w,z) is not a
        // solution of the LCP it was handed. DistOrientedBox3Cone3 believes
        // it and reports a pair whose "box point" is outside the box and
        // whose "frustum point" has height below hmin, with a distance below
        // the true minimum (about 0.538 for this configuration).
        //
        // The port preserves upstream behavior; update this test if the
        // upstream LCPSolver is ever made robust for degenerate problems.
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)],
            v(0.9824121509090933, 1.864152179336889, 0.8282953716840157));
        const cone = Cone.fromRayAngleMinMaxHeight(
            Ray.fromOriginDirection(v(0, 3, 0), v(0, 0, -1)),
            0.6731308200362891, 1.3644785852643415, 1.8644785852643424);
        const r = query.compute(box, cone);

        // The reported distance is still the separation of the reported pair.
        const d = sub(r.boxClosestPoint, r.coneClosestPoint);
        expect(Math.sqrt(dot(d, d))).toBeCloseTo(r.distance, 8);

        // But the "frustum point" is below the minimum height ...
        const cd = sub(r.coneClosestPoint, cone.ray.origin);
        const h = dot(cone.ray.direction, cd);
        expect(h).toBeCloseTo(1.1092690739354032, 9);
        expect(h).toBeLessThan(cone.getMinHeight());

        // ... and the "box point" is outside the box.
        expect(Math.abs(dot(box.axis[1], sub(r.boxClosestPoint, box.center))))
            .toBeGreaterThan(box.extent.values[1]);

        expect(r.distance).toBeCloseTo(0.5097726463860089, 9);
    });
});
