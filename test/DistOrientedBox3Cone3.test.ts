import { describe, expect, it } from 'vitest';
import { Cone } from '../src/Cone';
import { DistOrientedBox3Cone3 } from '../src/DistOrientedBox3Cone3';
import type { DistOrientedBox3Cone3Result }
    from '../src/DistOrientedBox3Cone3';
import { OrientedBox } from '../src/OrientedBox';
import { Ray } from '../src/Ray';
import { Vector, add, dot, mul, sub } from '../src/Vector';

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
