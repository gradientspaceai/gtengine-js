import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Cone } from '../src/Cone.js';
import { Ray } from '../src/Ray.js';
import { Vector, dot, normalize, sub } from '../src/Vector.js';
import {
    IntrAlignedBox3Cone3TI,
    defaultIntrAlignedBox3Cone3TIResult,
    intrAlignedBox3Cone3ComputeBoxHeightInterval,
    intrAlignedBox3Cone3ConeAxisIntersectsBox,
    intrAlignedBox3Cone3HasPointInsideCone
} from '../src/IntrAlignedBox3Cone3.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function box(x0: number, y0: number, z0: number, x1: number, y1: number,
    z1: number): AlignedBox {
    return AlignedBox.fromMinMax(v3(x0, y0, z0), v3(x1, y1, z1));
}

function makeRay(origin: Vector, direction: Vector): Ray {
    const ray = new Ray(3);
    ray.origin = origin.clone();
    ray.direction = direction.clone();
    normalize(ray.direction);
    return ray;
}

// A finite cone (frustum) with the given vertex, axis, angle and heights.
function frustum(vertex: Vector, axis: Vector, angle: number, hmin: number,
    hmax: number): Cone {
    return Cone.fromRayAngleMinMaxHeight(makeRay(vertex, axis), angle, hmin,
        hmax);
}

// An infinite cone with hmin = 0.
function infiniteCone(vertex: Vector, axis: Vector, angle: number): Cone {
    return Cone.fromRayAngle(makeRay(vertex, axis), angle);
}

// The "insideness" of the point X relative to the cone: the minimum of the
// three constraints that define the solid cone. The value is positive exactly
// when X is strictly inside the cone.
function insideness(cone: Cone, X: Vector): number {
    const d = sub(X, cone.ray.origin);
    const h = dot(cone.ray.direction, d);
    const radial = h - cone.cosAngle * Math.sqrt(dot(d, d));
    const above = h - cone.getMinHeight();
    const below = cone.isFinite() ? cone.getMaxHeight() - h : Number.MAX_VALUE;
    return Math.min(radial, above, below);
}

// The maximum insideness over a regular grid of samples in the box.
function maxInsideness(cone: Cone, b: AlignedBox, n: number): number {
    let best = -Number.MAX_VALUE;
    const lo = b.min.values;
    const hi = b.max.values;
    const p = v3(0, 0, 0);
    for (let i = 0; i <= n; ++i) {
        p.values[0] = lo[0] + (hi[0] - lo[0]) * i / n;
        for (let j = 0; j <= n; ++j) {
            p.values[1] = lo[1] + (hi[1] - lo[1]) * j / n;
            for (let k = 0; k <= n; ++k) {
                p.values[2] = lo[2] + (hi[2] - lo[2]) * k / n;
                const value = insideness(cone, p);
                if (value > best) {
                    best = value;
                }
            }
        }
    }
    return best;
}

describe('IntrAlignedBox3Cone3', () => {
    it('default result reports no intersection', () => {
        expect(defaultIntrAlignedBox3Cone3TIResult().intersect).toBe(false);
    });

    it('computes the box height interval relative to the cone axis', () => {
        const cone = infiniteCone(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 6);
        const interval = intrAlignedBox3Cone3ComputeBoxHeightInterval(
            box(-1, -1, 2, 1, 1, 5), cone);
        expect(interval.boxMinHeight).toBeCloseTo(2, 12);
        expect(interval.boxMaxHeight).toBeCloseTo(5, 12);

        // A diagonal axis projects the box corners onto the axis.
        const diagonal = infiniteCone(v3(0, 0, 0), v3(1, 1, 1), Math.PI / 6);
        const d = intrAlignedBox3Cone3ComputeBoxHeightInterval(
            box(0, 0, 0, 1, 1, 1), diagonal);
        expect(d.boxMinHeight).toBeCloseTo(0, 12);
        expect(d.boxMaxHeight).toBeCloseTo(Math.sqrt(3), 12);
    });

    it('detects when the cone axis meets the box', () => {
        const cone = infiniteCone(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 6);
        expect(intrAlignedBox3Cone3ConeAxisIntersectsBox(
            box(-1, -1, 3, 1, 1, 4), cone)).toBe(true);
        // The axis passes beside the box.
        expect(intrAlignedBox3Cone3ConeAxisIntersectsBox(
            box(5, 5, 3, 6, 6, 4), cone)).toBe(false);
        // The axis of a frustum stops before the box.
        const stub = frustum(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 6, 0, 1);
        expect(intrAlignedBox3Cone3ConeAxisIntersectsBox(
            box(-1, -1, 3, 1, 1, 4), stub)).toBe(false);
    });

    it('tests a segment for points strictly inside the cone', () => {
        // A cone with vertex at the origin, axis +z and half-angle pi/4. The
        // arguments are relative to the cone vertex.
        const cone = infiniteCone(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 4);
        // An endpoint inside.
        expect(intrAlignedBox3Cone3HasPointInsideCone(
            v3(0, 0, 1), v3(10, 0, 0), cone)).toBe(true);
        // Both endpoints outside but the segment interior passes inside.
        expect(intrAlignedBox3Cone3HasPointInsideCone(
            v3(-2, 0, 1.5), v3(2, 0, 1.5), cone)).toBe(true);
        // A segment entirely outside the cone.
        expect(intrAlignedBox3Cone3HasPointInsideCone(
            v3(-2, 0, -1), v3(2, 0, -1), cone)).toBe(false);
        // A segment on the cone boundary is not strictly inside.
        expect(intrAlignedBox3Cone3HasPointInsideCone(
            v3(0, 0, 0), v3(1, 0, 1), cone)).toBe(false);
    });

    it('accepts a box fully inside an infinite cone', () => {
        const query = new IntrAlignedBox3Cone3TI();
        const cone = infiniteCone(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 4);
        expect(query.test(box(-0.1, -0.1, 5, 0.1, 0.1, 6), cone).intersect)
            .toBe(true);
    });

    it('rejects a box that is separated from the cone', () => {
        const query = new IntrAlignedBox3Cone3TI();
        const cone = infiniteCone(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 6);
        // Below the cone vertex.
        expect(query.test(box(-1, -1, -5, 1, 1, -3), cone).intersect)
            .toBe(false);
        // Beside the cone, far from the axis.
        expect(query.test(box(20, 20, 0, 21, 21, 1), cone).intersect)
            .toBe(false);
    });

    it('reports no intersection for a box that only touches the cone', () => {
        const query = new IntrAlignedBox3Cone3TI();
        const cone = infiniteCone(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 4);
        // The box lies below the plane of the cone vertex and touches only
        // the vertex; the intersection has zero volume.
        expect(query.test(box(-1, -1, -1, 1, 1, 0), cone).intersect)
            .toBe(false);
        // The box lies above the maximum-height disk of a frustum and touches
        // it in a polygon of zero volume.
        const fin = frustum(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 4, 0, 2);
        expect(query.test(box(-1, -1, 2, 1, 1, 3), fin).intersect).toBe(false);
    });

    it('honors the minimum and maximum heights of a frustum', () => {
        const query = new IntrAlignedBox3Cone3TI();
        const cone = frustum(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 4, 2, 4);
        // Inside the height slab and near the axis.
        expect(query.test(box(-0.2, -0.2, 2.5, 0.2, 0.2, 3), cone).intersect)
            .toBe(true);
        // Near the axis but below the minimum height.
        expect(query.test(box(-0.2, -0.2, 0.5, 0.2, 0.2, 1), cone).intersect)
            .toBe(false);
        // Near the axis but above the maximum height.
        expect(query.test(box(-0.2, -0.2, 5, 0.2, 0.2, 6), cone).intersect)
            .toBe(false);
        // Straddling the minimum-height plane.
        expect(query.test(box(-0.2, -0.2, 1.5, 0.2, 0.2, 2.5), cone).intersect)
            .toBe(true);
        // An infinite truncated cone has no maximum height.
        const truncated = Cone.fromRayAngleMinHeight(
            makeRay(v3(0, 0, 0), v3(0, 0, 1)), Math.PI / 4, 2);
        expect(query.test(box(-0.2, -0.2, 500, 0.2, 0.2, 501), truncated)
            .intersect).toBe(true);
    });

    it('detects a box that meets only the flared side of the cone', () => {
        const query = new IntrAlignedBox3Cone3TI();
        // A wide cone with axis +z. The box sits high up and far off the
        // axis, so the axis misses the box but the cone surface does not.
        const cone = infiniteCone(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 3);
        const b = box(4, -0.5, 3, 6, 0.5, 4);
        expect(maxInsideness(cone, b, 12)).toBeGreaterThan(0);
        expect(query.test(b, cone).intersect).toBe(true);
    });

    it('is invariant under a translation of the box and the cone', () => {
        const query = new IntrAlignedBox3Cone3TI();
        const cone = frustum(v3(0, 0, 0), v3(1, 2, 3), 0.5, 1, 6);
        const b = box(0.5, 0.5, 1, 1.5, 2, 3);
        const base = query.test(b, cone).intersect;

        const t = v3(-7, 13, 2.5);
        const movedCone = frustum(
            v3(t.values[0], t.values[1], t.values[2]), v3(1, 2, 3), 0.5, 1, 6);
        const movedBox = AlignedBox.fromMinMax(
            v3(0.5 + t.values[0], 0.5 + t.values[1], 1 + t.values[2]),
            v3(1.5 + t.values[0], 2 + t.values[1], 3 + t.values[2]));
        expect(query.test(movedBox, movedCone).intersect).toBe(base);
    });

    it('gives the same answers when the query object is reused', () => {
        // The query caches per-call state (the clipped-box graph). Running a
        // sequence of queries on one object must agree with running each on a
        // fresh object.
        const cone = frustum(v3(0, 0, 0), v3(0, 0, 1), Math.PI / 4, 1, 3);
        const boxes = [
            box(-0.2, -0.2, 1.5, 0.2, 0.2, 2),   // fully inside the slab
            box(-3, -3, -3, 3, 3, 3),            // straddles both planes
            box(-0.2, -0.2, 1.2, 0.2, 0.2, 2.8), // fully inside the slab
            box(5, 5, 0, 6, 6, 10),              // outside the cone
            box(-1, -1, 0.5, 1, 1, 1.5),         // straddles the hmin plane
            box(-0.1, -0.1, 1.1, 0.1, 0.1, 2.9)  // fully inside the slab
        ];
        const shared = new IntrAlignedBox3Cone3TI();
        for (let pass = 0; pass < 3; ++pass) {
            for (const b of boxes) {
                const fresh = new IntrAlignedBox3Cone3TI().test(b, cone);
                expect(shared.test(b, cone).intersect).toBe(fresh.intersect);
            }
        }

        // The same check over randomized configurations. This is the
        // regression test for the upstream state-leak bug described in
        // boxFullyInConeSlab: with the upstream logic this loop reports
        // dozens of disagreements.
        let seed = 12345;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const reused = new IntrAlignedBox3Cone3TI();
        for (let trial = 0; trial < 4000; ++trial) {
            const trialCone = frustum(
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2),
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1),
                0.2 + 1.1 * rand(), rand(), 1 + 4 * rand());
            const c = v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3);
            const e = v3(0.2 + rand(), 0.2 + rand(), 0.2 + rand());
            const b = AlignedBox.fromMinMax(
                v3(c.values[0] - e.values[0], c.values[1] - e.values[1],
                    c.values[2] - e.values[2]),
                v3(c.values[0] + e.values[0], c.values[1] + e.values[1],
                    c.values[2] + e.values[2]));
            const fresh = new IntrAlignedBox3Cone3TI().test(b, trialCone);
            expect(reused.test(b, trialCone).intersect).toBe(fresh.intersect);
        }
    });

    it('agrees with dense sampling for randomized configurations', () => {
        // A deterministic linear congruential generator keeps the test
        // reproducible.
        let seed = 20260901;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const query = new IntrAlignedBox3Cone3TI();
        let numTrue = 0, numFalse = 0, numChecked = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const vertex = v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2);
            const axis = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (Math.sqrt(dot(axis, axis)) < 1e-3) {
                continue;
            }
            const angle = 0.2 + 1.1 * rand();
            const hmin = 2 * rand();
            const cone = rand() < 0.5
                ? infiniteCone(vertex, axis, angle)
                : frustum(vertex, axis, angle, hmin, hmin + 0.5 + 4 * rand());

            const c = v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3);
            const e = v3(0.2 + rand(), 0.2 + rand(), 0.2 + rand());
            const b = AlignedBox.fromMinMax(
                v3(c.values[0] - e.values[0], c.values[1] - e.values[1],
                    c.values[2] - e.values[2]),
                v3(c.values[0] + e.values[0], c.values[1] + e.values[1],
                    c.values[2] + e.values[2]));

            const best = maxInsideness(cone, b, 16);
            const result = query.test(b, cone).intersect;
            // A sampled point strictly inside the cone proves an overlap of
            // positive volume, so the query must report an intersection. The
            // grid maximum underestimates the true maximum, so only clearly
            // separated configurations are checked in the other direction.
            if (best > 1e-3) {
                expect(result).toBe(true);
                ++numTrue;
                ++numChecked;
            }
            else if (best < -0.2) {
                expect(result).toBe(false);
                ++numFalse;
                ++numChecked;
            }
        }
        // The sampling must exercise both outcomes.
        expect(numTrue).toBeGreaterThan(20);
        expect(numFalse).toBeGreaterThan(20);
        expect(numChecked).toBeGreaterThan(100);
    });
});
