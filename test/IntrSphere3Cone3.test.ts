import { describe, expect, it } from 'vitest';
import {
    IntrSphere3Cone3FI,
    IntrSphere3Cone3TI,
    defaultIntrSphere3Cone3FIResult,
    defaultIntrSphere3Cone3TIResult
} from '../src/IntrSphere3Cone3';
import { Cone } from '../src/Cone';
import { Hypersphere } from '../src/Hypersphere';
import { Ray } from '../src/Ray';
import { Vector, dot, length, normalize, sub } from '../src/Vector';

const V3 = (x: number, y: number, z: number) => Vector.fromArray([x, y, z]);

function unit(x: number, y: number, z: number): Vector {
    const v = V3(x, y, z);
    normalize(v);
    return v;
}

function makeCone(vertex: Vector, axis: Vector, angle: number,
    hmin: number, hmax: number | null): Cone {
    const ray = Ray.fromOriginDirection(vertex, axis);
    if (hmax === null) {
        return hmin > 0
            ? Cone.fromRayAngleMinHeight(ray, angle, hmin)
            : Cone.fromRayAngle(ray, angle);
    }
    return Cone.fromRayAngleMinMaxHeight(ray, angle, hmin, hmax);
}

const sphere = (c: Vector, r: number) => Hypersphere.fromCenterRadius(c, r);

// Distance from a point to a convex polygon in the plane, computed
// independently of the library. Returns 0 when the point is inside. The
// polygon must be convex and counterclockwise ordered.
function distPointConvexPolygon2(px: number, py: number,
    poly: readonly (readonly [number, number])[]): number {
    let inside = true;
    const n = poly.length;
    for (let i = 0; i < n; ++i) {
        const [ax, ay] = poly[i];
        const [bx, by] = poly[(i + 1) % n];
        if ((bx - ax) * (py - ay) - (by - ay) * (px - ax) < 0) {
            inside = false;
            break;
        }
    }
    if (inside) {
        return 0;
    }
    let best = Number.MAX_VALUE;
    for (let i = 0; i < n; ++i) {
        const [ax, ay] = poly[i];
        const [bx, by] = poly[(i + 1) % n];
        const ex = bx - ax, ey = by - ay;
        const denom = ex * ex + ey * ey;
        let t = denom > 0 ? ((px - ax) * ex + (py - ay) * ey) / denom : 0;
        t = Math.min(Math.max(t, 0), 1);
        const dx = px - (ax + t * ex), dy = py - (ay + t * ey);
        best = Math.min(best, Math.sqrt(dx * dx + dy * dy));
    }
    return best;
}

// The exact distance from a point to a solid cone (of possibly truncated
// height range). The solid cone is rotationally symmetric about its axis, so
// the 3D distance equals the 2D distance in the (h,s) half-plane, where h is
// the height along the axis and s >= 0 is the distance from the axis. The
// solid cone becomes the convex quadrilateral with corners (hmin,0),
// (hmax,0), (hmax,hmax*tan(A)), (hmin,hmin*tan(A)). This is an independent
// reference computation for the test-intersection query.
function distPointSolidCone(point: Vector, vertex: Vector, axis: Vector,
    angle: number, hmin: number, hmax: number): number {
    const diff = sub(point, vertex);
    const h = dot(diff, axis);
    const s = Math.sqrt(Math.max(dot(diff, diff) - h * h, 0));
    const tanA = Math.tan(angle);
    const poly: [number, number][] = [
        [hmin, 0], [hmax, 0], [hmax, hmax * tanA], [hmin, hmin * tanA]
    ];
    // Drop the degenerate corner for a non-truncated cone so that the polygon
    // stays a simple counterclockwise convex shape.
    const simple = hmin === 0 ? poly.slice(1) : poly;
    return distPointConvexPolygon2(h, s, simple);
}

describe('IntrSphere3Cone3 default results', () => {
    it('match the upstream default constructors', () => {
        expect(defaultIntrSphere3Cone3TIResult()).toEqual({ intersect: false });
        const fi = defaultIntrSphere3Cone3FIResult();
        expect(fi.intersect).toBe(false);
        expect(Array.from(fi.point.values)).toEqual([0, 0, 0]);
    });
});

describe('IntrSphere3Cone3TI known configurations', () => {
    const ti = new IntrSphere3Cone3TI();
    const axis = V3(0, 0, 1);
    const angle = Math.PI / 4;

    it('reports intersection when the sphere contains the cone vertex', () => {
        const cone = makeCone(V3(0, 0, 0), axis, angle, 0, null);
        expect(ti.test(sphere(V3(0, 0, -0.5), 1), cone).intersect).toBe(true);
    });

    it('rejects a sphere entirely behind the infinite cone vertex', () => {
        const cone = makeCone(V3(0, 0, 0), axis, angle, 0, null);
        expect(ti.test(sphere(V3(0, 0, -5), 1), cone).intersect).toBe(false);
    });

    it('reports intersection for a sphere on the cone axis', () => {
        const cone = makeCone(V3(0, 0, 0), axis, angle, 0, null);
        expect(ti.test(sphere(V3(0, 0, 4), 0.5), cone).intersect).toBe(true);
    });

    it('rejects a sphere lying outside the lateral surface', () => {
        // With a 45-degree half angle the cone is { s <= h }. The point
        // (10,0,2) is at signed distance (10-2)/sqrt(2) from the surface.
        const cone = makeCone(V3(0, 0, 0), axis, angle, 0, null);
        expect(ti.test(sphere(V3(10, 0, 2), 1), cone).intersect).toBe(false);
        // A radius large enough to reach the surface does intersect.
        expect(ti.test(sphere(V3(10, 0, 2), 8 / Math.SQRT2 + 1e-6), cone)
            .intersect).toBe(true);
    });

    it('handles a sphere just touching the lateral surface', () => {
        const cone = makeCone(V3(0, 0, 0), axis, angle, 0, null);
        const touching = 8 / Math.SQRT2;
        expect(ti.test(sphere(V3(10, 0, 2), touching * (1 + 1e-9)), cone)
            .intersect).toBe(true);
        expect(ti.test(sphere(V3(10, 0, 2), touching * (1 - 1e-9)), cone)
            .intersect).toBe(false);
    });

    it('respects the maximum height of a finite cone', () => {
        const cone = makeCone(V3(0, 0, 0), axis, angle, 0, 3);
        // On the axis beyond hmax by more than the radius.
        expect(ti.test(sphere(V3(0, 0, 5), 1), cone).intersect).toBe(false);
        expect(ti.test(sphere(V3(0, 0, 3.5), 1), cone).intersect).toBe(true);
    });

    it('respects the minimum height of a truncated cone', () => {
        const cone = makeCone(V3(0, 0, 0), axis, angle, 2, null);
        // Inside the untruncated cone but below hmin by more than the radius.
        expect(ti.test(sphere(V3(0, 0, 0.5), 0.25), cone).intersect).toBe(false);
        expect(ti.test(sphere(V3(0, 0, 1.5), 0.75), cone).intersect).toBe(true);
    });

    it('respects both heights of a cone frustum', () => {
        const cone = makeCone(V3(0, 0, 0), axis, angle, 2, 4);
        expect(ti.test(sphere(V3(0, 0, 3), 0.1), cone).intersect).toBe(true);
        expect(ti.test(sphere(V3(0, 0, 6), 1), cone).intersect).toBe(false);
        expect(ti.test(sphere(V3(0, 0, 0), 1), cone).intersect).toBe(false);
        // Near the frustum's outer rim at h = 4, s = 4.
        expect(ti.test(sphere(V3(4.5, 0, 4), 0.6), cone).intersect).toBe(true);
        expect(ti.test(sphere(V3(4.5, 0, 4), 0.4), cone).intersect).toBe(false);
    });

    it('is translation and rotation invariant', () => {
        const c0 = makeCone(V3(0, 0, 0), axis, angle, 1, 5);
        const s0 = sphere(V3(1.5, 0, 3), 0.75);
        const base = ti.test(s0, c0).intersect;
        // Rotate the whole configuration by 90 degrees about the y-axis and
        // translate; (x,y,z) -> (z,y,-x) + t.
        const t = V3(7, -3, 2);
        const map = (v: Vector) => V3(
            v.values[2] + t.values[0],
            v.values[1] + t.values[1],
            -v.values[0] + t.values[2]);
        const c1 = makeCone(map(V3(0, 0, 0)), V3(1, 0, 0), angle, 1, 5);
        const s1 = sphere(map(s0.center), s0.radius);
        expect(ti.test(s1, c1).intersect).toBe(base);
        expect(base).toBe(true);
    });
});

describe('IntrSphere3Cone3TI randomized cross-check', () => {
    it('agrees with an independent point-to-solid-cone distance', () => {
        const ti = new IntrSphere3Cone3TI();
        let seed = 987654321;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        let tested = 0, hits = 0;
        for (let k = 0; k < 3000; ++k) {
            const angle = 0.2 + 1.1 * rand();
            const hmin = k % 3 === 0 ? 0 : 2 * rand();
            const hmax = hmin + 0.5 + 4 * rand();
            const vertex = V3((rand() - 0.5) * 4, (rand() - 0.5) * 4,
                (rand() - 0.5) * 4);
            const axis = unit(rand() - 0.5, rand() - 0.5, rand() - 0.5);
            if (!Number.isFinite(axis.values[0])) {
                continue;
            }
            const cone = makeCone(vertex, axis, angle, hmin, hmax);
            const center = V3((rand() - 0.5) * 12, (rand() - 0.5) * 12,
                (rand() - 0.5) * 12);
            const radius = 0.1 + 3 * rand();
            const d = distPointSolidCone(center, vertex, axis, angle, hmin,
                hmax);
            // Skip near-tangential cases where floating-point round-off in
            // either computation can flip the answer.
            if (Math.abs(d - radius) < 1e-6) {
                continue;
            }
            const expected = d < radius;
            expect(ti.test(sphere(center, radius), cone).intersect)
                .toBe(expected);
            ++tested;
            if (expected) {
                ++hits;
            }
        }
        // Sanity: the sample covers both outcomes.
        expect(tested).toBeGreaterThan(2000);
        expect(hits).toBeGreaterThan(100);
        expect(hits).toBeLessThan(tested - 100);
    });

    it('agrees with the reference for infinite and truncated infinite cones',
        () => {
            const ti = new IntrSphere3Cone3TI();
            let seed = 24680;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) & 0x7fffffff;
                return seed / 0x7fffffff;
            };
            // A "large enough" height stands in for +infinity: every sphere
            // in the sample lies well within this height range.
            const bigH = 1e4;
            let tested = 0, hits = 0;
            for (let k = 0; k < 2000; ++k) {
                const angle = 0.2 + 1.1 * rand();
                const hmin = k % 2 === 0 ? 0 : 2 * rand();
                const vertex = V3((rand() - 0.5) * 4, (rand() - 0.5) * 4,
                    (rand() - 0.5) * 4);
                const axis = unit(rand() - 0.5, rand() - 0.5, rand() - 0.5);
                if (!Number.isFinite(axis.values[0])) {
                    continue;
                }
                const cone = makeCone(vertex, axis, angle, hmin, null);
                const center = V3((rand() - 0.5) * 12, (rand() - 0.5) * 12,
                    (rand() - 0.5) * 12);
                const radius = 0.1 + 3 * rand();
                const d = distPointSolidCone(center, vertex, axis, angle, hmin,
                    bigH);
                if (Math.abs(d - radius) < 1e-6) {
                    continue;
                }
                const expected = d < radius;
                expect(ti.test(sphere(center, radius), cone).intersect)
                    .toBe(expected);
                ++tested;
                if (expected) {
                    ++hits;
                }
            }
            expect(tested).toBeGreaterThan(1500);
            expect(hits).toBeGreaterThan(100);
            expect(hits).toBeLessThan(tested - 100);
        });
});

describe('IntrSphere3Cone3FI', () => {
    const fi = new IntrSphere3Cone3FI();
    const angle = Math.PI / 6;

    it('returns the cone vertex when the sphere contains it', () => {
        const cone = makeCone(V3(1, 2, 3), V3(0, 0, 1), angle, 0, null);
        const result = fi.find(sphere(V3(1, 2, 3.25), 1), cone);
        expect(result.intersect).toBe(true);
        expect(Array.from(result.point.values)).toEqual([1, 2, 3]);
    });

    it('returns the sphere center when it is inside the cone', () => {
        const cone = makeCone(V3(0, 0, 0), V3(0, 0, 1), angle, 0, null);
        const center = V3(0.1, 0, 5);
        const result = fi.find(sphere(center, 0.25), cone);
        expect(result.intersect).toBe(true);
        expect(Array.from(result.point.values)).toEqual([0.1, 0, 5]);
    });

    it('reports no intersection for a separated sphere', () => {
        const cone = makeCone(V3(0, 0, 0), V3(0, 0, 1), angle, 0, null);
        const result = fi.find(sphere(V3(20, 0, 1), 1), cone);
        expect(result.intersect).toBe(false);
    });

    it('reports no intersection for a sphere behind the vertex', () => {
        const cone = makeCone(V3(0, 0, 0), V3(0, 0, 1), angle, 0, null);
        expect(fi.find(sphere(V3(0, 0, -10), 1), cone).intersect).toBe(false);
    });

    it('returns a point on both the sphere and the cone surface', () => {
        // The cone vertex is NOT at the origin, which exercises the upstream
        // fix: upstream returns t * D instead of V + t * D.
        const vertex = V3(3, -1, 2);
        const axisDir = unit(1, 1, 1);
        const cone = makeCone(vertex, axisDir, angle, 0, null);
        // A sphere center outside the cone but close enough to reach it.
        const center = V3(3, -1, 8);
        const radius = 3.0;
        const result = fi.find(sphere(center, radius), cone);
        expect(result.intersect).toBe(true);

        const p = result.point;
        // On the sphere.
        expect(length(sub(p, center))).toBeCloseTo(radius, 8);
        // On the cone: Dot(D, P - V) = |P - V| * cos(angle) with a
        // nonnegative height.
        const pmv = sub(p, vertex);
        const h = dot(axisDir, pmv);
        expect(h).toBeGreaterThan(0);
        expect(h).toBeCloseTo(length(pmv) * Math.cos(angle), 9);
    });

    it('produces points on both surfaces over randomized configurations',
        () => {
            let seed = 777;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) & 0x7fffffff;
                return seed / 0x7fffffff;
            };
            let found = 0, missed = 0;
            for (let k = 0; k < 1500; ++k) {
                const angle2 = 0.2 + 1.0 * rand();
                const vertex = V3((rand() - 0.5) * 6, (rand() - 0.5) * 6,
                    (rand() - 0.5) * 6);
                const axisDir = unit(rand() - 0.5, rand() - 0.5, rand() - 0.5);
                if (!Number.isFinite(axisDir.values[0])) {
                    continue;
                }
                const cone = makeCone(vertex, axisDir, angle2, 0, null);
                const center = V3((rand() - 0.5) * 12, (rand() - 0.5) * 12,
                    (rand() - 0.5) * 12);
                const radius = 0.2 + 3 * rand();
                const result = fi.find(sphere(center, radius), cone);

                // The FI query uses the infinite cone, so cross-check the
                // outcome against the independent solid-cone distance.
                const d = distPointSolidCone(center, vertex, axisDir, angle2,
                    0, 1e4);
                if (Math.abs(d - radius) > 1e-6) {
                    expect(result.intersect).toBe(d < radius);
                }

                if (!result.intersect) {
                    ++missed;
                    continue;
                }
                ++found;

                const p = result.point;
                // The returned point must be inside (or on) both solids.
                expect(length(sub(p, center))).toBeLessThanOrEqual(
                    radius * (1 + 1e-8) + 1e-8);
                const pmv = sub(p, vertex);
                const h = dot(axisDir, pmv);
                expect(h).toBeGreaterThanOrEqual(-1e-8);
                const lp = length(pmv);
                expect(h).toBeGreaterThanOrEqual(
                    lp * Math.cos(angle2) - 1e-7 * Math.max(1, lp));
            }
            expect(found).toBeGreaterThan(100);
            expect(missed).toBeGreaterThan(100);
        });
});
