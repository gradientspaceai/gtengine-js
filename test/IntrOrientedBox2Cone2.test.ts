import { describe, it, expect } from 'vitest';
import { Cone } from '../src/Cone.js';
import { IntrRay2OrientedBox2TI } from '../src/IntrRay2OrientedBox2.js';
import {
    IntrOrientedBox2Cone2TI,
    defaultIntrOrientedBox2Cone2TIResult
} from '../src/IntrOrientedBox2Cone2.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import {
    check, fc, orientedBox as arbOrientedBox, positive, rotationFrame,
    seededRandom, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function obox(center: number[], angle: number, extent: number[]): OrientedBox {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center),
        [vec(c, s), vec(-s, c)], Vector.fromArray(extent));
}

function cone2(origin: number[], direction: number[], angle: number): Cone {
    const c = new Cone(2);
    const d = Vector.fromArray(direction);
    normalize(d);
    c.ray = Ray.fromOriginDirection(Vector.fromArray(origin), d);
    c.setAngle(angle);
    return c;
}

// An independent evaluation of the documented criterion in world coordinates:
// the cone axis ray meets the box, or some box vertex P satisfies
// Dot(D,P-V) > 0 and (Dot(D,P-V))^2 > |P-V|^2 * cosAngle^2.
function reference(box: OrientedBox, cone: Cone): boolean {
    const rbQuery = new IntrRay2OrientedBox2TI();
    if (rbQuery.test(cone.ray, box).intersect) {
        return true;
    }
    for (const P of box.getVertices()) {
        const diff = sub(P, cone.ray.origin);
        const num = dot(cone.ray.direction, diff);
        if (num > 0) {
            if (num * num > dot(diff, diff) * cone.cosAngle * cone.cosAngle) {
                return true;
            }
        }
    }
    return false;
}

const ti = new IntrOrientedBox2Cone2TI();

describe('IntrOrientedBox2Cone2', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrOrientedBox2Cone2TIResult().intersect).toBe(false);
    });

    it('detects a box straddling the cone axis', () => {
        const C = cone2([0, 0], [1, 0], Math.PI / 6);
        expect(ti.test(obox([5, 0], 0, [1, 1]), C).intersect).toBe(true);
    });

    it('rejects a box behind the cone apex', () => {
        const C = cone2([0, 0], [1, 0], Math.PI / 6);
        expect(ti.test(obox([-5, 0], 0, [1, 1]), C).intersect).toBe(false);
    });

    it('rejects a box outside the cone angle and accepts one inside', () => {
        // A 30-degree half-angle cone along +x. At x = 5 the cone half-width
        // is 5*tan(30 deg) = 2.887.
        const C = cone2([0, 0], [1, 0], Math.PI / 6);
        // A small box centered well above the cone.
        expect(ti.test(obox([5, 6], 0, [0.5, 0.5]), C).intersect).toBe(false);
        // A small box centered inside the cone.
        expect(ti.test(obox([5, 2], 0, [0.2, 0.2]), C).intersect).toBe(true);
    });

    it('does not report a box that only touches the cone boundary', () => {
        // The 45-degree cone boundary through the origin is the line y = x.
        // A box whose only cone-side corner is exactly on y = x touches but
        // does not overlap.
        const C = cone2([0, 0], [1, 0], Math.PI / 4);
        const B = obox([3, 5], 0, [1, 1]);  // corner (4, 4) is on y = x
        expect(ti.test(B, C).intersect).toBe(false);
        // Nudging the box towards the axis produces an overlap.
        expect(ti.test(obox([3, 4.99], 0, [1, 1]), C).intersect).toBe(true);
    });

    it('accounts for the box orientation', () => {
        const C = cone2([0, 0], [1, 0], Math.PI / 12);  // 15 degrees
        // An axis-aligned square well above the narrow cone misses it.
        expect(ti.test(obox([10, 4], 0, [1, 1]), C).intersect).toBe(false);
        // Rotating it 45 degrees pushes a corner down to y = 4 - sqrt(2),
        // which is inside the cone (10*tan(15 deg) = 2.679).
        expect(ti.test(obox([10, 4], Math.PI / 4, [1, 1]), C).intersect)
            .toBe(true);
    });

    it('agrees with a world-coordinate evaluation on random inputs', () => {
        let state = 30313;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const C = cone2([rand() * 2, rand() * 2],
                [rand(), rand() + 0.001],
                0.15 + Math.abs(rand()) * 1.2);
            const B = obox([rand() * 6, rand() * 6], rand() * Math.PI,
                [0.3 + Math.abs(rand()) * 2, 0.3 + Math.abs(rand()) * 2]);
            const actual = ti.test(B, C).intersect;
            expect(actual).toBe(reference(B, C));
            if (actual) {
                ++numHits;
            }
        }
        expect(numHits).toBeGreaterThan(50);
        expect(numHits).toBeLessThan(450);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrOrientedBox2Cone2.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrOrientedBox2Cone2 verification', () => {
    const q = new IntrOrientedBox2Cone2TI();

    type Pt = [number, number];

    // Sutherland-Hodgman clip of a convex polygon by the half-plane
    // Dot(n, X - o) >= 0.
    function clipHalfPlane(poly: Pt[], n: Pt, o: Pt): Pt[] {
        const out: Pt[] = [];
        const f = (p: Pt) => n[0] * (p[0] - o[0]) + n[1] * (p[1] - o[1]);
        for (let i = 0; i < poly.length; ++i) {
            const a = poly[i];
            const b = poly[(i + 1) % poly.length];
            const fa = f(a);
            const fb = f(b);
            if (fa >= 0) { out.push(a); }
            if ((fa > 0 && fb < 0) || (fa < 0 && fb > 0)) {
                const t = fa / (fa - fb);
                out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
            }
        }
        return out;
    }

    function polygonArea(poly: Pt[]): number {
        let s = 0;
        for (let i = 0; i < poly.length; ++i) {
            const a = poly[i];
            const b = poly[(i + 1) % poly.length];
            s += a[0] * b[1] - a[1] * b[0];
        }
        return Math.abs(s) / 2;
    }

    // The 2D solid cone with half-angle in (0, pi/2) is the wedge bounded by
    // the two rays at +-angle from the axis, that is, the intersection of two
    // half-planes through the cone vertex. This is an independent model of
    // the cone that does not use the query's quadratic inequality.
    function overlapArea(box: OrientedBox, cone: Cone): number {
        const c = box.center.values;
        const u0 = box.axis[0].values;
        const u1 = box.axis[1].values;
        const e = box.extent.values;
        let poly: Pt[] = [
            [c[0] - e[0] * u0[0] - e[1] * u1[0],
                c[1] - e[0] * u0[1] - e[1] * u1[1]],
            [c[0] + e[0] * u0[0] - e[1] * u1[0],
                c[1] + e[0] * u0[1] - e[1] * u1[1]],
            [c[0] + e[0] * u0[0] + e[1] * u1[0],
                c[1] + e[0] * u0[1] + e[1] * u1[1]],
            [c[0] - e[0] * u0[0] + e[1] * u1[0],
                c[1] - e[0] * u0[1] + e[1] * u1[1]]
        ];
        const V: Pt = [cone.ray.origin.values[0], cone.ray.origin.values[1]];
        const b = Math.atan2(cone.ray.direction.values[1],
            cone.ray.direction.values[0]);
        const a = cone.angle;
        const n0: Pt = [Math.cos(b + Math.PI / 2 - a),
            Math.sin(b + Math.PI / 2 - a)];
        const n1: Pt = [Math.cos(b - Math.PI / 2 + a),
            Math.sin(b - Math.PI / 2 + a)];
        poly = clipHalfPlane(poly, n0, V);
        if (poly.length === 0) { return 0; }
        poly = clipHalfPlane(poly, n1, V);
        return poly.length >= 3 ? polygonArea(poly) : 0;
    }

    it('agrees with exact clipping of the box against the cone wedge', () => {
        const rnd = seededRandom(0xc02e2);
        let numPositive = 0;
        let numSkipped = 0;
        for (let iter = 0; iter < 5000; ++iter) {
            const a = rnd() * Math.PI;
            const ca = Math.cos(a);
            const sa = Math.sin(a);
            const box = OrientedBox.fromCenterAxisExtent(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]),
                [Vector.fromArray([ca, sa]), Vector.fromArray([-sa, ca])],
                Vector.fromArray([0.2 + rnd() * 1.5, 0.2 + rnd() * 1.5]));
            const b = rnd() * 2 * Math.PI;
            const cone = Cone.fromRayAngle(Ray.fromOriginDirection(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]),
                Vector.fromArray([Math.cos(b), Math.sin(b)])),
                0.05 + rnd() * 1.45);

            const area = overlapArea(box, cone);
            // Skip the sliver band, where the two computations can disagree
            // through rounding of a near measure-zero contact.
            if (area > 1e-12 && area < 1e-6) { ++numSkipped; continue; }
            if (area > 0) { ++numPositive; }
            expect(q.test(box, cone).intersect).toBe(area > 1e-9);
        }
        expect(numPositive).toBeGreaterThan(1500);
        expect(numSkipped).toBeLessThan(20);
    }, 30000);

    it('reports a box strictly inside the cone', () => {
        check(fc.tuple(positive(1), positive(1)), ([e0, e1]) => {
            // The cone is the +x wedge of half-angle pi/4 at the origin; the
            // box sits far along +x, well inside.
            const cone = Cone.fromRayAngle(Ray.fromOriginDirection(
                Vector.zero(2), Vector.fromArray([1, 0])), Math.PI / 4);
            const box = OrientedBox.fromCenterAxisExtent(
                Vector.fromArray([100, 0]),
                [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
                Vector.fromArray([e0, e1]));
            expect(q.test(box, cone).intersect).toBe(true);
        });
    });

    it('does not report a box that only touches the cone boundary', () => {
        // The cone is the +x wedge of half-angle pi/4 at the origin, so its
        // upper boundary is the ray y = x. The box [1,3] x [3,5] lies in
        // y >= x and touches the boundary only at the corner (3,3).
        const cone = Cone.fromRayAngle(Ray.fromOriginDirection(
            Vector.zero(2), Vector.fromArray([1, 0])), Math.PI / 4);
        const touching = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([2, 4]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(q.test(touching, cone).intersect).toBe(false);
        // Nudging the box down puts the corner (3, 2.9) strictly inside.
        const inside = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([2, 3.9]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(q.test(inside, cone).intersect).toBe(true);
    });

    it('accepts a box that contains the cone vertex through the axis test', () => {
        // No box corner need be inside a narrow cone whose vertex is inside
        // the box; the ray-box quick-acceptance test covers this.
        check(fc.tuple(unitVector(2), positive(1)), ([d, s]) => {
            const box = OrientedBox.fromCenterAxisExtent(Vector.zero(2),
                [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
                Vector.fromArray([1, 1]));
            const cone = Cone.fromRayAngle(
                Ray.fromOriginDirection(Vector.zero(2), d), 1e-3 + s * 0.005);
            expect(q.test(box, cone).intersect).toBe(true);
        });
    });

    it('is equivariant under a rigid motion of box and cone', () => {
        check(fc.tuple(arbOrientedBox(2), rotationFrame(2), wellScaledVector(2),
            unitVector(2), positive(1.4, 0.05)),
            ([box, R, t, d, ang]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(R[0], v), dot(R[1], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const cone = Cone.fromRayAngle(
                    Ray.fromOriginDirection(Vector.zero(2), d), ang);
                const area = overlapArea(box, cone);
                // Skip the measure-zero contacts, where a rigid motion can
                // flip the strict comparisons.
                if (area < 1e-6) { return; }
                const box2 = OrientedBox.fromCenterAxisExtent(map(box.center),
                    box.axis.map(rot), box.extent);
                const cone2 = Cone.fromRayAngle(Ray.fromOriginDirection(
                    map(cone.ray.origin), rot(cone.ray.direction)), ang);
                expect(q.test(box2, cone2).intersect)
                    .toBe(q.test(box, cone).intersect);
            });
    });

    it('ignores the cone height range (the cone is treated as infinite)', () => {
        // The query never reads getMinHeight/getMaxHeight, so a truncated
        // cone gives the same answer as the infinite one. This pins the
        // upstream convention.
        const ray = Ray.fromOriginDirection(Vector.zero(2),
            Vector.fromArray([1, 0]));
        const infinite = Cone.fromRayAngle(ray, Math.PI / 6);
        const frustum = Cone.fromRayAngleMinMaxHeight(ray, Math.PI / 6, 1, 2);
        const far = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([100, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(q.test(far, infinite).intersect).toBe(true);
        expect(q.test(far, frustum).intersect).toBe(true);
    });

    it('uses only the box corners and the axis ray, as upstream documents', () => {
        // Reproduce upstream's algorithm independently (ray-box quick accept,
        // then the four corners tested with the quadratic inequality) and
        // check that the port matches it exactly, including the strict '>'.
        const rbQuery = new IntrRay2OrientedBox2TI();
        check(fc.tuple(arbOrientedBox(2), wellScaledVector(2), unitVector(2),
            positive(1.5, 0.05)),
            ([box, v, d, ang]) => {
                const cone = Cone.fromRayAngle(
                    Ray.fromOriginDirection(v, d), ang);
                let expected = rbQuery.test(cone.ray, box).intersect;
                if (!expected) {
                    for (let i1 = 0; i1 < 2 && !expected; ++i1) {
                        for (let i0 = 0; i0 < 2 && !expected; ++i0) {
                            const P = add(box.center, add(
                                mul((i0 * 2 - 1) * box.extent.values[0],
                                    box.axis[0]),
                                mul((i1 * 2 - 1) * box.extent.values[1],
                                    box.axis[1])));
                            const diff = sub(P, v);
                            const h = dot(d, diff);
                            if (h > 0
                                && h * h > dot(diff, diff) * cone.cosAngleSqr) {
                                expected = true;
                            }
                        }
                    }
                }
                expect(q.test(box, cone).intersect).toBe(expected);
            });
    });
});
