import { describe, it, expect } from 'vitest';
import { Cone } from '../src/Cone';
import { GTE_C_HALF_PI, GTE_C_QUARTER_PI, GTE_C_TWO_PI } from '../src/Constants';
import { Ray } from '../src/Ray';
import { Vector, dot, normalize, sub } from '../src/Vector';

function V(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function unit(x: number, y: number, z: number): Vector {
    const v = V(x, y, z);
    normalize(v);
    return v;
}

function makeRay(origin: Vector, direction: Vector): Ray {
    return Ray.fromOriginDirection(origin, direction);
}

// The solid-cone containment test from the class documentation.
function inCone(cone: Cone, point: Vector): boolean {
    const delta = sub(point, cone.ray.origin);
    const h = dot(cone.ray.direction, delta);
    return cone.heightInRange(h) && h * h >= dot(delta, delta) * cone.cosAngleSqr;
}

// The direct (square-root) containment test, for cross-checking.
function inConeDirect(cone: Cone, point: Vector): boolean {
    const delta = sub(point, cone.ray.origin);
    const h = dot(cone.ray.direction, delta);
    const len = Math.sqrt(dot(delta, delta));
    return cone.heightInRange(h) && h >= len * cone.cosAngle;
}

describe('Cone', () => {
    it('default-constructs an infinite cone with angle pi/4', () => {
        const cone = new Cone(3);
        expect(cone.ray.origin.values).toEqual([0, 0, 0]);
        expect(cone.ray.direction.values).toEqual([0, 0, 1]);
        expect(cone.angle).toBe(GTE_C_QUARTER_PI);
        expect(cone.getMinHeight()).toBe(0);
        expect(cone.getMaxHeight()).toBe(-1);
        expect(cone.isInfinite()).toBe(true);
        expect(cone.isFinite()).toBe(false);
        expect(cone.dimension).toBe(3);

        // A 2D cone (a wedge) has the axis (0,1).
        const wedge = new Cone(2);
        expect(wedge.ray.direction.values).toEqual([0, 1]);
    });

    it('computes the derived trigonometric quantities in setAngle', () => {
        const cone = new Cone(3);
        cone.setAngle(Math.PI / 6);
        expect(cone.angle).toBe(Math.PI / 6);
        expect(cone.cosAngle).toBeCloseTo(Math.sqrt(3) / 2, 15);
        expect(cone.sinAngle).toBeCloseTo(0.5, 15);
        expect(cone.tanAngle).toBeCloseTo(1 / Math.sqrt(3), 15);
        expect(cone.cosAngleSqr).toBeCloseTo(0.75, 15);
        expect(cone.sinAngleSqr).toBeCloseTo(0.25, 15);
        expect(cone.invSinAngle).toBeCloseTo(2, 15);
        expect(cone.cosAngleSqr + cone.sinAngleSqr).toBeCloseTo(1, 15);

        // The angle must be in (0,pi/2).
        expect(() => cone.setAngle(0)).toThrow();
        expect(() => cone.setAngle(-0.1)).toThrow();
        expect(() => cone.setAngle(GTE_C_HALF_PI)).toThrow();
        expect(() => cone.setAngle(2)).toThrow();
    });

    it('builds the four cone types with the -1 infinity convention', () => {
        const ray = makeRay(V(1, 2, 3), unit(0, 0, 1));

        // 1. Infinite cone.
        const infinite = Cone.fromRayAngle(ray, 0.5);
        expect(infinite.getMinHeight()).toBe(0);
        expect(infinite.getMaxHeight()).toBe(-1);
        expect(infinite.isInfinite()).toBe(true);
        expect(infinite.ray.origin.values).toEqual([1, 2, 3]);

        // 2. Infinite truncated cone.
        const truncated = Cone.fromRayAngleMinHeight(ray, 0.5, 2);
        expect(truncated.getMinHeight()).toBe(2);
        expect(truncated.isInfinite()).toBe(true);

        // 3. Finite cone.
        const finite = Cone.fromRayAngleMinMaxHeight(ray, 0.5, 0, 4);
        expect(finite.getMinHeight()).toBe(0);
        expect(finite.getMaxHeight()).toBe(4);
        expect(finite.isFinite()).toBe(true);

        // 4. Frustum of a cone.
        const frustum = Cone.fromRayAngleMinMaxHeight(ray, 0.5, 1, 4);
        expect(frustum.getMinHeight()).toBe(1);
        expect(frustum.getMaxHeight()).toBe(4);
        expect(frustum.isFinite()).toBe(true);

        // makeFiniteCone resets the minimum height to zero.
        const cone = frustum.clone();
        cone.makeFiniteCone(7);
        expect(cone.getMinHeight()).toBe(0);
        expect(cone.getMaxHeight()).toBe(7);

        // Invalid heights.
        expect(() => cone.makeInfiniteTruncatedCone(-1)).toThrow();
        expect(() => cone.makeFiniteCone(0)).toThrow();
        expect(() => cone.makeConeFrustum(2, 2)).toThrow();
        expect(() => cone.makeConeFrustum(-1, 3)).toThrow();

        // The ray is copied by the factories and by clone.
        const source = makeRay(V(0, 0, 0), unit(0, 0, 1));
        const copied = Cone.fromRayAngle(source, 0.5);
        source.origin.values[0] = 42;
        expect(copied.ray.origin.values[0]).toBe(0);
        const cloned = copied.clone();
        copied.ray.origin.values[1] = 9;
        expect(cloned.ray.origin.values[1]).toBe(0);
        expect(cloned.getMaxHeight()).toBe(copied.getMaxHeight());
    });

    it('answers the height-range queries for each cone type', () => {
        const ray = makeRay(V(0, 0, 0), unit(0, 0, 1));

        const infinite = Cone.fromRayAngle(ray, 0.5);
        expect(infinite.heightInRange(0)).toBe(true);
        expect(infinite.heightInRange(1e100)).toBe(true);
        expect(infinite.heightInRange(-1e-9)).toBe(false);
        expect(infinite.heightLessThanMin(-1)).toBe(true);
        expect(infinite.heightGreaterThanMax(1e100)).toBe(false);

        const frustum = Cone.fromRayAngleMinMaxHeight(ray, 0.5, 1, 4);
        expect(frustum.heightInRange(1)).toBe(true);
        expect(frustum.heightInRange(4)).toBe(true);
        expect(frustum.heightInRange(0.99)).toBe(false);
        expect(frustum.heightInRange(4.01)).toBe(false);
        expect(frustum.heightLessThanMin(0.99)).toBe(true);
        expect(frustum.heightLessThanMin(1)).toBe(false);
        expect(frustum.heightGreaterThanMax(4.01)).toBe(true);
        expect(frustum.heightGreaterThanMax(4)).toBe(false);
    });

    it('has a quadratic containment test that matches the direct one', () => {
        const cone = Cone.fromRayAngleMinMaxHeight(
            makeRay(V(1, -2, 0.5), unit(1, 2, -2)), Math.PI / 5, 0.5, 3);

        // Points on the axis at valid heights are inside.
        for (const h of [0.6, 1, 2, 2.9]) {
            const p = V(
                cone.ray.origin.values[0] + h * cone.ray.direction.values[0],
                cone.ray.origin.values[1] + h * cone.ray.direction.values[1],
                cone.ray.origin.values[2] + h * cone.ray.direction.values[2]);
            expect(inCone(cone, p)).toBe(true);
        }

        // The cone vertex is outside because the minimum height is positive.
        expect(inCone(cone, cone.ray.origin)).toBe(false);

        // The height-range endpoints are inclusive. Use an axis-aligned cone
        // so that the heights are computed exactly.
        const axisCone = Cone.fromRayAngleMinMaxHeight(
            makeRay(V(0, 0, 0), unit(0, 0, 1)), Math.PI / 5, 0.5, 3);
        expect(inCone(axisCone, V(0, 0, 0.5))).toBe(true);
        expect(inCone(axisCone, V(0, 0, 3))).toBe(true);
        expect(inCone(axisCone, V(0, 0, 0.25))).toBe(false);
        expect(inCone(axisCone, V(0, 0, 3.25))).toBe(false);
        // A point on the lateral surface at height 1: r = tan(A).
        expect(inCone(axisCone, V(axisCone.tanAngle * (1 - 1e-12), 0, 1)))
            .toBe(true);
        expect(inCone(axisCone, V(axisCone.tanAngle * (1 + 1e-9), 0, 1)))
            .toBe(false);

        // A point mirrored through the vertex lies in the negative cone: the
        // quadratic form alone accepts it, but the height test rejects it.
        const axisPoint = V(
            cone.ray.origin.values[0] - 2 * cone.ray.direction.values[0],
            cone.ray.origin.values[1] - 2 * cone.ray.direction.values[1],
            cone.ray.origin.values[2] - 2 * cone.ray.direction.values[2]);
        const delta = sub(axisPoint, cone.ray.origin);
        const h = dot(cone.ray.direction, delta);
        expect(h * h >= dot(delta, delta) * cone.cosAngleSqr).toBe(true);
        expect(cone.heightInRange(h)).toBe(false);
        expect(inCone(cone, axisPoint)).toBe(false);

        // Randomized cross-check of the quadratic and direct forms.
        let seed = 24680;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const R = (): number => 8 * rand() - 4;

        let numInside = 0;
        let numOutside = 0;
        const cones = [
            cone,
            Cone.fromRayAngle(makeRay(V(0, 0, 0), unit(0, 0, 1)), GTE_C_QUARTER_PI),
            Cone.fromRayAngleMinHeight(makeRay(V(-1, 0, 2), unit(-1, 1, 0)), 0.6, 1),
            Cone.fromRayAngleMinMaxHeight(makeRay(V(0, 1, 0), unit(0, -1, 0)), 1.2, 0, 5)
        ];
        for (let trial = 0; trial < 400; ++trial) {
            const p = V(R(), R(), R());
            for (const c of cones) {
                const a = inCone(c, p);
                const b = inConeDirect(c, p);
                expect(a).toBe(b);
                if (a) {
                    ++numInside;
                }
                else {
                    ++numOutside;
                }
            }
        }
        expect(numInside).toBeGreaterThan(20);
        expect(numOutside).toBeGreaterThan(20);

        // The 45-degree cone about the z-axis has the closed form
        // z >= sqrt(x^2+y^2).
        const c45 = cones[1];
        for (let trial = 0; trial < 200; ++trial) {
            const p = V(R(), R(), R());
            const expected = p.values[2] >= Math.hypot(p.values[0], p.values[1]);
            expect(inConeDirect(c45, p)).toBe(expected);
        }
    });

    it('supports the lexicographic comparisons', () => {
        const ray = makeRay(V(0, 0, 0), unit(0, 0, 1));
        const a = Cone.fromRayAngleMinMaxHeight(ray, 0.5, 1, 4);
        const b = a.clone();
        expect(a.equals(b)).toBe(true);
        expect(a.notEquals(b)).toBe(false);
        expect(a.lessThan(b)).toBe(false);
        expect(a.lessThanOrEqual(b)).toBe(true);
        expect(a.greaterThanOrEqual(b)).toBe(true);

        const wider = Cone.fromRayAngleMinMaxHeight(ray, 0.6, 1, 4);
        expect(a.lessThan(wider)).toBe(true);
        expect(wider.greaterThan(a)).toBe(true);
        expect(a.equals(wider)).toBe(false);

        const taller = Cone.fromRayAngleMinMaxHeight(ray, 0.5, 1, 5);
        expect(a.lessThan(taller)).toBe(true);

        const shifted = Cone.fromRayAngleMinMaxHeight(
            makeRay(V(1, 0, 0), unit(0, 0, 1)), 0.5, 1, 4);
        expect(a.lessThan(shifted)).toBe(true);
        expect(shifted.lessThan(a)).toBe(false);
    });

    it('creates an inscribed frustum mesh with the expected geometry', () => {
        const cone = Cone.fromRayAngleMinMaxHeight(
            makeRay(V(0, 0, 0), unit(0, 0, 1)), GTE_C_QUARTER_PI, 1, 2);
        const numMinVertices = 8;
        const { vertices, indices } = cone.createMesh(numMinVertices, true);

        const rMin = 1 * cone.tanAngle;
        const rMax = 2 * cone.tanAngle;
        // tNumExtra = 0.5*(1+rMax)/(1+rMin) - 1 = 0.5*3/2 - 1 < 0.
        const numExtra = 0;
        const numMaxVertices = 2 * numMinVertices * (1 + numExtra);

        expect(indices.length % 3).toBe(0);
        // The triangle soup has numMinVertices*(2*(numExtra+1)+1) lateral
        // triangles plus numMinVertices + numMaxVertices disk triangles.
        const numTriangles = numMinVertices * (2 * (numExtra + 1) + 1) +
            numMinVertices + numMaxVertices;
        expect(indices.length).toBe(3 * numTriangles);
        for (const index of indices) {
            expect(index).toBeGreaterThanOrEqual(0);
            expect(index).toBeLessThan(vertices.length);
        }
        // Every triangle references three distinct vertices (there are no
        // degenerate triangles in a proper frustum).
        for (let t = 0; t < numTriangles; ++t) {
            const tri = indices.slice(3 * t, 3 * t + 3);
            expect(new Set(tri).size).toBe(3);
        }

        // The unique vertices are the two disk centers plus the two polygon
        // rings.
        expect(vertices.length).toBe(2 + numMinVertices + numMaxVertices);
        const keys = new Set(vertices.map(v => v.values.join(',')));
        expect(keys.size).toBe(vertices.length);

        // The cone axis is the z-axis, so every vertex lies at height 1 or 2
        // and at radius 0, rMin or rMax.
        let numMinRing = 0;
        let numMaxRing = 0;
        let numCenters = 0;
        for (const v of vertices) {
            const radius = Math.hypot(v.values[0], v.values[1]);
            const z = v.values[2];
            expect([1, 2]).toContain(z);
            if (radius < 1e-12) {
                ++numCenters;
            }
            else if (z === 1) {
                expect(radius).toBeCloseTo(rMin, 12);
                ++numMinRing;
            }
            else {
                expect(radius).toBeCloseTo(rMax, 12);
                ++numMaxRing;
            }
        }
        expect(numCenters).toBe(2);
        expect(numMinRing).toBe(numMinVertices);
        expect(numMaxRing).toBe(numMaxVertices);

        // The lateral vertices satisfy the cone equation h = r/tan(A).
        for (const v of vertices) {
            const radius = Math.hypot(v.values[0], v.values[1]);
            if (radius > 1e-12) {
                expect(v.values[2] * cone.tanAngle).toBeCloseTo(radius, 12);
            }
        }
    });

    it('creates a circumscribed mesh and honors the cone frame', () => {
        // A circumscribed polygon has vertices at radius r/cos(theta/2).
        const numMinVertices = 6;
        const cone = Cone.fromRayAngleMinMaxHeight(
            makeRay(V(0, 0, 0), unit(0, 0, 1)), GTE_C_QUARTER_PI, 1, 2);
        const { vertices } = cone.createMesh(numMinVertices, false);
        const numMaxVertices = 2 * numMinVertices;
        const rMin = cone.tanAngle;
        const rMax = 2 * cone.tanAngle;
        const scaleMin = 1 / Math.cos(GTE_C_TWO_PI / (2 * numMinVertices));
        const scaleMax = 1 / Math.cos(GTE_C_TWO_PI / (2 * numMaxVertices));
        for (const v of vertices) {
            const radius = Math.hypot(v.values[0], v.values[1]);
            if (radius < 1e-12) {
                continue;
            }
            if (v.values[2] === 1) {
                expect(radius).toBeCloseTo(rMin * scaleMin, 12);
            }
            else {
                expect(radius).toBeCloseTo(rMax * scaleMax, 12);
            }
        }

        // A cone whose axis is not the z-axis: every vertex must have the
        // correct height along the axis and the correct distance from it.
        const oblique = Cone.fromRayAngleMinMaxHeight(
            makeRay(V(1, -2, 3), unit(1, 2, 2)), Math.PI / 6, 1, 3);
        const mesh = oblique.createMesh(7, true);
        const rMinO = 1 * oblique.tanAngle;
        const rMaxO = 3 * oblique.tanAngle;
        for (const v of mesh.vertices) {
            const delta = sub(v, oblique.ray.origin);
            const h = dot(oblique.ray.direction, delta);
            // The radius involves a cancellation, so it is accurate only to
            // about the square root of the machine epsilon.
            const radius = Math.sqrt(Math.max(dot(delta, delta) - h * h, 0));
            const isMinRing = (h < 2);
            expect(h).toBeCloseTo(isMinRing ? 1 : 3, 10);
            if (radius > 1e-6) {
                expect(radius).toBeCloseTo(isMinRing ? rMinO : rMaxO, 6);
            }
            else {
                // A disk center lies on the axis.
                expect(radius).toBeLessThan(1e-6);
            }
        }
    });

    it('creates a finite-cone mesh with an apex and rejects bad input', () => {
        // A finite cone (hMin = 0) degenerates the h-minimum disk to the apex.
        const cone = Cone.fromRayAngleMinMaxHeight(
            makeRay(V(0, 0, 0), unit(0, 0, 1)), GTE_C_QUARTER_PI, 0, 2);
        const numMinVertices = 8;
        const { vertices, indices } = cone.createMesh(numMinVertices, true);

        // tNumExtra = 0.5*(1+2)/(1+0) - 1 = 0.5, so numExtra = 1.
        const numExtra = 1;
        const numMaxVertices = 2 * numMinVertices * (1 + numExtra);
        // All h-minimum vertices collapse onto the apex, so the unique count
        // is one apex plus the maximum ring plus the h-maximum disk center.
        expect(vertices.length).toBe(2 + numMaxVertices);
        let numApex = 0;
        for (const v of vertices) {
            if (v.values[2] === 0) {
                expect(Math.hypot(v.values[0], v.values[1]))
                    .toBeLessThan(1e-15);
                ++numApex;
            }
        }
        expect(numApex).toBe(1);
        for (const index of indices) {
            expect(index).toBeLessThan(vertices.length);
        }

        // Meshes require finite, 3-dimensional cones.
        const infinite = Cone.fromRayAngle(
            makeRay(V(0, 0, 0), unit(0, 0, 1)), 0.5);
        expect(() => infinite.createMesh(8, true)).toThrow();
        expect(() => cone.createMesh(2, true)).toThrow();
        const cone2 = new Cone(2);
        cone2.makeFiniteCone(1);
        expect(() => cone2.createMesh(8, true)).toThrow();
    });
});
