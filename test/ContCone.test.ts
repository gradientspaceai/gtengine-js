import { describe, it, expect } from 'vitest';
import { inContainerCone } from '../src/ContCone.js';
import { Cone } from '../src/Cone.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, fc, rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray3(origin: Vector, direction: Vector): Ray {
    const d = direction.clone();
    normalize(d);
    return Ray.fromOriginDirection(origin, d);
}

// Independent containment test: the point is in the solid cone when its
// height along the axis is in range and the angle between (P-V) and the axis
// is at most the cone angle.
function referenceInCone(point: Vector, cone: Cone, angle: number): boolean {
    const diff = sub(point, cone.ray.origin);
    const h = dot(cone.ray.direction, diff);
    if (!cone.heightInRange(h)) {
        return false;
    }
    const lenSqr = dot(diff, diff);
    if (lenSqr === 0) {
        return true;
    }
    const cosTheta = h / Math.sqrt(lenSqr);
    return cosTheta >= Math.cos(angle);
}

describe('ContCone', () => {
    const angle = Math.PI / 6; // 30 degrees

    it('contains points on the axis of an infinite cone', () => {
        const cone = Cone.fromRayAngle(ray3(v3(0, 0, 0), v3(0, 0, 1)), angle);
        expect(inContainerCone(v3(0, 0, 0), cone)).toBe(true);
        expect(inContainerCone(v3(0, 0, 1), cone)).toBe(true);
        expect(inContainerCone(v3(0, 0, 1000), cone)).toBe(true);
        // Behind the vertex.
        expect(inContainerCone(v3(0, 0, -1), cone)).toBe(false);
    });

    it('matches the exact boundary of the cone', () => {
        const cone = Cone.fromRayAngle(ray3(v3(0, 0, 0), v3(0, 0, 1)), angle);
        const t = Math.tan(angle);
        // At height 1 the cone radius is tan(angle).
        expect(inContainerCone(v3(t * 0.999, 0, 1), cone)).toBe(true);
        expect(inContainerCone(v3(t * 1.001, 0, 1), cone)).toBe(false);
        // A point exactly on the lateral surface (up to round-off).
        expect(inContainerCone(v3(t, 0, 1), cone)).toBe(true);
    });

    it('honors the height range of a cone frustum', () => {
        const cone = Cone.fromRayAngleMinMaxHeight(
            ray3(v3(1, 2, 3), v3(0, 0, 1)), angle, 2, 5);
        expect(inContainerCone(v3(1, 2, 3), cone)).toBe(false); // h = 0
        expect(inContainerCone(v3(1, 2, 4), cone)).toBe(false); // h = 1 < min
        expect(inContainerCone(v3(1, 2, 5), cone)).toBe(true);  // h = 2
        expect(inContainerCone(v3(1, 2, 8), cone)).toBe(true);  // h = 5
        expect(inContainerCone(v3(1, 2, 8.5), cone)).toBe(false); // h > max
    });

    it('works in 2D', () => {
        const cone2 = Cone.fromRayAngle(
            Ray.fromOriginDirection(Vector.fromArray([0, 0]),
                Vector.fromArray([1, 0])), Math.PI / 4);
        expect(inContainerCone(Vector.fromArray([1, 0.5]), cone2)).toBe(true);
        expect(inContainerCone(Vector.fromArray([1, 1.5]), cone2)).toBe(false);
        expect(inContainerCone(Vector.fromArray([-1, 0]), cone2)).toBe(false);
    });

    it('rejects mismatched dimensions', () => {
        const cone = Cone.fromRayAngle(ray3(v3(0, 0, 0), v3(0, 0, 1)), angle);
        expect(() => inContainerCone(Vector.fromArray([0, 0]), cone)).toThrow();
    });

    it('agrees with the angle-based reference test on random inputs', () => {
        let seed = 20260902;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const axis = v3(1, -2, 3);
        const cone = Cone.fromRayAngleMinMaxHeight(
            ray3(v3(-1, 0.5, 2), axis), 0.4, 1, 6);
        let inside = 0;
        for (let i = 0; i < 3000; ++i) {
            const p = v3(8 * rand() - 4, 8 * rand() - 4, 10 * rand() - 2);
            const actual = inContainerCone(p, cone);
            if (actual) {
                ++inside;
            }
            // Skip points essentially on the lateral surface, where the two
            // formulations can disagree by round-off.
            const diff = sub(p, cone.ray.origin);
            const h = dot(cone.ray.direction, diff);
            if (Math.abs(h * h - cone.cosAngleSqr * dot(diff, diff)) < 1e-9) {
                continue;
            }
            expect(actual).toBe(referenceInCone(p, cone, 0.4));
        }
        // Sanity: the sample hits the cone sometimes but not always.
        expect(inside).toBeGreaterThan(0);
        expect(inside).toBeLessThan(3000);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContCone.h semantics.
// ---------------------------------------------------------------------------

describe('ContCone verification', () => {
    const coneArb = fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
        fc.double({ min: 0.05, max: 1.5, noNaN: true }),
        // The minimum height is bounded away from zero: fc.double({min: 0})
        // emits subnormals, and a cone whose height range starts at 1e-308
        // makes every height comparison below meaningless.
        fc.option(fc.tuple(fc.double({ min: 0.01, max: 3, noNaN: true }),
            fc.double({ min: 0.1, max: 4, noNaN: true })), { nil: undefined }))
        .map(([origin, dir, angle, heights]) => {
            const cone = Cone.fromRayAngle(
                Ray.fromOriginDirection(origin, dir), angle);
            if (heights !== undefined) {
                cone.makeConeFrustum(heights[0], heights[0] + heights[1]);
            }
            return { cone, angle };
        });

    // The quadratic-inequality test must agree with the angle test that the
    // upstream Cone.h documentation states it implements. Points near the
    // lateral surface or a height plane are skipped: the two forms differ
    // there only by rounding.
    it('agrees with the angle-and-height reference test', () => {
        check(fc.tuple(coneArb, wellScaledVector(3, -8, 8)),
            ([{ cone, angle }, p]: [{ cone: Cone, angle: number }, Vector]) => {
                const diff = sub(p, cone.ray.origin);
                const h = dot(cone.ray.direction, diff);
                const len = length(diff);
                if (len === 0) {
                    return;
                }
                const cosTheta = h / len;
                // Distance from the two boundary conditions.
                const angularSlack = Math.cos(angle) - cosTheta;
                const minSlack = h - cone.getMinHeight();
                const maxSlack = cone.isInfinite()
                    ? 1 : cone.getMaxHeight() - h;
                if (Math.abs(angularSlack) < 1e-9 || Math.abs(minSlack) < 1e-9
                    || Math.abs(maxSlack) < 1e-9) {
                    return;
                }
                expect(inContainerCone(p, cone))
                    .toBe(referenceInCone(p, cone, angle));
            });
    });

    // The h*h >= cosAngleSqr*|diff|^2 form is satisfied by the reflected cone
    // as well; upstream relies on HeightInRange to exclude it, because
    // minHeight is never negative. Pin that: reflecting a contained point
    // through the vertex must leave the cone.
    it('excludes the reflected cone', () => {
        check(fc.tuple(coneArb, wellScaledVector(3, -8, 8)),
            ([{ cone }, p]: [{ cone: Cone }, Vector]) => {
                if (!inContainerCone(p, cone)) {
                    return;
                }
                const reflected = sub(mul(2, cone.ray.origin), p);
                const h = dot(cone.ray.direction, sub(reflected, cone.ray.origin));
                if (h === 0) {
                    return;   // the vertex itself is in both halves
                }
                expect(inContainerCone(reflected, cone)).toBe(false);
            });
    });

    // Rigid motions: the answer depends only on the relative geometry.
    it('is equivariant under rigid motions', () => {
        check(fc.tuple(coneArb, wellScaledVector(3, -8, 8), rotationFrame(3),
            wellScaledVector(3)),
            ([{ cone }, p, frame, t]:
                [{ cone: Cone }, Vector, Vector[], Vector]) => {
                const rotate = (q: Vector): Vector =>
                    add(add(mul(q.get(0), frame[0]), mul(q.get(1), frame[1])),
                        mul(q.get(2), frame[2]));
                const xform = (q: Vector): Vector => add(rotate(q), t);
                const moved = cone.clone();
                moved.ray = Ray.fromOriginDirection(xform(cone.ray.origin),
                    rotate(cone.ray.direction));
                // Skip configurations within rounding distance of a boundary.
                const diff = sub(p, cone.ray.origin);
                const h = dot(cone.ray.direction, diff);
                const len = length(diff);
                if (len === 0) {
                    return;
                }
                const slack = h * h - cone.cosAngleSqr * len * len;
                if (Math.abs(slack) < 1e-9 * len * len) {
                    return;
                }
                expect(inContainerCone(xform(p), moved))
                    .toBe(inContainerCone(p, cone));
            });
    });

    // The cone vertex is always contained (h = 0 >= minHeight = 0 for an
    // untruncated cone, and 0 >= 0 for the quadratic inequality).
    it('contains the vertex of an untruncated cone', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
            fc.double({ min: 0.05, max: 1.5, noNaN: true })),
            ([origin, dir, angle]: [Vector, Vector, number]) => {
                const cone = Cone.fromRayAngle(
                    Ray.fromOriginDirection(origin, dir), angle);
                expect(inContainerCone(origin.clone(), cone)).toBe(true);
            });
    });

    // Points on the axis at an in-range height are contained for every angle.
    it('contains axis points at in-range heights', () => {
        check(fc.tuple(coneArb, fc.double({ min: 0, max: 3, noNaN: true })),
            ([{ cone }, s]: [{ cone: Cone }, number]) => {
                const h = cone.getMinHeight()
                    + s * (cone.isInfinite() ? 1
                        : (cone.getMaxHeight() - cone.getMinHeight()) / 3);
                if (!cone.heightInRange(h)) {
                    return;
                }
                const p = add(cone.ray.origin, mul(h, cone.ray.direction));
                // The height recovered from p can differ from h in the last
                // bits, which matters at the ends of the height range.
                const hp = dot(cone.ray.direction, sub(p, cone.ray.origin));
                if (!cone.heightInRange(hp)) {
                    return;
                }
                expect(inContainerCone(p, cone)).toBe(true);
            });
    });
});
