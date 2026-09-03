import { describe, it, expect } from 'vitest';
import { inContainerCone } from '../src/ContCone.js';
import { Cone } from '../src/Cone.js';
import { Ray } from '../src/Ray.js';
import { Vector, dot, normalize, sub } from '../src/Vector.js';

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
