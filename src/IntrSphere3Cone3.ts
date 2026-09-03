// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSphere3Cone3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection query is based on the document
// https://www.geometrictools.com/Documentation/IntersectionSphereCone.pdf
//
// The find-intersection query returns a single point in the set of
// intersection when that intersection is not empty.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// TIQuery and FIQuery specializations become IntrSphere3Cone3TI and
// IntrSphere3Cone3FI. The private DoQuery* helpers become module-private
// functions.

import type { Cone } from './Cone.js';
import type { Hypersphere } from './Hypersphere.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import { Vector, add, dot, length, mul, sub } from './Vector.js';
import { cross } from './Vector3.js';

// The result of IntrSphere3Cone3TI.test.
export interface IntrSphere3Cone3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSphere3Cone3TIResult(): IntrSphere3Cone3TIResult {
    return { intersect: false };
}

// The result of IntrSphere3Cone3FI.find.
export interface IntrSphere3Cone3FIResult {
    // If an intersection occurs, it is potentially an infinite set. If the
    // cone vertex is inside the sphere, 'point' is set to the cone vertex. If
    // the sphere center is inside the cone, 'point' is set to the sphere
    // center. Otherwise, 'point' is set to the cone point that is closest to
    // the cone vertex and inside the sphere.
    intersect: boolean;
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSphere3Cone3FIResult(): IntrSphere3Cone3FIResult {
    return { intersect: false, point: new Vector(3) };
}

// The cone is infinite (hmin = 0, hmax = +infinity).
function doQueryInfiniteCone(sphere: Hypersphere, cone: Cone): boolean {
    const U = sub(cone.ray.origin,
        mul(sphere.radius * cone.invSinAngle, cone.ray.direction));
    const CmU = sub(sphere.center, U);
    const AdCmU = dot(cone.ray.direction, CmU);
    if (AdCmU > 0) {
        const sqrLengthCmU = dot(CmU, CmU);
        if (AdCmU * AdCmU >= sqrLengthCmU * cone.cosAngleSqr) {
            const CmV = sub(sphere.center, cone.ray.origin);
            const AdCmV = dot(cone.ray.direction, CmV);
            if (AdCmV < -sphere.radius) {
                return false;
            }

            const rSinAngle = sphere.radius * cone.sinAngle;
            if (AdCmV >= -rSinAngle) {
                return true;
            }

            const sqrLengthCmV = dot(CmV, CmV);
            return sqrLengthCmV <= sphere.radius * sphere.radius;
        }
    }

    return false;
}

// The cone is truncated for h-minimum (hmin > 0, hmax = +infinity).
function doQueryInfiniteTruncatedCone(sphere: Hypersphere,
    cone: Cone): boolean {
    const U = sub(cone.ray.origin,
        mul(sphere.radius * cone.invSinAngle, cone.ray.direction));
    const CmU = sub(sphere.center, U);
    const AdCmU = dot(cone.ray.direction, CmU);
    if (AdCmU > 0) {
        const sqrLengthCmU = dot(CmU, CmU);
        if (AdCmU * AdCmU >= sqrLengthCmU * cone.cosAngleSqr) {
            const CmV = sub(sphere.center, cone.ray.origin);
            const AdCmV = dot(cone.ray.direction, CmV);
            const minHeight = cone.getMinHeight();
            if (AdCmV < minHeight - sphere.radius) {
                return false;
            }

            const rSinAngle = sphere.radius * cone.sinAngle;
            if (AdCmV >= -rSinAngle) {
                return true;
            }

            const D = sub(CmV, mul(minHeight, cone.ray.direction));
            const lengthAxD = length(cross(cone.ray.direction, D));
            const hminTanAngle = minHeight * cone.tanAngle;
            if (lengthAxD <= hminTanAngle) {
                return true;
            }

            const AdD = AdCmV - minHeight;
            const diff = lengthAxD - hminTanAngle;
            const sqrLengthCmK = AdD * AdD + diff * diff;
            return sqrLengthCmK <= sphere.radius * sphere.radius;
        }
    }

    return false;
}

// The cone is truncated for h-maximum (hmin = 0, hmax < +infinity).
function doQueryFiniteCone(sphere: Hypersphere, cone: Cone): boolean {
    const U = sub(cone.ray.origin,
        mul(sphere.radius * cone.invSinAngle, cone.ray.direction));
    const CmU = sub(sphere.center, U);
    const AdCmU = dot(cone.ray.direction, CmU);
    if (AdCmU > 0) {
        const sqrLengthCmU = dot(CmU, CmU);
        if (AdCmU * AdCmU >= sqrLengthCmU * cone.cosAngleSqr) {
            const CmV = sub(sphere.center, cone.ray.origin);
            const AdCmV = dot(cone.ray.direction, CmV);
            if (AdCmV < -sphere.radius) {
                return false;
            }

            const maxHeight = cone.getMaxHeight();
            if (AdCmV > maxHeight + sphere.radius) {
                return false;
            }

            const rSinAngle = sphere.radius * cone.sinAngle;
            if (AdCmV >= -rSinAngle) {
                if (AdCmV <= maxHeight - rSinAngle) {
                    return true;
                } else {
                    const barD = sub(CmV, mul(maxHeight, cone.ray.direction));
                    const lengthAxBarD = length(
                        cross(cone.ray.direction, barD));
                    const hmaxTanAngle = maxHeight * cone.tanAngle;
                    if (lengthAxBarD <= hmaxTanAngle) {
                        return true;
                    }

                    const AdBarD = AdCmV - maxHeight;
                    const diff = lengthAxBarD - hmaxTanAngle;
                    const sqrLengthCmBarK = AdBarD * AdBarD + diff * diff;
                    return sqrLengthCmBarK <= sphere.radius * sphere.radius;
                }
            } else {
                const sqrLengthCmV = dot(CmV, CmV);
                return sqrLengthCmV <= sphere.radius * sphere.radius;
            }
        }
    }

    return false;
}

// The cone is truncated for both h-minimum and h-maximum.
function doQueryConeFrustum(sphere: Hypersphere, cone: Cone): boolean {
    const U = sub(cone.ray.origin,
        mul(sphere.radius * cone.invSinAngle, cone.ray.direction));
    const CmU = sub(sphere.center, U);
    const AdCmU = dot(cone.ray.direction, CmU);
    if (AdCmU > 0) {
        const sqrLengthCmU = dot(CmU, CmU);
        if (AdCmU * AdCmU >= sqrLengthCmU * cone.cosAngleSqr) {
            const CmV = sub(sphere.center, cone.ray.origin);
            const AdCmV = dot(cone.ray.direction, CmV);
            const minHeight = cone.getMinHeight();
            if (AdCmV < minHeight - sphere.radius) {
                return false;
            }

            const maxHeight = cone.getMaxHeight();
            if (AdCmV > maxHeight + sphere.radius) {
                return false;
            }

            const rSinAngle = sphere.radius * cone.sinAngle;
            if (AdCmV >= minHeight - rSinAngle) {
                if (AdCmV <= maxHeight - rSinAngle) {
                    return true;
                } else {
                    const barD = sub(CmV, mul(maxHeight, cone.ray.direction));
                    const lengthAxBarD = length(
                        cross(cone.ray.direction, barD));
                    const hmaxTanAngle = maxHeight * cone.tanAngle;
                    if (lengthAxBarD <= hmaxTanAngle) {
                        return true;
                    }

                    const AdBarD = AdCmV - maxHeight;
                    const diff = lengthAxBarD - hmaxTanAngle;
                    const sqrLengthCmBarK = AdBarD * AdBarD + diff * diff;
                    return sqrLengthCmBarK <= sphere.radius * sphere.radius;
                }
            } else {
                const D = sub(CmV, mul(minHeight, cone.ray.direction));
                const lengthAxD = length(cross(cone.ray.direction, D));
                const hminTanAngle = minHeight * cone.tanAngle;
                if (lengthAxD <= hminTanAngle) {
                    return true;
                }

                const AdD = AdCmV - minHeight;
                const diff = lengthAxD - hminTanAngle;
                const sqrLengthCmK = AdD * AdD + diff * diff;
                return sqrLengthCmK <= sphere.radius * sphere.radius;
            }
        }
    }

    return false;
}

// Test-intersection query for a solid sphere and a solid cone.
export class IntrSphere3Cone3TI implements
    TIQuery<Hypersphere, Cone, IntrSphere3Cone3TIResult> {

    test(sphere: Hypersphere, cone: Cone): IntrSphere3Cone3TIResult {
        const result = defaultIntrSphere3Cone3TIResult();
        if (cone.getMinHeight() > 0) {
            if (cone.isFinite()) {
                result.intersect = doQueryConeFrustum(sphere, cone);
            } else {
                result.intersect = doQueryInfiniteTruncatedCone(sphere, cone);
            }
        } else {
            if (cone.isFinite()) {
                result.intersect = doQueryFiniteCone(sphere, cone);
            } else {
                result.intersect = doQueryInfiniteCone(sphere, cone);
            }
        }
        return result;
    }
}

// Find-intersection query for a solid sphere and an infinite solid cone. The
// query ignores the cone height range; it treats the cone as the infinite
// cone with the same vertex, axis and angle.
export class IntrSphere3Cone3FI implements
    FIQuery<Hypersphere, Cone, IntrSphere3Cone3FIResult> {

    find(sphere: Hypersphere, cone: Cone): IntrSphere3Cone3FIResult {
        const result = defaultIntrSphere3Cone3FIResult();

        // Test whether the cone vertex is inside the sphere.
        const diff = sub(sphere.center, cone.ray.origin);
        const rSqr = sphere.radius * sphere.radius;
        const lenSqr = dot(diff, diff);
        if (lenSqr <= rSqr) {
            // The cone vertex is inside the sphere, so the sphere and cone
            // intersect.
            result.intersect = true;
            result.point = cone.ray.origin.clone();
            return result;
        }

        // Test whether the sphere center is inside the cone.
        const dotAD = dot(diff, cone.ray.direction);
        const dotSqr = dotAD * dotAD;
        if (dotSqr >= lenSqr * cone.cosAngleSqr && dotAD > 0) {
            // The sphere center is inside the cone, so the sphere and cone
            // intersect.
            result.intersect = true;
            result.point = sphere.center.clone();
            return result;
        }

        // The sphere center is outside the cone. The problem now reduces to
        // computing an intersection between the circle and the ray in the
        // plane containing the cone vertex and spanned by the cone axis and
        // the vector from the cone vertex to the sphere center.
        //
        // The ray is parameterized by t * D + V with t >= 0, |D| = 1 and
        // Dot(A,D) = cos(angle). Also, D = e * A + f * (C - V). Substituting
        // the ray equation into the sphere equation yields
        // R^2 = |t * D + V - C|^2, so the quadratic for intersections is
        // t^2 - 2 * Dot(D, C - V) * t + |C - V|^2 - R^2 = 0. An intersection
        // occurs if and only if the discriminant is nonnegative. This test
        // becomes
        //     Dot(D, C - V)^2 >= Dot(C - V, C - V) - R^2
        // Note that if the right-hand side is nonpositive, then the
        // inequality is true (the sphere contains V). This is already ruled
        // out in the first block of code in this function.
        const uLen = Math.sqrt(Math.max(lenSqr - dotSqr, 0));
        const test = cone.cosAngle * dotAD + cone.sinAngle * uLen;
        const discr = test * test - lenSqr + rSqr;

        if (discr >= 0 && test >= 0) {
            // Compute the point of intersection closest to the cone vertex.
            result.intersect = true;
            const t = test - Math.sqrt(Math.max(discr, 0));
            const B = sub(diff, mul(dotAD, cone.ray.direction));
            const tmp = cone.sinAngle / uLen;
            const dir = add(mul(cone.cosAngle, cone.ray.direction), mul(tmp, B));
            // Upstream bug (fixed here): upstream assigns
            //   result.point = t * (cosAngle * A + tmp * B)
            // which is t * D, omitting the cone vertex V. The ray in the
            // reduction is X(t) = V + t * D, so the intersection point is
            // V + t * D. Upstream is correct only when V is the origin.
            result.point = add(cone.ray.origin, mul(t, dir));
        } else {
            result.intersect = false;
        }

        return result;
    }
}
