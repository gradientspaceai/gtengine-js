// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ChebyshevRatio.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Chebyshev ratio is f(t,A) = sin(t*A)/sin(A) for t in [0,1] and A in
// [0,pi). The implementation chebyshevRatio computes this function. The
// implementation chebyshevRatios computes the pair {f(1-t,A), f(t,A)}, which
// is useful for spherical linear interpolation.
//
// TODO (upstream): The evaluation for A near 0 or pi needs to be more robust.
// For A near 0, sin(t*A)/sin(A) has a removable singularity. Use the idea in
// RAEFGC for an approximation to remove the singularity. For A near pi, the
// singularity is not removable, so some approximation must be used such as
// those found in ChebyshevRatioEstimate.ts.

import { GTE_C_PI } from './Constants.js';
import { logError } from './Logger.js';

// The angle must be in [0,pi).
export function chebyshevRatio(t: number, angle: number): number {
    if (angle > 0) {
        if (angle < GTE_C_PI) {
            // The angle A is in (0,pi).
            return Math.sin(t * angle) / Math.sin(angle);
        }
    } else if (angle === 0) {
        // The angle A is 0. Using l'Hospital's rule,
        // lim_{A->0} sin(t*A)/sin(A) = lim_{A->0} t*cos(t*A)/cos(A) = t.
        return t;
    }

    // The angle A is not in [0,pi).
    logError('Invalid angle.');
}

// The angle extracted from cosAngle is in [0,pi).
export function chebyshevRatioUsingCosAngle(t: number, cosAngle: number): number {
    if (cosAngle < 1) {
        if (cosAngle > -1) {
            // The angle A is in (0,pi).
            const angle = Math.acos(cosAngle);
            return Math.sin(t * angle) / Math.sin(angle);
        } else {
            // The angle A is pi.
            logError('Invalid angle.');
        }
    } else {
        // The angle A is 0. Using l'Hospital's rule,
        // lim_{A->0} sin(t*A)/sin(A) = lim_{A->0} t*cos(t*A)/cos(A) = t.
        return t;
    }
}

// The angle must be in [0,pi). Although it is possible to compute
// invSin = 1/sin(angle) and perform two multiplications for f[0] and
// f[1], the resulting ratios typically do not match those from
// chebyshevRatio. Therefore, two divisions are performed in this
// function to ensure the resulting ratios are the same.
export function chebyshevRatios(t: number, angle: number): [number, number] {
    if (angle > 0) {
        if (angle < GTE_C_PI) {
            // The angle A is in (0,pi).
            const sinAngle = Math.sin(angle);
            return [
                Math.sin((1 - t) * angle) / sinAngle,
                Math.sin(t * angle) / sinAngle
            ];
        }
    } else if (angle === 0) {
        // The angle A is 0. Using l'Hospital's rule,
        // lim_{A->0} sin(t*A)/sin(A) = lim_{A->0} t*cos(t*A)/cos(A) = t.
        return [1 - t, t];
    }

    // The angle A is not in [0,pi).
    logError('Invalid angle.');
}

// The angle extracted from cosAngle is in [0,pi). Although it is possible
// to compute invSin = 1/sin(angle) and perform two multiplications for
// f[0] and f[1], the resulting ratios typically do not match those from
// chebyshevRatioUsingCosAngle. Therefore, two divisions are performed in
// this function to ensure the resulting ratios are the same.
export function chebyshevRatiosUsingCosAngle(t: number, cosAngle: number): [number, number] {
    if (cosAngle < 1) {
        if (cosAngle > -1) {
            // The angle A is in (0,pi).
            const angle = Math.acos(cosAngle);
            const sinAngle = Math.sin(angle);
            return [
                Math.sin((1 - t) * angle) / sinAngle,
                Math.sin(t * angle) / sinAngle
            ];
        } else {
            // The angle A is pi.
            logError('Invalid angle.');
        }
    } else {
        // The angle A is 0. Using l'Hospital's rule,
        // lim_{A->0} sin(t*A)/sin(A) = lim_{A->0} t*cos(t*A)/cos(A) = t.
        return [1 - t, t];
    }
}
