// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the arc to be a 1-dimensional object.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrLine2Arc2TI and
// IntrLine2Arc2FI. Upstream calls the single-argument Arc2::Contains, which
// the port names 'containsOnCircle'.

import type { Arc2 } from './Arc2.js';
import type { FIQuery } from './FIQuery.js';
import { IntrLine2Circle2FI } from './IntrLine2Circle2.js';
import { Hypersphere } from './Hypersphere.js';
import type { Line2 } from './Line.js';
import { Vector } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrLine2Arc2TI.test.
export interface IntrLine2Arc2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine2Arc2TIResult {
    return { intersect: false };
}

// The result of IntrLine2Arc2FI.find.
export interface IntrLine2Arc2FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine2Arc2FIResult(): IntrLine2Arc2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

// Test-intersection query for a line and an arc in 2D.
export class IntrLine2Arc2TI implements
    TIQuery<Line2, Arc2, IntrLine2Arc2TIResult> {

    test(line: Line2, arc: Arc2): IntrLine2Arc2TIResult {
        const result = defaultTIResult();
        const laQuery = new IntrLine2Arc2FI();
        const laResult = laQuery.find(line, arc);
        result.intersect = laResult.intersect;
        return result;
    }
}

// Find-intersection query for a line and an arc in 2D.
export class IntrLine2Arc2FI implements
    FIQuery<Line2, Arc2, IntrLine2Arc2FIResult> {

    find(line: Line2, arc: Arc2): IntrLine2Arc2FIResult {
        const result = defaultIntrLine2Arc2FIResult();

        const lcQuery = new IntrLine2Circle2FI();
        const circle = Hypersphere.fromCenterRadius(arc.center, arc.radius);
        const lcResult = lcQuery.find(line, circle);
        if (lcResult.intersect) {
            // Test whether line-circle intersections are on the arc.
            result.numIntersections = 0;
            for (let i = 0; i < lcResult.numIntersections; ++i) {
                if (arc.containsOnCircle(lcResult.point[i])) {
                    result.intersect = true;
                    result.parameter[result.numIntersections] =
                        lcResult.parameter[i];
                    result.point[result.numIntersections] =
                        lcResult.point[i].clone();
                    ++result.numIntersections;
                }
            }
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}
