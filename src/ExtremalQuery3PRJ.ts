// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ExtremalQuery3PRJ.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The extremal queries for convex objects are based on projection of the
// vertices onto the specified line. This is the brute-force O(n) query.
//
// Port notes:
// - 'GetExtremeVertices(direction, positiveDirection, negativeDirection)' has
//   output reference parameters upstream; the port returns the
//   'ExtremalQuery3Result' object declared by the base class.
// - The deleted copy constructor and copy assignment have no port; TypeScript
//   objects are reference types and are never implicitly copied.
// - 'std::numeric_limits<Real>::max()' is Number.MAX_VALUE, matching the
//   binary64 instantiation of the C++ template.
// - Upstream dereferences the shared_ptr vertex pool as
//   'vertexPool.get()->at(i)'; the port uses Polyhedron3.getVertices(), which
//   is the same dereference and throws on an invalid polyhedron (upstream has
//   undefined behavior there). The base-class constructor already requires a
//   valid polyhedron, so this cannot be reached through a constructed query.

import { ExtremalQuery3 } from './ExtremalQuery3';
import type { ExtremalQuery3Result } from './ExtremalQuery3';
import { Polyhedron3 } from './Polyhedron3';
import { Vector, dot, sub } from './Vector';

export class ExtremalQuery3PRJ extends ExtremalQuery3 {
    private mCentroid: Vector;

    // The caller must ensure that the input polyhedron is convex.
    constructor(polytope: Polyhedron3) {
        super(polytope);
        this.mCentroid = this.mPolytope.computeVertexAverage();
    }

    // Compute the extreme vertices in the specified direction and return the
    // indices of the vertices in the polyhedron vertex array. The projections
    // are taken relative to the centroid, which changes neither the argmin nor
    // the argmax: subtracting the centroid shifts every projection by the same
    // constant Dot(direction, centroid).
    getExtremeVertices(direction: Vector): ExtremalQuery3Result {
        let minValue = Number.MAX_VALUE;
        let maxValue = -minValue;
        let negativeDirection = -1;
        let positiveDirection = -1;

        const vertexPool = this.mPolytope.getVertices();
        for (const i of this.mPolytope.getUniqueIndices()) {
            const diff = sub(vertexPool[i], this.mCentroid);
            const d = dot(direction, diff);
            // The comparisons are strict, so among vertices with equal
            // projections the first one visited wins. The unique indices are
            // visited in ascending order (the std::set order upstream), so a
            // tie resolves to the smallest vertex index.
            if (d < minValue) {
                negativeDirection = i;
                minValue = d;
            }
            if (d > maxValue) {
                positiveDirection = i;
                maxValue = d;
            }
        }

        return {
            positiveDirection: positiveDirection,
            negativeDirection: negativeDirection
        };
    }
}
