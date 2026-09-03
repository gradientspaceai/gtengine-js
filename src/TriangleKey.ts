// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TriangleKey.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An ordered triangle has V[0] = min(v0, v1, v2). Choose (V[0], V[1], V[2])
// to be a permutation of (v0, v1, v2) so that the final storage is one of
// (v0, v1, v2), (v1, v2, v0) or (v2, v0, v1). The idea is that if v0
// corresponds to (1,0,0), v1 corresponds to (0,1,0), and v2 corresponds to
// (0,0,1), the ordering (v0, v1, v2) corresponds to the 3x3 identity matrix
// I; the rows are the specified 3-tuples. The permutation (V[0], V[1], V[2])
// induces a permutation of the rows of the identity matrix to form a
// permutation matrix P with det(P) = 1 = det(I).
//
// An unordered triangle stores a permutation of (v0, v1, v2) so that
// V[0] < V[1] < V[2].
//
// Port notes: upstream is the class template TriangleKey<Ordered> deriving
// from FeatureKey<3, Ordered>. Following the FeatureKey port, the
// compile-time Ordered becomes the first constructor argument and is stored
// in the readonly 'ordered' property of the base class. The two Initialize
// overloads selected by std::enable_if become two private methods dispatched
// on 'ordered'. The upstream default constructor initializes to invalid
// indices (-1, -1, -1); the port keeps that behavior when the vertex indices
// are omitted.

import { FeatureKey } from './FeatureKey.js';

export class TriangleKey extends FeatureKey {
    // With v0, v1 and v2 omitted, the key is initialized to invalid indices.
    constructor(ordered: boolean, v0?: number, v1?: number, v2?: number) {
        super(3, ordered);
        if (v0 === undefined || v1 === undefined || v2 === undefined) {
            this.V[0] = -1;
            this.V[1] = -1;
            this.V[2] = -1;
        } else if (ordered) {
            this.initializeOrdered(v0, v1, v2);
        } else {
            this.initializeUnordered(v0, v1, v2);
        }
    }

    private initializeOrdered(v0: number, v1: number, v2: number): void {
        if (v0 < v1) {
            if (v0 < v2) {
                // v0 is minimum
                this.V[0] = v0;
                this.V[1] = v1;
                this.V[2] = v2;
            } else {
                // v2 is minimum
                this.V[0] = v2;
                this.V[1] = v0;
                this.V[2] = v1;
            }
        } else {
            if (v1 < v2) {
                // v1 is minimum
                this.V[0] = v1;
                this.V[1] = v2;
                this.V[2] = v0;
            } else {
                // v2 is minimum
                this.V[0] = v2;
                this.V[1] = v0;
                this.V[2] = v1;
            }
        }
    }

    private initializeUnordered(v0: number, v1: number, v2: number): void {
        if (v0 < v1) {
            if (v0 < v2) {
                // v0 is minimum
                this.V[0] = v0;
                this.V[1] = Math.min(v1, v2);
                this.V[2] = Math.max(v1, v2);
            } else {
                // v2 is minimum
                this.V[0] = v2;
                this.V[1] = Math.min(v0, v1);
                this.V[2] = Math.max(v0, v1);
            }
        } else {
            if (v1 < v2) {
                // v1 is minimum
                this.V[0] = v1;
                this.V[1] = Math.min(v2, v0);
                this.V[2] = Math.max(v2, v0);
            } else {
                // v2 is minimum
                this.V[0] = v2;
                this.V[1] = Math.min(v0, v1);
                this.V[2] = Math.max(v0, v1);
            }
        }
    }
}
