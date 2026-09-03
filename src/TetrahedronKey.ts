// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TetrahedronKey.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An ordered tetrahedron has V[0] = min(v0, v1, v2, v3). Let {u1, u2, u3} be
// the set of inputs excluding the one assigned to V[0] and define
// V[1] = min(u1, u2, u3). Choose (V[1], V[2], V[3]) to be a permutation of
// (u1, u2, u3) so that the final storage is one of
//   (v0, v1, v2, v3), (v0, v2, v3, v1), (v0, v3, v1, v2)
//   (v1, v3, v2, v0), (v1, v2, v0, v3), (v1, v0, v3, v2)
//   (v2, v3, v0, v1), (v2, v0, v1, v3), (v2, v1, v3, v0)
//   (v3, v1, v0, v2), (v3, v0, v2, v1), (v3, v2, v1, v0)
// The idea is that if v0 corresponds to (1,0,0,0), v1 corresponds to
// (0,1,0,0), v2 corresponds to (0,0,1,0), and v3 corresponds to (0,0,0,1),
// the ordering (v0, v1, v2, v3) corresponds to the 4x4 identity matrix I;
// the rows are the specified 4-tuples. The permutation
// (V[0], V[1], V[2], V[3]) induces a permutation of the rows of the identity
// matrix to form a permutation matrix P with det(P) = 1 = det(I).
//
// An unordered tetrahedron stores a permutation of (v0, v1, v2, v3) so that
// V[0] < V[1] < V[2] < V[3].
//
// Port notes: upstream is the class template TetrahedronKey<Ordered>
// deriving from FeatureKey<4, Ordered>. Following the FeatureKey port, the
// compile-time Ordered becomes the first constructor argument and is stored
// in the readonly 'ordered' property of the base class. The Initialize and
// Permute overloads selected by std::enable_if become private methods
// dispatched on 'ordered'. The upstream default constructor initializes to
// invalid indices (-1, -1, -1, -1); the port keeps that behavior when the
// vertex indices are omitted. The function-local static of
// GetOppositeFace() becomes a module-level frozen table.

import { FeatureKey } from './FeatureKey.js';

// Indexing for the vertices of the triangle opposite a vertex. The triangle
// opposite vertex j is
//   <oppositeFace[j][0], oppositeFace[j][1], oppositeFace[j][2]>
// and is listed in counterclockwise order when viewed from outside the
// tetrahedron.
const sOppositeFace: ReadonlyArray<ReadonlyArray<number>> = Object.freeze([
    Object.freeze([1, 2, 3]),
    Object.freeze([0, 3, 2]),
    Object.freeze([0, 1, 3]),
    Object.freeze([0, 2, 1])
]);

export class TetrahedronKey extends FeatureKey {
    // With the vertex indices omitted, the key is initialized to invalid
    // indices.
    constructor(ordered: boolean, v0?: number, v1?: number, v2?: number,
        v3?: number) {
        super(4, ordered);
        if (v0 === undefined || v1 === undefined || v2 === undefined ||
            v3 === undefined) {
            this.V[0] = -1;
            this.V[1] = -1;
            this.V[2] = -1;
            this.V[3] = -1;
        } else if (ordered) {
            this.initializeOrdered(v0, v1, v2, v3);
        } else {
            this.initializeUnordered(v0, v1, v2, v3);
        }
    }

    // Indexing for the vertices of the triangle opposite a vertex; see the
    // comments on sOppositeFace above.
    static getOppositeFace(): ReadonlyArray<ReadonlyArray<number>> {
        return sOppositeFace;
    }

    private initializeOrdered(v0: number, v1: number, v2: number,
        v3: number): void {
        let imin = 0;
        this.V[0] = v0;
        if (v1 < this.V[0]) {
            this.V[0] = v1;
            imin = 1;
        }
        if (v2 < this.V[0]) {
            this.V[0] = v2;
            imin = 2;
        }
        if (v3 < this.V[0]) {
            this.V[0] = v3;
            imin = 3;
        }

        if (imin === 0) {
            this.permute(v1, v2, v3);
        } else if (imin === 1) {
            this.permute(v0, v3, v2);
        } else if (imin === 2) {
            this.permute(v0, v1, v3);
        } else {  // imin === 3
            this.permute(v0, v2, v1);
        }
    }

    private initializeUnordered(v0: number, v1: number, v2: number,
        v3: number): void {
        this.V[0] = v0;
        this.V[1] = v1;
        this.V[2] = v2;
        this.V[3] = v3;
        // std::sort on the 4 indices; the numeric comparator is required
        // because the default JavaScript sort compares as strings.
        this.V.sort((a, b) => a - b);
    }

    // Once V[0] is determined, create a permutation (V[1], V[2], V[3]) so
    // that (V[0], V[1], V[2], V[3]) is a permutation of (v0, v1, v2, v3)
    // that corresponds to the identity matrix as mentioned in the comments
    // at the beginning of this file.
    private permute(u0: number, u1: number, u2: number): void {
        if (u0 < u1) {
            if (u0 < u2) {
                // u0 is minimum
                this.V[1] = u0;
                this.V[2] = u1;
                this.V[3] = u2;
            } else {
                // u2 is minimum
                this.V[1] = u2;
                this.V[2] = u0;
                this.V[3] = u1;
            }
        } else {
            if (u1 < u2) {
                // u1 is minimum
                this.V[1] = u1;
                this.V[2] = u2;
                this.V[3] = u0;
            } else {
                // u2 is minimum
                this.V[1] = u2;
                this.V[2] = u0;
                this.V[3] = u1;
            }
        }
    }
}
