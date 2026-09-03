// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) EdgeKey.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An ordered edge has (V[0], V[1]) = (v0, v1). An unordered edge has
// (V[0], V[1]) = (min(V[0],V[1]), max(V[0],V[1])).
//
// Port notes: upstream is the class template EdgeKey<Ordered> deriving from
// FeatureKey<2, Ordered>. Following the FeatureKey port, the compile-time
// Ordered becomes the first constructor argument and is stored in the
// readonly 'ordered' property of the base class. The two Initialize
// overloads selected by std::enable_if become a single private method that
// branches on 'ordered'. The upstream default constructor initializes to
// invalid indices (-1, -1); the port keeps that behavior when the vertex
// indices are omitted.

import { FeatureKey } from './FeatureKey.js';

export class EdgeKey extends FeatureKey {
    // With v0 and v1 omitted, the key is initialized to invalid indices.
    constructor(ordered: boolean, v0?: number, v1?: number) {
        super(2, ordered);
        if (v0 === undefined || v1 === undefined) {
            this.V[0] = -1;
            this.V[1] = -1;
        } else {
            this.initialize(v0, v1);
        }
    }

    private initialize(v0: number, v1: number): void {
        if (this.ordered) {
            this.V[0] = v0;
            this.V[1] = v1;
        } else {
            if (v0 < v1) {
                // v0 is minimum
                this.V[0] = v0;
                this.V[1] = v1;
            } else {
                // v1 is minimum
                this.V[0] = v1;
                this.V[1] = v0;
            }
        }
    }
}
