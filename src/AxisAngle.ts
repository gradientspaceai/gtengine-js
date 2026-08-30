// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AxisAngle.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Axis-angle representation for N = 3 or N = 4. When N = 4, the axis must be
// a vector of the form (x,y,z,0) [affine representation of the 3-tuple
// direction].
//
// Port notes: upstream 'AxisAngle<N, Real>' has a compile-time dimension with
// 'static_assert(N == 3 || N == 4)'. The port checks the runtime size of the
// axis vector instead and throws for any other size. The default constructor
// (which upstream instantiates per N) produces the 3D zero axis with zero
// angle. As upstream copies the axis into the member, the constructor clones
// the input vector.

import { logAssert } from './Logger';
import { Vector } from './Vector';

export class AxisAngle {
    axis: Vector;
    angle: number;

    constructor(axis: Vector = Vector.zero(3), angle: number = 0) {
        logAssert(axis.size === 3 || axis.size === 4,
            'Dimension must be 3 or 4.');
        this.axis = axis.clone();
        this.angle = angle;
    }
}
