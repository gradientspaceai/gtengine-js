// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NURBSCircle.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithm for representing a circle as a NURBS curve is described in
//   https://www.geometrictools.com/Documentation/NURBSCircleSphere.pdf
// The implementations are related to the document as shown next.
//   NURBSQuarterCircleDegree2 implements equation (9)
//   NURBSQuarterCircleDegree4 implements equation (10)
//   NURBSHalfCircleDegree3 implements equation (12)
//   NURBSFullCircleDegree3 implements Section 2.3
//   NURBSCircularArcDegree2 implements Section 2.4
//
// Port notes: upstream has five class templates in this header; all five live
// in this file per the one-file-per-header rule. Each derives from the ported
// NURBSCurve, whose first constructor argument is the runtime dimension (2
// here). Upstream passes null pointers for the controls and weights and then
// fills the protected arrays; the port omits the optional arguments and does
// the same.

import { Arc2 } from './Arc2';
import { BasisFunctionInput, UniqueKnot } from './BasisFunction';
import { NURBSCurve } from './NURBSCurve';
import { Vector, div, dot, sub } from './Vector';
import { dotPerp, perp } from './Vector2';

// Local helper for the 2-tuple literals used by the control points.
function vector2(x: number, y: number): Vector {
    const v = new Vector(2);
    v.values[0] = x;
    v.values[1] = y;
    return v;
}

export class NURBSQuarterCircleDegree2 extends NURBSCurve {
    // Construction. The quarter circle is x^2 + y^2 = 1 for x >= 0 and
    // y >= 0. The direction of traversal is counterclockwise as u increases
    // from 0 to 1.
    constructor() {
        super(2, new BasisFunctionInput(3, 2));

        const sqrt2 = Math.sqrt(2);

        this.mWeights[0] = sqrt2;
        this.mWeights[1] = 1;
        this.mWeights[2] = sqrt2;

        this.mControls[0] = vector2(1, 0);
        this.mControls[1] = vector2(1, 1);
        this.mControls[2] = vector2(0, 1);
    }
}

export class NURBSQuarterCircleDegree4 extends NURBSCurve {
    // Construction. The quarter circle is x^2 + y^2 = 1 for x >= 0 and
    // y >= 0. The direction of traversal is counterclockwise as u increases
    // from 0 to 1.
    constructor() {
        super(2, new BasisFunctionInput(5, 4));

        const sqrt2 = Math.sqrt(2);

        this.mWeights[0] = 1;
        this.mWeights[1] = 1;
        this.mWeights[2] = 2 * sqrt2 / 3;
        this.mWeights[3] = 1;
        this.mWeights[4] = 1;

        const x1 = 1;
        const y1 = 0.5 / sqrt2;
        const x2 = 1 - sqrt2 / 8;
        this.mControls[0] = vector2(1, 0);
        this.mControls[1] = vector2(x1, y1);
        this.mControls[2] = vector2(x2, x2);
        this.mControls[3] = vector2(y1, x1);
        this.mControls[4] = vector2(0, 1);
    }
}

export class NURBSHalfCircleDegree3 extends NURBSCurve {
    // Construction. The half circle is x^2 + y^2 = 1 for y >= 0. The
    // direction of traversal is counterclockwise as u increases from 0 to 1,
    // from (1,0) through (0,1) to (-1,0).
    //
    // Upstream bug (documentation): the comment in NURBSCircle.h says the
    // half circle is the part with x >= 0. The control points
    // (1,0), (1,2), (-1,2), (-1,0) generate the upper half, y >= 0. Only the
    // comment is wrong; the construction is correct and is preserved.
    constructor() {
        super(2, new BasisFunctionInput(4, 3));

        const oneThird = 1 / 3;

        this.mWeights[0] = 1;
        this.mWeights[1] = oneThird;
        this.mWeights[2] = oneThird;
        this.mWeights[3] = 1;

        this.mControls[0] = vector2(1, 0);
        this.mControls[1] = vector2(1, 2);
        this.mControls[2] = vector2(-1, 2);
        this.mControls[3] = vector2(-1, 0);
    }
}

export class NURBSFullCircleDegree3 extends NURBSCurve {
    // Construction. The full circle is x^2 + y^2 = 1. The direction of
    // traversal is counterclockwise as u increases from 0 to 1.
    constructor() {
        super(2, NURBSFullCircleDegree3.createBasisFunctionInput());

        const oneThird = 1 / 3;

        this.mWeights[0] = 1;
        this.mWeights[1] = oneThird;
        this.mWeights[2] = oneThird;
        this.mWeights[3] = 1;
        this.mWeights[4] = oneThird;
        this.mWeights[5] = oneThird;
        this.mWeights[6] = 1;

        this.mControls[0] = vector2(1, 0);
        this.mControls[1] = vector2(1, 2);
        this.mControls[2] = vector2(-1, 2);
        this.mControls[3] = vector2(-1, 0);
        this.mControls[4] = vector2(-1, -2);
        this.mControls[5] = vector2(1, -2);
        this.mControls[6] = vector2(1, 0);
    }

    private static createBasisFunctionInput(): BasisFunctionInput {
        // We need knots (0,0,0,0,1/2,1/2,1/2,1,1,1,1).
        const input = new BasisFunctionInput();
        input.numControls = 7;
        input.degree = 3;
        input.uniform = true;
        input.periodic = false;
        input.numUniqueKnots = 3;
        input.uniqueKnots = [
            new UniqueKnot(0, 4),
            new UniqueKnot(0.5, 3),
            new UniqueKnot(1, 4)
        ];
        return input;
    }
}

export class NURBSCircularArcDegree2 extends NURBSCurve {
    // Construction from an arc. The arc endpoints are arc.end[0] and
    // arc.end[1], traversed counterclockwise about arc.center.
    constructor(arc: Arc2) {
        super(2, new BasisFunctionInput(3, 2));

        const p0 = div(sub(arc.end[0], arc.center), arc.radius);
        const p2 = div(sub(arc.end[1], arc.center), arc.radius);
        const p1 = div(perp(sub(p2, p0)), dotPerp(p0, p2));

        this.mWeights[0] = Math.sqrt(2 * (dot(p1, p1) - 1) / (1 - dot(p0, p2)));
        this.mWeights[1] = 1;
        this.mWeights[2] = this.mWeights[0];

        this.mControls[0] = vector2(
            arc.center.values[0] + arc.radius * p0.values[0],
            arc.center.values[1] + arc.radius * p0.values[1]);
        this.mControls[1] = vector2(
            arc.center.values[0] + arc.radius * p1.values[0],
            arc.center.values[1] + arc.radius * p1.values[1]);
        this.mControls[2] = vector2(
            arc.center.values[0] + arc.radius * p2.values[0],
            arc.center.values[1] + arc.radius * p2.values[1]);
    }
}
