// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OdeMidpoint.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Upstream is templated on TVector so that a solver can be created with
// Vector<N,Real> (dimension known at compile time) or GVector<Real>
// (dimension known at run time). The port has a single Vector class whose
// dimension is a run-time value, so the solvers are concrete classes derived
// from OdeSolver<Vector>.

import { OdeSolver, type OdeFunction } from './OdeSolver';
import { Vector, add, mul } from './Vector';

export class OdeMidpoint extends OdeSolver<Vector> {
    constructor(tDelta: number, F: OdeFunction<Vector>) {
        super(tDelta, F);
    }

    // Estimate x(t + tDelta) from x(t) using dx/dt = F(t,x). The estimate is
    // a new Vector, so xIn is never modified.
    update(tIn: number, xIn: Vector): { tOut: number; xOut: Vector } {
        // Compute the first step.
        const halfTDelta = 0.5 * this.mTDelta;
        let fVector = this.mFunction(tIn, xIn);
        const xTemp = add(xIn, mul(halfTDelta, fVector));

        // Compute the second step.
        const halfT = tIn + halfTDelta;
        fVector = this.mFunction(halfT, xTemp);
        const tOut = tIn + this.mTDelta;
        const xOut = add(xIn, mul(this.mTDelta, fVector));
        return { tOut, xOut };
    }
}
