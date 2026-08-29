// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OdeSolver.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The differential equation is dx/dt = F(t,x). The TVector type parameter
// allows you to create solvers with Vector when the dimension is known at
// construction time. (Upstream is templated on both Real and TVector; the
// port keeps TVector generic and maps Real to number.)

// The port of 'std::function<TVector(Real, TVector const&)>'.
export type OdeFunction<TVector> = (t: number, x: TVector) => TVector;

// Abstract base class.
export abstract class OdeSolver<TVector> {
    protected mTDelta: number;
    protected mFunction: OdeFunction<TVector>;

    protected constructor(tDelta: number, F: OdeFunction<TVector>) {
        this.mTDelta = tDelta;
        this.mFunction = F;
    }

    // Member access.
    setTDelta(tDelta: number): void {
        this.mTDelta = tDelta;
    }

    getTDelta(): number {
        return this.mTDelta;
    }

    // Estimate x(t + tDelta) from x(t) using dx/dt = F(t,x). Upstream
    // returns tOut and xOut through output references (and allows xIn and
    // xOut to be the same object); per PORTING.md the port returns them as
    // an object literal with named fields.
    abstract update(tIn: number, xIn: TVector): { tOut: number; xOut: TVector };
}
