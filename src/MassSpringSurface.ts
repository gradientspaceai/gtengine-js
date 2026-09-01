// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MassSpringSurface.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An R-by-C lattice of masses connected by springs. The state is advanced by
// the fourth-order Runge-Kutta solver of ParticleSystem.
//
// Port notes (matching the MassSpringVolume/MassSpringCurve precedents):
// - Upstream 'MassSpringSurface<N, Real>' becomes a runtime dimension:
//   'new MassSpringSurface(dimension, numRows, numCols, step)'.
// - Upstream overloads SetMass/SetPosition/SetVelocity (and the getters) on
//   (r, c) in addition to the linear-index versions inherited from
//   ParticleSystem. TypeScript cannot override a base method with a different
//   arity, so the lattice-indexed forms are named setMassAt, setPositionAt,
//   setVelocityAt, getMassAt, getPositionAt and getVelocityAt; the
//   linear-index forms of ParticleSystem are unchanged.
// - The upstream GetCoordinates(i, r, c) output-reference function returns the
//   lattice coordinates as an object.
// - The class is abstract in name only upstream (Acceleration is implemented
//   here), so the port is a concrete class; derive from it to supply nonzero
//   external forces by overriding externalAcceleration.

import { Vector, add, length as vectorLength, mul, sub } from './Vector';
import { ParticleSystem } from './ParticleSystem';

export class MassSpringSurface extends ParticleSystem {
    protected mNumRows: number;
    protected mNumCols: number;

    // The spring constants and resting lengths for the springs in the row and
    // column directions. Both arrays are of size R*C, with the entries for the
    // last row/column unused.
    protected mConstantR: number[];
    protected mLengthR: number[];
    protected mConstantC: number[];
    protected mLengthC: number[];

    // Construction. This class represents an RxC array of masses lying on a
    // surface and connected by an array of springs. The masses are indexed by
    // mass[r][c] for 0 <= r < R and 0 <= c < C. The mass at interior position
    // X[r][c] is connected by springs to the masses at positions X[r-1][c],
    // X[r+1][c], X[r][c-1] and X[r][c+1]. Boundary masses have springs
    // connecting them to the obvious neighbors (an "edge" mass has 3
    // neighbors, a "corner" mass has 2 neighbors). The masses are arranged in
    // row-major order: position[c+C*r] = X[r][c] for 0 <= r < R and
    // 0 <= c < C. The other arrays are stored similarly.
    constructor(dimension: number, numRows: number, numCols: number, step: number) {
        super(dimension, numRows * numCols, step);
        this.mNumRows = numRows;
        this.mNumCols = numCols;
        const numMasses = numRows * numCols;
        this.mConstantR = new Array<number>(numMasses).fill(0);
        this.mLengthR = new Array<number>(numMasses).fill(0);
        this.mConstantC = new Array<number>(numMasses).fill(0);
        this.mLengthC = new Array<number>(numMasses).fill(0);
    }

    // Member access.
    getNumRows(): number {
        return this.mNumRows;
    }

    getNumCols(): number {
        return this.mNumCols;
    }

    setMassAt(r: number, c: number, mass: number): void {
        this.setMass(this.getIndex(r, c), mass);
    }

    setPositionAt(r: number, c: number, position: Vector): void {
        this.setPosition(this.getIndex(r, c), position);
    }

    setVelocityAt(r: number, c: number, velocity: Vector): void {
        this.setVelocity(this.getIndex(r, c), velocity);
    }

    getMassAt(r: number, c: number): number {
        return this.getMass(this.getIndex(r, c));
    }

    // Upstream returns a 'Vector<N,Real> const&'. The port returns the live
    // internal object, so it must not be modified by the caller.
    getPositionAt(r: number, c: number): Vector {
        return this.getPosition(this.getIndex(r, c));
    }

    getVelocityAt(r: number, c: number): Vector {
        return this.getVelocity(this.getIndex(r, c));
    }

    // The interior mass at (r,c) has springs to the left, right, bottom and
    // top. Edge masses have only three neighbors and corner masses have only
    // two neighbors. The mass at (r,c) provides access to the springs
    // connecting to locations (r,c+1) and (r+1,c). Edge and corner masses
    // provide access to only a subset of these. The caller is responsible for
    // ensuring the validity of the (r,c) inputs.

    // to (r+1,c)
    setConstantR(r: number, c: number, constant: number): void {
        this.mConstantR[this.getIndex(r, c)] = constant;
    }

    // to (r+1,c)
    setLengthR(r: number, c: number, length: number): void {
        this.mLengthR[this.getIndex(r, c)] = length;
    }

    // to (r,c+1)
    setConstantC(r: number, c: number, constant: number): void {
        this.mConstantC[this.getIndex(r, c)] = constant;
    }

    // to (r,c+1)
    setLengthC(r: number, c: number, length: number): void {
        this.mLengthC[this.getIndex(r, c)] = length;
    }

    getConstantR(r: number, c: number): number {
        return this.mConstantR[this.getIndex(r, c)];
    }

    getLengthR(r: number, c: number): number {
        return this.mLengthR[this.getIndex(r, c)];
    }

    getConstantC(r: number, c: number): number {
        return this.mConstantC[this.getIndex(r, c)];
    }

    getLengthC(r: number, c: number): number {
        return this.mLengthC[this.getIndex(r, c)];
    }

    // The default external force is zero. Derive a class from this one to
    // provide nonzero external forces such as gravity, wind, friction, and so
    // on. This function is called by acceleration(...) to compute the impulse
    // F/m generated by the external force F.
    externalAcceleration(_i: number, _time: number,
        _position: readonly Vector[], _velocity: readonly Vector[]): Vector {
        return Vector.zero(this.mDimension);
    }

    // Callback for acceleration (the ODE solver uses x" = F/m) applied to
    // particle i. The positions and velocities are not necessarily mPosition
    // and mVelocity, because the ODE solver evaluates the impulse function at
    // intermediate positions.
    protected acceleration(i: number, time: number,
        position: readonly Vector[], velocity: readonly Vector[]): Vector {
        // Compute spring forces on position X[i]. The positions are not
        // necessarily mPosition, because the RK4 solver in ParticleSystem
        // evaluates the acceleration function at intermediate positions. The
        // edge and corner points of the surface of masses must be handled
        // separately, because each has fewer than four springs attached to it.

        let acceleration = this.externalAcceleration(i, time, position, velocity);
        let diff: Vector, force: Vector, ratio: number;
        let prev: number, next: number;

        const { r, c } = this.getCoordinates(i);

        if (r > 0) {
            prev = i - this.mNumCols;  // index to previous row-neighbor
            diff = sub(position[prev], position[i]);
            ratio = this.getLengthR(r - 1, c) / vectorLength(diff);
            force = mul(this.getConstantR(r - 1, c) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        if (r < this.mNumRows - 1) {
            next = i + this.mNumCols;  // index to next row-neighbor
            diff = sub(position[next], position[i]);
            ratio = this.getLengthR(r, c) / vectorLength(diff);
            force = mul(this.getConstantR(r, c) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        if (c > 0) {
            prev = i - 1;  // index to previous col-neighbor
            diff = sub(position[prev], position[i]);
            ratio = this.getLengthC(r, c - 1) / vectorLength(diff);
            force = mul(this.getConstantC(r, c - 1) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        if (c < this.mNumCols - 1) {
            next = i + 1;  // index to next col-neighbor
            diff = sub(position[next], position[i]);
            ratio = this.getLengthC(r, c) / vectorLength(diff);
            force = mul(this.getConstantC(r, c) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        return acceleration;
    }

    protected getIndex(r: number, c: number): number {
        return c + this.mNumCols * r;
    }

    // The port of the upstream GetCoordinates(i, r, c) output-reference
    // function; it returns the lattice coordinates as an object.
    protected getCoordinates(i: number): { r: number, c: number } {
        const c = i % this.mNumCols;
        const r = Math.trunc(i / this.mNumCols);
        return { r, c };
    }
}
