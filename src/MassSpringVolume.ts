// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MassSpringVolume.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An S-by-R-by-C lattice of masses connected by springs. The state is
// advanced by the fourth-order Runge-Kutta solver of ParticleSystem.
//
// Port notes:
// - Upstream 'MassSpringVolume<N, Real>' becomes a runtime dimension,
//   matching the ParticleSystem port: 'new MassSpringVolume(dimension,
//   numSlices, numRows, numCols, step)'.
// - Upstream overloads SetMass/SetPosition/SetVelocity (and the getters) on
//   (s, r, c) in addition to the linear-index versions inherited from
//   ParticleSystem. TypeScript cannot override a base method with a
//   different arity, so the lattice-indexed forms are named setMassAt,
//   setPositionAt, setVelocityAt, getMassAt, getPositionAt and
//   getVelocityAt; the linear-index forms of ParticleSystem are unchanged.
// - The class is abstract in name only upstream (Acceleration is implemented
//   here), so the port is a concrete class; derive from it to supply nonzero
//   external forces by overriding externalAcceleration.

import { Vector, add, length as vectorLength, mul, sub } from './Vector.js';
import { ParticleSystem } from './ParticleSystem.js';

export class MassSpringVolume extends ParticleSystem {
    protected mNumSlices: number;
    protected mNumRows: number;
    protected mNumCols: number;

    // The spring constants and resting lengths for the springs in the slice,
    // row and column directions. All arrays are of size S*R*C, with the
    // entries for the last slice/row/column unused.
    protected mConstantS: number[];
    protected mLengthS: number[];
    protected mConstantR: number[];
    protected mLengthR: number[];
    protected mConstantC: number[];
    protected mLengthC: number[];

    // Construction. This class represents an SxRxC array of masses lying in a
    // volume and connected by an array of springs. The masses are indexed by
    // mass[s][r][c] for 0 <= s < S, 0 <= r < R and 0 <= c < C. The mass at
    // interior position X[s][r][c] is connected by springs to the masses at
    // positions X[s][r-1][c], X[s][r+1][c], X[s][r][c-1], X[s][r][c+1],
    // X[s-1][r][c] and X[s+1][r][c]. Boundary masses have springs connecting
    // them to the obvious neighbors (a "face" mass has 5 neighbors, an "edge"
    // mass has 4 neighbors, a "corner" mass has 3 neighbors). The masses are
    // arranged in lexicographical order: position[c+C*(r+R*s)] = X[s][r][c].
    // The other arrays are stored similarly.
    constructor(dimension: number, numSlices: number, numRows: number,
        numCols: number, step: number) {
        super(dimension, numSlices * numRows * numCols, step);
        this.mNumSlices = numSlices;
        this.mNumRows = numRows;
        this.mNumCols = numCols;
        const numMasses = numSlices * numRows * numCols;
        this.mConstantS = new Array<number>(numMasses).fill(0);
        this.mLengthS = new Array<number>(numMasses).fill(0);
        this.mConstantR = new Array<number>(numMasses).fill(0);
        this.mLengthR = new Array<number>(numMasses).fill(0);
        this.mConstantC = new Array<number>(numMasses).fill(0);
        this.mLengthC = new Array<number>(numMasses).fill(0);
    }

    // Member access.
    getNumSlices(): number {
        return this.mNumSlices;
    }

    getNumRows(): number {
        return this.mNumRows;
    }

    getNumCols(): number {
        return this.mNumCols;
    }

    setMassAt(s: number, r: number, c: number, mass: number): void {
        this.setMass(this.getIndex(s, r, c), mass);
    }

    setPositionAt(s: number, r: number, c: number, position: Vector): void {
        this.setPosition(this.getIndex(s, r, c), position);
    }

    setVelocityAt(s: number, r: number, c: number, velocity: Vector): void {
        this.setVelocity(this.getIndex(s, r, c), velocity);
    }

    getMassAt(s: number, r: number, c: number): number {
        return this.getMass(this.getIndex(s, r, c));
    }

    // Upstream returns a 'Vector<N,Real> const&'. The port returns the live
    // internal object, so it must not be modified by the caller.
    getPositionAt(s: number, r: number, c: number): Vector {
        return this.getPosition(this.getIndex(s, r, c));
    }

    getVelocityAt(s: number, r: number, c: number): Vector {
        return this.getVelocity(this.getIndex(s, r, c));
    }

    // Each interior mass at (s,r,c) has 6 adjacent springs. Face masses have
    // only 5 neighbors, edge masses have only 4 neighbors, and corner masses
    // have only 3 neighbors. Each mass provides access to 3 adjacent springs
    // at (s,r,c+1), (s,r+1,c) and (s+1,r,c). The face, edge and corner masses
    // provide access to only an appropriate subset of these. The caller is
    // responsible for ensuring the validity of the (s,r,c) inputs.

    // to (s+1,r,c)
    setConstantS(s: number, r: number, c: number, constant: number): void {
        this.mConstantS[this.getIndex(s, r, c)] = constant;
    }

    // to (s+1,r,c)
    setLengthS(s: number, r: number, c: number, length: number): void {
        this.mLengthS[this.getIndex(s, r, c)] = length;
    }

    // to (s,r+1,c)
    setConstantR(s: number, r: number, c: number, constant: number): void {
        this.mConstantR[this.getIndex(s, r, c)] = constant;
    }

    // to (s,r+1,c)
    setLengthR(s: number, r: number, c: number, length: number): void {
        this.mLengthR[this.getIndex(s, r, c)] = length;
    }

    // to (s,r,c+1)
    setConstantC(s: number, r: number, c: number, constant: number): void {
        this.mConstantC[this.getIndex(s, r, c)] = constant;
    }

    // to (s,r,c+1)
    setLengthC(s: number, r: number, c: number, length: number): void {
        this.mLengthC[this.getIndex(s, r, c)] = length;
    }

    getConstantS(s: number, r: number, c: number): number {
        return this.mConstantS[this.getIndex(s, r, c)];
    }

    getLengthS(s: number, r: number, c: number): number {
        return this.mLengthS[this.getIndex(s, r, c)];
    }

    getConstantR(s: number, r: number, c: number): number {
        return this.mConstantR[this.getIndex(s, r, c)];
    }

    getLengthR(s: number, r: number, c: number): number {
        return this.mLengthR[this.getIndex(s, r, c)];
    }

    getConstantC(s: number, r: number, c: number): number {
        return this.mConstantC[this.getIndex(s, r, c)];
    }

    getLengthC(s: number, r: number, c: number): number {
        return this.mLengthC[this.getIndex(s, r, c)];
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
        // face, edge and corner points of the volume of masses must be
        // handled separately, because each has fewer than six springs
        // attached to it.

        let acceleration = this.externalAcceleration(i, time, position, velocity);
        let diff: Vector, force: Vector, ratio: number;
        let prev: number, next: number;

        const { s, r, c } = this.getCoordinates(i);

        if (s > 0) {
            prev = i - this.mNumRows * this.mNumCols;  // previous s-neighbor
            diff = sub(position[prev], position[i]);
            ratio = this.getLengthS(s - 1, r, c) / vectorLength(diff);
            force = mul(this.getConstantS(s - 1, r, c) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        if (s < this.mNumSlices - 1) {
            next = i + this.mNumRows * this.mNumCols;  // next s-neighbor
            diff = sub(position[next], position[i]);
            ratio = this.getLengthS(s, r, c) / vectorLength(diff);
            force = mul(this.getConstantS(s, r, c) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        if (r > 0) {
            prev = i - this.mNumCols;  // previous r-neighbor
            diff = sub(position[prev], position[i]);
            ratio = this.getLengthR(s, r - 1, c) / vectorLength(diff);
            force = mul(this.getConstantR(s, r - 1, c) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        if (r < this.mNumRows - 1) {
            next = i + this.mNumCols;  // next r-neighbor
            diff = sub(position[next], position[i]);
            ratio = this.getLengthR(s, r, c) / vectorLength(diff);
            force = mul(this.getConstantR(s, r, c) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        if (c > 0) {
            prev = i - 1;  // previous c-neighbor
            diff = sub(position[prev], position[i]);
            ratio = this.getLengthC(s, r, c - 1) / vectorLength(diff);
            force = mul(this.getConstantC(s, r, c - 1) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        if (c < this.mNumCols - 1) {
            next = i + 1;  // next c-neighbor
            diff = sub(position[next], position[i]);
            ratio = this.getLengthC(s, r, c) / vectorLength(diff);
            force = mul(this.getConstantC(s, r, c) * (1 - ratio), diff);
            acceleration = add(acceleration, mul(this.mInvMass[i], force));
        }

        return acceleration;
    }

    protected getIndex(s: number, r: number, c: number): number {
        return c + this.mNumCols * (r + this.mNumRows * s);
    }

    // The port of the upstream GetCoordinates(i, s, r, c) output-reference
    // function; it returns the lattice coordinates as an object.
    protected getCoordinates(i: number): { s: number, r: number, c: number } {
        const c = i % this.mNumCols;
        const j = (i - c) / this.mNumCols;
        const r = j % this.mNumRows;
        const s = Math.trunc(j / this.mNumRows);
        return { s, r, c };
    }
}
