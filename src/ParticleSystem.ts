// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ParticleSystem.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A system of particles, each with a mass, a position and a velocity. The
// state is advanced by a fourth-order Runge-Kutta solver applied to the
// second-order system x" = F/m, where the acceleration F/m is supplied by a
// derived class through the acceleration(...) callback.
//
// Port notes:
// - The upstream template parameter N (the dimension of the vectors) becomes
//   a runtime dimension passed to the constructor, matching the Vector port:
//   'ParticleSystem<N, Real>(numParticles, step)' -> 'new ParticleSystem(
//   dimension, numParticles, step)'.
// - The upstream protected nested struct 'Temporary' is exported as
//   'ParticleSystemTemporary' because src/index.ts star-exports every file.
// - C++ vector assignment copies; the port clones explicitly wherever
//   upstream relies on that value semantics, so that no two entries of the
//   position/velocity/temporary arrays ever alias one another.

import { Vector, add, mul } from './Vector';

// Temporary storage for the Runge-Kutta differential equation solver (the
// port of ParticleSystem::Temporary).
export interface ParticleSystemTemporary {
    d1: Vector;
    d2: Vector;
    d3: Vector;
    d4: Vector;
}

export abstract class ParticleSystem {
    protected mDimension: number;
    protected mNumParticles: number;
    protected mMass: number[];
    protected mInvMass: number[];
    protected mPosition: Vector[];
    protected mVelocity: Vector[];
    protected mStep: number;
    protected mHalfStep: number;
    protected mSixthStep: number;

    // Temporary storage for the Runge-Kutta differential equation solver.
    protected mPTmp: Vector[];
    protected mVTmp: Vector[];
    protected mPAllTmp: ParticleSystemTemporary[];
    protected mVAllTmp: ParticleSystemTemporary[];

    // Construction. If a particle is to be immovable, set its mass to
    // Number.MAX_VALUE (upstream std::numeric_limits<Real>::max()).
    //
    // The masses, positions and velocities are all zero-initialized. Note
    // that a mass of zero has an inverse mass of zero, so a particle whose
    // mass is never assigned is immovable.
    protected constructor(dimension: number, numParticles: number, step: number) {
        this.mDimension = dimension;
        this.mNumParticles = numParticles;
        this.mMass = new Array<number>(numParticles).fill(0);
        this.mInvMass = new Array<number>(numParticles).fill(0);
        this.mPosition = new Array<Vector>(numParticles);
        this.mVelocity = new Array<Vector>(numParticles);
        this.mStep = step;
        this.mHalfStep = step / 2;
        this.mSixthStep = step / 6;
        this.mPTmp = new Array<Vector>(numParticles);
        this.mVTmp = new Array<Vector>(numParticles);
        this.mPAllTmp = new Array<ParticleSystemTemporary>(numParticles);
        this.mVAllTmp = new Array<ParticleSystemTemporary>(numParticles);
        for (let i = 0; i < numParticles; ++i) {
            this.mPosition[i] = new Vector(dimension);
            this.mVelocity[i] = new Vector(dimension);
            this.mPTmp[i] = new Vector(dimension);
            this.mVTmp[i] = new Vector(dimension);
            this.mPAllTmp[i] = ParticleSystem.newTemporary(dimension);
            this.mVAllTmp[i] = ParticleSystem.newTemporary(dimension);
        }
    }

    private static newTemporary(dimension: number): ParticleSystemTemporary {
        return {
            d1: new Vector(dimension),
            d2: new Vector(dimension),
            d3: new Vector(dimension),
            d4: new Vector(dimension)
        };
    }

    // Member access.
    getDimension(): number {
        return this.mDimension;
    }

    getNumParticles(): number {
        return this.mNumParticles;
    }

    setMass(i: number, mass: number): void {
        if (0 < mass && mass < Number.MAX_VALUE) {
            this.mMass[i] = mass;
            this.mInvMass[i] = 1 / mass;
        }
        else {
            this.mMass[i] = Number.MAX_VALUE;
            this.mInvMass[i] = 0;
        }
    }

    setPosition(i: number, position: Vector): void {
        this.mPosition[i] = position.clone();
    }

    setVelocity(i: number, velocity: Vector): void {
        this.mVelocity[i] = velocity.clone();
    }

    setStep(step: number): void {
        this.mStep = step;
        this.mHalfStep = this.mStep / 2;
        this.mSixthStep = this.mStep / 6;
    }

    getMass(i: number): number {
        return this.mMass[i];
    }

    // Upstream returns a 'Vector<N,Real> const&'. The port returns the live
    // internal object, so it must not be modified by the caller; clone it if
    // the value has to outlive the next update(...) call.
    getPosition(i: number): Vector {
        return this.mPosition[i];
    }

    getVelocity(i: number): Vector {
        return this.mVelocity[i];
    }

    getStep(): number {
        return this.mStep;
    }

    // Update the particle positions based on current time and particle
    // state. The acceleration(...) function is called in this update for each
    // (movable) particle. Derived classes may override this to perform
    // pre-update and/or post-update semantics.
    update(time: number): void {
        // Runge-Kutta fourth-order solver.
        const halfTime = time + this.mHalfStep;
        const fullTime = time + this.mStep;

        // Compute the first step.
        let i: number;
        for (i = 0; i < this.mNumParticles; ++i) {
            if (this.mInvMass[i] > 0) {
                this.mPAllTmp[i].d1 = this.mVelocity[i].clone();
                this.mVAllTmp[i].d1 = this.acceleration(i, time, this.mPosition,
                    this.mVelocity);
            }
        }
        for (i = 0; i < this.mNumParticles; ++i) {
            if (this.mInvMass[i] > 0) {
                this.mPTmp[i] = add(this.mPosition[i],
                    mul(this.mHalfStep, this.mPAllTmp[i].d1));
                this.mVTmp[i] = add(this.mVelocity[i],
                    mul(this.mHalfStep, this.mVAllTmp[i].d1));
            }
            else {
                this.mPTmp[i] = this.mPosition[i].clone();
                this.mVTmp[i].makeZero();
            }
        }

        // Compute the second step.
        for (i = 0; i < this.mNumParticles; ++i) {
            if (this.mInvMass[i] > 0) {
                this.mPAllTmp[i].d2 = this.mVTmp[i].clone();
                this.mVAllTmp[i].d2 = this.acceleration(i, halfTime, this.mPTmp,
                    this.mVTmp);
            }
        }
        for (i = 0; i < this.mNumParticles; ++i) {
            if (this.mInvMass[i] > 0) {
                this.mPTmp[i] = add(this.mPosition[i],
                    mul(this.mHalfStep, this.mPAllTmp[i].d2));
                this.mVTmp[i] = add(this.mVelocity[i],
                    mul(this.mHalfStep, this.mVAllTmp[i].d2));
            }
            else {
                this.mPTmp[i] = this.mPosition[i].clone();
                this.mVTmp[i].makeZero();
            }
        }

        // Compute the third step.
        for (i = 0; i < this.mNumParticles; ++i) {
            if (this.mInvMass[i] > 0) {
                this.mPAllTmp[i].d3 = this.mVTmp[i].clone();
                this.mVAllTmp[i].d3 = this.acceleration(i, halfTime, this.mPTmp,
                    this.mVTmp);
            }
        }
        for (i = 0; i < this.mNumParticles; ++i) {
            if (this.mInvMass[i] > 0) {
                this.mPTmp[i] = add(this.mPosition[i],
                    mul(this.mStep, this.mPAllTmp[i].d3));
                this.mVTmp[i] = add(this.mVelocity[i],
                    mul(this.mStep, this.mVAllTmp[i].d3));
            }
            else {
                this.mPTmp[i] = this.mPosition[i].clone();
                this.mVTmp[i].makeZero();
            }
        }

        // Compute the fourth step.
        for (i = 0; i < this.mNumParticles; ++i) {
            if (this.mInvMass[i] > 0) {
                this.mPAllTmp[i].d4 = this.mVTmp[i].clone();
                this.mVAllTmp[i].d4 = this.acceleration(i, fullTime, this.mPTmp,
                    this.mVTmp);
            }
        }
        for (i = 0; i < this.mNumParticles; ++i) {
            if (this.mInvMass[i] > 0) {
                const p = this.mPAllTmp[i];
                this.mPosition[i] = add(this.mPosition[i], mul(this.mSixthStep,
                    add(add(p.d1, mul(2, add(p.d2, p.d3))), p.d4)));

                const v = this.mVAllTmp[i];
                this.mVelocity[i] = add(this.mVelocity[i], mul(this.mSixthStep,
                    add(add(v.d1, mul(2, add(v.d2, v.d3))), v.d4)));
            }
        }
    }

    // Callback for acceleration (the ODE solver uses x" = F/m) applied to
    // particle i. The positions and velocities are not necessarily mPosition
    // and mVelocity, because the ODE solver evaluates the impulse function at
    // intermediate positions.
    protected abstract acceleration(i: number, time: number,
        position: readonly Vector[], velocity: readonly Vector[]): Vector;
}
