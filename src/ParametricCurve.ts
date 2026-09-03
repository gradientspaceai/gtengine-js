// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ParametricCurve.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Abstract base class for a parameterized curve X(t), where t is the
// parameter in [tmin,tmax] and X is an N-tuple position. The first
// constructor form is for single-segment curves. The second constructor form
// is for multiple-segment curves. The times must be strictly increasing.
//
// Port notes (these set the precedent for the curve family: BezierCurve,
// BSplineCurve, NURBSCurve, ...):
// - Upstream 'template <int32_t N, typename Real>' becomes a runtime
//   dimension, matching the Vector port and ParametricSurface: the protected
//   constructor takes the dimension as its first argument.
// - The two protected C++ constructors become TypeScript constructor
//   overloads distinguished by the type of the third argument:
//   'new Derived(dimension, tmin, tmax)' and
//   'new Derived(dimension, numSegments, times)'.
// - 'operator bool' has no TS equivalent; it is ported as isConstructed().
// - The 'jet' output of evaluate() is a 'Vector[]' of at least SUP_ORDER
//   entries; use createJet() to allocate it (upstream sizes it from the
//   template parameter N via std::array).
// - The 'Evaluate(Real t, uint32_t order, Real* values)' overload is a C++
//   reinterpret_cast of a flat scalar buffer to Vector<N,Real>*; it has no
//   meaning in TS and is omitted.
// - SubdivideByTime/SubdivideByLength write into a caller-supplied pointer;
//   the port returns a freshly allocated 'Vector[]' instead.
// - mSegmentLength/mAccumulatedLength are 'mutable' in C++ (lazily filled by
//   the const member getLength). TS has no const methods, so they are plain
//   fields updated by getLength(), which preserves the lazy behavior.

import { Integration } from './Integration.js';
import { RootsBisection } from './RootsBisection.js';
import { Vector, length as vectorLength, normalize } from './Vector.js';

// The port of std::lower_bound on the strictly increasing time array: the
// index of the first entry >= value, or times.length when there is none.
function lowerBound(times: readonly number[], value: number): number {
    let first = 0;
    let count = times.length;
    while (count > 0) {
        const step = Math.floor(count / 2);
        const mid = first + step;
        if (times[mid] < value) {
            first = mid + 1;
            count -= step + 1;
        }
        else {
            count = step;
        }
    }
    return first;
}

export abstract class ParametricCurve {
    // The number of entries a 'jet' array must have: position, first
    // derivative, second derivative, third derivative.
    static readonly SUP_ORDER = 4;

    protected static readonly DEFAULT_ROMBERG_ORDER = 8;
    protected static readonly DEFAULT_MAX_BISECTIONS = 1024;

    protected mDimension: number;
    protected mTime: number[];
    protected mSegmentLength: number[];
    protected mAccumulatedLength: number[];
    protected mRombergOrder: number;
    protected mMaxBisections: number;
    protected mConstructed: boolean;

    // Single-segment curve on [tmin, tmax].
    protected constructor(dimension: number, tmin: number, tmax: number);

    // Multiple-segment curve; 'times' holds numSegments+1 strictly
    // increasing values.
    protected constructor(dimension: number, numSegments: number,
        times: readonly number[]);

    protected constructor(dimension: number, arg1: number,
        arg2: number | readonly number[]) {
        this.mDimension = dimension;
        this.mRombergOrder = ParametricCurve.DEFAULT_ROMBERG_ORDER;
        this.mMaxBisections = ParametricCurve.DEFAULT_MAX_BISECTIONS;
        this.mConstructed = false;

        if (typeof arg2 === 'number') {
            // ParametricCurve(tmin, tmax)
            this.mTime = [arg1, arg2];
            this.mSegmentLength = [0];
            this.mAccumulatedLength = [0];
        }
        else {
            // ParametricCurve(numSegments, times)
            const numSegments = arg1;
            this.mTime = new Array<number>(numSegments + 1).fill(0);
            for (let i = 0; i <= numSegments; ++i) {
                this.mTime[i] = arg2[i];
            }
            this.mSegmentLength = new Array<number>(numSegments).fill(0);
            this.mAccumulatedLength = new Array<number>(numSegments).fill(0);
        }
    }

    // To validate construction, create an object as shown:
    //     const curve = new DerivedClassCurve(parameters);
    //     if (!curve.isConstructed()) { /* constructor failed */ }
    isConstructed(): boolean {
        return this.mConstructed;
    }

    // Member access.
    getDimension(): number {
        return this.mDimension;
    }

    getTMin(): number {
        return this.mTime[0];
    }

    getTMax(): number {
        return this.mTime[this.mTime.length - 1];
    }

    getNumSegments(): number {
        return this.mSegmentLength.length;
    }

    // Upstream returns 'Real const*' into the internal array; the port
    // returns the live array, which must not be modified by the caller.
    getTimes(): readonly number[] {
        return this.mTime;
    }

    // This function applies only when the first constructor form is used
    // (two times rather than a sequence of three or more times).
    setTimeInterval(tmin: number, tmax: number): void {
        if (this.mTime.length === 2) {
            this.mTime[0] = tmin;
            this.mTime[1] = tmax;
        }
    }

    // Parameters used in getLength(...), getTotalLength() and getTime(...).

    // The default value is 8.
    setRombergOrder(order: number): void {
        this.mRombergOrder = Math.max(order, 1);
    }

    // The default value is 1024.
    setMaxBisections(maxBisections: number): void {
        this.mMaxBisections = Math.max(maxBisections, 1);
    }

    // Allocate storage for a jet of this curve: SUP_ORDER zero vectors of
    // the curve dimension. This is the port of upstream's
    // 'std::array<Vector<N, Real>, SUP_ORDER> jet{}'.
    createJet(): Vector[] {
        const jet = new Array<Vector>(ParametricCurve.SUP_ORDER);
        for (let i = 0; i < ParametricCurve.SUP_ORDER; ++i) {
            jet[i] = new Vector(this.mDimension);
        }
        return jet;
    }

    // Evaluation of the curve. The function supports derivative calculation
    // through order 3; that is, order <= 3 is required. If you want only the
    // position, pass in order of 0. If you want the position and first
    // derivative, pass in order of 1, and so on. The output array 'jet' must
    // have enough storage to support the maximum order. The values are
    // ordered as: position, first derivative, second derivative, third
    // derivative.
    abstract evaluate(t: number, order: number, jet: Vector[]): void;

    // Differential geometric quantities. Upstream's jet is an array of value
    // objects, so an implementation of Evaluate that writes a stored vector
    // into a jet slot writes a copy. TS arrays hold references, so these
    // accessors clone the jet entry before returning (and before
    // normalizing) to keep the C++ value semantics.
    getPosition(t: number): Vector {
        const jet = this.createJet();
        this.evaluate(t, 0, jet);
        return jet[0].clone();
    }

    getTangent(t: number): Vector {
        const jet = this.createJet();
        this.evaluate(t, 1, jet);
        const tangent = jet[1].clone();
        normalize(tangent);
        return tangent;
    }

    getSpeed(t: number): number {
        const jet = this.createJet();
        this.evaluate(t, 1, jet);
        return vectorLength(jet[1]);
    }

    getLength(t0: number, t1: number): number {
        const speed = (t: number): number => this.getSpeed(t);

        if (this.mSegmentLength[0] === 0) {
            // Lazy initialization of lengths of segments.
            const numSegments = this.mSegmentLength.length;
            let accumulated = 0;
            for (let i = 0, ip1 = 1; i < numSegments; ++i, ++ip1) {
                this.mSegmentLength[i] = Integration.romberg(this.mRombergOrder,
                    this.mTime[i], this.mTime[ip1], speed);
                accumulated += this.mSegmentLength[i];
                this.mAccumulatedLength[i] = accumulated;
            }
        }

        t0 = Math.max(t0, this.getTMin());
        t1 = Math.min(t1, this.getTMax());
        const index0 = lowerBound(this.mTime, t0);
        const index1 = lowerBound(this.mTime, t1);

        let length: number;
        if (index0 < index1) {
            length = 0;
            if (t0 < this.mTime[index0]) {
                length += Integration.romberg(this.mRombergOrder, t0,
                    this.mTime[index0], speed);
            }

            let isup: number;
            if (t1 < this.mTime[index1]) {
                length += Integration.romberg(this.mRombergOrder,
                    this.mTime[index1 - 1], t1, speed);
                isup = index1 - 1;
            }
            else {
                isup = index1;
            }
            for (let i = index0; i < isup; ++i) {
                length += this.mSegmentLength[i];
            }
        }
        else {
            length = Integration.romberg(this.mRombergOrder, t0, t1, speed);
        }
        return length;
    }

    getTotalLength(): number {
        const last = this.mAccumulatedLength[this.mAccumulatedLength.length - 1];
        if (last === 0) {
            // Lazy evaluation of the accumulated length array.
            return this.getLength(this.mTime[0], this.mTime[this.mTime.length - 1]);
        }

        return last;
    }

    // Inverse mapping of s = Length(t) given by t = Length^{-1}(s). The
    // inverse length function generally cannot be written in closed form, in
    // which case it is not directly computable. Instead, we can specify s and
    // estimate the root t for F(t) = Length(t) - s. The derivative is
    // F'(t) = Speed(t) >= 0, so F(t) is nondecreasing. To be robust, we use
    // bisection to locate the root, although it is possible to use a hybrid
    // of Newton's method and bisection. For details, see the document
    // https://www.geometrictools.com/Documentation/MovingAlongCurveSpecifiedSpeed.pdf
    getTime(length: number): number {
        const tmin = this.mTime[0];
        const tmax = this.mTime[this.mTime.length - 1];

        if (length > 0) {
            const totalLength = this.getTotalLength();
            if (length < totalLength) {
                const F = (t: number): number => {
                    return Integration.romberg(this.mRombergOrder, tmin, t,
                        (z: number) => this.getSpeed(z)) - length;
                };

                // We know that F(tmin) < 0 and F(tmax) > 0, which allows us
                // to use bisection. Rather than bisect the entire interval,
                // let's narrow it down with a reasonable initial guess.
                const ratio = length / totalLength;
                const omratio = 1 - ratio;
                let tmid = omratio * tmin + ratio * tmax;
                const fmid = F(tmid);
                if (fmid > 0) {
                    tmid = RootsBisection.find(F, tmin, tmid, -1, 1,
                        this.mMaxBisections).root;
                }
                else if (fmid < 0) {
                    tmid = RootsBisection.find(F, tmid, tmax, -1, 1,
                        this.mMaxBisections).root;
                }
                return tmid;
            }
            else {
                return tmax;
            }
        }
        else {
            return tmin;
        }
    }

    // Compute a subset of curve points according to the specified attribute.
    // The input 'numPoints' must be two or larger. Upstream writes into a
    // caller-supplied array; the port returns the points.
    subdivideByTime(numPoints: number): Vector[] {
        const tmin = this.mTime[0];
        const tmax = this.mTime[this.mTime.length - 1];
        const delta = (tmax - tmin) / (numPoints - 1);
        const points = new Array<Vector>(numPoints);
        for (let i = 0; i < numPoints; ++i) {
            const t = tmin + delta * i;
            points[i] = this.getPosition(t);
        }
        return points;
    }

    subdivideByLength(numPoints: number): Vector[] {
        const delta = this.getTotalLength() / (numPoints - 1);
        const points = new Array<Vector>(numPoints);
        for (let i = 0; i < numPoints; ++i) {
            const length = delta * i;
            const t = this.getTime(length);
            points[i] = this.getPosition(t);
        }
        return points;
    }
}
