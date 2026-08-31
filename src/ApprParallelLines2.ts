// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprParallelLines2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fit of two parallel lines to points that presumably are
// clustered on the lines. The algorithm is described in
//   https://www.geometrictools.com/Documentation/FitParallelLinesToPoints2D.pdf
//
// The two lines are C - r*U + s*V and C + r*U + s*V for s in R, where V is
// the unit-length common direction of the lines, U is the unit-length
// normal (perpendicular to V) and r is the half-distance between the lines.
// The direction is parameterized as V = (gamma,sigma) with
// gamma^2 + sigma^2 = 1, so that U = (-sigma,gamma). The moments Z[p][q] of
// the mean-centered samples determine an error function F(sigma,gamma) whose
// roots are the candidate minimizers.
//
// Port notes: upstream's 'Fit' has output reference parameters (C, V,
// radius); per PORTING.md they become the fields of the returned object
// literal. The Real-typed constants mR0..mR6 (which exist upstream only to
// avoid repeated construction costs for rational types) are numeric literals
// here.

import { Polynomial1 } from './Polynomial1';
import { RootsPolynomial } from './RootsPolynomial';
import { Vector, dot } from './Vector';

export interface ApprParallelLines2Result {
    // The point equidistant from the two lines, with only a U-component
    // relative to the mean of the samples.
    center: Vector;

    // The unit-length common direction (gamma,sigma) of the two lines.
    direction: Vector;

    // The half-distance between the two lines.
    radius: number;
}

// The port of C++ 'Polynomial1<Real>::operator[]' reads on an index that can
// exceed the polynomial degree. Upstream's arithmetic operators eliminate
// leading zero coefficients, so the extraction loops below can index past
// the end of the std::vector of coefficients (undefined behavior in C++)
// when the leading coefficients of a product happen to cancel. The port
// treats the missing coefficients as zero, which is the mathematically
// correct value.
function coefficient(p: Polynomial1, i: number): number {
    return i <= p.getDegree() ? p.get(i) : 0;
}

// The moments Z[p][q] = average(x^p * y^q) of the mean-centered samples.
class ZValues {
    Z20: number = 0; Z11: number = 0; Z02: number = 0;
    Z30: number = 0; Z21: number = 0; Z12: number = 0; Z03: number = 0;
    Z40: number = 0; Z31: number = 0; Z22: number = 0; Z13: number = 0;
    Z04: number = 0;

    constructor(P: readonly Vector[]) {
        const invN = 1 / P.length;
        for (const sample of P) {
            const x = sample.values[0];
            const y = sample.values[1];
            const xx = x * x;
            const xy = x * y;
            const yy = y * y;
            const xxx = xx * x;
            const xxy = xy * x;
            const xyy = xy * y;
            const yyy = yy * y;
            const xxxx = xxx * x;
            const xxxy = xxx * y;
            const xxyy = xx * yy;
            const xyyy = yyy * x;
            const yyyy = yyy * y;
            this.Z20 += xx;
            this.Z11 += xy;
            this.Z02 += yy;
            this.Z30 += xxx;
            this.Z21 += xxy;
            this.Z12 += xyy;
            this.Z03 += yyy;
            this.Z40 += xxxx;
            this.Z31 += xxxy;
            this.Z22 += xxyy;
            this.Z13 += xyyy;
            this.Z04 += yyyy;
        }
        this.Z20 *= invN;
        this.Z11 *= invN;
        this.Z02 *= invN;
        this.Z30 *= invN;
        this.Z21 *= invN;
        this.Z12 *= invN;
        this.Z03 *= invN;
        this.Z40 *= invN;
        this.Z31 *= invN;
        this.Z22 *= invN;
        this.Z13 *= invN;
        this.Z04 *= invN;
    }
}

// The running minimum of the error function and the parameters that attain
// it. Upstream passes these as five Real reference parameters to
// UpdateParameters.
interface Minimum {
    sigma: number;
    gamma: number;
    k: number;
    rSqr: number;
    error: number;
}

export class ApprParallelLines2 {
    fit(P: readonly Vector[], maxIterations: number): ApprParallelLines2Result {
        // Compute the average of the samples.
        const n = P.length;
        const invN = 1 / n;
        const PAdjust: Vector[] = new Array<Vector>(n);
        const A = new Vector(2);
        for (let i = 0; i < n; ++i) {
            PAdjust[i] = P[i].clone();
            A.values[0] += PAdjust[i].values[0];
            A.values[1] += PAdjust[i].values[1];
        }
        A.values[0] *= invN;
        A.values[1] *= invN;

        // Subtract the average from the samples so that the replacement
        // points have zero average.
        for (const sample of PAdjust) {
            sample.values[0] -= A.values[0];
            sample.values[1] -= A.values[1];
        }

        // Compute the Zpq terms.
        const data = new ZValues(PAdjust);

        // Compute F(sigma,gamma) = f0(sigma) + gamma * f1(sigma).
        const { f0, f1 } = ApprParallelLines2.computeF(data);
        const freduced0 = new Polynomial1(4);
        const freduced1 = new Polynomial1(3);
        for (let i = 0; i <= 4; ++i) {
            freduced0.set(i, coefficient(f0, 2 * i));
        }
        for (let i = 0; i <= 3; ++i) {
            freduced1.set(i, coefficient(f1, 2 * i + 1));
        }

        // Evaluate the error function at any (sigma,gamma). Choose (0,1) so
        // that we do not have to process a root sigma = 0 later.
        const minK = data.Z03 / (2 * data.Z02);
        const minKSqr = minK * minK;
        const minimum: Minimum = {
            sigma: 0,
            gamma: 1,
            k: minK,
            rSqr: minKSqr + data.Z02,
            error: data.Z04 - 4 * minK * data.Z03
                + (4 * minKSqr - data.Z02) * data.Z02
        };

        if (f1.notEquals(Polynomial1.fromCoefficients([0]))) {
            const sigmaSqrPoly = Polynomial1.fromCoefficients([0, 0, 1]);
            const f0Sqr = f0.mul(f0);
            const f1Sqr = f1.mul(f1);
            const h = sigmaSqrPoly.mul(f1Sqr).add(f0Sqr.sub(f1Sqr));
            const hreduced = new Polynomial1(8);
            for (let i = 0; i <= 8; ++i) {
                hreduced.set(i, coefficient(h, 2 * i));
            }

            const roots = RootsPolynomial.find(8, hreduced.getCoefficients(),
                maxIterations);
            for (const sigmaSqr of roots) {
                // Upstream tests only 'sigmaSqr > 0'. The direction
                // V = (gamma,sigma) is unit length, so a root with
                // sigmaSqr > 1 has gamma^2 = 1 - sigmaSqr < 0 and cannot
                // correspond to a direction. Such a root is a spurious
                // solution of h (which is F=0 after squaring away gamma) and
                // updateParameters can report a negative error for it, in
                // which case it displaces the true minimum and yields a
                // non-unit direction and a NaN radius. See the upstream bug
                // notes for this file.
                if (sigmaSqr > 0 && sigmaSqr <= 1) {
                    const sigma = Math.sqrt(sigmaSqr);
                    const gamma = -freduced0.evaluate(sigmaSqr)
                        / (sigma * freduced1.evaluate(sigmaSqr));
                    ApprParallelLines2.updateParameters(data, sigma, sigmaSqr,
                        gamma, minimum);
                }
            }
        }
        else {
            const hreduced = new Polynomial1(4);
            for (let i = 0; i <= 4; ++i) {
                hreduced.set(i, coefficient(f0, 2 * i));
            }

            const roots = RootsPolynomial.find(4, hreduced.getCoefficients(),
                maxIterations);
            for (const sigmaSqr of roots) {
                // See the comment on the same test above: a root with
                // sigmaSqr > 1 cannot come from a unit-length direction.
                if (sigmaSqr > 0 && sigmaSqr <= 1) {
                    const sigma = Math.sqrt(sigmaSqr);

                    // When f1 is identically zero, F(sigma,gamma) = f0(sigma)
                    // is independent of gamma, so gamma is determined only by
                    // the unit-length constraint gamma^2 + sigma^2 = 1 for the
                    // line direction V = (gamma,sigma). Upstream computes
                    // 'gamma = sqrt(sigma)', which does not satisfy the
                    // constraint; the port uses the constraint. See the
                    // upstream bug notes for this file.
                    let gamma = Math.sqrt(1 - sigmaSqr);
                    ApprParallelLines2.updateParameters(data, sigma, sigmaSqr,
                        gamma, minimum);

                    gamma = -gamma;
                    ApprParallelLines2.updateParameters(data, sigma, sigmaSqr,
                        gamma, minimum);
                }
            }
        }

        // Compute the minimizers V, C and radius. The center minK*U must
        // have A added to it because the inputs P had A subtracted from
        // them. The addition no longer guarantees that Dot(V,C) = 0, so the
        // V-component of A + minK*U is projected out so that the returned
        // center has only a U-component.
        const V = Vector.fromArray([minimum.gamma, minimum.sigma]);
        const C = Vector.fromArray([
            A.values[0] + minimum.k * -minimum.sigma,
            A.values[1] + minimum.k * minimum.gamma
        ]);
        const dotCV = dot(C, V);
        C.values[0] -= dotCV * V.values[0];
        C.values[1] -= dotCV * V.values[1];

        return { center: C, direction: V, radius: Math.sqrt(minimum.rSqr) };
    }

    // Given two polynomials A0+gamma*B0 and A1+gamma*B1, the product is
    // [A0*A1+(1-sigma^2)*B0*B1] + gamma*[A0*B1+B0*A1] = A2+gamma*B2, where
    // gamma^2 = 1-sigma^2.
    private static computeProduct(A0: Polynomial1, B0: Polynomial1,
        A1: Polynomial1, B1: Polynomial1): { A2: Polynomial1, B2: Polynomial1 } {
        const gammaSqr = Polynomial1.fromCoefficients([1, 0, -1]);
        return {
            A2: A0.mul(A1).add(gammaSqr.mul(B0).mul(B1)),
            B2: A0.mul(B1).add(B0.mul(A1))
        };
    }

    private static computeF(data: ZValues): { f0: Polynomial1, f1: Polynomial1 } {
        // Compute the apq and bpq terms, where
        // S[p][q](sigma,gamma) = a[p][q](sigma) + gamma * b[p][q](sigma) is
        // the moment average(u^p * v^q) with u = Dot(P,U) and v = Dot(P,V).
        const a11 = new Polynomial1(2);
        a11.set(0, data.Z11);
        a11.set(2, -2 * data.Z11);

        const b11 = new Polynomial1(1);
        b11.set(1, data.Z02 - data.Z20);

        const a20 = new Polynomial1(2);
        a20.set(0, data.Z02);
        a20.set(2, data.Z20 - data.Z02);

        const b20 = new Polynomial1(1);
        b20.set(1, -2 * data.Z11);

        const a30 = new Polynomial1(3);
        // Upstream sets a30[1] = -3, which is missing the factor Z12 that
        // the equivalent evaluation in updateParameters uses (A30 =
        // -sigma * (3*Z12 + (Z30 - 3*Z12) * sigma^2)). The port uses the
        // consistent coefficient. See the upstream bug notes for this file.
        a30.set(1, -3 * data.Z12);
        a30.set(3, 3 * data.Z12 - data.Z30);

        const b30 = new Polynomial1(2);
        b30.set(0, data.Z03);
        b30.set(2, 3 * data.Z21 - data.Z03);

        const a21 = new Polynomial1(3);
        a21.set(1, data.Z03 - 2 * data.Z21);
        a21.set(3, 3 * data.Z21 - data.Z03);

        const b21 = new Polynomial1(2);
        b21.set(0, data.Z12);
        b21.set(2, data.Z30 - 3 * data.Z12);

        const a40 = new Polynomial1(4);
        a40.set(0, data.Z04);
        a40.set(2, 6 * data.Z22 - 2 * data.Z04);
        a40.set(4, data.Z40 - 6 * data.Z22 + data.Z04);

        const b40 = new Polynomial1(3);
        b40.set(1, -4 * data.Z13);
        b40.set(3, 4 * (data.Z13 - data.Z31));

        const a31 = new Polynomial1(4);
        a31.set(0, data.Z13);
        a31.set(2, 3 * data.Z31 - 5 * data.Z13);
        a31.set(4, 4 * (data.Z13 - data.Z31));

        const b31 = new Polynomial1(3);
        b31.set(1, data.Z04 - 3 * data.Z22);
        b31.set(3, 6 * data.Z22 - data.Z40 - data.Z04);

        // Compute S20^2 = c0 + gamma*d0.
        const { A2: c0, B2: d0 } = ApprParallelLines2.computeProduct(a20, b20, a20, b20);

        // Compute S31 * S20^2 = c1 + gamma*d1.
        const { A2: c1, B2: d1 } = ApprParallelLines2.computeProduct(a31, b31, c0, d0);

        // Compute S21 * S20 = c2 + gamma*d2.
        const { A2: c2, B2: d2 } = ApprParallelLines2.computeProduct(a21, b21, a20, b20);

        // Compute S30 * (S21 * S20) = c3 + gamma*d3.
        const { A2: c3, B2: d3 } = ApprParallelLines2.computeProduct(a30, b30, c2, d2);

        // Compute S30 * S11 = c4 + gamma*d4.
        const { A2: c4, B2: d4 } = ApprParallelLines2.computeProduct(a30, b30, a11, b11);

        // Compute S30 * (S30 * S11) = c5 + gamma*d5.
        const { A2: c5, B2: d5 } = ApprParallelLines2.computeProduct(a30, b30, c4, d4);

        // Compute S20^2 * S11 = c6 + gamma*d6.
        const { A2: c6, B2: d6 } = ApprParallelLines2.computeProduct(c0, d0, a11, b11);

        // Compute S20 * (S20^2 * S11) = c7 + gamma*d7.
        const { A2: c7, B2: d7 } = ApprParallelLines2.computeProduct(a20, b20, c6, d6);

        // Compute F = 2*S31*S20^2 - 3*S30*S21*S20 + S30^2*S11 - 2*S20^3*S11
        // = f0 + gamma*f1, where f0 is even of degree 8 and f1 is odd of
        // degree 7.
        const f0 = c1.sub(c7).mul(2).sub(c3.mul(3)).add(c5);
        const f1 = d1.sub(d7).mul(2).sub(d3.mul(3)).add(d5);
        return { f0: f0, f1: f1 };
    }

    private static updateParameters(data: ZValues, sigma: number,
        sigmaSqr: number, gamma: number, minimum: Minimum): void {
        // Rather than evaluate apq(sigma) and bpq(sigma), the polynomials
        // are evaluated at sigmaSqr to avoid the rounding errors that are
        // inherent by computing s = sqrt(ssqr); ssqr = s * s;
        const A20 = data.Z02 + (data.Z20 - data.Z02) * sigmaSqr;
        const B20 = -2 * data.Z11 * sigma;
        const S20 = A20 + gamma * B20;
        const A30 = -sigma * (3 * data.Z12 + (data.Z30 - 3 * data.Z12) * sigmaSqr);
        const B30 = data.Z03 + (3 * data.Z21 - data.Z03) * sigmaSqr;
        const S30 = A30 + gamma * B30;
        const A40 = data.Z04 + ((6 * data.Z22 - 2 * data.Z04)
            + (data.Z40 - 6 * data.Z22 + data.Z04) * sigmaSqr) * sigmaSqr;
        const B40 = -4 * sigma * (data.Z13 + (data.Z31 - data.Z13) * sigmaSqr);
        const S40 = A40 + gamma * B40;
        const k = S30 / (2 * S20);
        const ksqr = k * k;
        const rsqr = ksqr + S20;
        const error = S40 - 4 * k * S30 + (4 * ksqr - S20) * S20;
        if (error < minimum.error) {
            minimum.sigma = sigma;
            minimum.gamma = gamma;
            minimum.k = k;
            minimum.rSqr = rsqr;
            minimum.error = error;
        }
    }
}
