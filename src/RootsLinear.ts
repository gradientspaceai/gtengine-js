// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsLinear.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the real-valued root of a linear polynomial with real-valued
// coefficients. The general linear polynomial is g(x) = g0 + g1 * x and the
// monic linear polynomial is m(x) = m0 + x.
//
// Port notes:
//   - Upstream RootsLinear<T> is instantiated with T = 'float', 'double' or
//     Rational = BSRational<UIntegerAP32>. The port has no templates, so the
//     'double' instantiation is exposed as RootsLinear.solve/solveMonic and
//     the rational instantiation as RootsLinear.solveRational and
//     RootsLinear.solveMonicRational. The shared algorithm lives in the
//     generic helpers solveLinearGeneric/solveMonicLinearGeneric, which are
//     parameterized by a RootsScalarOps<T, R> record. RootsQuadratic,
//     RootsCubic and RootsQuartic reuse the same mechanism, which is how they
//     recover the upstream 'RootsLinear<T>::Solve' calls with T still open.
//   - The upstream out-parameter 'PolynomialRoot<T>* roots' plus the size_t
//     count become a returned array of roots.
//   - PolynomialRoot<Rational> becomes the class PolynomialRootRational,
//     defined here because PolynomialRoot.ts (B114) ports only the
//     floating-point instantiation.

import { BSRational } from './BSRational';
import { PolynomialRoot } from './PolynomialRoot';

// The port of PolynomialRoot<Rational> from PolynomialRoot.h. As upstream,
// the comparison operators use only the root estimate x.
export class PolynomialRootRational {
    // x is the root estimate and m is the multiplicity of x. The object is
    // invalid when m is 0.
    x: BSRational;
    m: number;

    constructor(x?: BSRational, m: number = 0) {
        this.x = (x !== undefined ? x : BSRational.fromNumber(0));
        this.m = m;
    }

    equals(other: PolynomialRootRational): boolean {
        return this.x.equals(other.x);
    }

    lessThan(other: PolynomialRootRational): boolean {
        return this.x.lessThan(other.x);
    }
}

// The port-specific stand-in for the C++ template parameter T of the Roots*
// classes. T is the coefficient type ('number' or BSRational) and R is the
// matching root record type.
export interface RootsScalarOps<T, R extends { x: T, m: number }> {
    // A freshly created zero. This is a function rather than a constant
    // because BSRational objects are mutable.
    zero(): T;
    isZero(v: T): boolean;
    negate(v: T): T;
    div(a: T, b: T): T;
    lessThan(a: T, b: T): boolean;
    fromRational(r: BSRational): T;
    toRational(v: T): BSRational;
    makeRoot(x: T, m: number): R;
}

// The ops for the upstream T = 'double' instantiations.
export const rootsNumberOps: RootsScalarOps<number, PolynomialRoot> = {
    zero: () => 0,
    isZero: (v: number) => v === 0,
    negate: (v: number) => -v,
    div: (a: number, b: number) => a / b,
    lessThan: (a: number, b: number) => a < b,
    fromRational: (r: BSRational) => r.toNumber(),
    toRational: (v: number) => BSRational.fromNumber(v),
    makeRoot: (x: number, m: number) => new PolynomialRoot(x, m)
};

// The ops for the upstream T = Rational instantiations.
export const rootsRationalOps: RootsScalarOps<BSRational, PolynomialRootRational> = {
    zero: () => BSRational.fromNumber(0),
    isZero: (v: BSRational) => v.getSign() === 0,
    negate: (v: BSRational) => v.negated(),
    div: (a: BSRational, b: BSRational) => a.div(b),
    lessThan: (a: BSRational, b: BSRational) => a.lessThan(b),
    fromRational: (r: BSRational) => r.clone(),
    toRational: (v: BSRational) => v.clone(),
    makeRoot: (x: BSRational, m: number) => new PolynomialRootRational(x, m)
};

// The port of std::sort on the roots, which orders by the root estimate only.
export function sortRoots<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, roots: R[]): void {
    roots.sort((a, b) => (ops.lessThan(a.x, b.x) ? -1 : (ops.lessThan(b.x, a.x) ? +1 : 0)));
}

// The port of RootsLinear<T>::Solve(g0, g1, roots).
export function solveLinearGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, g0: T, g1: T): R[] {
    // Test whether the degree is smaller than 1.
    if (ops.isZero(g1)) {
        // The solution set is either all real-valued x (g0 = 0) or no
        // solution (g0 != 0). In either case, report no roots.
        return [];
    }

    // Test for zero-valued roots.
    if (ops.isZero(g0)) {
        return [ops.makeRoot(ops.zero(), 1)];
    }

    // At this time g0 and g1 are not zero.
    return [ops.makeRoot(ops.div(ops.negate(g0), g1), 1)];
}

// The port of RootsLinear<T>::Solve(m0, roots).
export function solveMonicLinearGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, m0: T): R[] {
    return [ops.makeRoot(ops.negate(m0), 1)];
}

export class RootsLinear {
    // Solve the general polynomial g0 + g1*x = 0.
    static solve(g0: number, g1: number): PolynomialRoot[] {
        return solveLinearGeneric(rootsNumberOps, g0, g1);
    }

    // Solve the monic polynomial m0 + x = 0.
    static solveMonic(m0: number): PolynomialRoot[] {
        return solveMonicLinearGeneric(rootsNumberOps, m0);
    }

    // Solve the general polynomial g0 + g1*x = 0 with rational coefficients.
    static solveRational(g0: BSRational, g1: BSRational): PolynomialRootRational[] {
        return solveLinearGeneric(rootsRationalOps, g0, g1);
    }

    // Solve the monic polynomial m0 + x = 0 with rational coefficients.
    static solveMonicRational(m0: BSRational): PolynomialRootRational[] {
        return solveMonicLinearGeneric(rootsRationalOps, m0);
    }
}
