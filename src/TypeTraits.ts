// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TypeTraits.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Upstream TypeTraits.h supplies compile-time traits (is_arbitrary_precision,
// has_division_operator) and SFINAE selectors used to choose between the
// floating-point and arbitrary-precision code paths of a template. TypeScript
// has no compile-time template dispatch, so the port makes the traits
// runtime predicates on values:
//
// - Plain 'number' (the port of float/double/long double) is not arbitrary
//   precision and supports division.
// - Arbitrary-precision types (BSNumber, BSRational, QFNumber, ported on top
//   of bigint) implement the ArbitraryPrecisionNumber marker interface below,
//   declaring 'isArbitraryPrecision: true' and whether they support division
//   (true for BSRational, false for BSNumber), exactly mirroring the trait
//   specializations upstream places in their headers.
//
// Dependent files that upstream splits into enable_if overloads are ported as
// a single function/class using these predicates for runtime dispatch, or,
// per PORTING.md, only the floating-point instantiation is ported when the
// exact-arithmetic path is not needed. The SFINAE helper aliases
// (TraitSelector, IsFPType, IsAPType, IsDivisionType, IsNotDivisionType) are
// compile-time machinery with no runtime equivalent and are intentionally
// omitted.

// Marker interface implemented by arbitrary-precision number types.
export interface ArbitraryPrecisionNumber {
    readonly isArbitraryPrecision: true;
    readonly hasDivisionOperator: boolean;
}

// The port of is_arbitrary_precision<T>::value: false for 'number', true for
// values implementing ArbitraryPrecisionNumber.
export function isArbitraryPrecision(x: unknown): x is ArbitraryPrecisionNumber {
    return typeof x === 'object' && x !== null
        && (x as ArbitraryPrecisionNumber).isArbitraryPrecision === true;
}

// The port of has_division_operator<T>::value: true for 'number' (float,
// double, long double upstream), and for arbitrary-precision values whatever
// the type declares (BSRational supports division, BSNumber does not).
export function hasDivisionOperator(x: unknown): boolean {
    if (typeof x === 'number') {
        return true;
    }
    if (isArbitraryPrecision(x)) {
        return x.hasDivisionOperator;
    }
    return false;
}
