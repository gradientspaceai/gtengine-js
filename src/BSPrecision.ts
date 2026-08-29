// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSPrecision.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Support for determining the number of bits of precision required to compute
// an expression using BSNumber or BSRational.
//
// Port notes: the nested enum BSPrecision::Type is exported as
// BSPrecisionType and the nested struct BSPrecision::Parameters as
// BSPrecisionParameters (global-export-uniqueness precedent). The C++ free
// operators become instance methods returning new objects:
//   operator+  -> add,  operator-  -> sub,  operator*  -> mul,
//   operator/  -> div,
//   operator== -> equal,       operator!= -> notEqual,
//   operator<  -> lessThan,    operator<= -> lessThanEqual,
//   operator>  -> greaterThan, operator>= -> greaterThanEqual
// The comparison methods do not compare BSPrecision objects; as upstream,
// they return the precision parameters required to evaluate the comparison
// exactly with BSNumber/BSRational operands.

export enum BSPrecisionType {
    IS_FLOAT,
    IS_DOUBLE,
    IS_INT32,
    IS_INT64,
    IS_UINT32,
    IS_UINT64
}

export class BSPrecisionParameters {
    minExponent: number;
    maxExponent: number;
    maxBits: number;
    maxWords: number;

    constructor(minExponent: number = 0, maxExponent: number = 0, maxBits: number = 0) {
        this.minExponent = minExponent;
        this.maxExponent = maxExponent;
        this.maxBits = maxBits;
        // The C++ default constructor sets maxWords to 0 and the 3-argument
        // constructor sets maxWords = GetMaxWords(); with maxBits = 0 the
        // formula also yields 0, so one path suffices.
        this.maxWords = this.getMaxWords();
    }

    getMaxWords(): number {
        return Math.trunc(this.maxBits / 32) + ((this.maxBits % 32) > 0 ? 1 : 0);
    }

    // C++ assignment (bsr = bsn) copies by value; TS objects alias, so
    // copies are made explicit.
    clone(): BSPrecisionParameters {
        return new BSPrecisionParameters(this.minExponent, this.maxExponent, this.maxBits);
    }
}

export class BSPrecision {
    bsn: BSPrecisionParameters;
    bsr: BSPrecisionParameters;

    constructor();
    constructor(type: BSPrecisionType);
    constructor(minExponent: number, maxExponent: number, maxBits: number);
    constructor(arg0?: BSPrecisionType | number, arg1?: number, arg2?: number) {
        if (arg0 === undefined) {
            // Default constructor: zero-valued parameters.
            this.bsn = new BSPrecisionParameters();
            this.bsr = new BSPrecisionParameters();
        } else if (arg1 === undefined) {
            // C++ BSPrecision(Type) value-initializes bsn and bsr, then
            // overwrites bsn in the switch and copies bsr = bsn.
            this.bsn = new BSPrecisionParameters();
            const type = arg0 as BSPrecisionType;
            switch (type) {
                case BSPrecisionType.IS_FLOAT:
                    this.bsn = new BSPrecisionParameters(-149, 127, 24);
                    break;
                case BSPrecisionType.IS_DOUBLE:
                    this.bsn = new BSPrecisionParameters(-1074, 1023, 53);
                    break;
                case BSPrecisionType.IS_INT32:
                    this.bsn = new BSPrecisionParameters(0, 30, 31);
                    break;
                case BSPrecisionType.IS_INT64:
                    this.bsn = new BSPrecisionParameters(0, 62, 63);
                    break;
                case BSPrecisionType.IS_UINT32:
                    this.bsn = new BSPrecisionParameters(0, 31, 32);
                    break;
                case BSPrecisionType.IS_UINT64:
                    this.bsn = new BSPrecisionParameters(0, 63, 64);
                    break;
            }
            this.bsr = this.bsn.clone();
        } else {
            this.bsn = new BSPrecisionParameters(arg0 as number, arg1, arg2 as number);
            this.bsr = new BSPrecisionParameters(arg0 as number, arg1, arg2 as number);
        }
    }

    // C++ operator+(bsp0, bsp1).
    add(bsp1: BSPrecision): BSPrecision {
        const bsp0 = this;
        const result = new BSPrecision();

        result.bsn.minExponent = Math.min(bsp0.bsn.minExponent, bsp1.bsn.minExponent);
        if (bsp0.bsn.maxExponent >= bsp1.bsn.maxExponent) {
            result.bsn.maxExponent = bsp0.bsn.maxExponent;
            if (bsp0.bsn.maxExponent - bsp0.bsn.maxBits + 1 <= bsp1.bsn.maxExponent) {
                ++result.bsn.maxExponent;
            }

            result.bsn.maxBits = bsp0.bsn.maxExponent - bsp1.bsn.minExponent + 1;
            if (result.bsn.maxBits <= bsp0.bsn.maxBits + bsp1.bsn.maxBits - 1) {
                ++result.bsn.maxBits;
            }
        } else {
            result.bsn.maxExponent = bsp1.bsn.maxExponent;
            if (bsp1.bsn.maxExponent - bsp1.bsn.maxBits + 1 <= bsp0.bsn.maxExponent) {
                ++result.bsn.maxExponent;
            }

            result.bsn.maxBits = bsp1.bsn.maxExponent - bsp0.bsn.minExponent + 1;
            if (result.bsn.maxBits <= bsp0.bsn.maxBits + bsp1.bsn.maxBits - 1) {
                ++result.bsn.maxBits;
            }
        }
        result.bsn.maxWords = result.bsn.getMaxWords();

        // Addition is n0/d0 + n1/d1 = (n0*d1 + n1*d0)/(d0*d1). The numerator
        // and denominator of a number are assumed to have the same
        // parameters, so for the addition, the numerator is used for the
        // parameter computations.

        // Compute the parameters for the multiplication.
        const mulMinExponent = bsp0.bsr.minExponent + bsp1.bsr.minExponent;
        const mulMaxExponent = bsp0.bsr.maxExponent + bsp1.bsr.maxExponent + 1;
        const mulMaxBits = bsp0.bsr.maxBits + bsp1.bsr.maxBits;

        // Compute the parameters for the addition. The numbers n0*d1 and
        // n1*d0 are in the same arbitrary-precision set.
        result.bsr.minExponent = mulMinExponent;
        result.bsr.maxExponent = mulMaxExponent + 1; // Always a carry-out.
        result.bsr.maxBits = mulMaxExponent - mulMinExponent + 1;
        if (result.bsr.maxBits <= 2 * mulMaxBits - 1) {
            ++result.bsr.maxBits;
        }
        result.bsr.maxWords = result.bsr.getMaxWords();

        return result;
    }

    // C++ operator-(bsp0, bsp1).
    sub(bsp1: BSPrecision): BSPrecision {
        return this.add(bsp1);
    }

    // C++ operator*(bsp0, bsp1).
    mul(bsp1: BSPrecision): BSPrecision {
        const bsp0 = this;
        const result = new BSPrecision();

        result.bsn.minExponent = bsp0.bsn.minExponent + bsp1.bsn.minExponent;
        result.bsn.maxExponent = bsp0.bsn.maxExponent + bsp1.bsn.maxExponent + 1;
        result.bsn.maxBits = bsp0.bsn.maxBits + bsp1.bsn.maxBits;
        result.bsn.maxWords = result.bsn.getMaxWords();

        // Multiplication is (n0/d0) * (n1/d1) = (n0 * n1) / (d0 * d1). The
        // parameters are the same as for numerator/denominator.
        result.bsr.minExponent = bsp0.bsr.minExponent + bsp1.bsr.minExponent;
        result.bsr.maxExponent = bsp0.bsr.maxExponent + bsp1.bsr.maxExponent + 1;
        result.bsr.maxBits = bsp0.bsr.maxBits + bsp1.bsr.maxBits;
        result.bsr.maxWords = result.bsr.getMaxWords();

        return result;
    }

    // C++ operator/(bsp0, bsp1).
    div(bsp1: BSPrecision): BSPrecision {
        const bsp0 = this;
        const result = new BSPrecision();

        // BSNumber does not support division, so result.bsn has all members
        // set to zero.

        // Division is (n0/d0) / (n1/d1) = (n0 * d1) / (n1 * d0). The
        // parameters are the same as for multiplication.
        result.bsr.minExponent = bsp0.bsr.minExponent + bsp1.bsr.minExponent;
        result.bsr.maxExponent = bsp0.bsr.maxExponent + bsp1.bsr.maxExponent + 1;
        result.bsr.maxBits = bsp0.bsr.maxBits + bsp1.bsr.maxBits;
        result.bsr.maxWords = result.bsr.getMaxWords();

        return result;
    }

    // Comparisons for BSNumber do not involve dynamic allocations, so the
    // results are the extremes of the inputs. Comparisons for BSRational
    // involve multiplications of numerators and denominators.

    // C++ operator==(bsp0, bsp1).
    equal(bsp1: BSPrecision): BSPrecision {
        const bsp0 = this;
        const result = new BSPrecision();

        result.bsn.minExponent = Math.min(bsp0.bsn.minExponent, bsp1.bsn.minExponent);
        result.bsn.maxExponent = Math.max(bsp0.bsn.maxExponent, bsp1.bsn.maxExponent);
        result.bsn.maxBits = Math.max(bsp0.bsn.maxBits, bsp1.bsn.maxBits);
        result.bsn.maxWords = result.bsn.getMaxWords();

        result.bsr.minExponent = bsp0.bsr.minExponent + bsp1.bsr.minExponent;
        result.bsr.maxExponent = bsp0.bsr.maxExponent + bsp1.bsr.maxExponent + 1;
        result.bsr.maxBits = bsp0.bsr.maxBits + bsp1.bsr.maxBits;
        result.bsr.maxWords = result.bsr.getMaxWords();

        return result;
    }

    // C++ operator!=(bsp0, bsp1).
    notEqual(bsp1: BSPrecision): BSPrecision {
        return this.equal(bsp1);
    }

    // C++ operator<(bsp0, bsp1).
    lessThan(bsp1: BSPrecision): BSPrecision {
        return this.equal(bsp1);
    }

    // C++ operator<=(bsp0, bsp1).
    lessThanEqual(bsp1: BSPrecision): BSPrecision {
        return this.equal(bsp1);
    }

    // C++ operator>(bsp0, bsp1).
    greaterThan(bsp1: BSPrecision): BSPrecision {
        return this.equal(bsp1);
    }

    // C++ operator>=(bsp0, bsp1).
    greaterThanEqual(bsp1: BSPrecision): BSPrecision {
        return this.equal(bsp1);
    }
}
