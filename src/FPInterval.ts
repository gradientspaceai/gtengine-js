// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) FPInterval.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

import { IEEEBinary64 } from './IEEEBinary.js';
import { logAssert } from './Logger.js';

// ===========================================================================
// PORT DEVIATION: directed rounding is emulated, not native.
// ===========================================================================
//
// Upstream FPInterval.h computes rigorous interval bounds by switching the
// hardware rounding mode around every arithmetic operation:
//
//     std::fesetround(FE_DOWNWARD);  w[0] = u + v;
//     std::fesetround(FE_UPWARD);    w[1] = u + v;
//
// JavaScript has no equivalent of <cfenv>: every double operation rounds to
// nearest-even and the mode cannot be changed. A direct transliteration would
// compute w[0] == w[1] == roundNearest(u + v), which is NOT an enclosure --
// the true real result can fall outside such an "interval", so every
// downstream guarantee of interval arithmetic would be lost.
//
// This port therefore emulates directed rounding by widening. Let x be the
// round-to-nearest double of an exact real value v. IEEE 754 guarantees
// |x - v| <= 0.5 ulp(x), so
//
//     roundDown(x) <= v <= roundUp(x)
//
// always holds when roundDown/roundUp step one representable value outward.
// The port implements them with the ported IEEEBinary64 next-down/next-up
// encodings, whose behavior at the infinities is exactly what enclosure
// needs:
//
//   roundDown(+inf) = +MAX_VALUE  (an overflowed upper computation still
//                                  bounds the true value from below)
//   roundDown(-inf) = -inf        (a true -inf bound stays -inf)
//   roundUp(-inf)   = -MAX_VALUE
//   roundUp(+inf)   = +inf
//
// The last two lines are the reason FPInterval does not reuse the nextUp and
// nextDown helpers of the sibling SWInterval port. Those reproduce
// std::nextafter(x, -/+max), which steps an infinite bound back to a finite
// one and so loses the enclosure for intervals with infinite endpoints; that
// is the defect reported upstream for SWInterval.h. Upstream FPInterval.h
// does not have it, because directed rounding leaves an infinity alone, and
// neither does this port. FPInterval produces intervals with infinite
// endpoints routinely (reciprocalDown, reciprocalUp and reals), and those
// endpoints then flow through mul, so getting this right matters here.
//
// Widening unconditionally would be correct but would lose the exact
// degenerate intervals that upstream produces whenever an operation happens
// to be exact (for example [0.5,0.5] + [0.25,0.25] is [0.75,0.75] with
// directed rounding, since 0.75 is representable). So each arithmetic
// primitive first runs a rigorous exactness test on the round-to-nearest
// result:
//
//   * addition/subtraction: Knuth's TwoSum computes the rounding error of
//     a + b exactly (for finite results); an error of 0 proves a + b is
//     representable, so no widening is needed.
//   * multiplication: Dekker's TwoProduct (with the classic guards against
//     overflow in the splitting and against underflow of the error terms)
//     computes the rounding error of a * b exactly; an error of 0 proves the
//     product is representable.
//   * division: q = a / b is exact iff the exact product q * b equals a,
//     which is decided with TwoProduct.
//
// When a guard is not satisfied the test answers "not proven exact" and the
// endpoint is widened. The enclosure is therefore conservative in every case;
// the exactness tests only ever make the result tighter, never wrong.
//
// Other port notes:
// - Only the 'double' instantiation is ported (the port maps C++ floating
//   point to IEEE binary64); the 'float' instantiation is not provided.
// - The API mirrors the sibling SWInterval port, whose upstream header has
//   the same structure. The C++ free operators become instance methods
//   returning new objects:
//     operator- (unary)                 -> negate
//     operator+/-/*//                   -> add/sub/mul/div (each accepts an
//                                          FPInterval or a scalar number,
//                                          porting the operator(interval, T)
//                                          overloads; s + v and s * v are
//                                          covered by v.add(s) and v.mul(s),
//                                          which produce identical endpoints)
//     operator-(T, interval)            -> static scalarSub
//     operator/(T, interval)            -> static scalarDiv
//   The compound assignments (+=, -=, *=, /=) are not ported; intervals are
//   immutable, so u = u.add(v) is the replacement. operator[] becomes get(i).
//   The static leaf-node operations Add/Sub/Mul/Div keep their upstream
//   2-argument and internal 4-argument overloads as add/sub/mul/div, and
//   Mul2, Reciprocal, ReciprocalDown, ReciprocalUp and Reals become mul2,
//   reciprocal, reciprocalDown, reciprocalUp and reals.
// - Upstream's ProductLowerBound/ProductUpperBound require the caller to have
//   selected FE_DOWNWARD/FE_UPWARD so that many bounds can be computed in one
//   batch without repeated FPU control-word changes. That batching has no
//   meaning here, so productLowerBound and productUpperBound apply the port's
//   own directed rounding internally and are directly usable.
// - GTE_THROW_ON_INVALID_INTERVAL is ported as the static flag
//   FPInterval.throwOnInvalid (default false, matching upstream's
//   commented-out #define), the same name the SWInterval port uses.

// Dekker's splitting constant, 2^27 + 1.
const SPLITTER = 134217729;

// 2^995. At or above this magnitude the Dekker split SPLITTER * a overflows.
const SPLIT_SAFE_MAX = 3.2311742677852644e+299;

// 2^-969. Below this magnitude the TwoProduct error terms can underflow, so
// the exactness test is not applied.
const TWO_PRODUCT_MIN = 5.293955920339377e-292;

// The FPInterval [e0, e1] must satisfy e0 <= e1. Set FPInterval.throwOnInvalid
// to trap invalid construction where e0 > e1.
export class FPInterval {
    // Convenient constants.
    static readonly zero: number = 0;
    static readonly one: number = 1;
    static readonly max: number = Number.MAX_VALUE;
    static readonly inf: number = Number.POSITIVE_INFINITY;

    // Port of the upstream macro GTE_THROW_ON_INVALID_INTERVAL.
    static throwOnInvalid = false;

    private mEndpoints: [number, number];

    // Construction. All such intervals are conceptually immutable once
    // created; the endpoints cannot be modified outside the arithmetic
    // operations. new FPInterval() creates [0, 0], new FPInterval(e) creates
    // the degenerate interval [e, e] and new FPInterval(e0, e1) creates
    // [e0, e1].
    constructor(e0: number = 0, e1?: number) {
        if (e1 === undefined) {
            this.mEndpoints = [e0, e0];
        } else {
            this.mEndpoints = [e0, e1];
            if (FPInterval.throwOnInvalid) {
                logAssert(e0 <= e1, 'Invalid FPInterval.');
            }
        }
    }

    // The port of the C++ constructor taking std::array<FPType, 2>.
    static fromEndpoints(endpoint: readonly [number, number]): FPInterval {
        return new FPInterval(endpoint[0], endpoint[1]);
    }

    // C++ copy construction/assignment copies by value; TS objects alias, so
    // copies are made explicit.
    clone(): FPInterval {
        return new FPInterval(this.mEndpoints[0], this.mEndpoints[1]);
    }

    // Member access (the port of operator[]). It is only possible to read the
    // endpoints.
    get(i: number): number {
        return this.mEndpoints[i];
    }

    getEndpoints(): [number, number] {
        return [this.mEndpoints[0], this.mEndpoints[1]];
    }

    // -----------------------------------------------------------------
    // Emulated directed rounding. See the notes at the top of the file.
    // -----------------------------------------------------------------

    // Given the round-to-nearest double x of an exact real value v, return a
    // double that is guaranteed to be <= v. The port's stand-in for
    // std::fesetround(FE_DOWNWARD).
    static roundDown(x: number): number {
        return new IEEEBinary64(IEEEBinary64.fromNumber(x).getNextDown()).number;
    }

    // Given the round-to-nearest double x of an exact real value v, return a
    // double that is guaranteed to be >= v. The port's stand-in for
    // std::fesetround(FE_UPWARD).
    static roundUp(x: number): number {
        return new IEEEBinary64(IEEEBinary64.fromNumber(x).getNextUp()).number;
    }

    // Knuth's TwoSum: for finite a, b and s = a + b (round to nearest), the
    // rounding error is computed exactly. Returns true when the error is
    // zero, that is, when s is exactly a + b.
    private static sumIsExact(a: number, b: number, s: number): boolean {
        if (!Number.isFinite(s)) {
            return false;
        }
        const bb = s - a;
        const err = (a - (s - bb)) + (b - bb);
        return err === 0;
    }

    // Dekker's TwoProduct: for p = a * b (round to nearest), the rounding
    // error is computed exactly provided the splitting does not overflow and
    // the error terms do not underflow. Returns true only when the product is
    // proven to be exactly representable.
    private static productIsExact(a: number, b: number, p: number): boolean {
        if (!Number.isFinite(p)) {
            return false;
        }
        if (a === 0 || b === 0) {
            // The product is exactly a signed zero.
            return true;
        }
        if (Math.abs(a) >= SPLIT_SAFE_MAX || Math.abs(b) >= SPLIT_SAFE_MAX
            || Math.abs(p) < TWO_PRODUCT_MIN) {
            // The guards for Dekker's algorithm are not satisfied, so the
            // error cannot be computed reliably. Report "not proven exact".
            return false;
        }

        let c = SPLITTER * a;
        const ah = c - (c - a);
        const al = a - ah;
        c = SPLITTER * b;
        const bh = c - (c - b);
        const bl = b - bh;
        const err = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
        return err === 0;
    }

    // q = a / b is exact iff the exact real product q * b equals a. The
    // product is checked with TwoProduct: the rounded product must equal a
    // and must itself be exact.
    private static quotientIsExact(a: number, b: number, q: number): boolean {
        if (!Number.isFinite(q)) {
            return false;
        }
        if (a === 0) {
            // b is nonzero at every call site, so q is a signed zero and the
            // quotient is exact.
            return true;
        }
        const p = q * b;
        return p === a && FPInterval.productIsExact(q, b, p);
    }

    // The eight directed primitives used by every operation below.
    private static addDown(a: number, b: number): number {
        const s = a + b;
        return FPInterval.sumIsExact(a, b, s) ? s : FPInterval.roundDown(s);
    }

    private static addUp(a: number, b: number): number {
        const s = a + b;
        return FPInterval.sumIsExact(a, b, s) ? s : FPInterval.roundUp(s);
    }

    private static subDown(a: number, b: number): number {
        const d = a - b;
        return FPInterval.sumIsExact(a, -b, d) ? d : FPInterval.roundDown(d);
    }

    private static subUp(a: number, b: number): number {
        const d = a - b;
        return FPInterval.sumIsExact(a, -b, d) ? d : FPInterval.roundUp(d);
    }

    private static mulDown(a: number, b: number): number {
        const p = a * b;
        return FPInterval.productIsExact(a, b, p) ? p : FPInterval.roundDown(p);
    }

    private static mulUp(a: number, b: number): number {
        const p = a * b;
        return FPInterval.productIsExact(a, b, p) ? p : FPInterval.roundUp(p);
    }

    private static divDown(a: number, b: number): number {
        const q = a / b;
        return FPInterval.quotientIsExact(a, b, q) ? q : FPInterval.roundDown(q);
    }

    private static divUp(a: number, b: number): number {
        const q = a / b;
        return FPInterval.quotientIsExact(a, b, q) ? q : FPInterval.roundUp(q);
    }

    // -----------------------------------------------------------------
    // Arithmetic operations to compute intervals at the leaf nodes of an
    // expression tree (the 2-argument overloads, where u and v are the raw
    // floating-point variables of the expression). The 4-argument overloads
    // are the upstream internal-use functions consumed by the instance
    // operator methods at the interior nodes of the expression tree.
    // -----------------------------------------------------------------

    static add(u: number, v: number): FPInterval;
    static add(u0: number, u1: number, v0: number, v1: number): FPInterval;
    static add(a0: number, a1: number, a2?: number, a3?: number): FPInterval {
        if (a2 === undefined || a3 === undefined) {
            return new FPInterval(
                FPInterval.addDown(a0, a1), FPInterval.addUp(a0, a1));
        }
        return new FPInterval(
            FPInterval.addDown(a0, a2), FPInterval.addUp(a1, a3));
    }

    static sub(u: number, v: number): FPInterval;
    static sub(u0: number, u1: number, v0: number, v1: number): FPInterval;
    static sub(a0: number, a1: number, a2?: number, a3?: number): FPInterval {
        if (a2 === undefined || a3 === undefined) {
            return new FPInterval(
                FPInterval.subDown(a0, a1), FPInterval.subUp(a0, a1));
        }
        return new FPInterval(
            FPInterval.subDown(a0, a3), FPInterval.subUp(a1, a2));
    }

    static mul(u: number, v: number): FPInterval;
    static mul(u0: number, u1: number, v0: number, v1: number): FPInterval;
    static mul(a0: number, a1: number, a2?: number, a3?: number): FPInterval {
        if (a2 === undefined || a3 === undefined) {
            return new FPInterval(
                FPInterval.mulDown(a0, a1), FPInterval.mulUp(a0, a1));
        }
        return new FPInterval(
            FPInterval.mulDown(a0, a2), FPInterval.mulUp(a1, a3));
    }

    static mul2(u0: number, u1: number, v0: number, v1: number): FPInterval {
        const u0mv1 = FPInterval.mulDown(u0, v1);
        const u1mv0 = FPInterval.mulDown(u1, v0);
        const u0mv0 = FPInterval.mulUp(u0, v0);
        const u1mv1 = FPInterval.mulUp(u1, v1);
        return new FPInterval(Math.min(u0mv1, u1mv0), Math.max(u0mv0, u1mv1));
    }

    static div(u: number, v: number): FPInterval;
    static div(u0: number, u1: number, v0: number, v1: number): FPInterval;
    static div(a0: number, a1: number, a2?: number, a3?: number): FPInterval {
        if (a2 === undefined || a3 === undefined) {
            if (a1 !== FPInterval.zero) {
                return new FPInterval(
                    FPInterval.divDown(a0, a1), FPInterval.divUp(a0, a1));
            }
            // Division by zero does not lead to a determinate FPInterval.
            // Return the entire set of real numbers.
            return FPInterval.reals();
        }
        return new FPInterval(
            FPInterval.divDown(a0, a3), FPInterval.divUp(a1, a2));
    }

    static reciprocal(v0: number, v1: number): FPInterval {
        return new FPInterval(
            FPInterval.divDown(FPInterval.one, v1),
            FPInterval.divUp(FPInterval.one, v0));
    }

    static reciprocalDown(v: number): FPInterval {
        const recpv = FPInterval.divDown(FPInterval.one, v);
        return new FPInterval(recpv, FPInterval.inf);
    }

    static reciprocalUp(v: number): FPInterval {
        const recpv = FPInterval.divUp(FPInterval.one, v);
        return new FPInterval(-FPInterval.inf, recpv);
    }

    static reals(): FPInterval {
        return new FPInterval(-FPInterval.inf, FPInterval.inf);
    }

    // This function computes the lower bound on the product of the two
    // intervals given by their endpoint pairs u and v.
    //
    // UPSTREAM BUG (preserved): in the branch where both u and v straddle
    // zero, upstream returns u[0]*v[0], which is a product of two negative
    // numbers and therefore positive; the true lower bound is
    // min(u[0]*v[1], u[1]*v[0]), the value that Mul2 (and therefore
    // operator*) computes. For u = v = [-1, 1] upstream returns +1 while the
    // set of products is [-1, 1]. The corresponding branch of
    // ProductUpperBound returns u[1]*v[1] instead of
    // max(u[0]*v[0], u[1]*v[1]), which is too small for, say,
    // u = v = [-4, 1] (it returns 1 rather than 16). Nothing in GTE calls
    // either function, so the port preserves the upstream behavior rather
    // than silently diverging from it; use mul (operator*) for a correct
    // product interval.
    static productLowerBound(u: readonly number[], v: readonly number[]): number {
        let w0: number;
        if (u[0] >= FPInterval.zero) {
            if (v[0] >= FPInterval.zero) {
                w0 = FPInterval.mulDown(u[0], v[0]);
            } else if (v[1] <= FPInterval.zero) {
                w0 = FPInterval.mulDown(u[1], v[0]);
            } else {
                w0 = FPInterval.mulDown(u[1], v[0]);
            }
        } else if (u[1] <= FPInterval.zero) {
            if (v[0] >= FPInterval.zero) {
                w0 = FPInterval.mulDown(u[0], v[1]);
            } else if (v[1] <= FPInterval.zero) {
                w0 = FPInterval.mulDown(u[1], v[1]);
            } else {
                w0 = FPInterval.mulDown(u[0], v[1]);
            }
        } else {
            if (v[0] >= FPInterval.zero) {
                w0 = FPInterval.mulDown(u[0], v[1]);
            } else if (v[1] <= FPInterval.zero) {
                w0 = FPInterval.mulDown(u[1], v[0]);
            } else {
                w0 = FPInterval.mulDown(u[0], v[0]);
            }
        }
        return w0;
    }

    // This function computes the upper bound on the product of the two
    // intervals given by their endpoint pairs u and v.
    static productUpperBound(u: readonly number[], v: readonly number[]): number {
        let w1: number;
        if (u[0] >= FPInterval.zero) {
            if (v[0] >= FPInterval.zero) {
                w1 = FPInterval.mulUp(u[1], v[1]);
            } else if (v[1] <= FPInterval.zero) {
                w1 = FPInterval.mulUp(u[0], v[1]);
            } else {
                w1 = FPInterval.mulUp(u[1], v[1]);
            }
        } else if (u[1] <= FPInterval.zero) {
            if (v[0] >= FPInterval.zero) {
                w1 = FPInterval.mulUp(u[1], v[0]);
            } else if (v[1] <= FPInterval.zero) {
                w1 = FPInterval.mulUp(u[0], v[0]);
            } else {
                w1 = FPInterval.mulUp(u[0], v[0]);
            }
        } else {
            if (v[0] >= FPInterval.zero) {
                w1 = FPInterval.mulUp(u[1], v[1]);
            } else if (v[1] <= FPInterval.zero) {
                w1 = FPInterval.mulUp(u[0], v[0]);
            } else {
                w1 = FPInterval.mulUp(u[1], v[1]);
            }
        }
        return w1;
    }

    // Unary operator-. Negation of [e0, e1] produces [-e1, -e0]. This
    // operation needs to be supported in the sense of negating a "number" in
    // an arithmetic expression.
    negate(): FPInterval {
        return new FPInterval(-this.mEndpoints[1], -this.mEndpoints[0]);
    }

    // Addition operations.
    add(v: FPInterval | number): FPInterval {
        const u = this.mEndpoints;
        if (typeof v === 'number') {
            return FPInterval.add(u[0], u[1], v, v);
        }
        return FPInterval.add(u[0], u[1], v.mEndpoints[0], v.mEndpoints[1]);
    }

    // Subtraction operations.
    sub(v: FPInterval | number): FPInterval {
        const u = this.mEndpoints;
        if (typeof v === 'number') {
            return FPInterval.sub(u[0], u[1], v, v);
        }
        return FPInterval.sub(u[0], u[1], v.mEndpoints[0], v.mEndpoints[1]);
    }

    // The port of C++ operator-(FPType u, FPInterval const& v).
    static scalarSub(u: number, v: FPInterval): FPInterval {
        return FPInterval.sub(u, u, v.mEndpoints[0], v.mEndpoints[1]);
    }

    // Multiplication operations.
    mul(v: FPInterval | number): FPInterval {
        const u = this.mEndpoints;
        if (typeof v === 'number') {
            if (v >= FPInterval.zero) {
                return FPInterval.mul(u[0], u[1], v, v);
            } else {
                return FPInterval.mul(u[1], u[0], v, v);
            }
        }
        const w = v.mEndpoints;
        if (u[0] >= FPInterval.zero) {
            if (w[0] >= FPInterval.zero) {
                return FPInterval.mul(u[0], u[1], w[0], w[1]);
            } else if (w[1] <= FPInterval.zero) {
                return FPInterval.mul(u[1], u[0], w[0], w[1]);
            } else { // w[0] < 0 < w[1]
                return FPInterval.mul(u[1], u[1], w[0], w[1]);
            }
        } else if (u[1] <= FPInterval.zero) {
            if (w[0] >= FPInterval.zero) {
                return FPInterval.mul(u[0], u[1], w[1], w[0]);
            } else if (w[1] <= FPInterval.zero) {
                return FPInterval.mul(u[1], u[0], w[1], w[0]);
            } else { // w[0] < 0 < w[1]
                return FPInterval.mul(u[0], u[0], w[1], w[0]);
            }
        } else { // u[0] < 0 < u[1]
            if (w[0] >= FPInterval.zero) {
                return FPInterval.mul(u[0], u[1], w[1], w[1]);
            } else if (w[1] <= FPInterval.zero) {
                return FPInterval.mul(u[1], u[0], w[0], w[0]);
            } else { // w[0] < 0 < w[1]
                return FPInterval.mul2(u[0], u[1], w[0], w[1]);
            }
        }
    }

    // Division operations. If the divisor FPInterval is [v0, v1] with
    // v0 < 0 < v1, then the returned FPInterval is (-inf, +inf) instead of
    // Union((-inf, 1/v0), (1/v1, +inf)). An application should try to avoid
    // this case by branching based on [v0, 0] and [0, v1].
    div(v: FPInterval | number): FPInterval {
        const u = this.mEndpoints;
        if (typeof v === 'number') {
            if (v > FPInterval.zero) {
                return FPInterval.div(u[0], u[1], v, v);
            } else if (v < FPInterval.zero) {
                return FPInterval.div(u[1], u[0], v, v);
            } else { // v = 0
                return FPInterval.reals();
            }
        }
        const w = v.mEndpoints;
        if (w[0] > FPInterval.zero || w[1] < FPInterval.zero) {
            return this.mul(FPInterval.reciprocal(w[0], w[1]));
        } else {
            if (w[0] === FPInterval.zero) {
                return this.mul(FPInterval.reciprocalDown(w[1]));
            } else if (w[1] === FPInterval.zero) {
                return this.mul(FPInterval.reciprocalUp(w[0]));
            } else { // w[0] < 0 < w[1]
                return FPInterval.reals();
            }
        }
    }

    // The port of C++ operator/(FPType u, FPInterval const& v). The scalar u
    // is multiplied by the reciprocal interval; scalar multiplication is
    // commutative, so reciprocal(...).mul(u) produces the same endpoints as
    // the upstream u * Reciprocal(...).
    static scalarDiv(u: number, v: FPInterval): FPInterval {
        const w = v.mEndpoints;
        if (w[0] > FPInterval.zero || w[1] < FPInterval.zero) {
            return FPInterval.reciprocal(w[0], w[1]).mul(u);
        } else {
            if (w[0] === FPInterval.zero) {
                return FPInterval.reciprocalDown(w[1]).mul(u);
            } else if (w[1] === FPInterval.zero) {
                return FPInterval.reciprocalUp(w[0]).mul(u);
            } else { // w[0] < 0 < w[1]
                return FPInterval.reals();
            }
        }
    }
}
