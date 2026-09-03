import { describe, it, expect } from 'vitest';
import { check, fc, seededRandom } from './helpers/arbitraries.js';
import { BSPrecision, BSPrecisionParameters, BSPrecisionType } from '../src/BSPrecision.js';

function expectParams(
    p: BSPrecisionParameters,
    minExponent: number, maxExponent: number, maxBits: number, maxWords: number): void {
    expect(p.minExponent).toBe(minExponent);
    expect(p.maxExponent).toBe(maxExponent);
    expect(p.maxBits).toBe(maxBits);
    expect(p.maxWords).toBe(maxWords);
}

describe('BSPrecisionParameters', () => {
    it('default-constructs to all zeros', () => {
        expectParams(new BSPrecisionParameters(), 0, 0, 0, 0);
    });

    it('computes maxWords = ceil(maxBits / 32)', () => {
        expect(new BSPrecisionParameters(0, 0, 1).maxWords).toBe(1);
        expect(new BSPrecisionParameters(0, 0, 31).maxWords).toBe(1);
        expect(new BSPrecisionParameters(0, 0, 32).maxWords).toBe(1);
        expect(new BSPrecisionParameters(0, 0, 33).maxWords).toBe(2);
        expect(new BSPrecisionParameters(0, 0, 64).maxWords).toBe(2);
        expect(new BSPrecisionParameters(0, 0, 65).maxWords).toBe(3);
    });

    it('clone produces an independent copy', () => {
        const p = new BSPrecisionParameters(-149, 127, 24);
        const q = p.clone();
        expectParams(q, -149, 127, 24, 1);
        q.maxBits = 999;
        expect(p.maxBits).toBe(24);
    });
});

describe('BSPrecision type constructor', () => {
    it('produces the upstream parameters for each native type', () => {
        const cases: Array<[BSPrecisionType, number, number, number, number]> = [
            [BSPrecisionType.IS_FLOAT, -149, 127, 24, 1],
            [BSPrecisionType.IS_DOUBLE, -1074, 1023, 53, 2],
            [BSPrecisionType.IS_INT32, 0, 30, 31, 1],
            [BSPrecisionType.IS_INT64, 0, 62, 63, 2],
            [BSPrecisionType.IS_UINT32, 0, 31, 32, 1],
            [BSPrecisionType.IS_UINT64, 0, 63, 64, 2]
        ];
        for (const [type, minE, maxE, bits, words] of cases) {
            const bsp = new BSPrecision(type);
            expectParams(bsp.bsn, minE, maxE, bits, words);
            expectParams(bsp.bsr, minE, maxE, bits, words);
        }
    });

    it('bsr is a copy of bsn, not an alias (C++ value semantics)', () => {
        const bsp = new BSPrecision(BSPrecisionType.IS_FLOAT);
        bsp.bsn.maxBits = 999;
        expect(bsp.bsr.maxBits).toBe(24);
    });

    it('default constructor zeros both parameter sets', () => {
        const bsp = new BSPrecision();
        expectParams(bsp.bsn, 0, 0, 0, 0);
        expectParams(bsp.bsr, 0, 0, 0, 0);
    });

    it('explicit (minExponent, maxExponent, maxBits) constructor', () => {
        const bsp = new BSPrecision(-10, 10, 20);
        expectParams(bsp.bsn, -10, 10, 20, 1);
        expectParams(bsp.bsr, -10, 10, 20, 1);
    });
});

describe('BSPrecision add/sub', () => {
    it('float + float: BSNumber requires 277 bits (9 words)', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const sum = f.add(f);
        // maxBits = 127 - (-149) + 1 = 277, the well-known bound for the
        // exact sum of two arbitrary floats.
        expectParams(sum.bsn, -149, 128, 277, 9);
        // Rational path treats the sum as (n0*d1 + n1*d0)/(d0*d1).
        expectParams(sum.bsr, -298, 256, 554, 18);
    });

    it('double + double: BSNumber requires 2098 bits (66 words)', () => {
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        const sum = d.add(d);
        // maxBits = 1023 - (-1074) + 1 = 2098, the well-known bound for the
        // exact sum of two arbitrary doubles.
        expectParams(sum.bsn, -1074, 1024, 2098, 66);
        expectParams(sum.bsr, -2148, 2048, 4196, 132);
    });

    it('int32 + int32 requires 32 bits', () => {
        const i = new BSPrecision(BSPrecisionType.IS_INT32);
        const sum = i.add(i);
        expectParams(sum.bsn, 0, 31, 32, 1);
    });

    it('mixed-type addition is symmetric (both branches agree)', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        const fd = f.add(d);
        const df = d.add(f);
        expect(fd.bsn).toEqual(df.bsn);
        expect(fd.bsr).toEqual(df.bsr);
        // maxBits = 1023 - (-149) + 1 = 1173; no carry-out increment since
        // 1023 - 53 + 1 = 971 > 127.
        expectParams(fd.bsn, -1074, 1023, 1173, 37);
    });

    it('sub is identical to add (worst case is the same)', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        expect(f.sub(d)).toEqual(f.add(d));
        expect(d.sub(d)).toEqual(d.add(d));
    });
});

describe('BSPrecision mul', () => {
    it('float * float: BSNumber requires 48 bits (2 words)', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const prod = f.mul(f);
        expectParams(prod.bsn, -298, 255, 48, 2);
        // Product of rationals multiplies numerators/denominators
        // independently, so bsr matches bsn here.
        expectParams(prod.bsr, -298, 255, 48, 2);
    });

    it('double * double: BSNumber requires 106 bits (4 words)', () => {
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        const prod = d.mul(d);
        expectParams(prod.bsn, -2148, 2047, 106, 4);
        expectParams(prod.bsr, -2148, 2047, 106, 4);
    });
});

describe('BSPrecision div', () => {
    it('BSNumber does not support division: bsn is zeroed', () => {
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        const quot = d.div(d);
        expectParams(quot.bsn, 0, 0, 0, 0);
    });

    it('BSRational division has the same parameters as multiplication', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        expect(f.div(d).bsr).toEqual(f.mul(d).bsr);
        expect(d.div(d).bsr).toEqual(d.mul(d).bsr);
    });
});

describe('BSPrecision comparisons', () => {
    it('equal takes extremes for BSNumber and multiplies for BSRational', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        const cmp = f.equal(d);
        expectParams(cmp.bsn, -1074, 1023, 53, 2);
        // n0*d1 vs n1*d0: same parameters as multiplication.
        expectParams(cmp.bsr, -1223, 1151, 77, 3);
    });

    it('all comparison operators produce the same requirements', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        const expected = f.equal(d);
        expect(f.notEqual(d)).toEqual(expected);
        expect(f.lessThan(d)).toEqual(expected);
        expect(f.lessThanEqual(d)).toEqual(expected);
        expect(f.greaterThan(d)).toEqual(expected);
        expect(f.greaterThanEqual(d)).toEqual(expected);
    });
});

// Cross-checks against the worst-case BSNumber precision requirements
// documented in the upstream PrimalQuery2.h and PrimalQuery3.h comments
// ("Choice of N for UIntegerFP32<N>", compute type BSNumber). The
// expressions below mirror the arithmetic in those queries exactly.
// (The BSRational N values documented there were produced by an older
// BSPrecision with different rational formulas and are not asserted.)
describe('BSPrecision documented query requirements (BSNumber)', () => {
    // PrimalQuery2::ToLine: det = (u-u)*(u-u) - (u-u)*(u-u).
    function toLine(u: BSPrecision): BSPrecision {
        const s = u.sub(u);
        const p = s.mul(s);
        return p.sub(p);
    }

    // PrimalQuery3::ToPlane: 3x3 determinant of difference vectors,
    // det = x0*c0 + x1*c1 + x2*c2 with ci = yjzk - ykzj.
    function toPlane(u: BSPrecision): BSPrecision {
        const s = u.sub(u);
        const p = s.mul(s);
        const c = p.sub(p);
        const t = s.mul(c);
        const term = t.add(t);
        return term.add(t);
    }

    // PrimalQuery2::ToCircumcircle: zi = (u+u)*(u-u) + (u+u)*(u-u),
    // ci = yj*zk - yk*zj, det = x0*c0 + x1*c1 + x2*c2. The sums u+u and
    // differences u-u have identical parameters.
    function toCircumcircle(u: BSPrecision): BSPrecision {
        const s = u.sub(u);
        const t = s.mul(s);
        const z = t.add(t);
        const yz = s.mul(z);
        const c = yz.sub(yz);
        const xc = s.mul(c);
        const term = xc.add(xc);
        return term.add(xc);
    }

    // PrimalQuery3::ToCircumsphere: wi = sum of three products of
    // sums/differences, ai = xi*yj - xj*yi, bi = zi*wj - zj*wi,
    // det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0.
    function toCircumsphere(u: BSPrecision): BSPrecision {
        const s = u.sub(u);
        const t = s.mul(s);
        const tt = t.add(t);
        const w = tt.add(t);
        const xy = s.mul(s);
        const a = xy.sub(xy);
        const zw = s.mul(w);
        const b = zw.sub(zw);
        const ab = a.mul(b);
        const term0 = ab.sub(ab);
        const term1 = term0.add(ab);
        const term2 = term1.add(ab);
        const term3 = term2.sub(ab);
        return term3.add(ab);
    }

    const float = new BSPrecision(BSPrecisionType.IS_FLOAT);
    const double = new BSPrecision(BSPrecisionType.IS_DOUBLE);

    it('ToLine: float needs N = 18 words, double needs N = 132 words', () => {
        expect(toLine(float).bsn.maxWords).toBe(18);
        expect(toLine(double).bsn.maxWords).toBe(132);
    });

    it('ToPlane: float needs N = 27 words, double needs N = 197 words', () => {
        expect(toPlane(float).bsn.maxWords).toBe(27);
        expect(toPlane(double).bsn.maxWords).toBe(197);
    });

    it('ToCircumcircle: float needs N = 35 words, double needs N = 263 words', () => {
        expect(toCircumcircle(float).bsn.maxWords).toBe(35);
        expect(toCircumcircle(double).bsn.maxWords).toBe(263);
    });

    it('ToCircumsphere: float needs N = 44 words, double needs N = 329 words', () => {
        expect(toCircumsphere(float).bsn.maxWords).toBe(44);
        expect(toCircumsphere(double).bsn.maxWords).toBe(329);
    });

    it('ToLine intermediate bit counts match hand derivation', () => {
        // float: s = (-149, 128, 277), p = (-298, 257, 554),
        // det = (-298, 258, 557) -> 18 words.
        const s = float.sub(float);
        expectParams(s.bsn, -149, 128, 277, 9);
        const p = s.mul(s);
        expectParams(p.bsn, -298, 257, 554, 18);
        const det = p.sub(p);
        expectParams(det.bsn, -298, 258, 557, 18);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V05): property-based cross-checks of the upstream bounds
// against exact bigint arithmetic.
//
// Model of the arbitrary-precision set described by Parameters
// (minExponent, maxExponent, maxBits): the representable nonzero values are
//   v = s * M * 2^p,  s in {-1,+1},  M a positive odd integer,
//   p >= minExponent, bitlen(M) <= maxBits, p + bitlen(M) - 1 <= maxExponent.
// That is, minExponent bounds the trailing set bit, maxExponent bounds the
// leading set bit and maxBits bounds their span. BSPrecision's job is to
// produce a set that contains the exact result of an operation.
// ---------------------------------------------------------------------------
describe('BSPrecision verification', () => {
    interface PSet { minE: number; maxE: number; maxBits: number; }

    const namedSets: PSet[] = [
        { minE: -149, maxE: 127, maxBits: 24 },      // IS_FLOAT
        { minE: -1074, maxE: 1023, maxBits: 53 },    // IS_DOUBLE
        { minE: 0, maxE: 30, maxBits: 31 },          // IS_INT32
        { minE: 0, maxE: 62, maxBits: 63 },          // IS_INT64
        { minE: 0, maxE: 31, maxBits: 32 },          // IS_UINT32
        { minE: 0, maxE: 63, maxBits: 64 }           // IS_UINT64
    ];

    // Small synthetic sets keep the bigint fuzzing fast while exercising the
    // same branch structure as the native-type sets.
    const smallSet: fc.Arbitrary<PSet> = fc.tuple(
        fc.integer({ min: -30, max: 10 }),
        fc.integer({ min: 1, max: 24 }),
        fc.integer({ min: 0, max: 20 })
    ).map(([minE, bits, span]) => ({ minE, maxE: minE + bits - 1 + span, maxBits: bits }));

    const anySet: fc.Arbitrary<PSet> = fc.oneof(smallSet, fc.constantFrom(...namedSets));

    const toBSP = (s: PSet): BSPrecision => new BSPrecision(s.minE, s.maxE, s.maxBits);

    function bitLength(n: bigint): number {
        let m = n < 0n ? -n : n;
        let b = 0;
        while (m > 0n) { m >>= 1n; ++b; }
        return b;
    }

    /** Reduce M * 2^p to the canonical form with M odd. Returns null for zero. */
    function normalize(M: bigint, p: number): { M: bigint; p: number } | null {
        if (M === 0n) { return null; }
        let m = M;
        let e = p;
        while ((m & 1n) === 0n) { m >>= 1n; ++e; }
        return { M: m, p: e };
    }

    /** Draw a uniformly random member of the set using a seeded PRNG. */
    function sampleSet(s: PSet, rand: () => number): { M: bigint; p: number } {
        const bl = 1 + Math.floor(rand() * s.maxBits);
        let M = 1n << BigInt(bl - 1);
        for (let i = 1; i < bl - 1; ++i) {
            if (rand() < 0.5) { M |= 1n << BigInt(i); }
        }
        if (bl > 1) { M |= 1n; }   // odd mantissa
        const lo = s.minE;
        const hi = s.maxE - bl + 1;
        const p = lo + Math.floor(rand() * (hi - lo + 1));
        return { M: rand() < 0.5 ? -M : M, p };
    }

    function inSet(v: { M: bigint; p: number }, p: BSPrecisionParameters): boolean {
        const bits = bitLength(v.M);
        return v.p >= p.minExponent
            && v.p + bits - 1 <= p.maxExponent
            && bits <= p.maxBits;
    }

    it('maxWords is ceil(maxBits/32) on every produced result', () => {
        check(fc.tuple(anySet, anySet), ([a, b]) => {
            const x = toBSP(a);
            const y = toBSP(b);
            for (const r of [x.add(y), x.sub(y), x.mul(y), x.div(y), x.equal(y)]) {
                for (const p of [r.bsn, r.bsr]) {
                    expect(p.maxWords).toBe(Math.ceil(p.maxBits / 32));
                    expect(p.maxWords).toBe(p.getMaxWords());
                }
            }
            return true;
        });
    });

    it('mul and equal are symmetric in their operands', () => {
        check(fc.tuple(anySet, anySet), ([a, b]) => {
            const x = toBSP(a);
            const y = toBSP(b);
            expect(x.mul(y)).toEqual(y.mul(x));
            expect(x.equal(y)).toEqual(y.equal(x));
            return true;
        });
    });

    // add is symmetric whenever the two operands have distinct maxExponent:
    // both orderings then select the same "larger" operand. Ties are the
    // upstream asymmetry pinned further below.
    it('add is symmetric when the operand maxExponents differ', () => {
        check(fc.tuple(anySet, anySet).filter(([a, b]) => a.maxE !== b.maxE), ([a, b]) => {
            expect(toBSP(a).add(toBSP(b))).toEqual(toBSP(b).add(toBSP(a)));
            return true;
        });
    });

    it('sub aliases add and every comparison aliases equal', () => {
        check(fc.tuple(anySet, anySet), ([a, b]) => {
            const x = toBSP(a);
            const y = toBSP(b);
            expect(x.sub(y)).toEqual(x.add(y));
            const e = x.equal(y);
            expect(x.notEqual(y)).toEqual(e);
            expect(x.lessThan(y)).toEqual(e);
            expect(x.lessThanEqual(y)).toEqual(e);
            expect(x.greaterThan(y)).toEqual(e);
            expect(x.greaterThanEqual(y)).toEqual(e);
            return true;
        });
    });

    it('div zeroes bsn (BSNumber has no division) and reuses the mul bsr', () => {
        check(fc.tuple(anySet, anySet), ([a, b]) => {
            const q = toBSP(a).div(toBSP(b));
            expect(q.bsn).toEqual(new BSPrecisionParameters());
            expect(q.bsr).toEqual(toBSP(a).mul(toBSP(b)).bsr);
            return true;
        });
    });

    it('results never alias their operands or each other', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        const r = f.add(d);
        expect(r.bsn).not.toBe(f.bsn);
        expect(r.bsn).not.toBe(d.bsn);
        expect(r.bsn).not.toBe(r.bsr);
        r.bsn.maxBits = -1;
        expect(f.bsn.maxBits).toBe(24);
        expect(d.bsn.maxBits).toBe(53);
        expect(r.bsr.maxBits).not.toBe(-1);
    });

    // Exact cross-check: the product of any two members of the operand sets
    // lies in the mul bsn set. The multiplication bound is exact in general.
    it('mul bsn contains the exact bigint product for arbitrary operand sets', () => {
        const rand = seededRandom(0x5a17c0de);
        check(fc.tuple(anySet, anySet), ([a, b]) => {
            const r = toBSP(a).mul(toBSP(b)).bsn;
            for (let k = 0; k < 12; ++k) {
                const x = sampleSet(a, rand);
                const y = sampleSet(b, rand);
                const prod = normalize(x.M * y.M, x.p + y.p);
                expect(prod).not.toBeNull();
                expect(inSet(prod as { M: bigint; p: number }, r)).toBe(true);
            }
            return true;
        }, 60);
    });

    // The addition bound is only valid when both operands come from the SAME
    // set (the documented usage: both operands are the same native type, or
    // the same intermediate expression type). See the heterogeneous case
    // below for the upstream limitation.
    it('add bsn contains the exact bigint sum for identical operand sets', () => {
        const rand = seededRandom(0x0add5eed);
        check(anySet, a => {
            const r = toBSP(a).add(toBSP(a)).bsn;
            for (let k = 0; k < 12; ++k) {
                const x = sampleSet(a, rand);
                const y = sampleSet(a, rand);
                const p = Math.min(x.p, y.p);
                const M = (x.M << BigInt(x.p - p)) + (y.M << BigInt(y.p - p));
                const sum = normalize(M, p);
                if (sum === null) { continue; }   // exact cancellation
                expect(inSet(sum, r)).toBe(true);
            }
            return true;
        }, 60);
    });

    // Upstream limitation, preserved by the port: operator+ computes maxBits
    // from the maxExponent of one operand and the minExponent of the OTHER,
    // so when the operand with the larger maxExponent also has the smaller
    // minExponent the bit count is under-estimated. Pinned here so a future
    // "cleanup" of the formula has to be a deliberate deviation.
    it('add bsn under-estimates maxBits for heterogeneous operand sets', () => {
        const f = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const d = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        const r = d.add(f);
        expect(r.bsn.minExponent).toBe(-1074);
        expect(r.bsn.maxExponent).toBe(1023);
        expect(r.bsn.maxBits).toBe(1173);   // 1023 - (-149) + 1

        // A concrete exact sum needing more than 1173 bits: the smallest
        // positive double (2^-1074) plus a float with leading exponent 127.
        const x = { M: 1n, p: -1074 };   // in the IS_DOUBLE set
        const y = { M: 1n, p: 127 };     // in the IS_FLOAT set
        const sum = normalize((y.M << BigInt(y.p - x.p)) + x.M, x.p) as { M: bigint; p: number };
        expect(bitLength(sum.M)).toBe(1202);   // 127 + 1074 + 1
        expect(bitLength(sum.M)).toBeGreaterThan(r.bsn.maxBits);
        expect(inSet(sum, r.bsn)).toBe(false);
    });

    // Same upstream defect seen from the other side: when the two operands
    // tie on maxExponent, operator+ still pairs one operand's maxExponent
    // with the other's minExponent, so a.add(b) != b.add(a).
    it('add is asymmetric when the operand maxExponents tie', () => {
        const a = new BSPrecision(-17, 20, 24);
        const b = new BSPrecision(-5, 20, 12);
        expect(a.add(b).bsn.maxBits).toBe(27);   // 20 - (-5) + 1, then +1 for overlap
        expect(b.add(a).bsn.maxBits).toBe(38);   // 20 - (-17) + 1, no overlap bump
    });
});
