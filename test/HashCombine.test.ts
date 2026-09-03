import { describe, it, expect } from 'vitest';
import { check, finite, seededRandom, fc } from './helpers/arbitraries.js';
import { hashCombine, hashValue, hashValueWithSeed } from '../src/HashCombine.js';

describe('HashCombine', () => {
    it('hashCombine(0, 0) applies the magic constant', () => {
        // hashNumber(0) == 0, so the result is 0 ^ (0 + 0x9e3779b9 + 0 + 0).
        expect(hashCombine(0, 0)).toBe(0x9e3779b9);
    });

    it('produces unsigned 32-bit integers', () => {
        const values = [0, 1, -1, 0.5, 1e300, -1e-300, Math.PI, Number.MAX_VALUE];
        let seed = 0;
        for (const v of values) {
            seed = hashCombine(seed, v);
            expect(Number.isInteger(seed)).toBe(true);
            expect(seed).toBeGreaterThanOrEqual(0);
            expect(seed).toBeLessThanOrEqual(0xFFFFFFFF);
        }
    });

    it('is deterministic', () => {
        expect(hashValue(1, 2, 3)).toBe(hashValue(1, 2, 3));
        expect(hashCombine(12345, Math.E)).toBe(hashCombine(12345, Math.E));
    });

    it('is order sensitive', () => {
        expect(hashValue(1, 2)).not.toBe(hashValue(2, 1));
    });

    it('distinguishes nearby values', () => {
        expect(hashValue(1)).not.toBe(hashValue(1 + Number.EPSILON));
        expect(hashValue(0)).not.toBe(hashValue(Number.MIN_VALUE));
        expect(hashValue(2)).not.toBe(hashValue(-2));
    });

    it('hashes +0 and -0 identically (they compare equal)', () => {
        expect(hashValue(0)).toBe(hashValue(-0));
    });

    it('hashValue folds values like repeated hashCombine from seed 0', () => {
        const expected = hashCombine(hashCombine(hashCombine(0, 1.5), -2.25), 1e9);
        expect(hashValue(1.5, -2.25, 1e9)).toBe(expected);
    });

    it('hashValueWithSeed starts from the given seed', () => {
        const seed = 0xDEADBEEF;
        expect(hashValueWithSeed(seed)).toBe(seed);
        expect(hashValueWithSeed(seed, 7)).toBe(hashCombine(seed, 7));
        expect(hashValueWithSeed(seed, 7, 8)).toBe(hashCombine(hashCombine(seed, 7), 8));
    });

    it('hashValue of no arguments is the zero seed', () => {
        expect(hashValue()).toBe(0);
    });
});

describe('HashCombine verification', () => {
    // Independent reference implementation of the upstream formula
    //   seed ^= hash(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2)
    // evaluated in exact 32-bit unsigned arithmetic with BigInt, using a
    // separately written DataView-based fold of the binary64 bits as the
    // stand-in for std::hash<double>.
    const M32 = (1n << 32n) - 1n;

    function refHashNumber(value: number): bigint {
        const dv = new DataView(new ArrayBuffer(8));
        dv.setFloat64(0, value === 0 ? 0 : value, true);
        const lo = BigInt(dv.getUint32(0, true));
        const hi = BigInt(dv.getUint32(4, true));
        return lo ^ hi;
    }

    function refCombine(seed: bigint, value: number): bigint {
        const s = seed & M32;
        const term = (refHashNumber(value) + 0x9e3779b9n + ((s << 6n) & M32) + (s >> 2n)) & M32;
        return (s ^ term) & M32;
    }

    it('hashCombine matches the exact 32-bit reference formula', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 0xFFFFFFFF }), finite(-1e6, 1e6)),
            ([seed, value]) =>
                BigInt(hashCombine(seed, value)) === refCombine(BigInt(seed), value));
    });

    it('hashCombine matches the reference on extreme values too', () => {
        const specials = [0, -0, 1, -1, Number.MIN_VALUE, Number.MAX_VALUE,
            Number.EPSILON, 1e300, -1e-300, Math.PI, NaN, Infinity, -Infinity];
        check(fc.tuple(fc.integer({ min: 0, max: 0xFFFFFFFF }),
            fc.constantFrom(...specials)),
            ([seed, value]) =>
                BigInt(hashCombine(seed, value)) === refCombine(BigInt(seed), value));
    });

    it('hashValue folds left over the argument list from seed 0', () => {
        check(fc.array(finite(-1e6, 1e6), { minLength: 0, maxLength: 8 }), values => {
            let seed = 0;
            for (const v of values) { seed = hashCombine(seed, v); }
            return hashValue(...values) === seed;
        });
    });

    it('hashValueWithSeed folds left from the given seed', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 0xFFFFFFFF }),
            fc.array(finite(-1e6, 1e6), { minLength: 0, maxLength: 6 })),
            ([seed, values]) => {
                let expected = seed >>> 0;
                for (const v of values) { expected = hashCombine(expected, v); }
                return hashValueWithSeed(seed, ...values) === expected;
            });
    });

    it('every result is an unsigned 32-bit integer', () => {
        check(fc.array(finite(-1e9, 1e9), { minLength: 0, maxLength: 8 }), values => {
            const h = hashValue(...values);
            return Number.isInteger(h) && h >= 0 && h <= 0xFFFFFFFF;
        });
    });

    it('is deterministic across calls', () => {
        check(fc.array(finite(-1e6, 1e6), { minLength: 1, maxLength: 6 }), values =>
            hashValue(...values) === hashValue(...values));
    });

    it('+0 and -0 hash alike (they compare equal, as in C++)', () => {
        check(fc.integer({ min: 0, max: 0xFFFFFFFF }), seed =>
            hashCombine(seed, 0) === hashCombine(seed, -0));
    });

    it('distinct values usually produce distinct hashes (collision sanity)', () => {
        // Not a cryptographic claim; catches a hash that ignores its input.
        const rand = seededRandom(0x9E3779B9);
        const seen = new Set<number>();
        for (let k = 0; k < 2000; ++k) {
            seen.add(hashValue(rand() * 1e6, rand() * 1e-6));
        }
        expect(seen.size).toBeGreaterThan(1900);
    });
});
