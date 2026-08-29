import { describe, it, expect } from 'vitest';
import { hashCombine, hashValue, hashValueWithSeed } from '../src/HashCombine';

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
