import { describe, it, expect } from 'vitest';
import {
    chebyshevRatio, chebyshevRatioUsingCosAngle,
    chebyshevRatios, chebyshevRatiosUsingCosAngle
} from '../src/ChebyshevRatio.js';
import { GTE_C_PI } from '../src/Constants.js';

describe('chebyshevRatio', () => {
    it('computes sin(t*A)/sin(A) for angles in (0,pi)', () => {
        for (const angle of [0.1, 0.5, 1, Math.PI / 2, 2, 3, 3.14]) {
            for (const t of [0, 0.25, 0.5, 0.75, 1]) {
                expect(chebyshevRatio(t, angle))
                    .toBe(Math.sin(t * angle) / Math.sin(angle));
            }
        }
    });

    it('returns t for angle 0 (removable singularity, l\'Hospital)', () => {
        for (const t of [0, 0.25, 0.5, 0.75, 1]) {
            expect(chebyshevRatio(t, 0)).toBe(t);
        }
    });

    it('approaches t continuously as the angle approaches 0', () => {
        for (const t of [0.25, 0.5, 0.75]) {
            expect(Math.abs(chebyshevRatio(t, 1e-8) - t)).toBeLessThan(1e-15);
        }
    });

    it('is exact at the endpoints of t', () => {
        for (const angle of [0.5, 1.5, 3]) {
            expect(chebyshevRatio(0, angle)).toBe(0);
            expect(chebyshevRatio(1, angle)).toBe(1);
        }
    });

    it('throws for angles outside [0,pi)', () => {
        expect(() => chebyshevRatio(0.5, -0.1)).toThrow('Invalid angle.');
        expect(() => chebyshevRatio(0.5, GTE_C_PI)).toThrow('Invalid angle.');
        expect(() => chebyshevRatio(0.5, 4)).toThrow('Invalid angle.');
    });
});

describe('chebyshevRatioUsingCosAngle', () => {
    it('agrees with chebyshevRatio for the angle acos(cosAngle)', () => {
        for (const cosAngle of [-0.99, -0.5, 0, 0.5, 0.99]) {
            const angle = Math.acos(cosAngle);
            for (const t of [0, 0.25, 0.5, 0.75, 1]) {
                expect(chebyshevRatioUsingCosAngle(t, cosAngle))
                    .toBe(Math.sin(t * angle) / Math.sin(angle));
            }
        }
    });

    it('returns t for cosAngle = 1 (angle 0)', () => {
        for (const t of [0, 0.3, 1]) {
            expect(chebyshevRatioUsingCosAngle(t, 1)).toBe(t);
            // Values beyond 1 also take the angle-0 branch, as upstream.
            expect(chebyshevRatioUsingCosAngle(t, 1.5)).toBe(t);
        }
    });

    it('throws for cosAngle <= -1 (angle pi)', () => {
        expect(() => chebyshevRatioUsingCosAngle(0.5, -1)).toThrow('Invalid angle.');
        expect(() => chebyshevRatioUsingCosAngle(0.5, -2)).toThrow('Invalid angle.');
    });
});

describe('chebyshevRatios', () => {
    it('returns the pair {f(1-t,A), f(t,A)}', () => {
        for (const angle of [0.25, 1, 2.5]) {
            for (const t of [0, 0.25, 0.5, 0.9, 1]) {
                const [f0, f1] = chebyshevRatios(t, angle);
                expect(f0).toBe(chebyshevRatio(1 - t, angle));
                expect(f1).toBe(chebyshevRatio(t, angle));
            }
        }
    });

    it('returns {1-t, t} for angle 0', () => {
        expect(chebyshevRatios(0.25, 0)).toEqual([0.75, 0.25]);
        expect(chebyshevRatios(1, 0)).toEqual([0, 1]);
    });

    it('gives slerp-style barycentric weights whose sum is 1/cos(A/2)', () => {
        // f(1/2,A) + f(1/2,A) = 2*sin(A/2)/sin(A) = 1/cos(A/2), which
        // approaches 1 as the angle approaches 0.
        for (const angle of [0.01, 0.1, 1]) {
            const [f0, f1] = chebyshevRatios(0.5, angle);
            expect(f0).toBe(f1);
            expect(Math.abs(f0 + f1 - 1 / Math.cos(angle / 2))).toBeLessThan(1e-14);
        }
    });

    it('throws for angles outside [0,pi)', () => {
        expect(() => chebyshevRatios(0.5, -1)).toThrow('Invalid angle.');
        expect(() => chebyshevRatios(0.5, GTE_C_PI)).toThrow('Invalid angle.');
    });
});

describe('chebyshevRatiosUsingCosAngle', () => {
    it('agrees with chebyshevRatios for the angle acos(cosAngle)', () => {
        for (const cosAngle of [-0.9, 0, 0.9]) {
            const angle = Math.acos(cosAngle);
            for (const t of [0, 0.25, 0.5, 1]) {
                expect(chebyshevRatiosUsingCosAngle(t, cosAngle))
                    .toEqual(chebyshevRatios(t, angle));
            }
        }
    });

    it('returns {1-t, t} for cosAngle >= 1 (angle 0)', () => {
        expect(chebyshevRatiosUsingCosAngle(0.25, 1)).toEqual([0.75, 0.25]);
        expect(chebyshevRatiosUsingCosAngle(0.25, 2)).toEqual([0.75, 0.25]);
    });

    it('throws for cosAngle <= -1 (angle pi)', () => {
        expect(() => chebyshevRatiosUsingCosAngle(0.5, -1)).toThrow('Invalid angle.');
    });
});
