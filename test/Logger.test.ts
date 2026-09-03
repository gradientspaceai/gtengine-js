import { describe, it, expect } from 'vitest';
import { logAssert, logError } from '../src/Logger.js';

describe('Logger', () => {
    it('logAssert does nothing when the condition is true', () => {
        expect(() => logAssert(true, 'should not throw')).not.toThrow();
    });

    it('logAssert throws Error with the message when the condition is false', () => {
        expect(() => logAssert(false, 'Invalid maximum iterations.'))
            .toThrowError(new Error('Invalid maximum iterations.'));
    });

    it('logAssert narrows types (asserts condition)', () => {
        const maybe: number | null = Math.random() < 2 ? 3 : null;
        logAssert(maybe !== null, 'value must exist');
        // After the assertion, TypeScript treats 'maybe' as number.
        const doubled: number = maybe * 2;
        expect(doubled).toBe(6);
    });

    it('logError always throws Error with the message', () => {
        expect(() => logError('The input is invalid.'))
            .toThrowError(new Error('The input is invalid.'));
    });

    it('thrown values are instances of Error', () => {
        try {
            logError('boom');
            expect.unreachable();
        } catch (e) {
            expect(e).toBeInstanceOf(Error);
            expect((e as Error).message).toBe('boom');
        }
    });
});
