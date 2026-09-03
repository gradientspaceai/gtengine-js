import { describe, it, expect } from 'vitest';
import { check, fc } from './helpers/arbitraries.js';
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

describe('Logger verification', () => {
    it('logAssert throws exactly when the condition is false', () => {
        check(fc.tuple(fc.boolean(), fc.string()), ([cond, msg]) => {
            let threw = false;
            try {
                logAssert(cond, msg);
            } catch {
                threw = true;
            }
            return threw === !cond;
        });
    });

    it('the thrown message is the caller message verbatim (no file/line prefix)', () => {
        check(fc.string(), msg => {
            try {
                logError(msg);
            } catch (e) {
                return e instanceof Error && e.message === msg;
            }
            return false;
        });
    });

    it('logAssert(false, msg) and logError(msg) throw identical messages', () => {
        check(fc.string(), msg => {
            let m0 = '', m1 = '';
            try { logAssert(false, msg); } catch (e) { m0 = (e as Error).message; }
            try { logError(msg); } catch (e) { m1 = (e as Error).message; }
            return m0 === m1 && m0 === msg;
        });
    });

    it('logError never returns (upstream GTE_ERROR always throws)', () => {
        check(fc.string(), msg => {
            let returned = false;
            try {
                logError(msg);
                returned = true;
            } catch {
                // expected
            }
            return !returned;
        });
    });
});
