import { describe, it, expect } from 'vitest';
import { dot, length } from '../src/Vector.js';
import { check, orthonormalFrame, unitVector, triangle, segment, alignedBox, expectClose } from './helpers/arbitraries.js';

describe('test helpers', () => {
    it('unitVector has length 1', () => {
        check(unitVector(3), v => { expectClose(length(v), 1, 1e-12); });
    });
    it('orthonormalFrame is orthonormal', () => {
        check(orthonormalFrame(3), f => {
            for (let i = 0; i < 3; ++i) {
                for (let j = 0; j < 3; ++j) {
                    expectClose(dot(f[i], f[j]), i === j ? 1 : 0, 1e-10);
                }
            }
        });
    });
    it('triangle, segment and alignedBox generators produce valid primitives', () => {
        check(triangle(2), t => { expect(t.v.length).toBe(3); });
        check(segment(3), s => { expect(s.p.length).toBe(2); });
        check(alignedBox(3), b => {
            for (let i = 0; i < 3; ++i) { expect(b.min.get(i)).toBeLessThanOrEqual(b.max.get(i)); }
        });
    });
});
