import { describe, it, expect } from 'vitest';
import { EulerAngles, EulerResult } from '../src/EulerAngles.js';

describe('EulerResult', () => {
    it('has the four upstream enumerators in declaration order', () => {
        expect(EulerResult.INVALID).toBe(0);
        expect(EulerResult.UNIQUE).toBe(1);
        expect(EulerResult.NOT_UNIQUE_SUM).toBe(2);
        expect(EulerResult.NOT_UNIQUE_DIF).toBe(3);
    });
});

describe('EulerAngles construction', () => {
    it('default constructor zeroes members and is INVALID', () => {
        const ea = new EulerAngles();
        expect(ea.axis).toEqual([0, 0, 0]);
        expect(ea.angle).toEqual([0, 0, 0]);
        expect(ea.result).toBe(EulerResult.INVALID);
    });

    it('6-argument constructor stores axes/angles and is UNIQUE', () => {
        const ea = new EulerAngles(2, 0, 1, 0.1, -0.2, 0.3);
        expect(ea.axis).toEqual([2, 0, 1]);
        expect(ea.angle).toEqual([0.1, -0.2, 0.3]);
        expect(ea.result).toBe(EulerResult.UNIQUE);
    });

    it('members are mutable, matching the upstream public fields', () => {
        const ea = new EulerAngles();
        ea.axis = [0, 1, 0];
        ea.angle = [Math.PI, 0, -Math.PI];
        ea.result = EulerResult.NOT_UNIQUE_SUM;
        expect(ea.axis).toEqual([0, 1, 0]);
        expect(ea.angle[0]).toBe(Math.PI);
        expect(ea.result).toBe(EulerResult.NOT_UNIQUE_SUM);
    });
});
