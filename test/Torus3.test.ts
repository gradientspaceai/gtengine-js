import { describe, it, expect } from 'vitest';
import { Torus3 } from '../src/Torus3';
import { Vector, dot, sub, length } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function expectClose(v: Vector, expected: readonly number[],
    digits = 12): void {
    expect(v.size).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(v.get(i)).toBeCloseTo(expected[i], digits);
    }
}

describe('Torus3 construction', () => {
    it('the default constructor uses the documented values', () => {
        const torus = new Torus3();
        expect(torus.center.values).toEqual([0, 0, 0]);
        expect(torus.direction0.values).toEqual([1, 0, 0]);
        expect(torus.direction1.values).toEqual([0, 1, 0]);
        expect(torus.normal.values).toEqual([0, 0, 1]);
        expect(torus.radius0).toBe(2);
        expect(torus.radius1).toBe(1);
    });

    it('fromCenterFrameRadii copies the vectors', () => {
        const center = v3(1, 2, 3);
        const d0 = v3(0, 1, 0);
        const d1 = v3(0, 0, 1);
        const n = v3(1, 0, 0);
        const torus = Torus3.fromCenterFrameRadii(center, d0, d1, n, 5, 2);
        center.set(0, 99);
        d0.set(1, 99);
        expect(torus.center.values).toEqual([1, 2, 3]);
        expect(torus.direction0.values).toEqual([0, 1, 0]);
        expect(torus.radius0).toBe(5);
        expect(torus.radius1).toBe(2);
    });

    it('rejects vectors that are not 3D', () => {
        expect(() => Torus3.fromCenterFrameRadii(Vector.fromArray([0, 0]),
            v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1), 2, 1)).toThrow();
    });

    it('clone is a deep copy', () => {
        const torus = new Torus3();
        const copy = torus.clone();
        copy.center.set(0, 5);
        copy.radius1 = 0.5;
        expect(torus.center.values).toEqual([0, 0, 0]);
        expect(torus.radius1).toBe(1);
    });
});

describe('Torus3 evaluate', () => {
    const torus = new Torus3(); // r0 = 2, r1 = 1, standard frame at origin.

    it('maxOrder 0 gives only the position', () => {
        const jet = torus.evaluate(0, 0, 0);
        expect(jet.length).toBe(1);
        expectClose(jet[0], [3, 0, 0]);
    });

    it('gives the hand-computed positions at the axis parameters', () => {
        expectClose(torus.evaluate(Math.PI / 2, 0, 0)[0], [0, 3, 0]);
        expectClose(torus.evaluate(0, Math.PI / 2, 0)[0], [2, 0, 1]);
        expectClose(torus.evaluate(0, Math.PI, 0)[0], [1, 0, 0]);
        expectClose(torus.evaluate(Math.PI, Math.PI, 0)[0], [-1, 0, 0]);
        expectClose(torus.evaluate(Math.PI / 2, -Math.PI / 2, 0)[0],
            [0, 2, -1]);
    });

    it('maxOrder 1 gives the hand-computed first derivatives', () => {
        // At (u,v) = (0,0): X = (3,0,0), dX/du = 3*(0,1,0), dX/dv = (0,0,1).
        const jet = torus.evaluate(0, 0, 1);
        expect(jet.length).toBe(3);
        expectClose(jet[0], [3, 0, 0]);
        expectClose(jet[1], [0, 3, 0]);
        expectClose(jet[2], [0, 0, 1]);
    });

    it('maxOrder 2 gives the hand-computed second derivatives', () => {
        // At (u,v) = (0,0): d2X/du2 = -(3,0,0), d2X/dudv = 0,
        // d2X/dv2 = -(1,0,0).
        const jet = torus.evaluate(0, 0, 2);
        expect(jet.length).toBe(6);
        expectClose(jet[3], [-3, 0, 0]);
        expectClose(jet[4], [0, 0, 0]);
        expectClose(jet[5], [-1, 0, 0]);
    });

    it('the derivatives match central finite differences', () => {
        const general = Torus3.fromCenterFrameRadii(v3(1, -2, 3),
            v3(0, 1, 0), v3(0, 0, 1), v3(1, 0, 0), 4, 1.5);
        const h = 1e-5;
        const u = 0.7;
        const v = -1.3;
        const jet = general.evaluate(u, v, 2);

        const pos = (a: number, b: number): Vector =>
            general.evaluate(a, b, 0)[0];

        const dU = sub(pos(u + h, v), pos(u - h, v));
        const dV = sub(pos(u, v + h), pos(u, v - h));
        for (let i = 0; i < 3; ++i) {
            expect(jet[1].get(i)).toBeCloseTo(dU.get(i) / (2 * h), 7);
            expect(jet[2].get(i)).toBeCloseTo(dV.get(i) / (2 * h), 7);
        }

        // Second derivatives from the first-derivative jets.
        const jetUp = general.evaluate(u + h, v, 1);
        const jetUm = general.evaluate(u - h, v, 1);
        const jetVp = general.evaluate(u, v + h, 1);
        const jetVm = general.evaluate(u, v - h, 1);
        for (let i = 0; i < 3; ++i) {
            expect(jet[3].get(i)).toBeCloseTo(
                (jetUp[1].get(i) - jetUm[1].get(i)) / (2 * h), 6);
            expect(jet[4].get(i)).toBeCloseTo(
                (jetVp[1].get(i) - jetVm[1].get(i)) / (2 * h), 6);
            expect(jet[5].get(i)).toBeCloseTo(
                (jetVp[2].get(i) - jetVm[2].get(i)) / (2 * h), 6);
        }
    });

    it('evaluated points satisfy the implicit torus equation', () => {
        const general = Torus3.fromCenterFrameRadii(v3(1, -2, 3),
            v3(0, 1, 0), v3(0, 0, 1), v3(1, 0, 0), 4, 1.5);
        const r0 = general.radius0;
        const r1 = general.radius1;
        for (let i = 0; i < 8; ++i) {
            for (let j = 0; j < 8; ++j) {
                const u = (i / 8) * 2 * Math.PI;
                const v = (j / 8) * 2 * Math.PI;
                const p = general.evaluate(u, v, 0)[0];
                const delta = sub(p, general.center);
                const d2 = dot(delta, delta);
                const nd = dot(general.normal, delta);
                const term = d2 + r0 * r0 - r1 * r1;
                const value = term * term - 4 * r0 * r0 * (d2 - nd * nd);
                expect(value).toBeCloseTo(0, 9);
            }
        }
    });

    it('the tube center is at distance r0 and the point at r1 from it', () => {
        const torusG = Torus3.fromCenterFrameRadii(v3(0, 0, 0), v3(1, 0, 0),
            v3(0, 1, 0), v3(0, 0, 1), 5, 2);
        for (let i = 0; i < 8; ++i) {
            const u = (i / 8) * 2 * Math.PI;
            const tubeCenter = v3(5 * Math.cos(u), 5 * Math.sin(u), 0);
            expect(length(tubeCenter)).toBeCloseTo(5, 12);
            for (let j = 0; j < 8; ++j) {
                const v = (j / 8) * 2 * Math.PI;
                const p = torusG.evaluate(u, v, 0)[0];
                expect(length(sub(p, tubeCenter))).toBeCloseTo(2, 12);
            }
        }
    });
});

describe('Torus3 getParameters', () => {
    it('inverts evaluate for parameters in (-pi, pi)', () => {
        const torus = Torus3.fromCenterFrameRadii(v3(1, -2, 3), v3(0, 1, 0),
            v3(0, 0, 1), v3(1, 0, 0), 4, 1.5);
        const samples = [-3, -2, -1, -0.25, 0, 0.4, 1.1, 2.5, 3];
        for (const u of samples) {
            for (const v of samples) {
                const p = torus.evaluate(u, v, 0)[0];
                const params = torus.getParameters(p);
                expect(params.u).toBeCloseTo(u, 10);
                expect(params.v).toBeCloseTo(v, 10);
            }
        }
    });

    it('returns values in (-pi, pi] and wraps parameters outside it', () => {
        const torus = new Torus3();
        const p = torus.evaluate(3 * Math.PI / 2, 0, 0)[0];
        const params = torus.getParameters(p);
        expect(params.u).toBeCloseTo(-Math.PI / 2, 12);
        expect(params.v).toBeCloseTo(0, 12);
    });

    it('gives the hand-computed parameters of a known point', () => {
        // Default torus: (0,3,0) is at u = pi/2, v = 0.
        const torus = new Torus3();
        const params = torus.getParameters(v3(0, 3, 0));
        expect(params.u).toBeCloseTo(Math.PI / 2, 12);
        expect(params.v).toBeCloseTo(0, 12);

        // (2,0,1) is at u = 0, v = pi/2.
        const params2 = torus.getParameters(v3(2, 0, 1));
        expect(params2.u).toBe(0);
        expect(params2.v).toBeCloseTo(Math.PI / 2, 12);
    });
});

describe('Torus3 comparisons', () => {
    const base = new Torus3();

    it('equals compares every member', () => {
        expect(base.equals(new Torus3())).toBe(true);
        expect(base.notEquals(new Torus3())).toBe(false);

        const other = base.clone();
        other.radius1 = 0.5;
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);
    });

    it('lessThan follows the upstream member order', () => {
        const mutations: Array<(t: Torus3) => void> = [
            t => { t.center = v3(-1, 0, 0); },
            t => { t.direction0 = v3(0.5, 0, 0); },
            t => { t.direction1 = v3(0, 0.5, 0); },
            t => { t.normal = v3(0, 0, 0.5); },
            t => { t.radius0 = 1; },
            t => { t.radius1 = 0.5; }
        ];
        for (const mutate of mutations) {
            const smaller = base.clone();
            mutate(smaller);
            expect(smaller.lessThan(base)).toBe(true);
            expect(base.lessThan(smaller)).toBe(false);
            expect(base.greaterThan(smaller)).toBe(true);
        }
    });

    it('the derived comparisons are consistent', () => {
        const smaller = base.clone();
        smaller.radius1 = 0.5;
        expect(smaller.lessThanOrEqual(base)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(base.greaterThanOrEqual(smaller)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});
