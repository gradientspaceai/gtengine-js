import { describe, it, expect } from 'vitest';
import { Frustum3 } from '../src/Frustum3.js';
import { Vector, sub, dot } from '../src/Vector.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, positive,
    rotationFrame, vector } from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

describe('Frustum3 construction', () => {
    it('the default constructor uses the documented values', () => {
        const frustum = new Frustum3();
        expect(frustum.origin.values).toEqual([0, 0, 0]);
        expect(frustum.dVector.values).toEqual([0, 0, 1]);
        expect(frustum.uVector.values).toEqual([0, 1, 0]);
        expect(frustum.rVector.values).toEqual([1, 0, 0]);
        expect(frustum.dMin).toBe(1);
        expect(frustum.dMax).toBe(2);
        expect(frustum.uBound).toBe(1);
        expect(frustum.rBound).toBe(1);
    });

    it('the default constructor calls update()', () => {
        const frustum = new Frustum3();
        expect(frustum.getDRatio()).toBe(2);
        expect(frustum.getMTwoUF()).toBe(-4);
        expect(frustum.getMTwoRF()).toBe(-4);
    });

    it('fromParameters copies the vectors and updates the derived values', () => {
        const origin = v3(1, 2, 3);
        const d = v3(0, 0, 1);
        const u = v3(0, 1, 0);
        const r = v3(1, 0, 0);
        const frustum = Frustum3.fromParameters(origin, d, u, r, 2, 10, 3, 4);
        origin.set(0, 99);
        d.set(2, 99);
        expect(frustum.origin.values).toEqual([1, 2, 3]);
        expect(frustum.dVector.values).toEqual([0, 0, 1]);
        expect(frustum.getDRatio()).toBe(5);
        expect(frustum.getMTwoUF()).toBe(-60);
        expect(frustum.getMTwoRF()).toBe(-80);
    });

    it('rejects vectors that are not 3D', () => {
        expect(() => Frustum3.fromParameters(Vector.fromArray([0, 0]),
            v3(0, 0, 1), v3(0, 1, 0), v3(1, 0, 0), 1, 2, 1, 1)).toThrow();
    });

    it('clone is a deep copy', () => {
        const frustum = new Frustum3();
        const copy = frustum.clone();
        copy.origin.set(0, 5);
        copy.dMax = 8;
        copy.update();
        expect(frustum.origin.values).toEqual([0, 0, 0]);
        expect(frustum.getDRatio()).toBe(2);
        expect(copy.getDRatio()).toBe(8);
    });
});

describe('Frustum3 update', () => {
    it('recomputes the derived values from the current members', () => {
        const frustum = new Frustum3();
        frustum.dMin = 2;
        frustum.dMax = 6;
        frustum.uBound = 3;
        frustum.rBound = 5;
        // The derived values are stale until update() is called.
        expect(frustum.getDRatio()).toBe(2);
        frustum.update();
        expect(frustum.getDRatio()).toBe(3);
        expect(frustum.getMTwoUF()).toBe(-36);
        expect(frustum.getMTwoRF()).toBe(-60);
    });
});

describe('Frustum3 computeVertices', () => {
    it('gives the hand-computed vertices of the default frustum', () => {
        const vertex = new Frustum3().computeVertices();
        expect(vertex.length).toBe(8);
        expect(vertex[0].values).toEqual([-1, -1, 1]);
        expect(vertex[1].values).toEqual([1, -1, 1]);
        expect(vertex[2].values).toEqual([1, 1, 1]);
        expect(vertex[3].values).toEqual([-1, 1, 1]);
        expect(vertex[4].values).toEqual([-2, -2, 2]);
        expect(vertex[5].values).toEqual([2, -2, 2]);
        expect(vertex[6].values).toEqual([2, 2, 2]);
        expect(vertex[7].values).toEqual([-2, 2, 2]);
    });

    it('translates with the origin', () => {
        const frustum = Frustum3.fromParameters(v3(10, 20, 30), v3(0, 0, 1),
            v3(0, 1, 0), v3(1, 0, 0), 1, 2, 1, 1);
        const vertex = frustum.computeVertices();
        expect(vertex[0].values).toEqual([9, 19, 31]);
        expect(vertex[6].values).toEqual([12, 22, 32]);
    });

    it('the vertices satisfy the frustum definition for a general frustum', () => {
        // E + n*D + s0*u*U + s1*r*R for the near plane; the far plane scales
        // the U and R offsets by f/n and uses f*D.
        const n = 2;
        const f = 5;
        const u = 3;
        const r = 4;
        const frustum = Frustum3.fromParameters(v3(1, -1, 2), v3(0, 1, 0),
            v3(0, 0, 1), v3(1, 0, 0), n, f, u, r);
        const vertex = frustum.computeVertices();
        const ratio = f / n;
        for (let i = 0; i < 4; ++i) {
            const near = sub(vertex[i], frustum.origin);
            expect(dot(near, frustum.dVector)).toBeCloseTo(n, 12);
            expect(Math.abs(dot(near, frustum.uVector))).toBeCloseTo(u, 12);
            expect(Math.abs(dot(near, frustum.rVector))).toBeCloseTo(r, 12);

            const far = sub(vertex[i + 4], frustum.origin);
            expect(dot(far, frustum.dVector)).toBeCloseTo(f, 12);
            expect(Math.abs(dot(far, frustum.uVector)))
                .toBeCloseTo(ratio * u, 12);
            expect(Math.abs(dot(far, frustum.rVector)))
                .toBeCloseTo(ratio * r, 12);
        }
    });

    it('does not alias the frustum vectors', () => {
        const frustum = new Frustum3();
        const vertex = frustum.computeVertices();
        vertex[0].set(0, 99);
        expect(frustum.origin.values).toEqual([0, 0, 0]);
        expect(frustum.dVector.values).toEqual([0, 0, 1]);
        expect(frustum.computeVertices()[0].values).toEqual([-1, -1, 1]);
    });
});

describe('Frustum3 comparisons', () => {
    const base = new Frustum3();

    it('equals compares every public member', () => {
        expect(base.equals(new Frustum3())).toBe(true);
        expect(base.notEquals(new Frustum3())).toBe(false);

        const other = base.clone();
        other.rBound = 3;
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);
    });

    it('lessThan follows the upstream member order', () => {
        const fields: Array<(f: Frustum3) => void> = [
            f => { f.origin = v3(-1, 0, 0); },
            f => { f.dVector = v3(0, 0, 0.5); },
            f => { f.uVector = v3(0, 0.5, 0); },
            f => { f.rVector = v3(0.5, 0, 0); },
            f => { f.dMin = 0.5; },
            f => { f.dMax = 1.5; },
            f => { f.uBound = 0.5; },
            f => { f.rBound = 0.5; }
        ];
        for (const mutate of fields) {
            const smaller = base.clone();
            mutate(smaller);
            expect(smaller.lessThan(base)).toBe(true);
            expect(base.lessThan(smaller)).toBe(false);
            expect(base.greaterThan(smaller)).toBe(true);
        }
    });

    it('the derived comparisons are consistent', () => {
        const bigger = base.clone();
        bigger.rBound = 2;
        expect(base.lessThanOrEqual(bigger)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(bigger.greaterThanOrEqual(base)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});

describe('Frustum3 verification', () => {
    const frustum = () => fc.tuple(vector(3, -5, 5), rotationFrame(3),
        positive(3, 0.1), positive(5), positive(3), positive(3))
        .map(([e, frame, n, extra, u, r]) => Frustum3.fromParameters(e,
            frame[0], frame[1], frame[2], n, n + extra, u, r));
    const key = (f: Frustum3) => [...f.origin.values, ...f.dVector.values,
        ...f.uVector.values, ...f.rVector.values, f.dMin, f.dMax, f.uBound,
        f.rBound];

    it('update recomputes the derived quantities from the members', () => {
        check(frustum(), f => {
            expectClose(f.getDRatio(), f.dMax / f.dMin, 1e-12, 1e-12);
            expect(f.getMTwoUF()).toBe(-2 * f.uBound * f.dMax);
            expect(f.getMTwoRF()).toBe(-2 * f.rBound * f.dMax);
        });
    });

    it('the derived quantities are stale until update is called', () => {
        // Upstream documents this: the accessors are valid only after
        // Update() has been called following a member change.
        check(fc.tuple(frustum(), positive(5)), ([f, newMax]) => {
            const before = f.getDRatio();
            f.dMax = f.dMin + newMax;
            expect(f.getDRatio()).toBe(before);
            f.update();
            expectClose(f.getDRatio(), f.dMax / f.dMin, 1e-12, 1e-12);
            expect(f.getMTwoUF()).toBe(-2 * f.uBound * f.dMax);
        });
    });

    it('computeVertices puts four corners on each of the two planes', () => {
        check(frustum(), f => {
            const vertex = f.computeVertices();
            expect(vertex.length).toBe(8);
            for (let i = 0; i < 4; ++i) {
                const near = sub(vertex[i], f.origin);
                const far = sub(vertex[i + 4], f.origin);
                // Distances along D are the near and far extents.
                expectClose(dot(f.dVector, near), f.dMin, 1e-9, 1e-9);
                expectClose(dot(f.dVector, far), f.dMax, 1e-9, 1e-9);
                // The U and R offsets scale by f/n between the planes.
                const scale = f.getDRatio();
                expectClose(dot(f.uVector, far),
                    scale * dot(f.uVector, near), 1e-9, 1e-9);
                expectClose(dot(f.rVector, far),
                    scale * dot(f.rVector, near), 1e-9, 1e-9);
                expect(Math.abs(dot(f.uVector, near)) - f.uBound)
                    .toBeLessThan(1e-9);
                expect(Math.abs(dot(f.rVector, near)) - f.rBound)
                    .toBeLessThan(1e-9);
            }
        });
    });

    it('the near-plane corners follow the (-u,-r), (-u,+r), (+u,+r), (+u,-r) '
        + 'order', () => {
            check(frustum(), f => {
                const signs = [[-1, -1], [-1, 1], [1, 1], [1, -1]];
                const vertex = f.computeVertices();
                for (let i = 0; i < 4; ++i) {
                    const near = sub(vertex[i], f.origin);
                    expectClose(dot(f.uVector, near),
                        signs[i][0] * f.uBound, 1e-9, 1e-9);
                    expectClose(dot(f.rVector, near),
                        signs[i][1] * f.rBound, 1e-9, 1e-9);
                }
            });
        });

    it('the comparisons follow the upstream member order', () => {
        check(fc.tuple(frustum(), frustum()), ([a, b]) => {
            const cmp = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(cmp < 0);
            expect(a.greaterThan(b)).toBe(cmp > 0);
            expect(a.lessThanOrEqual(b)).toBe(cmp <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(cmp >= 0);
            expect(a.equals(b)).toBe(cmp === 0);
        });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(frustum(), { minLength: 4, maxLength: 5 }), fs => {
            expectStrictWeakOrder(fs, (x, y) => x.lessThan(y));
        }, 30);
    });

    it('equals is element equality, so a NaN member breaks self-equality',
        () => {
            const f = Frustum3.fromParameters(Vector.fromArray([NaN, 0, 0]),
                Vector.fromArray([0, 0, 1]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([1, 0, 0]), 1, 2, 1, 1);
            expect(f.equals(f)).toBe(false);
            expect(f.lessThan(f)).toBe(false);
            expect(f.lessThanOrEqual(f)).toBe(true);
        });

    it('the factory and clone copy every vector and the derived values', () => {
        check(frustum(), f => {
            const cloned = f.clone();
            f.origin.set(0, 999);
            expect(cloned.origin.get(0)).not.toBe(999);
            expect(cloned.getDRatio()).toBe(f.getDRatio());
            expect(cloned.getMTwoUF()).toBe(f.getMTwoUF());
            expect(cloned.getMTwoRF()).toBe(f.getMTwoRF());
        }, 50);
    });
});
