import { describe, it, expect } from 'vitest';
import {
    getContainerOrientedBox3,
    inContainerOrientedBox3,
    mergeContainersOrientedBox3
} from '../src/ContOrientedBox3';
import { OrientedBox, type OrientedBox3 } from '../src/OrientedBox';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A right-handed orthonormal frame from a rotation about the given axis.
function frame(axis: Vector, angle: number): Vector[] {
    const w = axis.clone();
    normalize(w);
    // Rodrigues rotation of the standard basis about w.
    const rot = (u: Vector) => {
        const c = Math.cos(angle), s = Math.sin(angle);
        return add(add(mul(c, u), mul(s, cross(w, u))),
            mul((1 - c) * dot(w, u), w));
    };
    return [rot(v3(1, 0, 0)), rot(v3(0, 1, 0)), rot(v3(0, 0, 1))];
}

function box3(center: Vector, axis: Vector[], e: Vector): OrientedBox3 {
    return OrientedBox.fromCenterAxisExtent(center, axis, e);
}

// Independent containment check with a tolerance.
function containedWithin(point: Vector, box: OrientedBox3,
    tol: number): boolean {
    const diff = sub(point, box.center);
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(dot(diff, box.axis[i])) > box.extent.values[i] + tol) {
            return false;
        }
    }
    return true;
}

let seed = 24681357;
function rand(): number {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
}

describe('ContOrientedBox3', () => {
    describe('inContainerOrientedBox3', () => {
        it('matches an axis-aligned box', () => {
            const box = box3(v3(1, 2, 3),
                [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)], v3(2, 3, 4));
            expect(inContainerOrientedBox3(v3(1, 2, 3), box)).toBe(true);
            expect(inContainerOrientedBox3(v3(3, 5, 7), box)).toBe(true);
            expect(inContainerOrientedBox3(v3(-1, -1, -1), box)).toBe(true);
            expect(inContainerOrientedBox3(v3(3.001, 2, 3), box)).toBe(false);
            expect(inContainerOrientedBox3(v3(1, 2, -1.001), box)).toBe(false);
        });

        it('matches a rotated box', () => {
            const axis = frame(v3(1, 1, 1), 0.9);
            const box = box3(v3(-2, 0, 5), axis, v3(1, 2, 0.5));
            for (const vertex of box.getVertices()) {
                // The corners themselves are on the boundary, where round-off
                // can push a projection just past the extent; test points
                // slightly inside and slightly outside instead.
                const inward = add(box.center,
                    mul(0.999, sub(vertex, box.center)));
                expect(inContainerOrientedBox3(inward, box)).toBe(true);
                const outward = add(box.center,
                    mul(1.001, sub(vertex, box.center)));
                expect(inContainerOrientedBox3(outward, box)).toBe(false);
            }
        });

        it('rejects mismatched dimensions', () => {
            const box = new OrientedBox(3);
            expect(() => inContainerOrientedBox3(Vector.fromArray([0, 0]),
                box)).toThrow();
        });
    });

    describe('getContainerOrientedBox3', () => {
        it('recovers an axis-aligned box from its corner points', () => {
            const points: Vector[] = [];
            for (let i = 0; i < 8; ++i) {
                points.push(v3(
                    (i & 1) ? 4 : 0, (i & 2) ? 6 : 0, (i & 4) ? 2 : 0));
            }
            const box = getContainerOrientedBox3(points);
            expect(box).not.toBeNull();
            expect(box!.center.values[0]).toBeCloseTo(2, 10);
            expect(box!.center.values[1]).toBeCloseTo(3, 10);
            expect(box!.center.values[2]).toBeCloseTo(1, 10);
            const sorted = [...box!.extent.values].sort((a, b) => a - b);
            expect(sorted[0]).toBeCloseTo(1, 10);
            expect(sorted[1]).toBeCloseTo(2, 10);
            expect(sorted[2]).toBeCloseTo(3, 10);
        });

        it('throws on an empty point set', () => {
            expect(() => getContainerOrientedBox3([])).toThrow();
        });

        it('rejects non-3D points', () => {
            expect(() => getContainerOrientedBox3([Vector.fromArray([0, 0])]))
                .toThrow();
        });

        it('contains every input point for random clouds', () => {
            for (let trial = 0; trial < 40; ++trial) {
                const axis = frame(v3(2 * rand() - 1, 2 * rand() - 1,
                    2 * rand() - 1 + 0.1), 2 * Math.PI * rand());
                const c = v3(10 * rand() - 5, 10 * rand() - 5, 10 * rand() - 5);
                const e = v3(0.5 + 3 * rand(), 0.5 + 3 * rand(),
                    0.5 + 3 * rand());
                const points: Vector[] = [];
                for (let i = 0; i < 60; ++i) {
                    let p = c.clone();
                    for (let j = 0; j < 3; ++j) {
                        p = add(p, mul(e.values[j] * (2 * rand() - 1),
                            axis[j]));
                    }
                    points.push(p);
                }
                const box = getContainerOrientedBox3(points);
                expect(box).not.toBeNull();
                for (const p of points) {
                    expect(containedWithin(p, box!, 1e-9)).toBe(true);
                }
                // The box axes are orthonormal.
                for (let j = 0; j < 3; ++j) {
                    expect(length(box!.axis[j])).toBeCloseTo(1, 10);
                }
                expect(dot(box!.axis[0], box!.axis[1])).toBeCloseTo(0, 10);
                expect(dot(box!.axis[0], box!.axis[2])).toBeCloseTo(0, 10);
                expect(dot(box!.axis[1], box!.axis[2])).toBeCloseTo(0, 10);
            }
        });
    });

    describe('mergeContainersOrientedBox3', () => {
        it('reproduces the input when both boxes are the same', () => {
            const axis = frame(v3(0, 0, 1), 0.3);
            const box = box3(v3(1, -1, 2), axis, v3(2, 3, 1));
            const m = mergeContainersOrientedBox3(box, box);
            for (const vertex of box.getVertices()) {
                expect(containedWithin(vertex, m, 1e-9)).toBe(true);
            }
            // The merged volume matches the input volume.
            const vol = 8 * m.extent.values[0] * m.extent.values[1]
                * m.extent.values[2];
            expect(vol).toBeCloseTo(8 * 2 * 3 * 1, 8);
        });

        it('contains both input boxes (randomized)', () => {
            for (let trial = 0; trial < 60; ++trial) {
                const a0 = frame(v3(2 * rand() - 1, 2 * rand() - 1,
                    2 * rand() - 1 + 0.1), 2 * Math.PI * rand());
                const a1 = frame(v3(2 * rand() - 1, 2 * rand() - 1,
                    2 * rand() - 1 + 0.1), 2 * Math.PI * rand());
                const b0 = box3(v3(6 * rand() - 3, 6 * rand() - 3,
                    6 * rand() - 3), a0,
                    v3(0.5 + rand(), 0.5 + rand(), 0.5 + rand()));
                const b1 = box3(v3(6 * rand() - 3, 6 * rand() - 3,
                    6 * rand() - 3), a1,
                    v3(0.5 + rand(), 0.5 + rand(), 0.5 + rand()));
                const m = mergeContainersOrientedBox3(b0, b1);

                for (let j = 0; j < 3; ++j) {
                    expect(length(m.axis[j])).toBeCloseTo(1, 9);
                    expect(m.extent.values[j]).toBeGreaterThanOrEqual(0);
                }
                expect(dot(m.axis[0], m.axis[1])).toBeCloseTo(0, 9);
                expect(dot(m.axis[0], m.axis[2])).toBeCloseTo(0, 9);
                expect(dot(m.axis[1], m.axis[2])).toBeCloseTo(0, 9);

                for (const box of [b0, b1]) {
                    for (const vertex of box.getVertices()) {
                        expect(containedWithin(vertex, m, 1e-9)).toBe(true);
                    }
                }
            }
        });

        it('rejects mismatched dimensions', () => {
            expect(() => mergeContainersOrientedBox3(new OrientedBox(3),
                new OrientedBox(2))).toThrow();
        });
    });
});
