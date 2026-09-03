import { describe, it, expect } from 'vitest';
import {
    getContainerEllipsoid3,
    inContainerEllipsoid3,
    mergeContainersEllipsoid3
} from '../src/ContEllipsoid3.js';
import { Hyperellipsoid, type Ellipsoid3 } from '../src/Hyperellipsoid.js';
import { Line } from '../src/Line.js';
import { projectEllipsoid3 } from '../src/Projection.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A right-handed orthonormal frame: the standard basis rotated about 'axis'.
function frame(axis: Vector, angle: number): Vector[] {
    const w = axis.clone();
    normalize(w);
    const rot = (u: Vector) => {
        const c = Math.cos(angle), s = Math.sin(angle);
        return add(add(mul(c, u), mul(s, cross(w, u))),
            mul((1 - c) * dot(w, u), w));
    };
    return [rot(v3(1, 0, 0)), rot(v3(0, 1, 0)), rot(v3(0, 0, 1))];
}

function ellipsoid3(center: Vector, axis: Vector[], e: Vector): Ellipsoid3 {
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, e);
}

// Independent evaluation of the ellipsoid quadratic form: <= 1 inside.
function quadratic(p: Vector, e: Ellipsoid3): number {
    const d = sub(p, e.center);
    let sum = 0;
    for (let j = 0; j < 3; ++j) {
        const t = dot(d, e.axis[j]) / e.extent.values[j];
        sum += t * t;
    }
    return sum;
}

let seed = 777333111;
function rand(): number {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
}

describe('ContEllipsoid3', () => {
    describe('inContainerEllipsoid3', () => {
        it('matches the unit sphere', () => {
            const e = new Hyperellipsoid(3);
            expect(inContainerEllipsoid3(v3(0, 0, 0), e)).toBe(true);
            expect(inContainerEllipsoid3(v3(1, 0, 0), e)).toBe(true);
            expect(inContainerEllipsoid3(v3(0, 0, -1), e)).toBe(true);
            expect(inContainerEllipsoid3(v3(1.0001, 0, 0), e)).toBe(false);
            // sqrt(3) * 0.577 < 1 < sqrt(3) * 0.578.
            expect(inContainerEllipsoid3(v3(0.577, 0.577, 0.577), e))
                .toBe(true);
            expect(inContainerEllipsoid3(v3(0.578, 0.578, 0.578), e))
                .toBe(false);
        });

        it('matches a rotated, translated, anisotropic ellipsoid', () => {
            const e = ellipsoid3(v3(2, -3, 1), frame(v3(1, 2, -1), 0.8),
                v3(3, 1, 0.5));
            for (let j = 0; j < 3; ++j) {
                // Just inside the axis endpoint (exactly on the boundary the
                // comparison is at the mercy of round-off).
                const tip = add(e.center,
                    mul(e.extent.values[j] * 0.999, e.axis[j]));
                expect(inContainerEllipsoid3(tip, e)).toBe(true);
                const beyond = add(e.center,
                    mul(e.extent.values[j] * 1.001, e.axis[j]));
                expect(inContainerEllipsoid3(beyond, e)).toBe(false);
            }
        });

        it('rejects mismatched dimensions', () => {
            expect(() => inContainerEllipsoid3(Vector.fromArray([0, 0]),
                new Hyperellipsoid(3))).toThrow();
        });
    });

    describe('getContainerEllipsoid3', () => {
        it('contains every input point and is tight on at least one', () => {
            const points: Vector[] = [];
            for (let i = 0; i < 8; ++i) {
                points.push(v3((i & 1) ? 4 : -4, (i & 2) ? 2 : -2,
                    (i & 4) ? 1 : -1));
            }
            points.push(v3(0, 0, 0), v3(1, 0.5, 0.25));
            const e = getContainerEllipsoid3(points);
            expect(e).not.toBeNull();
            for (const p of points) {
                expect(quadratic(p, e!)).toBeLessThanOrEqual(1 + 1e-9);
            }
            expect(Math.max(...points.map(p => quadratic(p, e!))))
                .toBeCloseTo(1, 9);
        });

        it('produces a sphere for a spherically symmetric set', () => {
            const points: Vector[] = [];
            const r = 3;
            for (const s of [-1, 1]) {
                points.push(v3(s * r, 0, 0), v3(0, s * r, 0), v3(0, 0, s * r));
            }
            const e = getContainerEllipsoid3(points);
            expect(e).not.toBeNull();
            for (let j = 0; j < 3; ++j) {
                expect(e!.center.values[j]).toBeCloseTo(0, 10);
                expect(e!.extent.values[j]).toBeCloseTo(r, 8);
            }
        });

        it('yields an infinite extent for coplanar points (ill-conditioned)',
            () => {
                const points = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0),
                    v3(1, 1, 0), v3(0.3, 0.7, 0)];
                const e = getContainerEllipsoid3(points);
                expect(e).not.toBeNull();
                expect([...e!.extent.values].some(x => !Number.isFinite(x)))
                    .toBe(true);
            });

        it('throws on an empty point set', () => {
            expect(() => getContainerEllipsoid3([])).toThrow();
        });

        it('rejects non-3D points', () => {
            expect(() => getContainerEllipsoid3([Vector.fromArray([0, 0])]))
                .toThrow();
        });

        it('contains every input point for random clouds', () => {
            for (let trial = 0; trial < 30; ++trial) {
                const axis = frame(v3(2 * rand() - 1, 2 * rand() - 1,
                    2 * rand() - 1 + 0.1), 2 * Math.PI * rand());
                const c = v3(8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4);
                const ext = [0.5 + 3 * rand(), 0.5 + 3 * rand(),
                    0.5 + 3 * rand()];
                const points: Vector[] = [];
                for (let i = 0; i < 60; ++i) {
                    let p = c.clone();
                    for (let j = 0; j < 3; ++j) {
                        p = add(p, mul(ext[j] * (2 * rand() - 1), axis[j]));
                    }
                    points.push(p);
                }
                const e = getContainerEllipsoid3(points);
                expect(e).not.toBeNull();
                for (const p of points) {
                    expect(quadratic(p, e!)).toBeLessThanOrEqual(1 + 1e-9);
                }
                for (let j = 0; j < 3; ++j) {
                    expect(length(e!.axis[j])).toBeCloseTo(1, 10);
                }
                expect(dot(e!.axis[0], e!.axis[1])).toBeCloseTo(0, 10);
                expect(dot(e!.axis[0], e!.axis[2])).toBeCloseTo(0, 10);
                expect(dot(e!.axis[1], e!.axis[2])).toBeCloseTo(0, 10);
            }
        });
    });

    describe('mergeContainersEllipsoid3', () => {
        it('reproduces the input when both ellipsoids are the same', () => {
            const e = ellipsoid3(v3(1, 2, -1), frame(v3(0, 0, 1), 0.4),
                v3(3, 2, 1));
            const m = mergeContainersEllipsoid3(e, e);
            for (let j = 0; j < 3; ++j) {
                expect(m.center.values[j]).toBeCloseTo(e.center.values[j], 8);
            }
            const sorted = [...m.extent.values].sort((a, b) => a - b);
            expect(sorted[0]).toBeCloseTo(1, 8);
            expect(sorted[1]).toBeCloseTo(2, 8);
            expect(sorted[2]).toBeCloseTo(3, 8);
        });

        it('spans the projection intervals of both inputs (randomized)', () => {
            // The upstream algorithm builds the bounding box of the two
            // projected intervals along the averaged axes and uses its
            // half-widths as the merged extents.
            for (let trial = 0; trial < 40; ++trial) {
                const e0 = ellipsoid3(
                    v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3),
                    frame(v3(2 * rand() - 1, 2 * rand() - 1,
                        2 * rand() - 1 + 0.1), 2 * Math.PI * rand()),
                    v3(0.5 + 2 * rand(), 0.5 + 2 * rand(), 0.5 + 2 * rand()));
                const e1 = ellipsoid3(
                    v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3),
                    frame(v3(2 * rand() - 1, 2 * rand() - 1,
                        2 * rand() - 1 + 0.1), 2 * Math.PI * rand()),
                    v3(0.5 + 2 * rand(), 0.5 + 2 * rand(), 0.5 + 2 * rand()));
                const m = mergeContainersEllipsoid3(e0, e1);

                for (let j = 0; j < 3; ++j) {
                    expect(length(m.axis[j])).toBeCloseTo(1, 9);
                }
                expect(dot(m.axis[0], m.axis[1])).toBeCloseTo(0, 9);
                expect(dot(m.axis[0], m.axis[2])).toBeCloseTo(0, 9);
                expect(dot(m.axis[1], m.axis[2])).toBeCloseTo(0, 9);

                for (let j = 0; j < 3; ++j) {
                    const line = Line.fromOriginDirection(m.center, m.axis[j]);
                    const p0 = projectEllipsoid3(e0, line);
                    const p1 = projectEllipsoid3(e1, line);
                    const lo = Math.min(p0.smin, p1.smin);
                    const hi = Math.max(p0.smax, p1.smax);
                    expect(0.5 * (lo + hi)).toBeCloseTo(0, 8);
                    expect(m.extent.values[j]).toBeCloseTo(0.5 * (hi - lo), 8);
                }
            }
        });

        it('rejects mismatched dimensions', () => {
            expect(() => mergeContainersEllipsoid3(new Hyperellipsoid(3),
                new Hyperellipsoid(2))).toThrow();
        });
    });
});
