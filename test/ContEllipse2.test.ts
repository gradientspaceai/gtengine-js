import { describe, it, expect } from 'vitest';
import {
    getContainerEllipse2,
    inContainerEllipse2,
    mergeContainersEllipse2
} from '../src/ContEllipse2.js';
import { Hyperellipsoid, type Ellipse2 } from '../src/Hyperellipsoid.js';
import { Line } from '../src/Line.js';
import { projectEllipse2 } from '../src/Projection.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';

function v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function ellipse2(cx: number, cy: number, angle: number, e0: number,
    e1: number): Ellipse2 {
    const c = Math.cos(angle), s = Math.sin(angle);
    return Hyperellipsoid.fromCenterAxisExtent(v(cx, cy),
        [v(c, s), v(-s, c)], v(e0, e1));
}

// Independent evaluation of the ellipse quadratic form: <= 1 inside.
function quadratic(p: Vector, e: Ellipse2): number {
    const d = sub(p, e.center);
    const t0 = dot(d, e.axis[0]) / e.extent.values[0];
    const t1 = dot(d, e.axis[1]) / e.extent.values[1];
    return t0 * t0 + t1 * t1;
}

let seed = 555000111;
function rand(): number {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
}

describe('ContEllipse2', () => {
    describe('inContainerEllipse2', () => {
        it('matches the axis-aligned unit circle', () => {
            const e = ellipse2(0, 0, 0, 1, 1);
            expect(inContainerEllipse2(v(0, 0), e)).toBe(true);
            expect(inContainerEllipse2(v(1, 0), e)).toBe(true);
            expect(inContainerEllipse2(v(0, -1), e)).toBe(true);
            expect(inContainerEllipse2(v(1.0001, 0), e)).toBe(false);
            expect(inContainerEllipse2(v(0.71, 0.71), e)).toBe(false);
            expect(inContainerEllipse2(v(0.7, 0.7), e)).toBe(true);
        });

        it('matches a rotated, translated, elongated ellipse', () => {
            const e = ellipse2(3, -2, Math.PI / 3, 4, 1);
            // Endpoints of the axes are on the boundary.
            for (let j = 0; j < 2; ++j) {
                const tip = add(e.center, mul(e.extent.values[j], e.axis[j]));
                expect(inContainerEllipse2(tip, e)).toBe(true);
                const beyond = add(e.center,
                    mul(e.extent.values[j] * 1.001, e.axis[j]));
                expect(inContainerEllipse2(beyond, e)).toBe(false);
            }
        });

        it('rejects mismatched dimensions', () => {
            const e = ellipse2(0, 0, 0, 1, 1);
            expect(() => inContainerEllipse2(Vector.fromArray([0, 0, 0]), e))
                .toThrow();
        });
    });

    describe('getContainerEllipse2', () => {
        it('contains every input point', () => {
            const points = [v(0, 0), v(4, 0), v(4, 1), v(0, 1), v(2, 0.5),
                v(1, 0.9), v(3, 0.2)];
            const e = getContainerEllipse2(points);
            expect(e).not.toBeNull();
            for (const p of points) {
                expect(quadratic(p, e!)).toBeLessThanOrEqual(1 + 1e-9);
            }
            // At least one point is on the boundary (the growth is tight).
            const maxQ = Math.max(...points.map(p => quadratic(p, e!)));
            expect(maxQ).toBeCloseTo(1, 9);
        });

        it('produces a circle for a symmetric point set', () => {
            const points: Vector[] = [];
            for (let i = 0; i < 16; ++i) {
                const a = 2 * Math.PI * i / 16;
                points.push(v(5 + 2 * Math.cos(a), -1 + 2 * Math.sin(a)));
            }
            const e = getContainerEllipse2(points);
            expect(e).not.toBeNull();
            expect(e!.center.values[0]).toBeCloseTo(5, 10);
            expect(e!.center.values[1]).toBeCloseTo(-1, 10);
            expect(e!.extent.values[0]).toBeCloseTo(2, 8);
            expect(e!.extent.values[1]).toBeCloseTo(2, 8);
        });

        it('yields infinite extents for collinear points (ill-conditioned)',
            () => {
                const points = [v(0, 0), v(1, 1), v(2, 2), v(3, 3)];
                const e = getContainerEllipse2(points);
                expect(e).not.toBeNull();
                // The degenerate direction has a zero eigenvalue, so the
                // extent is infinite. Upstream has the same behavior; the
                // header documents the construction as ill-conditioned here.
                const extents = [e!.extent.values[0], e!.extent.values[1]];
                expect(extents.some(x => !Number.isFinite(x))).toBe(true);
            });

        it('throws on an empty point set', () => {
            expect(() => getContainerEllipse2([])).toThrow();
        });

        it('rejects non-2D points', () => {
            expect(() => getContainerEllipse2([Vector.fromArray([0, 0, 0])]))
                .toThrow();
        });

        it('contains every input point for random clouds', () => {
            for (let trial = 0; trial < 40; ++trial) {
                const ca = 2 * Math.PI * rand();
                const cs = Math.cos(ca), sn = Math.sin(ca);
                const a = 1 + 4 * rand();
                const b = 1 + 4 * rand();
                const points: Vector[] = [];
                for (let i = 0; i < 50; ++i) {
                    const x = a * (2 * rand() - 1);
                    const y = b * (2 * rand() - 1);
                    points.push(v(cs * x - sn * y + 7, sn * x + cs * y - 3));
                }
                const e = getContainerEllipse2(points);
                expect(e).not.toBeNull();
                for (const p of points) {
                    expect(quadratic(p, e!)).toBeLessThanOrEqual(1 + 1e-9);
                }
                // The axes are orthonormal.
                expect(length(e!.axis[0])).toBeCloseTo(1, 10);
                expect(length(e!.axis[1])).toBeCloseTo(1, 10);
                expect(dot(e!.axis[0], e!.axis[1])).toBeCloseTo(0, 10);
            }
        });
    });

    describe('mergeContainersEllipse2', () => {
        it('reproduces the input when both ellipses are the same', () => {
            const e = ellipse2(1, 2, 0.4, 3, 1);
            const m = mergeContainersEllipse2(e, e);
            expect(m.center.values[0]).toBeCloseTo(1, 9);
            expect(m.center.values[1]).toBeCloseTo(2, 9);
            expect(m.extent.values[0]).toBeCloseTo(3, 9);
            expect(m.extent.values[1]).toBeCloseTo(1, 9);
            expect(Math.abs(dot(m.axis[0], e.axis[0]))).toBeCloseTo(1, 9);
        });

        it('spans the projection intervals of both inputs', () => {
            // The upstream algorithm builds the bounding box of the two
            // projected intervals along the averaged axes and uses its
            // half-widths as the merged extents.
            for (let trial = 0; trial < 40; ++trial) {
                const e0 = ellipse2(10 * rand() - 5, 10 * rand() - 5,
                    Math.PI * rand(), 0.5 + 3 * rand(), 0.5 + 3 * rand());
                const e1 = ellipse2(10 * rand() - 5, 10 * rand() - 5,
                    Math.PI * rand(), 0.5 + 3 * rand(), 0.5 + 3 * rand());
                const m = mergeContainersEllipse2(e0, e1);

                expect(length(m.axis[0])).toBeCloseTo(1, 10);
                expect(length(m.axis[1])).toBeCloseTo(1, 10);
                expect(dot(m.axis[0], m.axis[1])).toBeCloseTo(0, 10);

                for (let j = 0; j < 2; ++j) {
                    const line = Line.fromOriginDirection(m.center, m.axis[j]);
                    const p0 = projectEllipse2(e0, line);
                    const p1 = projectEllipse2(e1, line);
                    const lo = Math.min(p0.smin, p1.smin);
                    const hi = Math.max(p0.smax, p1.smax);
                    // The merged center is the midpoint of the union interval
                    // and the extent is its half-width.
                    expect(0.5 * (lo + hi)).toBeCloseTo(0, 9);
                    expect(m.extent.values[j]).toBeCloseTo(0.5 * (hi - lo), 9);
                }
            }
        });

        it('rejects mismatched dimensions', () => {
            const e2 = ellipse2(0, 0, 0, 1, 1);
            const e3 = new Hyperellipsoid(3);
            expect(() => mergeContainersEllipse2(e2, e3)).toThrow();
        });
    });
});
