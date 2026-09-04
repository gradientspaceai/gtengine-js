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
import {
    check, expectClose, fc, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContEllipse2.h semantics.
// ---------------------------------------------------------------------------

describe('ContEllipse2 verification', () => {
    // Anisotropic lattice cloud so the fitted frame is unambiguous.
    const baseGrid: Vector[] = [];
    for (let i = -3; i <= 3; ++i) {
        for (let j = -1; j <= 1; ++j) {
            baseGrid.push(v(i, j));
        }
    }

    // Q(X) = sum_j D[j]*Dot(U[j],X-C)^2, the quadratic form the construction
    // is based on. inContainer is |standardized| <= 1, that is Q(X) <= V^2.
    const quadratic = (ellipse: Ellipse2, p: Vector): number => {
        const diff = sub(p, ellipse.center);
        let sum = 0;
        for (let j = 0; j < 2; ++j) {
            const s = dot(diff, ellipse.axis[j]) / ellipse.extent.values[j];
            sum += s * s;
        }
        return sum;
    };

    // The design claim: the fitted ellipse contains every input point, and it
    // is tight, i.e. at least one point is on the boundary (the point that
    // attained maxValue).
    it('contains every input point and is tight on one of them', () => {
        check(fc.array(wellScaledVector(2, -6, 6), { minLength: 2, maxLength: 12 }),
            (points: Vector[]) => {
                const ellipse = getContainerEllipse2(points);
                if (ellipse === null
                    || !ellipse.extent.values.every(e => Number.isFinite(e) && e > 0)) {
                    // Collinear or coincident input: upstream's zero
                    // eigenvalue gives an infinite extent (issue #292).
                    return;
                }
                // The construction makes max_i Q(X_i) exactly 1 in exact
                // arithmetic, so the extremal point sits on the boundary and
                // the strict '<= 1' of InContainer can reject it by one ulp.
                // That is upstream behavior; the property is that Q is at
                // most 1 and that every point strictly inside is accepted.
                let maxQ = 0;
                for (const p of points) {
                    const q = quadratic(ellipse, p);
                    maxQ = Math.max(maxQ, q);
                    expect(q).toBeLessThanOrEqual(1 + 1e-9);
                    if (q < 1 - 1e-9) {
                        expect(inContainerEllipse2(p, ellipse)).toBe(true);
                    }
                }
                expectClose(maxQ, 1, 1e-9, 1e-9);
            });
    });

    // Rigid motions: the extents are invariant and the center follows.
    it('is equivariant under rigid motions', () => {
        const reference = getContainerEllipse2(baseGrid)!;
        check(fc.tuple(rotationFrame(2), wellScaledVector(2)),
            ([frame, t]: [Vector[], Vector]) => {
                const xform = (p: Vector): Vector =>
                    add(add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])), t);
                const moved = getContainerEllipse2(baseGrid.map(xform))!;
                const want = [reference.extent.get(0), reference.extent.get(1)]
                    .sort((a, b) => a - b);
                const got = [moved.extent.get(0), moved.extent.get(1)]
                    .sort((a, b) => a - b);
                expectClose(got[0], want[0], 1e-9, 1e-9);
                expectClose(got[1], want[1], 1e-9, 1e-9);
                expect(length(sub(moved.center, xform(reference.center))))
                    .toBeLessThanOrEqual(1e-8);
            });
    });

    // Upstream issue #292 item 1: the comment says nonpositive eigenvalues are
    // adjusted, but only strictly negative ones are. A zero eigenvalue
    // (collinear input) divides by zero, the extent becomes infinite and
    // InContainer then accepts every point. Preserved; pinned here.
    it('collinear input yields an infinite extent that accepts everything (#292)',
        () => {
            const points = [v(-2, 0), v(-1, 0), v(0, 0), v(1, 0), v(2, 0)];
            const ellipse = getContainerEllipse2(points)!;
            const extents = [ellipse.extent.get(0), ellipse.extent.get(1)];
            expect(extents.some(e => !Number.isFinite(e))).toBe(true);
            expect(inContainerEllipse2(v(0, 1e6), ellipse)).toBe(true);
        });

    // Upstream issue #292 item 2: the merged ellipse is inscribed in the box
    // of the projected intervals, so it does not contain its inputs. The
    // counterexample from the issue, pinned so the deviation stays visible.
    it('merge does not contain its inputs (#292)', () => {
        const unit = (cx: number): Ellipse2 => Hyperellipsoid.fromCenterAxisExtent(
            v(cx, 0), [v(1, 0), v(0, 1)], v(1, 1));
        const merge = mergeContainersEllipse2(unit(-1), unit(1));
        expectClose(merge.extent.get(0), 2, 1e-12, 1e-12);
        expectClose(merge.extent.get(1), 1, 1e-12, 1e-12);
        // (-1,1) is on the first input circle but outside the merge.
        expect(inContainerEllipse2(v(-1, 0), unit(-1))).toBe(true);
        expect(inContainerEllipse2(v(-1, 1), unit(-1))).toBe(true);
        expect(inContainerEllipse2(v(-1, 1), merge)).toBe(false);
    });

    // What the merge does guarantee: its extents span the projected intervals
    // of both inputs along the merged axes.
    it('merge spans the projection intervals of both inputs', () => {
        const ellipseArb = fc.tuple(wellScaledVector(2, -5, 5), rotationFrame(2),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0.2, max: 3, noNaN: true }))
            .map(([c, frame, e0, e1]) => Hyperellipsoid.fromCenterAxisExtent(
                c, frame, Vector.fromArray([e0, e1])));
        check(fc.tuple(ellipseArb, ellipseArb),
            ([e0, e1]: [Ellipse2, Ellipse2]) => {
                const merge = mergeContainersEllipse2(e0, e1);
                expectClose(dot(merge.axis[0], merge.axis[1]), 0, 1e-12, 1e-12);
                for (let j = 0; j < 2; ++j) {
                    const line = Line.fromOriginDirection(merge.center,
                        merge.axis[j]);
                    const p0 = projectEllipse2(e0, line);
                    const p1 = projectEllipse2(e1, line);
                    const lo = Math.min(p0.smin, p1.smin);
                    const hi = Math.max(p0.smax, p1.smax);
                    // The merged center is the midpoint of that interval and
                    // the extent is its half-width.
                    expectClose(lo + hi, 0, 1e-8, 1e-8);
                    expectClose(merge.extent.get(j), 0.5 * (hi - lo),
                        1e-8, 1e-8);
                }
            });
    });
});
