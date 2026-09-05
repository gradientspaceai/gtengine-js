import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { IntrLine3Cylinder3FI } from '../src/IntrLine3Cylinder3.js';
import {
    intrLine3Cylinder3FIDoQuery,
    defaultIntrLine3Cylinder3FIResult
} from '../src/IntrLine3Cylinder3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function line(p: number[], d: number[]): Line {
    const dir = vec(d);
    normalize(dir);
    return Line.fromOriginDirection(vec(p), dir);
}

function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const dir = vec(direction);
    normalize(dir);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(vec(origin), dir), radius, height);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function insideCylinder(c: Cylinder3, x: Vector): boolean {
    const diff = sub(x, c.axis.origin);
    const z = dot(diff, c.axis.direction);
    if (Math.abs(z) > 0.5 * c.height) {
        return false;
    }
    const radial = sub(diff, mul(z, c.axis.direction));
    return dot(radial, radial) <= c.radius * c.radius;
}

describe('IntrLine3Cylinder3', () => {
    const fi = new IntrLine3Cylinder3FI();

    // Cylinder about the z-axis, radius 1, height 4 (z in [-2,2]).
    const cyl = cylinder([0, 0, 0], [0, 0, 1], 1, 4);

    it('finds the chord of a line perpendicular to the axis', () => {
        const l = line([0, 0, 0], [1, 0, 0]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a line parallel to the cylinder axis', () => {
        const l = line([0.5, 0, -10], [0, 0, 1]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(8, 12);
        expect(result.parameter[1]).toBeCloseTo(12, 12);
        expect(result.point[0].values[2]).toBeCloseTo(-2, 12);
        expect(result.point[1].values[2]).toBeCloseTo(2, 12);
    });

    it('reports no intersection for a parallel line outside the radius', () => {
        const l = line([2, 0, -10], [0, 0, 1]);
        expect(fi.find(l, cyl).intersect).toBe(false);
    });

    it('reports no intersection for a perpendicular line beyond the end disk', () => {
        const l = line([0, 0, 3], [1, 0, 0]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('finds the chord of an oblique line that exits through the wall', () => {
        const l = line([0, 0, 0], [1, 0, 1]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-Math.SQRT2, 10);
        expect(result.parameter[1]).toBeCloseTo(Math.SQRT2, 10);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 10);
        expect(result.point[1].values[0]).toBeCloseTo(1, 10);
    });

    it('finds an oblique line entering and leaving through the end disks', () => {
        const l = line([0, 0, -10], [0.01, 0, 1]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[2]).toBeCloseTo(-2, 10);
        expect(result.point[1].values[2]).toBeCloseTo(2, 10);
    });

    it('reports a single point for a line tangent to the wall', () => {
        const l = line([1, -10, 0], [0, 1, 0]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(10, 10);
        expect(result.parameter[1]).toBeCloseTo(result.parameter[0], 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 10);
    });

    it('handles a zero-height cylinder as a disk', () => {
        const flat = cylinder([0, 0, 0], [0, 0, 1], 1, 0);
        const l = line([0.5, 0, -3], [0, 0, 1]);
        const result = fi.find(l, flat);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
    });

    it('does not support an infinite cylinder', () => {
        // Cylinder3 stores height = -1 for the infinite state and this
        // upstream query has no infinite branch, so it computes with
        // halfHeight = -0.5 and finds nothing. The behavior is preserved.
        const infinite = cylinder([0, 0, 0], [0, 0, 1], 1, 1);
        infinite.makeInfiniteCylinder();
        expect(infinite.isInfinite()).toBe(true);
        const l = line([0, 0, 0], [1, 0, 0]);
        expect(fi.find(l, infinite).intersect).toBe(false);
    });

    it('agrees with dense sampling along the line', () => {
        const rand = makeRandom(112358);
        const cyl2 = cylinder([0.5, -1, 0.25], [1, 2, -1], 0.9, 3);
        for (let trial = 0; trial < 120; ++trial) {
            const l = line(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(l, cyl2);

            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = -10 + (20 * k) / n;
                const x = add(l.origin, mul(t, l.direction));
                if (insideCylinder(cyl2, x)) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(3e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(3e-3);
            }

            if (result.intersect) {
                // Every reported point lies in the closed cylinder.
                for (let i = 0; i < 2; ++i) {
                    const p = add(l.origin, mul(result.parameter[i], l.direction));
                    const diff = sub(p, cyl2.axis.origin);
                    const z = dot(diff, cyl2.axis.direction);
                    const radial = sub(diff, mul(z, cyl2.axis.direction));
                    expect(Math.abs(z)).toBeLessThanOrEqual(1.5 + 1e-9);
                    expect(Math.sqrt(dot(radial, radial)))
                        .toBeLessThanOrEqual(0.9 + 1e-9);
                }
            }
        }
    });
});

describe('intrLine3Cylinder3FIDoQuery', () => {
    const c = cylinder([0, 0, 0], [0, 0, 1], 1, 2);

    it('matches the class query but does not compute points', () => {
        const l = line([-5, 0, 0], [1, 0, 0]);
        const result = defaultIntrLine3Cylinder3FIResult();
        intrLine3Cylinder3FIDoQuery(l.origin, l.direction, c, result);
        const expected = new IntrLine3Cylinder3FI().find(l, c);
        expect(result.intersect).toBe(expected.intersect);
        expect(result.numIntersections).toBe(expected.numIntersections);
        expect(result.parameter[0]).toBeCloseTo(expected.parameter[0], 12);
        expect(result.parameter[1]).toBeCloseTo(expected.parameter[1], 12);
        // DoQuery leaves 'point' at its default value.
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([0, 0, 0]);
    });

    it('reports no intersection for a line missing the cylinder', () => {
        const result = defaultIntrLine3Cylinder3FIResult();
        intrLine3Cylinder3FIDoQuery(vec([-5, 3, 0]), vec([1, 0, 0]), c,
            result);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, expectClose,
    expectVectorClose, seededRandom
} from './helpers/arbitraries.js';
import { length as vlength } from '../src/Vector.js';

const lineCylinder = fc.tuple(wellScaledVector(3, -6, 6), unitVector(3),
    wellScaledVector(3, -4, 4), unitVector(3),
    fc.double({ min: 0.2, max: 3, noNaN: true }),
    fc.double({ min: 0.5, max: 5, noNaN: true }))
    .map(([o, d, c, w, radius, height]) => ({
        line: Line.fromOriginDirection(o, d),
        cylinder: Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(c, w), radius, height)
    }));

// The cylinder coordinates (radial distance, axial offset) of a point.
function cylCoordsV31(cyl: Cylinder3, p: Vector): { radial: number, z: number } {
    const diff = sub(p, cyl.axis.origin);
    const z = dot(diff, cyl.axis.direction);
    const radial = vlength(sub(diff, mul(z, cyl.axis.direction)));
    return { radial, z };
}

function insideCylinderV31(cyl: Cylinder3, p: Vector, tol: number): boolean {
    const { radial, z } = cylCoordsV31(cyl, p);
    return radial <= cyl.radius - tol
        && Math.abs(z) <= 0.5 * cyl.height - tol;
}

describe('IntrLine3Cylinder3 verification', () => {
    const fiq = new IntrLine3Cylinder3FI();

    it('the reported points are on the line and on the cylinder', () => {
        check(lineCylinder, ({ line: l, cylinder: c }) => {
            const f = fiq.find(l, c);
            if (!f.intersect) {
                return;
            }
            expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
            for (let i = 0; i < 2; ++i) {
                expectVectorClose(f.point[i],
                    add(l.origin, mul(f.parameter[i], l.direction)), 0, 0);
                const { radial, z } = cylCoordsV31(c, f.point[i]);
                const half = 0.5 * c.height;
                // The point is on the wall (radial == r, |z| <= h/2) or on an
                // end disk (|z| == h/2, radial <= r). The parameters carry a
                // square root or a division, so allow a small residual.
                const onWall = Math.abs(radial - c.radius) < 1e-7
                    && Math.abs(z) <= half + 1e-7;
                const onCap = Math.abs(Math.abs(z) - half) < 1e-7
                    && radial <= c.radius + 1e-7;
                expect(onWall || onCap).toBe(true);
                expect(Number.isNaN(f.parameter[i])).toBe(false);
            }
        });
    });

    it('a fine sweep of the line agrees with the reported interval', () => {
        const rnd = seededRandom(0x51d0b7e4);
        for (let trial = 0; trial < 150; ++trial) {
            const w = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (vlength(w) < 0.3) {
                continue;
            }
            normalize(w);
            const c = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(Vector.fromArray([rnd() * 4 - 2,
                    rnd() * 4 - 2, rnd() * 4 - 2]), w),
                0.3 + rnd() * 1.5, 0.5 + rnd() * 3);
            const d = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (vlength(d) < 0.3) {
                continue;
            }
            normalize(d);
            const l = Line.fromOriginDirection(
                Vector.fromArray([rnd() * 8 - 4, rnd() * 8 - 4,
                    rnd() * 8 - 4]), d);
            const f = fiq.find(l, c);
            for (let k = 0; k <= 600; ++k) {
                const t = -12 + (24 * k) / 600;
                if (insideCylinderV31(c, add(l.origin, mul(t, d)), 1e-5)) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-6);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-6);
                }
            }
        }
    }, 30000);

    it('a line along the axis is clipped by the end disks', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0.5, max: 5, noNaN: true }), fc.boolean()),
            ([o, w, radius, height, flip]) => {
                const c = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(o, w), radius, height);
                const dir = flip ? mul(-1, w) : w;
                const l = Line.fromOriginDirection(o, dir);
                const f = fiq.find(l, c);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(2);
                expectClose(f.parameter[0], -0.5 * height, 1e-9, 1e-9);
                expectClose(f.parameter[1], 0.5 * height, 1e-9, 1e-9);
            });
    });

    it('a line parallel to the axis but outside the wall misses', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0.5, max: 5, noNaN: true }),
            fc.double({ min: 1.05, max: 3, noNaN: true })),
            ([o, R, radius, height, scale]) => {
                const c = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(o, R[0]), radius, height);
                const l = Line.fromOriginDirection(
                    add(o, mul(scale * radius, R[1])), R[0]);
                expect(fiq.find(l, c).intersect).toBe(false);
            });
    });

    it('a line perpendicular to the axis through the center hits at +-r', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0.5, max: 5, noNaN: true }),
            fc.double({ min: -0.45, max: 0.45, noNaN: true })),
            ([o, R, radius, height, zfrac]) => {
                const c = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(o, R[0]), radius, height);
                const origin = add(o, mul(zfrac * height, R[0]));
                const l = Line.fromOriginDirection(origin, R[1]);
                const f = fiq.find(l, c);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(2);
                expectClose(f.parameter[0], -radius, 1e-8, 1e-9);
                expectClose(f.parameter[1], radius, 1e-8, 1e-9);
            });
    });

    it('a line through a sampled interior point always hits', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0.5, max: 5, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }),
            fc.double({ min: 0, max: 0.9, noNaN: true }),
            fc.double({ min: 0, max: 6.28, noNaN: true }), unitVector(3)),
            ([o, R, radius, height, zf, rf, theta, d]) => {
                const c = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(o, R[0]), radius, height);
                const target = add(o, add(mul(zf * 0.5 * height, R[0]),
                    add(mul(rf * radius * Math.cos(theta), R[1]),
                        mul(rf * radius * Math.sin(theta), R[2]))));
                const l = Line.fromOriginDirection(target, d);
                const f = fiq.find(l, c);
                expect(f.intersect).toBe(true);
                expect(f.parameter[0]).toBeLessThanOrEqual(1e-9);
                expect(f.parameter[1]).toBeGreaterThanOrEqual(-1e-9);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineCylinder, rotationFrame(3),
            wellScaledVector(3, -4, 4)),
            ([{ line: l, cylinder: c }, R, tr]) => {
                const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                    add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
                const xf = (p: Vector): Vector => add(tr, rot(p));
                const l2 = Line.fromOriginDirection(xf(l.origin),
                    rot(l.direction));
                const c2 = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(xf(c.axis.origin),
                        rot(c.axis.direction)), c.radius, c.height);
                const f1 = fiq.find(l, c);
                const f2 = fiq.find(l2, c2);
                const graze = (f: typeof f1): boolean => f.intersect
                    && f.parameter[1] - f.parameter[0] < 1e-4;
                if (graze(f1) || graze(f2)) {
                    return;
                }
                expect(f1.intersect).toBe(f2.intersect);
                if (!f1.intersect) {
                    return;
                }
                expectClose(f1.parameter[0], f2.parameter[0], 1e-6, 1e-7);
                expectClose(f1.parameter[1], f2.parameter[1], 1e-6, 1e-7);
            });
    });

    it('the exported DoQuery reproduces the class result', () => {
        check(lineCylinder, ({ line: l, cylinder: c }) => {
            const res = defaultIntrLine3Cylinder3FIResult();
            intrLine3Cylinder3FIDoQuery(l.origin, l.direction, c, res);
            const f = fiq.find(l, c);
            expect(res.intersect).toBe(f.intersect);
            expect(res.numIntersections).toBe(f.numIntersections);
            expect(res.parameter[0]).toBe(f.parameter[0]);
            expect(res.parameter[1]).toBe(f.parameter[1]);
            expect(res.point[0].values).toEqual([0, 0, 0]);
        });
    });

    it('an infinite cylinder (height = -1) reports no intersection', () => {
        // Upstream has no infinite-cylinder branch; it reads height directly,
        // so the sentinel gives a negative half-height and every query misses.
        // Pinned so the behaviour is a deliberate, documented limitation.
        const c = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(Vector.fromArray([0, 0, 0]),
                Vector.fromArray([0, 0, 1])), 1, 2);
        c.makeInfiniteCylinder();
        expect(c.isInfinite()).toBe(true);
        const l = Line.fromOriginDirection(Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]));
        expect(fiq.find(l, c).intersect).toBe(false);
    });
});
