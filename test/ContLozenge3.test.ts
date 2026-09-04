import { describe, it, expect } from 'vitest';
import {
    getContainerLozenge3,
    inContainerLozenge3
} from '../src/ContLozenge3.js';
import { Lozenge3 } from '../src/Lozenge3.js';
import { Rectangle } from '../src/Rectangle.js';
import { DistPointRectangle } from '../src/DistPointRectangle.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, fc, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function distanceToRectangle(point: Vector, rectangle: Rectangle): number {
    return new DistPointRectangle().compute(point, rectangle).distance;
}

function makeRandom(seed: number): () => number {
    let s = seed;
    return () => {
        s = (s * 1103515245 + 12345) % 2147483648;
        return s / 2147483648;
    };
}

describe('getContainerLozenge3', () => {
    it('contains every input point (random slab-like clouds)', () => {
        const rand = makeRandom(31337);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 50; ++i) {
                points.push(v(
                    6 * (rand() - 0.5),
                    4 * (rand() - 0.5),
                    1 * (rand() - 0.5)));
            }
            const lozenge = getContainerLozenge3(points);
            expect(lozenge.radius).toBeGreaterThan(0);
            for (const p of points) {
                expect(distanceToRectangle(p, lozenge.rectangle))
                    .toBeLessThanOrEqual(lozenge.radius + 1e-9);
            }
        }
    });

    it('contains every input point of a random isotropic cloud', () => {
        const rand = makeRandom(90210);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 40; ++i) {
                points.push(v(
                    4 * (rand() - 0.5),
                    4 * (rand() - 0.5),
                    4 * (rand() - 0.5)));
            }
            const lozenge = getContainerLozenge3(points);
            for (const p of points) {
                // Extreme points land on the lozenge boundary, where the
                // strict upstream test can fail by an ulp; allow rounding.
                expect(distanceToRectangle(p, lozenge.rectangle))
                    .toBeLessThanOrEqual(lozenge.radius + 1e-12);
            }
            // Points pulled toward the center are strictly inside.
            for (const p of points) {
                const shrunk = v(
                    lozenge.rectangle.center.values[0]
                        + 0.5 * (p.values[0] - lozenge.rectangle.center.values[0]),
                    lozenge.rectangle.center.values[1]
                        + 0.5 * (p.values[1] - lozenge.rectangle.center.values[1]),
                    lozenge.rectangle.center.values[2]
                        + 0.5 * (p.values[2] - lozenge.rectangle.center.values[2]));
                expect(inContainerLozenge3(shrunk, lozenge)).toBe(true);
            }
        }
    });

    it('fits a planar point set with a near-zero radius', () => {
        // All points in the z = 0 plane, so the thin direction has zero
        // spread and the lozenge radius is 0.
        const points: Vector[] = [];
        for (let i = -3; i <= 3; ++i) {
            for (let j = -2; j <= 2; ++j) {
                points.push(v(i, j, 0));
            }
        }
        const lozenge = getContainerLozenge3(points);
        expect(lozenge.radius).toBeCloseTo(0, 10);
        // The rectangle plane is z = 0.
        expect(lozenge.rectangle.center.values[2]).toBeCloseTo(0, 10);
        expect(Math.abs(lozenge.rectangle.axis[0].values[2])).toBeCloseTo(0, 10);
        expect(Math.abs(lozenge.rectangle.axis[1].values[2])).toBeCloseTo(0, 10);
        // The rectangle covers the 6 x 4 extent of the data.
        const extents = [
            lozenge.rectangle.extent.values[0],
            lozenge.rectangle.extent.values[1]
        ].sort((a, b) => a - b);
        expect(extents[0]).toBeCloseTo(2, 8);
        expect(extents[1]).toBeCloseTo(3, 8);
        for (const p of points) {
            expect(distanceToRectangle(p, lozenge.rectangle))
                .toBeLessThanOrEqual(lozenge.radius + 1e-9);
        }
    });

    it('degenerates to a sphere-like lozenge for coincident points', () => {
        const p = v(-2, 3, 1);
        const lozenge = getContainerLozenge3([p, p.clone(), p.clone(), p.clone()]);
        expect(lozenge.radius).toBeCloseTo(0, 12);
        expect(lozenge.rectangle.extent.values[0]).toBeCloseTo(0, 12);
        expect(lozenge.rectangle.extent.values[1]).toBeCloseTo(0, 12);
        expect(lozenge.rectangle.center.values[0]).toBeCloseTo(-2, 12);
        expect(lozenge.rectangle.center.values[1]).toBeCloseTo(3, 12);
        expect(lozenge.rectangle.center.values[2]).toBeCloseTo(1, 12);
    });

    it('handles collinear points (a capsule-shaped lozenge)', () => {
        const points = [v(0, 0, 0), v(1, 0, 0), v(2, 0, 0), v(3, 0, 0)];
        const lozenge = getContainerLozenge3(points);
        expect(lozenge.radius).toBeCloseTo(0, 10);
        const extents = [
            lozenge.rectangle.extent.values[0],
            lozenge.rectangle.extent.values[1]
        ].sort((a, b) => a - b);
        expect(extents[0]).toBeCloseTo(0, 10);
        expect(extents[1]).toBeCloseTo(1.5, 10);
        for (const p of points) {
            expect(distanceToRectangle(p, lozenge.rectangle))
                .toBeLessThanOrEqual(1e-9);
        }
    });

    it('keeps the rectangle axes orthonormal', () => {
        const rand = makeRandom(1717);
        const points: Vector[] = [];
        for (let i = 0; i < 30; ++i) {
            points.push(v(5 * (rand() - 0.5), 3 * (rand() - 0.5),
                1.5 * (rand() - 0.5)));
        }
        const lozenge = getContainerLozenge3(points);
        const a0 = lozenge.rectangle.axis[0];
        const a1 = lozenge.rectangle.axis[1];
        expect(dot(a0, a0)).toBeCloseTo(1, 10);
        expect(dot(a1, a1)).toBeCloseTo(1, 10);
        expect(dot(a0, a1)).toBeCloseTo(0, 10);
    });

    it('is translation covariant', () => {
        const rand = makeRandom(606);
        const base: Vector[] = [];
        for (let i = 0; i < 25; ++i) {
            base.push(v(4 * (rand() - 0.5), 3 * (rand() - 0.5),
                (rand() - 0.5)));
        }
        const shifted = base.map(p => v(p.values[0] + 10, p.values[1] - 5,
            p.values[2] + 2));
        const l0 = getContainerLozenge3(base);
        const l1 = getContainerLozenge3(shifted);
        expect(l1.radius).toBeCloseTo(l0.radius, 8);
        expect(l1.rectangle.center.values[0])
            .toBeCloseTo(l0.rectangle.center.values[0] + 10, 8);
        expect(l1.rectangle.center.values[1])
            .toBeCloseTo(l0.rectangle.center.values[1] - 5, 8);
        expect(l1.rectangle.center.values[2])
            .toBeCloseTo(l0.rectangle.center.values[2] + 2, 8);
    });

    it('throws for an empty point set and for non-3D points', () => {
        expect(() => getContainerLozenge3([])).toThrow();
        expect(() => getContainerLozenge3([Vector.fromArray([1, 2])])).toThrow();
    });
});

describe('inContainerLozenge3', () => {
    const lozenge = Lozenge3.fromRectangleRadius(
        Rectangle.fromCenterAxisExtent(v(0, 0, 0),
            [v(1, 0, 0), v(0, 1, 0)], Vector.fromArray([2, 1])), 1);

    it('accepts points inside the slab over the rectangle', () => {
        expect(inContainerLozenge3(v(0, 0, 0), lozenge)).toBe(true);
        expect(inContainerLozenge3(v(2, 1, 0), lozenge)).toBe(true);
        expect(inContainerLozenge3(v(0, 0, 1), lozenge)).toBe(true);
        expect(inContainerLozenge3(v(0, 0, -1), lozenge)).toBe(true);
    });

    it('accepts points in the rounded side and corner regions', () => {
        expect(inContainerLozenge3(v(3, 0, 0), lozenge)).toBe(true);
        expect(inContainerLozenge3(v(0, 2, 0), lozenge)).toBe(true);
        expect(inContainerLozenge3(
            v(2 + Math.SQRT1_2, 1 + Math.SQRT1_2, 0), lozenge)).toBe(true);
    });

    it('rejects points beyond the radius', () => {
        expect(inContainerLozenge3(v(0, 0, 1.0001), lozenge)).toBe(false);
        expect(inContainerLozenge3(v(3.0001, 0, 0), lozenge)).toBe(false);
        expect(inContainerLozenge3(v(3, 2, 0), lozenge)).toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContLozenge3.h semantics, including a pin for the
// corner-centering fix (upstream issue #174).
// ---------------------------------------------------------------------------

describe('ContLozenge3 verification', () => {
    // Anisotropic slab clouds: 7 x 5 lattice in the plane, thin in z, then
    // rigidly moved. The three principal variances (about 4, 2 and 0.02) are
    // well separated, so the fitted frame is not ambiguous and equivariance
    // properties hold to a tight tolerance.
    const baseGrid: Vector[] = [];
    for (let i = -3; i <= 3; ++i) {
        for (let j = -2; j <= 2; ++j) {
            baseGrid.push(v(i, j, 0.05 * ((i * 7 + j * 3) % 5 - 2)));
        }
    }

    const rigid = (frame: Vector[], t: Vector) => (p: Vector): Vector =>
        add(add(add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])),
            mul(p.get(2), frame[2])), t);

    // The design claim of GetContainer, and the point of the #174 fix: every
    // input point is within 'radius' of the rectangle.
    it('the fitted lozenge contains every input point', () => {
        check(fc.array(wellScaledVector(3, -6, 6), { minLength: 1, maxLength: 14 }),
            (points: Vector[]) => {
                const lozenge = getContainerLozenge3(points);
                for (const p of points) {
                    expect(distanceToRectangle(p, lozenge.rectangle))
                        .toBeLessThanOrEqual(lozenge.radius + 1e-9);
                }
            });
    });

    // Same claim on the anisotropic slab clouds, under a random rigid motion.
    it('contains every input point of a rigidly moved slab cloud', () => {
        check(fc.tuple(rotationFrame(3), wellScaledVector(3)),
            ([frame, t]: [Vector[], Vector]) => {
                const points = baseGrid.map(rigid(frame, t));
                const lozenge = getContainerLozenge3(points);
                for (const p of points) {
                    expect(distanceToRectangle(p, lozenge.rectangle))
                        .toBeLessThanOrEqual(lozenge.radius + 1e-9);
                }
            });
    });

    // Regression pin for upstream issue #174: the rectangle is a *centered*
    // primitive, so the center must be the midpoint of the fitted parameter
    // interval, not its lower corner. Upstream returns center (-3,-2,0) here
    // (a corner of the data), which leaves the point (3,2,0) three units
    // outside the lozenge.
    it('centers the rectangle on the data, not on a corner (#174)', () => {
        const points: Vector[] = [];
        for (let i = -3; i <= 3; ++i) {
            for (let j = -2; j <= 2; ++j) {
                points.push(v(i, j, 0));
            }
        }
        const lozenge = getContainerLozenge3(points);
        expect(lozenge.radius).toBeLessThanOrEqual(1e-9);
        for (let d = 0; d < 3; ++d) {
            expectClose(lozenge.rectangle.center.get(d), 0, 1e-9, 1e-9);
        }
        // The corner-centered upstream result would put (3,2,0) at distance 3.
        for (const p of points) {
            expect(distanceToRectangle(p, lozenge.rectangle))
                .toBeLessThanOrEqual(lozenge.radius + 1e-9);
        }
    });

    // Rigid motions: the radius and the two extents are invariant and the
    // rectangle center follows the motion.
    it('is equivariant under rigid motions', () => {
        const reference = getContainerLozenge3(baseGrid);
        check(fc.tuple(rotationFrame(3), wellScaledVector(3)),
            ([frame, t]: [Vector[], Vector]) => {
                const xform = rigid(frame, t);
                const moved = getContainerLozenge3(baseGrid.map(xform));
                expectClose(moved.radius, reference.radius, 1e-9, 1e-9);
                const want = [reference.rectangle.extent.get(0),
                    reference.rectangle.extent.get(1)].sort((a, b) => a - b);
                const got = [moved.rectangle.extent.get(0),
                    moved.rectangle.extent.get(1)].sort((a, b) => a - b);
                expectClose(got[0], want[0], 1e-9, 1e-9);
                expectClose(got[1], want[1], 1e-9, 1e-9);
                const wantCenter = xform(reference.rectangle.center);
                expect(length(sub(moved.rectangle.center, wantCenter)))
                    .toBeLessThanOrEqual(1e-8);
            });
    });

    // The rectangle axes come from the fitted box axes, so they stay
    // orthonormal for every cloud.
    it('keeps the rectangle axes orthonormal for random clouds', () => {
        check(fc.array(wellScaledVector(3, -6, 6), { minLength: 3, maxLength: 14 }),
            (points: Vector[]) => {
                const r = getContainerLozenge3(points).rectangle;
                expectClose(dot(r.axis[0], r.axis[0]), 1, 1e-9, 1e-9);
                expectClose(dot(r.axis[1], r.axis[1]), 1, 1e-9, 1e-9);
                expectClose(dot(r.axis[0], r.axis[1]), 0, 1e-9, 1e-9);
                expect(r.extent.get(0)).toBeGreaterThanOrEqual(0);
                expect(r.extent.get(1)).toBeGreaterThanOrEqual(0);
            });
    });

    // inContainerLozenge3 is distance-to-rectangle <= radius; cross-check
    // against the DCP query, skipping near-boundary points.
    it('inContainer agrees with the point-rectangle distance', () => {
        check(fc.tuple(fc.array(wellScaledVector(3, -6, 6),
            { minLength: 3, maxLength: 10 }), wellScaledVector(3, -10, 10)),
            ([points, q]: [Vector[], Vector]) => {
                const lozenge = getContainerLozenge3(points);
                const d = distanceToRectangle(q, lozenge.rectangle);
                if (Math.abs(d - lozenge.radius) < 1e-9) {
                    return;
                }
                expect(inContainerLozenge3(q, lozenge)).toBe(d < lozenge.radius);
            });
    });
});
