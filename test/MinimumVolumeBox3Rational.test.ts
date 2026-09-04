import { describe, it, expect } from 'vitest';
import { Vector, dot, normalize } from '../src/Vector.js';
import { cross, computeOrthogonalComplement3, dotCross } from '../src/Vector3.js';
import { type OrientedBox3 } from '../src/OrientedBox.js';
import { ConvexHull3 } from '../src/ConvexHull3.js';
import { MinimumAreaBox2 } from '../src/MinimumAreaBox2.js';
import { MinimumVolumeBox3FloatingPoint } from '../src/MinimumVolumeBox3FloatingPoint.js';
import {
    MinimumVolumeBox3Rational,
    type MinimumVolumeBox3RationalResult
} from '../src/MinimumVolumeBox3Rational.js';
import {
    check, expectClose, fc, latticeVector
} from './helpers/arbitraries.js';

// The exact pipeline is much slower than the floating-point one, so the point
// sets are kept small and the number of trials low. A generous per-test
// timeout is used for the tests that run the full pipeline several times.
const SLOW = 30000;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function V(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function sub3(u: Vector, v: Vector): Vector {
    return V(u.values[0] - v.values[0], u.values[1] - v.values[1],
        u.values[2] - v.values[2]);
}

function add3(u: Vector, v: Vector): Vector {
    return V(u.values[0] + v.values[0], u.values[1] + v.values[1],
        u.values[2] + v.values[2]);
}

function scale3(s: number, u: Vector): Vector {
    return V(s * u.values[0], s * u.values[1], s * u.values[2]);
}

// A simple linear congruential generator so the "random" tests are
// deterministic and reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The 8 corners of the box with the given center, orthonormal axes and
// extents.
function boxCorners(center: Vector, axis: readonly Vector[],
    extent: readonly number[]): Vector[] {
    const corners: Vector[] = [];
    for (let i = 0; i < 8; ++i) {
        let p = center.clone();
        for (let d = 0; d < 3; ++d) {
            const sign = ((i >> d) & 1) === 0 ? -1 : +1;
            p = add3(p, scale3(sign * extent[d], axis[d]));
        }
        corners.push(p);
    }
    return corners;
}

// Every point must be inside the box. The support projections are exact, so
// only the conversion of the box to floating point introduces error.
function expectContainsAll(box: OrientedBox3, points: readonly Vector[],
    epsilon: number): void {
    for (const p of points) {
        const diff = sub3(p, box.center);
        for (let i = 0; i < 3; ++i) {
            const d = Math.abs(dot(diff, box.axis[i]));
            expect(d).toBeLessThanOrEqual(box.extent.values[i] + epsilon);
        }
    }
}

// The box axes must be an orthonormal set.
function expectOrthonormalAxes(box: OrientedBox3, epsilon: number): void {
    for (let i = 0; i < 3; ++i) {
        expect(dot(box.axis[i], box.axis[i])).toBeCloseTo(1, 10);
        for (let j = i + 1; j < 3; ++j) {
            expect(Math.abs(dot(box.axis[i], box.axis[j]))).toBeLessThan(epsilon);
        }
    }
}

// The volume of the convex hull of the points; a lower bound for the volume of
// any containing box.
function convexHullVolume(points: readonly Vector[]): number {
    const ch3 = new ConvexHull3();
    ch3.compute(points);
    if (ch3.getDimension() < 3) {
        return 0;
    }
    const hull = ch3.getHull();
    let volume = 0;
    for (let i = 0; i < hull.length; i += 3) {
        volume += dotCross(points[hull[i]], points[hull[i + 1]], points[hull[i + 2]]);
    }
    return volume / 6;
}

// An independent, brute-force upper bound on the minimum-volume box: for each
// face of the convex hull, build the box whose first axis is the face normal
// and whose other two axes come from the minimum-area rectangle of the points
// projected onto the face plane.
function bruteForceFaceAlignedVolume(points: readonly Vector[]): number {
    const ch3 = new ConvexHull3();
    ch3.compute(points);
    expect(ch3.getDimension()).toBe(3);
    const hull = ch3.getHull();

    let best = Number.MAX_VALUE;
    for (let i = 0; i < hull.length; i += 3) {
        const p0 = points[hull[i]];
        const p1 = points[hull[i + 1]];
        const p2 = points[hull[i + 2]];
        const normal = cross(sub3(p1, p0), sub3(p2, p0));
        if (normalize(normal) === 0) {
            continue;
        }

        const basis: Vector[] = [normal, Vector.zero(3), Vector.zero(3)];
        computeOrthogonalComplement3(1, basis);

        let hmin = Number.MAX_VALUE, hmax = -Number.MAX_VALUE;
        const projection: Vector[] = [];
        for (const p of points) {
            const diff = sub3(p, p0);
            const h = dot(basis[0], diff);
            hmin = Math.min(hmin, h);
            hmax = Math.max(hmax, h);
            projection.push(Vector.fromArray([dot(basis[1], diff), dot(basis[2], diff)]));
        }

        const rectangle = new MinimumAreaBox2().compute(projection);
        const area = 4 * rectangle.extent.values[0] * rectangle.extent.values[1];
        best = Math.min(best, area * (hmax - hmin));
    }
    return best;
}

// ---------------------------------------------------------------------------
// Known point sets
// ---------------------------------------------------------------------------

const unitCube: Vector[] = [
    V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(1, 1, 0),
    V(0, 0, 1), V(1, 0, 1), V(0, 1, 1), V(1, 1, 1)
];

// A right-handed orthonormal frame that is not axis aligned. The vectors
// (2,1,-2) and (1,2,2) are orthogonal and have length 3, so the frame has
// exactly representable components.
const rotatedFrame: Vector[] = (() => {
    const u0 = scale3(1 / 3, V(2, 1, -2));
    const u1 = scale3(1 / 3, V(1, 2, 2));
    const u2 = cross(u0, u1);
    return [u0, u1, u2];
})();

describe('MinimumVolumeBox3Rational', () => {
    describe('construction', () => {
        it('accepts and reports a thread count that the port ignores', () => {
            expect(new MinimumVolumeBox3Rational().getNumThreads()).toBe(0);
            expect(new MinimumVolumeBox3Rational(4).getNumThreads()).toBe(4);
        });

        it('produces the same result whether or not threads are requested', () => {
            const r0 = new MinimumVolumeBox3Rational(0).compute(unitCube, 3);
            const r4 = new MinimumVolumeBox3Rational(4).compute(unitCube, 3);
            expect(r4.volume).toBe(r0.volume);
            expect(r4.box.center.values).toEqual(r0.box.center.values);
            expect(r4.box.extent.values).toEqual(r0.box.extent.values);
        }, SLOW);
    });

    describe('invalid arguments', () => {
        it('requires at least four points and lgMaxSample >= 2', () => {
            const mvb = new MinimumVolumeBox3Rational();
            expect(() => mvb.compute(unitCube.slice(0, 3), 3)).toThrow('Invalid argument.');
            expect(() => mvb.compute(unitCube, 1)).toThrow('Invalid argument.');
        });

        it('requires a valid triangle mesh in the indexed query', () => {
            const mvb = new MinimumVolumeBox3Rational();
            const tetra = [V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)];
            const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
            expect(() => mvb.computeHull(tetra, indices.slice(0, 9), 3))
                .toThrow('Invalid argument.');
            expect(() => mvb.computeHull(tetra, indices.concat([0]), 3))
                .toThrow('Invalid argument.');
            expect(() => mvb.computeHull(tetra.slice(0, 3), indices, 3))
                .toThrow('Invalid argument.');
        });
    });

    describe('axis-aligned boxes', () => {
        it('recovers the unit cube with no rounding error at all', () => {
            const result = new MinimumVolumeBox3Rational().compute(unitCube, 4);
            expect(result.dimension).toBe(3);
            // The exact pipeline produces exactly representable numbers here,
            // so the comparisons are exact rather than approximate.
            expect(result.volume).toBe(1);
            expect(result.box.center.values).toEqual([0.5, 0.5, 0.5]);
            expect(result.box.extent.values).toEqual([0.5, 0.5, 0.5]);
            expectOrthonormalAxes(result.box, 1e-15);
            expectContainsAll(result.box, unitCube, 0);
        }, SLOW);

        it('recovers a non-cubic axis-aligned box', () => {
            const center = V(-3, 7, 2.5);
            const axis = [V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)];
            const extent = [2, 0.5, 3];
            const points = boxCorners(center, axis, extent);
            // Add interior points; they must not change the box.
            const rnd = makeRandom(7);
            for (let i = 0; i < 6; ++i) {
                points.push(add3(center, V(
                    (2 * rnd() - 1) * extent[0],
                    (2 * rnd() - 1) * extent[1],
                    (2 * rnd() - 1) * extent[2])));
            }

            const result = new MinimumVolumeBox3Rational().compute(points, 3);
            expect(result.dimension).toBe(3);
            expect(result.volume).toBeCloseTo(8 * extent[0] * extent[1] * extent[2], 10);
            const sorted = [...result.box.extent.values].sort((a, b) => a - b);
            const expected = [...extent].sort((a, b) => a - b);
            for (let i = 0; i < 3; ++i) {
                expect(sorted[i]).toBeCloseTo(expected[i], 10);
            }
            expectContainsAll(result.box, points, 1e-10);
        }, SLOW);

        it('scales the volume by the cube of a uniform scale factor', () => {
            const mvb = new MinimumVolumeBox3Rational();
            const v1 = mvb.compute(unitCube, 3).volume;
            const v2 = mvb.compute(unitCube.map(p => scale3(2, p)), 3).volume;
            expect(v1).toBe(1);
            expect(v2).toBe(8);
        }, SLOW);
    });

    describe('rotated boxes', () => {
        it('recovers the analytic minimum volume of a rotated box', () => {
            const center = V(1, -2, 0.5);
            const extent = [1.5, 1, 0.25];
            const points = boxCorners(center, rotatedFrame, extent);
            const result = new MinimumVolumeBox3Rational().compute(points, 4);
            expect(result.dimension).toBe(3);

            const analytic = 8 * extent[0] * extent[1] * extent[2];
            expect(result.volume).toBeCloseTo(analytic, 10);
            expect(result.box.center.values[0]).toBeCloseTo(center.values[0], 10);
            expect(result.box.center.values[1]).toBeCloseTo(center.values[1], 10);
            expect(result.box.center.values[2]).toBeCloseTo(center.values[2], 10);

            const sorted = [...result.box.extent.values].sort((a, b) => a - b);
            const expected = [...extent].sort((a, b) => a - b);
            for (let i = 0; i < 3; ++i) {
                expect(sorted[i]).toBeCloseTo(expected[i], 10);
            }
            // Each recovered axis must be parallel to one of the source axes.
            for (let i = 0; i < 3; ++i) {
                const best = Math.max(...rotatedFrame.map(
                    a => Math.abs(dot(a, result.box.axis[i]))));
                expect(best).toBeCloseTo(1, 8);
            }
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, points, 1e-10);
        }, SLOW);
    });

    describe('tetrahedra', () => {
        it('computes a box for the indexed-hull query', () => {
            const vertices = [V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)];
            // Counterclockwise when viewed from outside the tetrahedron.
            const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
            const result = new MinimumVolumeBox3Rational()
                .computeHull(vertices, indices, 4);

            // The axis-aligned box of these vertices is the unit cube, so the
            // minimum volume box is no larger than 1 and no smaller than the
            // tetrahedron volume 1/6.
            expect(result.volume).toBeGreaterThan(1 / 6);
            expect(result.volume).toBeLessThanOrEqual(1);
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, vertices, 1e-12);
        }, SLOW);

        it('agrees with the point-set query for the same tetrahedron', () => {
            const vertices = [V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)];
            const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
            const fromPoints = new MinimumVolumeBox3Rational().compute(vertices, 4);
            const fromHull = new MinimumVolumeBox3Rational()
                .computeHull(vertices, indices, 4);
            expect(fromPoints.dimension).toBe(3);
            // Both queries run the same exact pipeline on the same polytope.
            expect(fromPoints.volume).toBe(fromHull.volume);
        }, SLOW);

        it('handles a regular tetrahedron', () => {
            const vertices = [V(1, 1, 1), V(1, -1, -1), V(-1, 1, -1), V(-1, -1, 1)];
            const result = new MinimumVolumeBox3Rational().compute(vertices, 4);
            expect(result.dimension).toBe(3);
            // The tetrahedron has volume 8/3. Its minimum-volume box is the
            // cube spanned by its edge midpoints, of volume 8.
            expect(result.volume).toBeGreaterThan(8 / 3);
            expect(result.volume).toBeLessThanOrEqual(8 + 1e-9);
            expectContainsAll(result.box, vertices, 1e-10);
        }, SLOW);
    });

    describe('degenerate inputs', () => {
        it('reports dimension 0 for coincident points', () => {
            const p = V(2, -1, 4);
            const points = [p.clone(), p.clone(), p.clone(), p.clone(), p.clone()];
            const result = new MinimumVolumeBox3Rational().compute(points, 3);
            expect(result.dimension).toBe(0);
            expect(result.volume).toBe(0);
            expect(result.box.center.values).toEqual([2, -1, 4]);
            expect(result.box.extent.values).toEqual([0, 0, 0]);
            expect(result.box.axis[0].values).toEqual([1, 0, 0]);
            expect(result.box.axis[1].values).toEqual([0, 1, 0]);
            expect(result.box.axis[2].values).toEqual([0, 0, 1]);
        });

        it('reports dimension 1 for collinear points', () => {
            const origin = V(1, 1, 1);
            const direction = V(2, -1, 2);  // length 3
            const points = [0, 0.25, 0.5, 0.75, 1].map(
                t => add3(origin, scale3(t, direction)));
            const result = new MinimumVolumeBox3Rational().compute(points, 3);
            expect(result.dimension).toBe(1);
            expect(result.volume).toBe(0);
            expect(result.box.extent.values[0]).toBeCloseTo(1.5, 12);
            expect(result.box.extent.values[1]).toBe(0);
            expect(result.box.extent.values[2]).toBe(0);
            expect(result.box.center.values[0]).toBeCloseTo(2, 12);
            expect(result.box.center.values[1]).toBeCloseTo(0.5, 12);
            expect(result.box.center.values[2]).toBeCloseTo(2, 12);
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, points, 1e-12);
        });

        it('reports dimension 2 for coplanar points with a quadrilateral hull', () => {
            // A unit square in the plane z = x + y.
            const points = [
                V(0, 0, 0), V(1, 0, 1), V(1, 1, 2), V(0, 1, 1),
                V(0.5, 0.5, 1)
            ];
            const result = new MinimumVolumeBox3Rational().compute(points, 3);
            expect(result.dimension).toBe(2);
            expect(result.volume).toBe(0);
            const zeroCount = [...result.box.extent.values].filter(e => e === 0).length;
            expect(zeroCount).toBe(1);
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, points, 1e-12);
        });

        it('reports dimension 2 for coplanar points with a triangular hull', () => {
            // This exercises the upstream Newell-loop bug shared by both MVB3
            // pipelines: for a triangular hull the upstream sum of cross
            // products is exactly zero, which produces a degenerate basis. The
            // port closes the loop, so the basis and the box are well defined.
            const points = [
                V(0, 0, 0), V(2, 0, 2), V(0, 2, 2),
                V(0.5, 0.5, 1), V(0.25, 1, 1.25)
            ];
            const result = new MinimumVolumeBox3Rational().compute(points, 3);
            expect(result.dimension).toBe(2);
            expect(result.volume).toBe(0);
            for (const e of result.box.extent.values) {
                expect(Number.isFinite(e)).toBe(true);
            }
            for (const axis of result.box.axis) {
                for (const value of axis.values) {
                    expect(Number.isFinite(value)).toBe(true);
                }
            }
            expect(result.box.extent.values[2]).toBe(0);
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, points, 1e-12);
        });

        it('handles a polytope with coplanar hull faces (a box with a split face)',
            () => {
                // The unit cube plus points in the interior of the z = 1 face.
                // The hull triangulation of that face has coplanar triangles,
                // which exercises removeCoplanarTriangleAdjacencies. The exact
                // coplanarity and colinearity tests must remove them all, so
                // the answer is again the unit cube.
                const points = unitCube.concat([
                    V(0.5, 0, 1), V(0, 0.5, 1), V(1, 0.5, 1), V(0.5, 1, 1),
                    V(0.5, 0.5, 1)
                ]);
                const result = new MinimumVolumeBox3Rational().compute(points, 3);
                expect(result.dimension).toBe(3);
                expect(result.volume).toBe(1);
                expectContainsAll(result.box, points, 0);
            }, SLOW);
    });

    describe('randomized cross-checks', () => {
        it('bounds the hull and is no larger than the best face-aligned box', () => {
            const rnd = makeRandom(20260127);
            for (let trial = 0; trial < 3; ++trial) {
                const points: Vector[] = [];
                for (let i = 0; i < 12; ++i) {
                    points.push(V(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2));
                }

                const result: MinimumVolumeBox3RationalResult =
                    new MinimumVolumeBox3Rational().compute(points, 3);
                expect(result.dimension).toBe(3);

                expectContainsAll(result.box, points, 1e-9);
                expectOrthonormalAxes(result.box, 1e-10);

                const hullVolume = convexHullVolume(points);
                const faceAligned = bruteForceFaceAlignedVolume(points);
                expect(result.volume).toBeGreaterThan(hullVolume * (1 - 1e-9));
                expect(result.volume).toBeLessThanOrEqual(faceAligned * (1 + 1e-9));

                // The reported volume agrees with the product of the extents.
                const e = result.box.extent.values;
                expect(result.volume).toBeCloseTo(8 * e[0] * e[1] * e[2], 8);
            }
        }, SLOW);

        it('agrees closely with the floating-point pipeline', () => {
            // The two pipelines are not expected to agree exactly. The
            // floating-point pipeline interpolates unit-length triangle
            // normals whereas the exact pipeline interpolates the unnormalized
            // normals, so the two sample different points of the same level
            // curve. Both are approximations of the same minimum, so the
            // volumes are close but either may be the smaller one.
            const rnd = makeRandom(4242);
            for (let trial = 0; trial < 3; ++trial) {
                const points: Vector[] = [];
                for (let i = 0; i < 12; ++i) {
                    points.push(V(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2));
                }
                const exact = new MinimumVolumeBox3Rational().compute(points, 3);
                const approx = new MinimumVolumeBox3FloatingPoint().compute(points, 3);
                expect(exact.dimension).toBe(approx.dimension);
                expect(Math.abs(exact.volume / approx.volume - 1)).toBeLessThan(0.05);
                expectContainsAll(exact.box, points, 1e-9);
            }
        }, SLOW);

        it('agrees with the floating-point pipeline for a rotated box', () => {
            // For a point set that is the corner set of a box, both pipelines
            // find the box itself.
            const points = boxCorners(V(1, -2, 0.5), rotatedFrame, [1.5, 1, 0.25]);
            const exact = new MinimumVolumeBox3Rational().compute(points, 3);
            const approx = new MinimumVolumeBox3FloatingPoint().compute(points, 3);
            expect(exact.volume).toBeCloseTo(approx.volume, 10);
            expect(exact.volume).toBeCloseTo(3, 10);
        }, SLOW);

        it('does not get worse as the number of samples increases', () => {
            // The sample sets are dyadic subdivisions, so the samples used for
            // lgMaxSample = k are a subset of those used for k + 1 and the
            // minimum volume is nonincreasing in k. With exact arithmetic the
            // comparison is exact.
            const rnd = makeRandom(98765);
            for (let trial = 0; trial < 2; ++trial) {
                const points: Vector[] = [];
                for (let i = 0; i < 10; ++i) {
                    points.push(V(rnd(), rnd() * 2, rnd() * 0.5));
                }
                const mvb = new MinimumVolumeBox3Rational();
                const coarse = mvb.compute(points, 2);
                const fine = mvb.compute(points, 4);
                expect(fine.volume).toBeLessThanOrEqual(coarse.volume * (1 + 1e-12));
                expectContainsAll(fine.box, points, 1e-9);
            }
        }, SLOW);

        it('is invariant under a rigid motion with exact rational coefficients', () => {
            const rnd = makeRandom(555);
            const points: Vector[] = [];
            for (let i = 0; i < 10; ++i) {
                points.push(V(rnd() * 2, rnd() * 3, rnd()));
            }
            const translation = V(100, -50, 25);
            const moved = points.map(p => add3(V(
                dot(p, rotatedFrame[0]),
                dot(p, rotatedFrame[1]),
                dot(p, rotatedFrame[2])), translation));

            const mvb = new MinimumVolumeBox3Rational();
            const v0 = mvb.compute(points, 3).volume;
            const v1 = mvb.compute(moved, 3).volume;
            expect(v1).toBeCloseTo(v0, 8);
        }, SLOW);

        it('reuses the query object across data sets', () => {
            const mvb = new MinimumVolumeBox3Rational();
            const first = mvb.compute(unitCube, 3);
            const scaled = unitCube.map(p => scale3(2, p));
            const second = mvb.compute(scaled, 3);
            const again = mvb.compute(unitCube, 3);
            expect(first.volume).toBe(1);
            expect(second.volume).toBe(8);
            expect(again.volume).toBe(first.volume);
        }, SLOW);
    });
});

// ---------------------------------------------------------------------------
// Verification (V12): property-based checks.
//
// The generators place the points on an integer lattice, so the exact
// (BSNumber / BSRational) pipeline sees exact input and its only inexactness
// is the final conversion of the box to binary64. Point counts are kept small
// because every candidate volume is an exact rational computation; the heavy
// properties carry an explicit timeout.
// ---------------------------------------------------------------------------

// A lattice point cloud whose convex hull is 3-dimensional.
const latticeCloud = (count: number, spread = 5): fc.Arbitrary<Vector[]> =>
    fc.array(latticeVector(3, -spread, spread),
        { minLength: count, maxLength: count })
        .filter(points => {
            const ch3 = new ConvexHull3();
            ch3.compute(points);
            return ch3.getDimension() === 3;
        });

// A lattice point cloud lying exactly in the plane z = a*x + b*y. 'hullSize'
// selects a triangular planar hull (the upstream Newell-loop case) or a
// quadrilateral one.
const coplanarLatticeCloud = (hullSize: 'triangle' | 'polygon'):
    fc.Arbitrary<Vector[]> =>
    fc.tuple(fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 }),
        fc.array(fc.tuple(fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: -4, max: 4 })), { minLength: 2, maxLength: 4 }))
        .map(([a, b, extra]) => {
            const lift = (x: number, y: number): Vector => V(x, y, a * x + b * y);
            const points = hullSize === 'triangle'
                ? [lift(0, 0), lift(6, 0), lift(0, 6)]
                : [lift(0, 0), lift(6, 0), lift(6, 6), lift(0, 6)];
            for (const [x, y] of extra) {
                points.push(lift(1 + (x + 4) % 3, 1 + (y + 4) % 3));
            }
            return points;
        });

// The volume of the axis-aligned bounding box of the points.
function alignedBoxVolume(points: readonly Vector[]): number {
    const lo = [Infinity, Infinity, Infinity];
    const hi = [-Infinity, -Infinity, -Infinity];
    for (const p of points) {
        for (let j = 0; j < 3; ++j) {
            lo[j] = Math.min(lo[j], p.values[j]);
            hi[j] = Math.max(hi[j], p.values[j]);
        }
    }
    return (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
}

describe('MinimumVolumeBox3Rational verification', () => {
    it('returns a box that contains the points and is no larger than the '
        + 'aligned or face-aligned candidates', () => {
        check(latticeCloud(8), points => {
            const result = new MinimumVolumeBox3Rational().compute(points, 2);
            expect(result.dimension).toBe(3);
            expectOrthonormalAxes(result.box, 1e-10);

            // The whole search is exact; the only rounding is the final
            // conversion of the rational box to binary64, so a tolerance of
            // 1e-9 on coordinates of magnitude at most 5 is generous.
            expectContainsAll(result.box, points, 1e-9);

            const e = result.box.extent.values;
            expectClose(result.volume, 8 * e[0] * e[1] * e[2], 1e-8, 1e-9);

            const faceAligned = bruteForceFaceAlignedVolume(points);
            expect(result.volume).toBeLessThanOrEqual(faceAligned * (1 + 1e-8));
            expect(result.volume)
                .toBeLessThanOrEqual(alignedBoxVolume(points) * (1 + 1e-8));
            expect(result.volume)
                .toBeGreaterThanOrEqual(convexHullVolume(points) * (1 - 1e-8));
        }, 40);
    }, SLOW);

    it('agrees with the floating-point pipeline on lattice input', () => {
        check(latticeCloud(8), points => {
            const lgMaxSample = 4;
            const exact = new MinimumVolumeBox3Rational()
                .compute(points, lgMaxSample);
            const approx = new MinimumVolumeBox3FloatingPoint()
                .compute(points, lgMaxSample);
            expect(approx.dimension).toBe(exact.dimension);

            // Both results must be valid containing boxes that are no worse
            // than the best hull-face-aligned candidate. These are the
            // assertions with content; the two pipelines are not required to
            // return the same box.
            expectContainsAll(approx.box, points, 1e-9);
            expectContainsAll(exact.box, points, 1e-9);
            const faceAligned = bruteForceFaceAlignedVolume(points);
            expect(approx.volume).toBeLessThanOrEqual(faceAligned * (1 + 1e-8));
            expect(exact.volume).toBeLessThanOrEqual(faceAligned * (1 + 1e-8));

            // The two pipelines search the same family of level curves but
            // with different parameterizations: the floating-point pipeline
            // interpolates unit-length face normals while the exact pipeline
            // interpolates the unnormalized integer normals, because
            // normalization is not an exact operation. The sampled points
            // therefore differ and the two minima differ by the accuracy of
            // the sampling search rather than by the accuracy of the
            // arithmetic. Over 120 random lattice clouds of 8 points the
            // relative gap has median 0.09%, 95th percentile 1.0% and maximum
            // 3.5% at lgMaxSample = 4; the bound below is loose enough to
            // cover the tail and still catches a pipeline that diverges.
            const gap = Math.abs(approx.volume - exact.volume)
                / Math.max(1, exact.volume);
            expect(gap).toBeLessThan(0.15);
        }, 30);
    }, SLOW);

    it('scales the volume by the cube of a uniform scale', () => {
        check(fc.tuple(latticeCloud(7), fc.integer({ min: 2, max: 4 })),
            ([points, s]) => {
                const mvb = new MinimumVolumeBox3Rational();
                const v0 = mvb.compute(points, 2).volume;
                const v1 = mvb.compute(points.map(
                    p => V(s * p.values[0], s * p.values[1], s * p.values[2])),
                2).volume;
                // The scaled problem is the same problem in exact arithmetic,
                // so only the final conversion to binary64 differs.
                expectClose(v1, s * s * s * v0, 1e-9, 1e-12);
            }, 20);
    }, SLOW);

    it('does not get worse when the sample count increases', () => {
        check(latticeCloud(7), points => {
            const mvb = new MinimumVolumeBox3Rational();
            const coarse = mvb.compute(points, 2).volume;
            const fine = mvb.compute(points, 3).volume;
            // The dyadic samples for lgMaxSample = 2 are a subset of those for
            // lgMaxSample = 3, and the comparison of candidate volumes is
            // exact, so the minimum cannot increase.
            expect(fine).toBeLessThanOrEqual(coarse * (1 + 1e-12));
        }, 20);
    }, SLOW);

    it('classifies coplanar point sets as dimension 2 with a finite basis', () => {
        // The triangular-hull case is the upstream Newell-loop defect: the
        // loop starts at i1 = 1, so for a 3-gon the sum of cross products is
        // Cross(P2,P1) + Cross(P1,P2) = 0 and the plane basis is degenerate.
        check(fc.oneof(coplanarLatticeCloud('triangle'),
            coplanarLatticeCloud('polygon')), points => {
            const result = new MinimumVolumeBox3Rational().compute(points, 2);
            expect(result.dimension).toBe(2);
            expect(result.volume).toBe(0);
            expect(result.box.extent.values[2]).toBe(0);
            expect(result.box.extent.values[0]).toBeGreaterThan(0);
            expect(result.box.extent.values[1]).toBeGreaterThan(0);
            for (const value of result.box.center.values) {
                expect(Number.isFinite(value)).toBe(true);
            }
            expectOrthonormalAxes(result.box, 1e-10);
            expectContainsAll(result.box, points, 1e-9);
        }, 40);
    }, SLOW);

    it('computes a nonzero Newell normal for a triangular coplanar hull', () => {
        check(fc.tuple(fc.integer({ min: -3, max: 3 }),
            fc.integer({ min: -3, max: 3 })), ([a, b]) => {
            const lift = (x: number, y: number): Vector => V(x, y, a * x + b * y);
            const points = [lift(0, 0), lift(4, 0), lift(0, 4), lift(1, 1)];
            const result = new MinimumVolumeBox3Rational().compute(points, 2);
            expect(result.dimension).toBe(2);
            const n = result.box.axis[2];
            expect(dot(n, n)).toBeCloseTo(1, 12);
            const exact = V(a, b, -1);
            normalize(exact);
            expect(Math.abs(dot(n, exact))).toBeCloseTo(1, 12);
        }, 49);
    }, SLOW);

    it('classifies collinear and coincident point sets', () => {
        check(fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -4, 4),
            fc.array(fc.integer({ min: -6, max: 6 }),
                { minLength: 4, maxLength: 7 })),
        ([origin, direction, ts]) => {
            const d = direction.values;
            if (d[0] === 0 && d[1] === 0 && d[2] === 0) {
                const points = ts.map(() => origin.clone());
                const result = new MinimumVolumeBox3Rational().compute(points, 2);
                expect(result.dimension).toBe(0);
                expect(result.volume).toBe(0);
                expect(result.box.extent.values).toEqual([0, 0, 0]);
                expect(result.box.center.values).toEqual([...origin.values]);
                return;
            }
            const points = ts.map(t => V(
                origin.values[0] + t * d[0],
                origin.values[1] + t * d[1],
                origin.values[2] + t * d[2]));
            const distinct = new Set(points.map(p => p.values.join(',')));
            const result = new MinimumVolumeBox3Rational().compute(points, 2);
            expect(result.dimension).toBe(distinct.size === 1 ? 0 : 1);
            expect(result.volume).toBe(0);
            expect(result.box.extent.values[1]).toBe(0);
            expect(result.box.extent.values[2]).toBe(0);
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, points, 1e-9);
        }, 40);
    }, SLOW);
});
