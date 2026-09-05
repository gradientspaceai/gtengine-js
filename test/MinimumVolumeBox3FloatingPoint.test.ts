import { describe, it, expect } from 'vitest';
import { Vector, dot, normalize } from '../src/Vector.js';
import { cross, computeOrthogonalComplement3, dotCross } from '../src/Vector3.js';
import { type OrientedBox3 } from '../src/OrientedBox.js';
import { ConvexHull3 } from '../src/ConvexHull3.js';
import { MinimumAreaBox2 } from '../src/MinimumAreaBox2.js';
import {
    MinimumVolumeBox3FloatingPoint,
    type MinimumVolumeBox3FloatingPointResult
} from '../src/MinimumVolumeBox3FloatingPoint.js';
import {
    check, expectClose, fc, latticeVector, rotationFrame
} from './helpers/arbitraries.js';

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

// Every point must be inside the box, allowing for floating-point rounding.
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

// The volume of the convex hull of the points, computed by summing the signed
// volumes of the tetrahedra formed by the origin and the hull triangles. This
// is a lower bound for the volume of any containing box.
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
// projected onto the face plane. The true minimum-volume box is no larger than
// the smallest of these.
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
        const volume = area * (hmax - hmin);
        best = Math.min(best, volume);
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

// A right-handed orthonormal frame that is not axis aligned.
const rotatedFrame: Vector[] = (() => {
    const u0 = V(2, 1, -2);
    normalize(u0);
    const u1 = V(1, 2, 2);
    // Gram-Schmidt (the vectors above are already orthogonal, so this is a
    // normalization).
    const proj = dot(u1, u0);
    const u1p = sub3(u1, scale3(proj, u0));
    normalize(u1p);
    const u2 = cross(u0, u1p);
    return [u0, u1p, u2];
})();

describe('MinimumVolumeBox3FloatingPoint', () => {
    describe('construction', () => {
        it('accepts and reports a thread count that the port ignores', () => {
            expect(new MinimumVolumeBox3FloatingPoint().getNumThreads()).toBe(0);
            expect(new MinimumVolumeBox3FloatingPoint(4).getNumThreads()).toBe(4);
        });

        it('produces the same result whether or not threads are requested', () => {
            const r0 = new MinimumVolumeBox3FloatingPoint(0).compute(unitCube, 3);
            const r4 = new MinimumVolumeBox3FloatingPoint(4).compute(unitCube, 3);
            expect(r4.volume).toBe(r0.volume);
            expect(r4.box.center.values).toEqual(r0.box.center.values);
            expect(r4.box.extent.values).toEqual(r0.box.extent.values);
        });
    });

    describe('invalid arguments', () => {
        it('requires at least four points and lgMaxSample >= 2', () => {
            const mvb = new MinimumVolumeBox3FloatingPoint();
            expect(() => mvb.compute(unitCube.slice(0, 3), 3)).toThrow('Invalid argument.');
            expect(() => mvb.compute(unitCube, 1)).toThrow('Invalid argument.');
        });

        it('requires a valid triangle mesh in the indexed query', () => {
            const mvb = new MinimumVolumeBox3FloatingPoint();
            const tetra = [V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)];
            const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
            // Too few indices.
            expect(() => mvb.computeHull(tetra, indices.slice(0, 9), 3))
                .toThrow('Invalid argument.');
            // Not a multiple of three.
            expect(() => mvb.computeHull(tetra, indices.concat([0]), 3))
                .toThrow('Invalid argument.');
            // Too few vertices.
            expect(() => mvb.computeHull(tetra.slice(0, 3), indices, 3))
                .toThrow('Invalid argument.');
        });
    });

    describe('axis-aligned boxes', () => {
        it('recovers the unit cube exactly', () => {
            const result = new MinimumVolumeBox3FloatingPoint().compute(unitCube, 4);
            expect(result.dimension).toBe(3);
            expect(result.volume).toBeCloseTo(1, 12);
            expect(result.box.center.values[0]).toBeCloseTo(0.5, 12);
            expect(result.box.center.values[1]).toBeCloseTo(0.5, 12);
            expect(result.box.center.values[2]).toBeCloseTo(0.5, 12);
            const extents = [...result.box.extent.values].sort((a, b) => a - b);
            for (const e of extents) {
                expect(e).toBeCloseTo(0.5, 12);
            }
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, unitCube, 1e-12);
        });

        it('recovers a non-cubic axis-aligned box', () => {
            const center = V(-3, 7, 2.5);
            const axis = [V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)];
            const extent = [2, 0.5, 3];
            const points = boxCorners(center, axis, extent);
            // Add interior points; they must not change the box.
            const rnd = makeRandom(7);
            for (let i = 0; i < 8; ++i) {
                points.push(add3(center, V(
                    (2 * rnd() - 1) * extent[0],
                    (2 * rnd() - 1) * extent[1],
                    (2 * rnd() - 1) * extent[2])));
            }

            const result = new MinimumVolumeBox3FloatingPoint().compute(points, 4);
            expect(result.dimension).toBe(3);
            expect(result.volume).toBeCloseTo(8 * extent[0] * extent[1] * extent[2], 10);
            const sorted = [...result.box.extent.values].sort((a, b) => a - b);
            const expected = [...extent].sort((a, b) => a - b);
            for (let i = 0; i < 3; ++i) {
                expect(sorted[i]).toBeCloseTo(expected[i], 10);
            }
            expectContainsAll(result.box, points, 1e-10);
        });
    });

    describe('rotated boxes', () => {
        it('recovers the analytic minimum volume of a rotated box', () => {
            const center = V(1, -2, 0.5);
            const extent = [1.5, 1, 0.25];
            const points = boxCorners(center, rotatedFrame, extent);
            const result = new MinimumVolumeBox3FloatingPoint().compute(points, 5);
            expect(result.dimension).toBe(3);

            const analytic = 8 * extent[0] * extent[1] * extent[2];
            expect(result.volume).toBeCloseTo(analytic, 8);
            expect(result.box.center.values[0]).toBeCloseTo(center.values[0], 8);
            expect(result.box.center.values[1]).toBeCloseTo(center.values[1], 8);
            expect(result.box.center.values[2]).toBeCloseTo(center.values[2], 8);

            // The recovered axes must match the source frame up to permutation
            // and sign; pair them by extent.
            const sorted = [...result.box.extent.values].sort((a, b) => a - b);
            const expected = [...extent].sort((a, b) => a - b);
            for (let i = 0; i < 3; ++i) {
                expect(sorted[i]).toBeCloseTo(expected[i], 8);
            }
            // Each recovered axis must be parallel to one of the source axes.
            for (let i = 0; i < 3; ++i) {
                const best = Math.max(...rotatedFrame.map(
                    a => Math.abs(dot(a, result.box.axis[i]))));
                expect(best).toBeCloseTo(1, 6);
            }
            expectOrthonormalAxes(result.box, 1e-10);
            expectContainsAll(result.box, points, 1e-8);
        });

        it('recovers a rotated box that also contains interior points', () => {
            const center = V(0, 0, 0);
            const extent = [3, 2, 1];
            const points = boxCorners(center, rotatedFrame, extent);
            const rnd = makeRandom(31);
            for (let i = 0; i < 10; ++i) {
                let p = center.clone();
                for (let d = 0; d < 3; ++d) {
                    p = add3(p, scale3((2 * rnd() - 1) * extent[d], rotatedFrame[d]));
                }
                points.push(p);
            }
            const result = new MinimumVolumeBox3FloatingPoint().compute(points, 4);
            expect(result.volume).toBeCloseTo(8 * 3 * 2 * 1, 6);
            expectContainsAll(result.box, points, 1e-8);
        });
    });

    describe('tetrahedra', () => {
        it('computes a box for the indexed-hull query', () => {
            const vertices = [V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)];
            // Counterclockwise when viewed from outside the tetrahedron.
            const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
            const result = new MinimumVolumeBox3FloatingPoint()
                .computeHull(vertices, indices, 4);

            // The axis-aligned box of these vertices is the unit cube, so the
            // minimum volume box is no larger than 1 and no smaller than the
            // tetrahedron volume 1/6.
            expect(result.volume).toBeGreaterThan(1 / 6);
            expect(result.volume).toBeLessThanOrEqual(1 + 1e-12);
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, vertices, 1e-12);
        });

        it('agrees with the point-set query for the same tetrahedron', () => {
            const vertices = [V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)];
            const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
            const fromPoints = new MinimumVolumeBox3FloatingPoint().compute(vertices, 4);
            const fromHull = new MinimumVolumeBox3FloatingPoint()
                .computeHull(vertices, indices, 4);
            expect(fromPoints.dimension).toBe(3);
            expect(fromPoints.volume).toBeCloseTo(fromHull.volume, 12);
        });

        it('handles a regular tetrahedron', () => {
            const vertices = [V(1, 1, 1), V(1, -1, -1), V(-1, 1, -1), V(-1, -1, 1)];
            const result = new MinimumVolumeBox3FloatingPoint().compute(vertices, 5);
            expect(result.dimension).toBe(3);
            // The tetrahedron has volume 8/3 and its axis-aligned bounding box
            // has volume 8. The minimum-volume box of a regular tetrahedron is
            // the cube spanned by its edge midpoints, of volume 8.
            expect(result.volume).toBeGreaterThan(8 / 3);
            expect(result.volume).toBeLessThanOrEqual(8 + 1e-9);
            expectContainsAll(result.box, vertices, 1e-10);
        });
    });

    describe('degenerate inputs', () => {
        it('reports dimension 0 for coincident points', () => {
            const p = V(2, -1, 4);
            const points = [p.clone(), p.clone(), p.clone(), p.clone(), p.clone()];
            const result = new MinimumVolumeBox3FloatingPoint().compute(points, 3);
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
            const result = new MinimumVolumeBox3FloatingPoint().compute(points, 3);
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
            const result = new MinimumVolumeBox3FloatingPoint().compute(points, 3);
            expect(result.dimension).toBe(2);
            expect(result.volume).toBe(0);
            const zeroCount = [...result.box.extent.values].filter(e => e === 0).length;
            expect(zeroCount).toBe(1);
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, points, 1e-12);
        });

        it('reports dimension 2 for coplanar points with a triangular hull', () => {
            // This is the upstream Newell-loop bug: for a triangular hull the
            // upstream sum of cross products is exactly zero, which produces a
            // degenerate basis. The port closes the loop, so the basis and the
            // box are well defined.
            const points = [
                V(0, 0, 0), V(2, 0, 2), V(0, 2, 2),
                V(0.5, 0.5, 1), V(0.25, 1, 1.25)
            ];
            const result = new MinimumVolumeBox3FloatingPoint().compute(points, 3);
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
                // which exercises removeCoplanarTriangleAdjacencies.
                const points = unitCube.concat([
                    V(0.5, 0, 1), V(0, 0.5, 1), V(1, 0.5, 1), V(0.5, 1, 1),
                    V(0.5, 0.5, 1)
                ]);
                const result = new MinimumVolumeBox3FloatingPoint().compute(points, 4);
                expect(result.dimension).toBe(3);
                expect(result.volume).toBeCloseTo(1, 10);
                expectContainsAll(result.box, points, 1e-10);
            });
    });

    describe('randomized cross-checks', () => {
        it('bounds the hull and is no larger than the best face-aligned box', () => {
            const rnd = makeRandom(20260127);
            for (let trial = 0; trial < 6; ++trial) {
                const points: Vector[] = [];
                for (let i = 0; i < 12; ++i) {
                    points.push(V(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2));
                }

                const result: MinimumVolumeBox3FloatingPointResult =
                    new MinimumVolumeBox3FloatingPoint().compute(points, 3);
                expect(result.dimension).toBe(3);

                // The box contains every input point.
                expectContainsAll(result.box, points, 1e-9);
                expectOrthonormalAxes(result.box, 1e-10);

                // The volume is bounded below by the hull volume and above by
                // the smallest face-aligned candidate box.
                const hullVolume = convexHullVolume(points);
                const faceAligned = bruteForceFaceAlignedVolume(points);
                expect(result.volume).toBeGreaterThan(hullVolume * (1 - 1e-9));
                expect(result.volume).toBeLessThanOrEqual(faceAligned * (1 + 1e-9));

                // The reported volume agrees with the product of the extents.
                const e = result.box.extent.values;
                expect(result.volume).toBeCloseTo(8 * e[0] * e[1] * e[2], 8);
            }
        });

        it('does not get worse as the number of samples increases', () => {
            // The sample sets are dyadic subdivisions, so the samples used for
            // lgMaxSample = k are a subset of those used for k + 1 and the
            // minimum volume is nonincreasing in k.
            const rnd = makeRandom(98765);
            for (let trial = 0; trial < 3; ++trial) {
                const points: Vector[] = [];
                for (let i = 0; i < 10; ++i) {
                    points.push(V(rnd(), rnd() * 2, rnd() * 0.5));
                }
                const mvb = new MinimumVolumeBox3FloatingPoint();
                const coarse = mvb.compute(points, 2);
                const fine = mvb.compute(points, 4);
                expect(fine.volume).toBeLessThanOrEqual(coarse.volume * (1 + 1e-9));
                expectContainsAll(fine.box, points, 1e-9);
            }
        });

        it('is invariant under rigid motions of the point set', () => {
            const rnd = makeRandom(555);
            const points: Vector[] = [];
            for (let i = 0; i < 12; ++i) {
                points.push(V(rnd() * 2, rnd() * 3, rnd()));
            }
            const translation = V(100, -50, 25);
            const moved = points.map(p => {
                const q = V(
                    dot(p, rotatedFrame[0]),
                    dot(p, rotatedFrame[1]),
                    dot(p, rotatedFrame[2]));
                return add3(q, translation);
            });

            const mvb = new MinimumVolumeBox3FloatingPoint();
            const v0 = mvb.compute(points, 3).volume;
            const v1 = mvb.compute(moved, 3).volume;
            expect(v1).toBeCloseTo(v0, 8);
        });

        it('reuses the query object across data sets', () => {
            const mvb = new MinimumVolumeBox3FloatingPoint();
            const first = mvb.compute(unitCube, 3);
            const scaled = unitCube.map(p => scale3(2, p));
            const second = mvb.compute(scaled, 3);
            const again = mvb.compute(unitCube, 3);
            expect(first.volume).toBeCloseTo(1, 10);
            expect(second.volume).toBeCloseTo(8, 10);
            expect(again.volume).toBeCloseTo(first.volume, 12);
        });
    });
});

// ---------------------------------------------------------------------------
// Verification (V12): property-based checks.
//
// The generators place the points on an integer lattice. Small integers are
// exact in binary64, so the exact convex-hull predicates and the exact
// (BSNumber) coplanarity tests inside the query see clean input, and the
// remaining error is only that of the floating-point axis search. Compare with
// the generator-hygiene note in VERIFYING.md: the default vector() generator
// emits subnormal components, which underflow in the covariance-like products
// this algorithm computes.
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

// A lattice point cloud of 'count' + 1 points: 'count' lattice points plus
// the midpoint of two of them. Halves of small integers are exact in
// binary64, so the extra point lies exactly on the segment between the two.
// When that segment is on the hull, the midpoint is a hull vertex that is
// only weakly extreme.
const cloudWithMidpoint = (count: number): fc.Arbitrary<Vector[]> =>
    fc.tuple(latticeCloud(count), fc.nat(), fc.nat()).map(([points, i, j]) => {
        const p = points[i % count];
        const q = points[(i + 1 + j % (count - 1)) % count];
        return points.concat([V(
            0.5 * (p.values[0] + q.values[0]),
            0.5 * (p.values[1] + q.values[1]),
            0.5 * (p.values[2] + q.values[2]))]);
    });

// A lattice point cloud that lies exactly in a plane through the origin with
// the integer normal (a, b, -1): the points are (x, y, a*x + b*y), all of them
// exact. 'hullSize' controls whether the planar hull is a triangle (the
// upstream Newell-loop case) or a larger polygon.
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
            // The extra points are strictly inside the hull, so they do not
            // change it: they are convex combinations with small weights.
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

describe('MinimumVolumeBox3FloatingPoint verification', () => {
    it('contains the points for a nearly degenerate bilinear form', () => {
        // Regression test for the upstream ComputeVolume defect described in
        // src/MinimumVolumeBox3FloatingPoint.ts. For this hull the edge pair
        // (5, 11) has d = f00*f11 - f10*f01 = -1e-7, which is cancellation
        // noise, so the endpoint sample of MinimizerVariableS builds its axis
        // from q0 = -4e-17 and q1 = -7e-18 -- a random direction. Upstream
        // then assumes a vertex of the edge realizes the minimum along that
        // axis, which is false, and reports the box volume 166.5979... whose
        // box misses the point (-1, 0, 5) by 0.957. The exact pipeline reports
        // 180 and the fixed floating-point pipeline agrees with it.
        const points = [
            V(0, 0, 0), V(-2, -4, 1), V(0, -2, 5),
            V(-1, 0, 5), V(0, 3, -5), V(0, 3, 0), V(-4, -1, 4)
        ];
        for (const lgMaxSample of [2, 3, 4]) {
            const result = new MinimumVolumeBox3FloatingPoint()
                .compute(points, lgMaxSample);
            expect(result.dimension).toBe(3);
            expectContainsAll(result.box, points, 1e-9);
            expectClose(result.volume, 180, 1e-9, 1e-12);
        }

        // The same point set with the origin repeated: the duplicate vertices
        // are removed before the box search, so the answer is unchanged.
        const withDuplicates = [
            V(0, 0, 0), V(-2, -4, 1), V(0, 0, 0), V(0, -2, 5), V(0, 0, 0),
            V(-1, 0, 5), V(0, 3, -5), V(0, 3, 0), V(-4, -1, 4)
        ];
        const duplicated = new MinimumVolumeBox3FloatingPoint()
            .compute(withDuplicates, 3);
        expectContainsAll(duplicated.box, withDuplicates, 1e-9);
        expectClose(duplicated.volume, 180, 1e-9, 1e-12);
    });

    it('contains the points for a rotated lattice cloud whose hull has a '
        + 'nearly flat face', () => {
        // Regression test for the upstream getExtreme defect described in
        // src/MinimumVolumeBox3FloatingPoint.ts (issue #426). In the lattice
        // cloud below the point (3, 1, 1) is the midpoint of the segment from
        // (2, -1, 3) to (4, 3, -1) and lies in a flat hull face, so the exact
        // tests of removeCoplanarTriangleAdjacencies drop it from the
        // adjacency graph. Rotating the cloud by 0.1 degrees about the y-axis
        // makes the face flat only to within rounding: the vertex survives in
        // the graph, and for a search direction along that face's normal every
        // neighbor of the vertex ties with it in double precision. The hill
        // climb then stalls on the resulting plateau, the candidate whose
        // second axis is that normal gets the extent 8.5e-17 instead of 3.90,
        // its volume of 5.7e-15 wins the minimization, and the query returns a
        // degenerate box that misses the cloud by 3.90 units. With the plateau
        // traversal the rotated cloud reproduces the unrotated answer to
        // machine precision.
        const points = [
            V(3, 1, 1), V(0, 0, 0), V(0, -1, -3), V(2, -1, 3),
            V(-2, 0, 1), V(0, 0, 0), V(0, -4, 0), V(4, 3, -1)
        ];
        const frame = [
            V(0.9999983994110458, 0, -0.0017891828711411968),
            V(0, 1, 0),
            V(0.0017891828711411968, 0, 0.9999983994110458)
        ];
        const moved = points.map(p => V(
            dot(p, frame[0]), dot(p, frame[1]), dot(p, frame[2])));

        for (const lgMaxSample of [2, 3, 4]) {
            const unrotated = new MinimumVolumeBox3FloatingPoint()
                .compute(points, lgMaxSample);
            const rotated = new MinimumVolumeBox3FloatingPoint()
                .compute(moved, lgMaxSample);
            expect(rotated.dimension).toBe(3);
            expectOrthonormalAxes(rotated.box, 1e-10);
            expectContainsAll(rotated.box, moved, 1e-9);

            // No extent is degenerate: the cloud is 3-dimensional and its
            // smallest containing box has side lengths of at least 3.
            for (const extent of rotated.box.extent.values) {
                expect(extent).toBeGreaterThan(1);
            }

            // A rotation cannot change the volume of the minimum-volume box,
            // and here the candidate search follows the same path on both
            // clouds, so the two volumes agree to rounding.
            expectClose(rotated.volume, unrotated.volume, 1e-9, 1e-12);
        }
    });

    it('returns a box that contains the points and is no larger than the '
        + 'aligned or face-aligned candidates', () => {
        check(latticeCloud(9), points => {
            const result = new MinimumVolumeBox3FloatingPoint()
                .compute(points, 3);
            expect(result.dimension).toBe(3);
            expectOrthonormalAxes(result.box, 1e-10);

            // The scale of the coordinates is at most 5, so an absolute
            // tolerance of 1e-9 is far above the accumulated rounding error of
            // the dot products and well below any real containment failure.
            expectContainsAll(result.box, points, 1e-9);

            // The reported volume is the product of the box side lengths.
            const e = result.box.extent.values;
            expectClose(result.volume, 8 * e[0] * e[1] * e[2], 1e-8, 1e-9);

            // O'Rourke's theorem: the minimum-volume box has a face flush with
            // a hull face, so the search must do at least as well as the best
            // face-aligned candidate box and as the axis-aligned box.
            const faceAligned = bruteForceFaceAlignedVolume(points);
            expect(result.volume).toBeLessThanOrEqual(faceAligned * (1 + 1e-8));
            expect(result.volume)
                .toBeLessThanOrEqual(alignedBoxVolume(points) * (1 + 1e-8));

            // The box contains the hull, so its volume is at least the hull
            // volume.
            expect(result.volume)
                .toBeGreaterThanOrEqual(convexHullVolume(points) * (1 - 1e-8));
        }, 40);
    }, 30000);

    it('is invariant under a rigid motion of the point set', () => {
        check(fc.tuple(latticeCloud(8), rotationFrame(3), latticeVector(3, -9, 9)),
            ([points, frame, shift]) => {
                const moved = points.map(p => V(
                    dot(p, frame[0]) + shift.values[0],
                    dot(p, frame[1]) + shift.values[1],
                    dot(p, frame[2]) + shift.values[2]));
                const mvb = new MinimumVolumeBox3FloatingPoint();
                const v0 = mvb.compute(points, 3).volume;
                const v1 = mvb.compute(moved, 3).volume;
                // The minimizer is iterative and path dependent: rotating the
                // input changes which dyadic sample of each level curve is the
                // best, so the two volumes agree only to the accuracy of the
                // search. With lgMaxSample = 3 discrepancies of ~1.7% were
                // observed on 8-point lattice clouds; a mis-ported candidate
                // or axis formula shows up as O(1) or as a non-containing box
                // (covered by the containment property above).
                expectClose(v1, v0, 1e-6, 5e-2);
            }, 25);
    }, 30000);

    it('contains the points for a deterministic sweep of rotated clouds with '
        + 'a weakly extreme hull vertex', () => {
        // The property below samples the same space, but the plateau stall of
        // issue #426 needs the degenerate candidate to also win the
        // minimization, which happens for roughly 0.3% of these clouds -- too
        // rare for 25 fast-check runs to be a dependable guard. This sweep is
        // a fixed pseudo-random stream of 300 clouds; the pre-fix code returns
        // a box that misses the cloud by 3.77 on the cloud of trial 160.
        const rnd = makeRandom(20260905);
        let tested = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 7; ++i) {
                points.push(V(Math.round(rnd() * 10 - 5),
                    Math.round(rnd() * 10 - 5), Math.round(rnd() * 10 - 5)));
            }
            const ch3 = new ConvexHull3();
            ch3.compute(points);
            if (ch3.getDimension() !== 3) {
                continue;
            }

            // Append the midpoint of an edge of a hull triangle. Halves of
            // these integers are exact in binary64, so the new point lies
            // exactly on the hull and is only weakly extreme.
            const hull = ch3.getHull();
            const t = 3 * Math.floor(rnd() * (hull.length / 3));
            const k = Math.floor(rnd() * 3);
            const p = points[hull[t + k]], q = points[hull[t + (k + 1) % 3]];
            points.push(V(0.5 * (p.values[0] + q.values[0]),
                0.5 * (p.values[1] + q.values[1]),
                0.5 * (p.values[2] + q.values[2])));

            const a = rnd() * 2 * Math.PI, b = rnd() * 2 * Math.PI;
            const c = rnd() * 2 * Math.PI;
            const ca = Math.cos(a), sa = Math.sin(a);
            const cb = Math.cos(b), sb = Math.sin(b);
            const cc = Math.cos(c), sc = Math.sin(c);
            const frame = [
                V(ca * cb, sa * cb, -sb),
                V(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc, cb * sc),
                V(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc, cb * cc)
            ];
            const moved = points.map(v => V(
                dot(v, frame[0]), dot(v, frame[1]), dot(v, frame[2])));

            const result = new MinimumVolumeBox3FloatingPoint()
                .compute(moved, [2, 3, 5][trial % 3]);
            expect(result.dimension).toBe(3);
            expectContainsAll(result.box, moved, 1e-9);
            ++tested;
        }
        expect(tested).toBeGreaterThan(250);
    }, 60000);

    it('contains the points for a rotated cloud whose plateau is flat only to '
        + 'within rounding', () => {
        // A second cloud from the sweep above, kept explicitly because it
        // pins the tolerance of the plateau traversal rather than the
        // traversal itself. The dot products across the flat face differ by a
        // few units in the last place instead of being equal, so a traversal
        // that expanded only through exact ties would still stall: the pre-fix
        // code and a tie-only traversal both report the extent 1.5e-16 along
        // the first axis, a volume of 1.7e-14 and a box that misses the cloud
        // by 6.53. The points are a lattice cloud plus the midpoint of a hull
        // edge, rotated by a pseudo-random frame.
        const points = [
            V(-2.5592634747427425, -1.8167041808335136, 0.3869837027443625),
            V(3.984271911243738, 2.4157801038470605, 1.513137081407649),
            V(-3.2845609605519765, -3.024000308835922, 0.2590008273695822),
            V(0.30841092652503344, 4.5135544318554, 2.5559164875008027),
            V(4.691911868787187, -1.5090305348123525, 3.2724287403300125),
            V(-0.41688655928420043, 3.3062583038529927, 2.4279336121260227),
            V(0.4424730944016727, -1.4907419708617227, -4.072088645528473),
            V(0.01279326755873611, 0.9077581664956351, -0.8220775167012255)
        ];
        const hullVolume = convexHullVolume(points);
        const faceAligned = bruteForceFaceAlignedVolume(points);
        for (const lgMaxSample of [2, 3, 5]) {
            const result = new MinimumVolumeBox3FloatingPoint()
                .compute(points, lgMaxSample);
            expect(result.dimension).toBe(3);
            expectOrthonormalAxes(result.box, 1e-10);
            expectContainsAll(result.box, points, 1e-9);
            expect(result.volume).toBeGreaterThan(hullVolume * (1 - 1e-8));
            expect(result.volume).toBeLessThanOrEqual(faceAligned * (1 + 1e-8));
        }
    });

    it('returns a box that contains the points after a rigid motion of a '
        + 'cloud with a weakly extreme hull vertex', () => {
        // The generator appends the midpoint of two of the points, which is
        // exact in binary64 for these coordinates. When the segment lies on
        // the hull, the midpoint is a hull vertex that is only weakly extreme:
        // it sits in the relative interior of a hull edge or face. Those are
        // the vertices removeCoplanarTriangleAdjacencies drops with its exact
        // colinearity and coplanarity tests, and the vertices whose
        // neighborhood becomes a rounding-level plateau for the hill climb
        // once the rigid motion makes the degeneracy inexact (issue #426).
        check(fc.tuple(cloudWithMidpoint(7), rotationFrame(3),
            latticeVector(3, -9, 9)),
            ([points, frame, shift]) => {
                const moved = points.map(p => V(
                    dot(p, frame[0]) + shift.values[0],
                    dot(p, frame[1]) + shift.values[1],
                    dot(p, frame[2]) + shift.values[2]));
                const result = new MinimumVolumeBox3FloatingPoint()
                    .compute(moved, 3);
                expect(result.dimension).toBe(3);
                expectContainsAll(result.box, moved, 1e-9);
                const e = result.box.extent.values;
                expect(8 * e[0] * e[1] * e[2]).toBeGreaterThan(1 / 6);
            }, 25);
    }, 60000);

    it('returns a box that contains the points after a rigid motion', () => {
        // The companion of the invariance property above and the property
        // that exposed issue #426. A rigid motion turns the exact
        // degeneracies of a lattice cloud -- coplanar hull faces, hull
        // vertices in the interior of a face or of an edge -- into
        // near-degeneracies, which the exact tests of
        // removeCoplanarTriangleAdjacencies no longer detect. Containment is
        // the invariant that a plateau stall in the hill climb breaks: the
        // stalled axis gets a rounding-level extent and the box collapses.
        check(fc.tuple(latticeCloud(8), rotationFrame(3), latticeVector(3, -9, 9)),
            ([points, frame, shift]) => {
                const moved = points.map(p => V(
                    dot(p, frame[0]) + shift.values[0],
                    dot(p, frame[1]) + shift.values[1],
                    dot(p, frame[2]) + shift.values[2]));
                const result = new MinimumVolumeBox3FloatingPoint()
                    .compute(moved, 3);
                expect(result.dimension).toBe(3);
                expectOrthonormalAxes(result.box, 1e-10);
                // The coordinates are bounded by 5 + 9 = 14 in magnitude, so
                // 1e-9 is far above the rounding error of the dot products.
                expectContainsAll(result.box, moved, 1e-9);
                // The box is not degenerate: it contains the hull, whose
                // volume is a positive integer multiple of 1/6 and therefore
                // at least 1/6.
                const e = result.box.extent.values;
                expect(8 * e[0] * e[1] * e[2]).toBeGreaterThan(1 / 6);
            }, 25);
    }, 60000);

    it('scales the volume by the cube of a uniform scale', () => {
        check(fc.tuple(latticeCloud(8), fc.integer({ min: 2, max: 5 })),
            ([points, s]) => {
                const mvb = new MinimumVolumeBox3FloatingPoint();
                const v0 = mvb.compute(points, 3).volume;
                const v1 = mvb.compute(points.map(
                    p => V(s * p.values[0], s * p.values[1], s * p.values[2])),
                3).volume;
                // A power-of-two scale would be exact; a general integer scale
                // only rescales the same floating-point search, so allow a
                // relative tolerance.
                expectClose(v1, s * s * s * v0, 1e-9, 1e-9);
            }, 25);
    }, 30000);

    it('does not get worse when the sample count increases', () => {
        check(latticeCloud(8), points => {
            const mvb = new MinimumVolumeBox3FloatingPoint();
            const coarse = mvb.compute(points, 2).volume;
            const fine = mvb.compute(points, 3).volume;
            // The dyadic sample set for lgMaxSample = 2 is a subset of the one
            // for lgMaxSample = 3, so the minimum cannot increase except by
            // rounding.
            expect(fine).toBeLessThanOrEqual(coarse * (1 + 1e-9));
        }, 25);
    }, 30000);

    it('classifies coplanar point sets as dimension 2 with a finite basis', () => {
        // The triangular-hull case is the one for which the upstream Newell
        // loop sums to the zero vector (issue: the loop starts at i1 = 1, so
        // it computes Cross(P2,P1) + Cross(P1,P2) = 0 for a 3-gon). A zero
        // normal makes ComputeOrthogonalComplement produce a degenerate basis
        // and every output value NaN. The port closes the loop.
        check(fc.oneof(coplanarLatticeCloud('triangle'),
            coplanarLatticeCloud('polygon')), points => {
            const result = new MinimumVolumeBox3FloatingPoint()
                .compute(points, 3);
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
        }, 60);
    }, 30000);

    it('computes a nonzero Newell normal for a triangular coplanar hull', () => {
        // The direct statement of the fixed defect: the third box axis is the
        // plane normal, so it must be a unit vector parallel to (a, b, -1).
        check(fc.tuple(fc.integer({ min: -3, max: 3 }),
            fc.integer({ min: -3, max: 3 })), ([a, b]) => {
            const lift = (x: number, y: number): Vector => V(x, y, a * x + b * y);
            const points = [lift(0, 0), lift(4, 0), lift(0, 4), lift(1, 1)];
            const result = new MinimumVolumeBox3FloatingPoint()
                .compute(points, 3);
            expect(result.dimension).toBe(2);
            const n = result.box.axis[2];
            expectClose(dot(n, n), 1, 1e-12, 1e-12);
            // The plane is z = a*x + b*y, whose normal is parallel to
            // (a, b, -1).
            const exact = V(a, b, -1);
            normalize(exact);
            expectClose(Math.abs(dot(n, exact)), 1, 1e-12, 1e-12);
        }, 49);
    }, 30000);

    it('classifies collinear and coincident point sets', () => {
        check(fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -4, 4),
            fc.array(fc.integer({ min: -6, max: 6 }),
                { minLength: 4, maxLength: 7 })),
        ([origin, direction, ts]) => {
            const d = direction.values;
            if (d[0] === 0 && d[1] === 0 && d[2] === 0) {
                // Coincident points.
                const points = ts.map(() => origin.clone());
                const result = new MinimumVolumeBox3FloatingPoint()
                    .compute(points, 3);
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
            const result = new MinimumVolumeBox3FloatingPoint()
                .compute(points, 3);
            expect(result.dimension).toBe(distinct.size === 1 ? 0 : 1);
            expect(result.volume).toBe(0);
            expect(result.box.extent.values[1]).toBe(0);
            expect(result.box.extent.values[2]).toBe(0);
            expectOrthonormalAxes(result.box, 1e-12);
            expectContainsAll(result.box, points, 1e-9);
        }, 60);
    }, 30000);
});
