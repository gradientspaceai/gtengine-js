import { describe, it, expect } from 'vitest';
import { Vector, dot, normalize } from '../src/Vector';
import { cross, computeOrthogonalComplement3, dotCross } from '../src/Vector3';
import { type OrientedBox3 } from '../src/OrientedBox';
import { ConvexHull3 } from '../src/ConvexHull3';
import { MinimumAreaBox2 } from '../src/MinimumAreaBox2';
import {
    MinimumVolumeBox3FloatingPoint,
    type MinimumVolumeBox3FloatingPointResult
} from '../src/MinimumVolumeBox3FloatingPoint';

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
