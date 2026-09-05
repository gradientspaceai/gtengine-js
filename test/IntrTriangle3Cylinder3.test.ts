import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import { Line } from '../src/Line.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import {
    IntrTriangle3Cylinder3TI
} from '../src/IntrTriangle3Cylinder3.js';
import { cross } from '../src/Vector3.js';
import {
    check, fc, positive, rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function tri(a: Vector, b: Vector, c: Vector): Triangle {
    return Triangle.fromVertices(a, b, c);
}

function cylinder(c: Vector, d: Vector, r: number, h: number): Cylinder3 {
    const dir = d.clone();
    normalize(dir);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(c, dir), r, h);
}

// The standard test cylinder: center at the origin, axis (0,0,1), radius 1,
// height 2 (so the solid is x^2 + y^2 <= 1 and |z| <= 1).
function unitCylinder(): Cylinder3 {
    return cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
}

function test(t: Triangle, c: Cylinder3): boolean {
    return new IntrTriangle3Cylinder3TI().test(t, c).intersect;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// An independent exact distance from a point to the SOLID finite cylinder.
// In cylinder coordinates the distance is
//   sqrt(max(rho - r, 0)^2 + max(|z| - h/2, 0)^2)
// where rho is the distance from the axis.
function distanceToSolidCylinder(p: Vector, cyl: Cylinder3): number {
    const diff = sub(p, cyl.axis.origin);
    const z = dot(diff, cyl.axis.direction);
    const radial = sub(diff, mul(cyl.axis.direction, z));
    const rho = Math.sqrt(dot(radial, radial));
    const dr = Math.max(rho - cyl.radius, 0);
    const dz = Math.max(Math.abs(z) - 0.5 * cyl.height, 0);
    return Math.sqrt(dr * dr + dz * dz);
}

// The minimum of distanceToSolidCylinder over a barycentric grid on the
// triangle. Because the distance function is 1-Lipschitz, the true minimum
// distance is at least (this value - the grid resolution).
function gridMinDistance(t: Triangle, cyl: Cylinder3, n: number): number {
    let minDistance = Number.MAX_VALUE;
    const e1 = sub(t.v[1], t.v[0]);
    const e2 = sub(t.v[2], t.v[0]);
    for (let i = 0; i <= n; ++i) {
        for (let j = 0; i + j <= n; ++j) {
            const p = add(t.v[0],
                add(mul(e1, i / n), mul(e2, j / n)));
            const d = distanceToSolidCylinder(p, cyl);
            if (d < minDistance) {
                minDistance = d;
            }
        }
    }
    return minDistance;
}

// The maximum distance between adjacent grid samples, an upper bound on how
// much gridMinDistance can exceed the true minimum distance.
function gridResolution(t: Triangle, n: number): number {
    const e1 = sub(t.v[1], t.v[0]);
    const e2 = sub(t.v[2], t.v[0]);
    const e3 = sub(t.v[2], t.v[1]);
    const lengths = [e1, e2, e3].map(e => Math.sqrt(dot(e, e)));
    return Math.max(...lengths) / n;
}

describe('IntrTriangle3Cylinder3TI', () => {
    it('detects a triangle that contains the cylinder axis', () => {
        // The triangle spans the axis in the plane y = 0.
        const t = tri(v3(-3, 0, 0), v3(3, 0, 0), v3(0, 0, 3));
        expect(test(t, unitCylinder())).toBe(true);
    });

    it('handles a triangle in the plane of the top disk', () => {
        const cyl = unitCylinder();

        // The triangle contains the axis point (0,0,1) of the top disk.
        const inside = tri(v3(-2, -2, 1), v3(2, -2, 1), v3(0, 2, 1));
        expect(test(inside, cyl)).toBe(true);

        // The triangle is in the top-disk plane and touches the disk
        // boundary at (1,0,1) only.
        const touching = tri(v3(1, 0, 1), v3(3, -1, 1), v3(3, 1, 1));
        expect(test(touching, cyl)).toBe(true);

        // The triangle is in the top-disk plane but misses the disk.
        const outside = tri(v3(1.25, 0, 1), v3(3, -1, 1), v3(3, 1, 1));
        expect(test(outside, cyl)).toBe(false);
    });

    it('rejects triangles strictly outside the cylinder slab', () => {
        const cyl = unitCylinder();
        const above = tri(v3(-1, -1, 1.5), v3(1, -1, 2), v3(0, 1, 3));
        expect(test(above, cyl)).toBe(false);
        const below = tri(v3(-1, -1, -1.5), v3(1, -1, -2), v3(0, 1, -3));
        expect(test(below, cyl)).toBe(false);
    });

    it('handles a triangle parallel to the axis and tangent to the wall',
        () => {
            const cyl = unitCylinder();

            // The plane x = 1 is tangent to the cylinder wall at (1,0,z).
            const tangent = tri(v3(1, -2, 0), v3(1, 2, 0), v3(1, 0, 0.5));
            expect(test(tangent, cyl)).toBe(true);

            // The plane x = 1.001 misses the cylinder.
            const separated = tri(v3(1.001, -2, 0), v3(1.001, 2, 0),
                v3(1.001, 0, 0.5));
            expect(test(separated, cyl)).toBe(false);

            // A plane inside the wall cuts the cylinder.
            const cutting = tri(v3(0.5, -2, 0), v3(0.5, 2, 0),
                v3(0.5, 0, 0.5));
            expect(test(cutting, cyl)).toBe(true);
        });

    it('detects vertex-only and edge-only contact', () => {
        const cyl = unitCylinder();

        // Only the vertex (1,0,0) touches the cylinder wall.
        const vertexContact = tri(v3(1, 0, 0), v3(3, 0, 0), v3(3, 2, 0));
        expect(test(vertexContact, cyl)).toBe(true);

        // Moving the vertex outward separates the triangle.
        const vertexMissed = tri(v3(1.0001, 0, 0), v3(3, 0, 0), v3(3, 2, 0));
        expect(test(vertexMissed, cyl)).toBe(false);

        // The whole edge x = 1, z = 0 is tangent to the disk at (1,0,0).
        const edgeContact = tri(v3(1, -1, 0), v3(1, 1, 0), v3(3, 0, 0));
        expect(test(edgeContact, cyl)).toBe(true);

        // The triangle touches only at the top rim point (1,0,1).
        const rimContact = tri(v3(1, 0, 1), v3(3, 0, 3), v3(3, 2, 3));
        expect(test(rimContact, cyl)).toBe(true);
    });

    it('detects a triangle fully inside the cylinder', () => {
        const t = tri(v3(0, 0, 0), v3(0.5, 0, 0.2), v3(0, 0.5, -0.3));
        expect(test(t, unitCylinder())).toBe(true);
    });

    it('rejects a clearly disjoint triangle', () => {
        const t = tri(v3(10, 10, 10), v3(12, 10, 10), v3(10, 12, 11));
        expect(test(t, unitCylinder())).toBe(false);
    });

    it('handles triangles that must be clipped by both disk planes', () => {
        const cyl = unitCylinder();

        // The triangle straddles both planes and passes through the axis.
        const through = tri(v3(0, 0, -5), v3(0.2, 0, 5), v3(0, 0.2, 0));
        expect(test(through, cyl)).toBe(true);

        // The triangle straddles both planes but is far from the axis, so
        // the clipped polygon does not overlap the disk.
        const missing = tri(v3(5, 0, -5), v3(5, 0, 5), v3(5, 5, 0));
        expect(test(missing, cyl)).toBe(false);

        // A single vertex is below the slab, the other two are inside it
        // (case 4e), with the clipped quadrilateral covering the axis.
        const oneBelow = tri(v3(0, 0, -3), v3(0.5, 0, 0.5),
            v3(-0.5, 0.2, 0.5));
        expect(test(oneBelow, cyl)).toBe(true);

        // A single vertex is above the slab (cases 3d/3e).
        const oneAbove = tri(v3(0, 0, 0.5), v3(0.5, 0, 3), v3(-0.5, 0.2, 3));
        expect(test(oneAbove, cyl)).toBe(true);
    });

    it('handles the touching-plane special cases (1a, 1b, 2a, 2b)', () => {
        const cyl = unitCylinder();

        // Case 1a: exactly one vertex lies in the bottom-disk plane and the
        // other two are strictly below it.
        const case1aHit = tri(v3(0.5, 0, -1), v3(0, 0, -4), v3(1, 1, -4));
        expect(test(case1aHit, cyl)).toBe(true);
        const case1aMiss = tri(v3(2, 0, -1), v3(0, 0, -4), v3(1, 1, -4));
        expect(test(case1aMiss, cyl)).toBe(false);

        // Case 2a: two vertices lie in the bottom-disk plane.
        const case2aHit = tri(v3(-2, 0, -1), v3(2, 0, -1), v3(0, 0, -4));
        expect(test(case2aHit, cyl)).toBe(true);
        const case2aMiss = tri(v3(-2, 3, -1), v3(2, 3, -1), v3(0, 0, -4));
        expect(test(case2aMiss, cyl)).toBe(false);

        // Case 1b: exactly one vertex lies in the top-disk plane.
        const case1bHit = tri(v3(0.5, 0, 1), v3(0, 0, 4), v3(1, 1, 4));
        expect(test(case1bHit, cyl)).toBe(true);
        const case1bMiss = tri(v3(2, 0, 1), v3(0, 0, 4), v3(1, 1, 4));
        expect(test(case1bMiss, cyl)).toBe(false);

        // Case 2b: two vertices lie in the top-disk plane.
        const case2bHit = tri(v3(-2, 0, 1), v3(2, 0, 1), v3(0, 0, 4));
        expect(test(case2bHit, cyl)).toBe(true);
        const case2bMiss = tri(v3(-2, 3, 1), v3(2, 3, 1), v3(0, 0, 4));
        expect(test(case2bMiss, cyl)).toBe(false);
    });

    it('is invariant to the order of the triangle vertices', () => {
        const cyl = unitCylinder();
        const a = v3(0.3, -2, -3);
        const b = v3(-1.5, 0.4, 4);
        const c = v3(2.5, 1.5, 0.25);
        const expected = test(tri(a, b, c), cyl);
        expect(test(tri(a, c, b), cyl)).toBe(expected);
        expect(test(tri(b, a, c), cyl)).toBe(expected);
        expect(test(tri(b, c, a), cyl)).toBe(expected);
        expect(test(tri(c, a, b), cyl)).toBe(expected);
        expect(test(tri(c, b, a), cyl)).toBe(expected);
        expect(expected).toBe(true);
    });

    it('handles degenerate triangles', () => {
        const cyl = unitCylinder();

        // A triangle degenerated to a single point.
        expect(test(tri(v3(0.2, 0.2, 0.2), v3(0.2, 0.2, 0.2),
            v3(0.2, 0.2, 0.2)), cyl)).toBe(true);
        expect(test(tri(v3(5, 5, 5), v3(5, 5, 5), v3(5, 5, 5)),
            cyl)).toBe(false);

        // A triangle degenerated to a segment that pierces the cylinder.
        expect(test(tri(v3(-3, 0, 0), v3(3, 0, 0), v3(-3, 0, 0)),
            cyl)).toBe(true);

        // A triangle degenerated to a segment that misses the cylinder.
        expect(test(tri(v3(-3, 2, 0), v3(3, 2, 0), v3(-3, 2, 0)),
            cyl)).toBe(false);

        // A degenerate segment PARALLEL to the axis and far from it. The
        // upstream containment test in DiskOverlapsPolygon reports that a
        // polygon degenerated to a single point contains the origin (all of
        // the DotPerp values are zero, so both the positive and negative
        // counts are zero), so upstream reports a false positive here. The
        // port instead falls through to the edge tests, which classify the
        // degenerate polygon correctly. See the V32 verification block below.
        expect(test(tri(v3(4, 0, -4), v3(4, 0, 4), v3(4, 0, -4)),
            cyl)).toBe(false);

        // A degenerate segment crossing the slab far from the axis but not
        // parallel to it (the projection is a segment, not a point).
        expect(test(tri(v3(4, 0, -4), v3(4, 2, 4), v3(4, 0, -4)),
            cyl)).toBe(false);
    });

    it('handles a zero-radius cylinder (a segment)', () => {
        const cyl = cylinder(v3(0, 0, 0), v3(0, 0, 1), 0, 2);

        // The triangle contains the axis segment.
        expect(test(tri(v3(-1, 0, -2), v3(1, 0, -2), v3(0, 0, 2)),
            cyl)).toBe(true);

        // The triangle crosses the axis line inside the slab.
        expect(test(tri(v3(-1, -1, 0), v3(1, -1, 0), v3(0, 1, 0)),
            cyl)).toBe(true);

        // The triangle is offset from the axis.
        expect(test(tri(v3(1, -1, 0), v3(3, -1, 0), v3(2, 1, 0)),
            cyl)).toBe(false);
    });

    it('handles a zero-height cylinder (a disk)', () => {
        const cyl = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 0);

        // The triangle crosses the plane z = 0 inside the disk.
        expect(test(tri(v3(0, 0, -1), v3(0.5, 0, 1), v3(-0.5, 0.2, 1)),
            cyl)).toBe(true);

        // The triangle crosses z = 0 outside the disk.
        expect(test(tri(v3(3, 0, -1), v3(3.5, 0, 1), v3(2.5, 0.2, 1)),
            cyl)).toBe(false);

        // The triangle lies in the disk plane and overlaps the disk.
        expect(test(tri(v3(-2, -2, 0), v3(2, -2, 0), v3(0, 2, 0)),
            cyl)).toBe(true);

        // The triangle does not reach the plane z = 0.
        expect(test(tri(v3(0, 0, 1), v3(1, 0, 2), v3(0, 1, 2)),
            cyl)).toBe(false);
    });

    it('rejects infinite cylinders and non-3D triangles', () => {
        const infinite = unitCylinder();
        infinite.makeInfiniteCylinder();
        const t = tri(v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0));
        expect(() => test(t, infinite)).toThrow(
            'Infinite cylinders are not yet supported.');

        const t2 = Triangle.fromVertices(Vector.fromArray([0, 0]),
            Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
        expect(() => test(t2, unitCylinder())).toThrow('mismatched sizes');
    });

    it('agrees with a dense sampling of the triangle (randomized)', () => {
        const rand = makeRandom(20260901);
        const n = 40;
        let numCases = 0;
        let numDecided = 0;
        let numFalseNegatives = 0;
        let numFalsePositives = 0;

        for (let trial = 0; trial < 200; ++trial) {
            const t = tri(
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2),
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2),
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2));
            const cyl = cylinder(
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1),
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1),
                0.2 + rand(), 0.4 + 2 * rand());
            ++numCases;

            const intersect = test(t, cyl);
            const gridMin = gridMinDistance(t, cyl, n);
            const resolution = gridResolution(t, n);

            if (gridMin === 0) {
                // A sampled triangle point is in the solid cylinder, so the
                // query must report an intersection.
                ++numDecided;
                if (!intersect) {
                    ++numFalseNegatives;
                }
            } else if (gridMin > 2 * resolution) {
                // The true distance between the triangle and the cylinder is
                // positive, so the query must report no intersection.
                ++numDecided;
                if (intersect) {
                    ++numFalsePositives;
                }
            }
        }

        expect(numCases).toBe(200);
        expect(numFalseNegatives).toBe(0);
        expect(numFalsePositives).toBe(0);
        expect(numDecided).toBeGreaterThan(150);
    });

    it('agrees with a dense sampling of the cylinder (randomized)', () => {
        // Sample points of the solid cylinder. For each sample, the closest
        // point of the solid triangle is computed with DistPointTriangle. If
        // that triangle point is in the solid cylinder, the triangle and the
        // cylinder must intersect.
        const rand = makeRandom(987654321);
        const dpt = new DistPointTriangle();
        let numProved = 0;
        let numMismatches = 0;

        for (let trial = 0; trial < 30; ++trial) {
            const t = tri(
                v3(3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5),
                v3(3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5),
                v3(3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5));
            const cyl = cylinder(
                v3(rand() - 0.5, rand() - 0.5, rand() - 0.5),
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1),
                0.3 + rand(), 0.5 + 2 * rand());
            const intersect = test(t, cyl);

            // Build the cylinder frame and sample the solid cylinder.
            const basis: Vector[] = [cyl.axis.direction.clone(),
                new Vector(3), new Vector(3)];
            computeOrthogonalComplement3(1, basis);
            let proved = false;
            for (let iz = 0; iz <= 6 && !proved; ++iz) {
                const zc = cyl.height * (iz / 6 - 0.5);
                for (let ir = 0; ir <= 4 && !proved; ++ir) {
                    const rho = cyl.radius * ir / 4;
                    for (let ia = 0; ia < 16 && !proved; ++ia) {
                        const angle = 2 * Math.PI * ia / 16;
                        const p = add(cyl.axis.origin,
                            add(mul(basis[0], zc),
                                add(mul(basis[1], rho * Math.cos(angle)),
                                    mul(basis[2],
                                        rho * Math.sin(angle)))));
                        const closest = dpt.compute(p, t).closest[1];
                        if (distanceToSolidCylinder(closest, cyl) === 0) {
                            proved = true;
                        }
                    }
                }
            }

            if (proved) {
                ++numProved;
                if (!intersect) {
                    ++numMismatches;
                }
            }
        }

        expect(numMismatches).toBe(0);
        expect(numProved).toBeGreaterThan(5);
    });

    it('is invariant under rigid motions (randomized)', () => {
        const rand = makeRandom(13579);
        let numMismatches = 0;
        let numIntersecting = 0;

        for (let trial = 0; trial < 500; ++trial) {
            const t = tri(
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2),
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2),
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2));
            const cyl = cylinder(
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1),
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1),
                0.2 + rand(), 0.4 + 2 * rand());
            const expected = test(t, cyl);
            if (expected) {
                ++numIntersecting;
            }

            // A random rotation (a right-handed orthonormal basis) and a
            // random translation.
            const frame: Vector[] = [
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1),
                new Vector(3), new Vector(3)
            ];
            normalize(frame[0]);
            computeOrthogonalComplement3(1, frame);
            const translation = v3(6 * rand() - 3, 6 * rand() - 3,
                6 * rand() - 3);
            const transform = (p: Vector): Vector => add(translation,
                add(mul(frame[0], p.values[0]),
                    add(mul(frame[1], p.values[1]),
                        mul(frame[2], p.values[2]))));
            const rotate = (p: Vector): Vector =>
                add(mul(frame[0], p.values[0]),
                    add(mul(frame[1], p.values[1]),
                        mul(frame[2], p.values[2])));

            const tXform = tri(transform(t.v[0]), transform(t.v[1]),
                transform(t.v[2]));
            const cylXform = cylinder(transform(cyl.axis.origin),
                rotate(cyl.axis.direction), cyl.radius, cyl.height);

            // Rigid motions can flip the answer only for configurations
            // within roundoff of tangency; require agreement away from the
            // boundary, which the sampled configurations are.
            if (test(tXform, cylXform) !== expected) {
                ++numMismatches;
            }
        }

        expect(numIntersecting).toBeGreaterThan(50);
        expect(numMismatches).toBe(0);
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrTriangle3Cylinder3.h.
// ---------------------------------------------------------------------------

// The signed cylinder function: negative strictly inside the solid, zero on
// the boundary, positive outside.
function signedCylinder(p: Vector, cyl: Cylinder3): number {
    const diff = sub(p, cyl.axis.origin);
    const z = dot(diff, cyl.axis.direction);
    const radial = sub(diff, mul(cyl.axis.direction, z));
    const rho = Math.sqrt(dot(radial, radial));
    return Math.max(rho - cyl.radius, Math.abs(z) - 0.5 * cyl.height);
}

function gridMinSigned(t: Triangle, cyl: Cylinder3, n: number): number {
    let best = Number.MAX_VALUE;
    const e1 = sub(t.v[1], t.v[0]);
    const e2 = sub(t.v[2], t.v[0]);
    for (let i = 0; i <= n; ++i) {
        for (let j = 0; i + j <= n; ++j) {
            const s = signedCylinder(
                add(t.v[0], add(mul(e1, i / n), mul(e2, j / n))), cyl);
            if (s < best) {
                best = s;
            }
        }
    }
    return best;
}

describe('IntrTriangle3Cylinder3 verification', () => {
    const ti = new IntrTriangle3Cylinder3TI();

    const arbTriangle = fc.tuple(wellScaledVector(3, -3, 3),
        wellScaledVector(3, -3, 3), wellScaledVector(3, -3, 3))
        .filter(([a, b, c]) => {
            const e0 = sub(b, a), e1 = sub(c, a);
            const n = cross(e0, e1);
            return dot(n, n) > 0.25;
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));
    const arbCylinder = fc.tuple(wellScaledVector(3, -2, 2), unitVector(3),
        positive(2, 0.25), positive(4, 0.5))
        .map(([o, d, r, h]) => Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(o, d), r, h));
    const arbPair = fc.tuple(arbTriangle, arbCylinder);

    it('reports an intersection whenever a sampled triangle point is well'
        + ' inside the solid cylinder', () => {
        check(arbPair, ([t, cyl]) => {
            if (gridMinSigned(t, cyl, 24) < -1e-6) {
                expect(ti.test(t, cyl).intersect).toBe(true);
            }
        }, 100);
    });

    // The containment test of DiskOverlapsPolygon is a set of sign tests on
    // DotPerp values. When the plane of the triangle nearly contains the
    // cylinder axis, the projection of the triangle onto the plane
    // perpendicular to the axis is a nearly degenerate polygon, those values
    // are pure round-off, and the sign tests are unreliable: if they all come
    // out nonnegative the query reports that the polygon contains the axis
    // and returns a false positive. The exactly degenerate case is handled
    // (see the regression test below), but the nearly degenerate one is a
    // conditioning limitation of the upstream algorithm, so the properties
    // that predict a 'false' answer exclude it.
    function axisTransversal(t: Triangle, cyl: Cylinder3): boolean {
        const n = cross(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0]));
        const len = Math.sqrt(dot(n, n));
        return len > 0
            && Math.abs(dot(n, cyl.axis.direction)) / len > 1e-2;
    }

    it('reports no intersection when the whole triangle is provably outside',
        () => {
            // distanceToSolidCylinder is 1-Lipschitz, so the true minimum over
            // the triangle is at least gridMinDistance - gridResolution.
            check(arbPair, ([t, cyl]) => {
                const n = 24;
                if (axisTransversal(t, cyl)
                    && gridMinDistance(t, cyl, n) - gridResolution(t, n)
                        > 1e-6) {
                    expect(ti.test(t, cyl).intersect).toBe(false);
                }
            }, 100);
        });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(arbPair, rotationFrame(3), wellScaledVector(3, -3, 3)),
            ([[t, cyl], frame, shift]) => {
                const rot = (v: Vector): Vector => Vector.fromArray([
                    dot(frame[0], v), dot(frame[1], v), dot(frame[2], v)]);
                const map = (v: Vector): Vector => add(rot(v), shift);
                const t2 = Triangle.fromVertices(map(t.v[0]), map(t.v[1]),
                    map(t.v[2]));
                const cyl2 = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(map(cyl.axis.origin),
                        rot(cyl.axis.direction)), cyl.radius, cyl.height);
                // The query classifies with exact comparisons, so a tangential
                // configuration can flip. Assert only when the answer is
                // robust: a sampled point is well inside, or the triangle is
                // provably outside.
                const n = 24;
                if (gridMinSigned(t, cyl, n) < -1e-4) {
                    expect(ti.test(t2, cyl2).intersect).toBe(true);
                }
                else if (axisTransversal(t, cyl)
                    && gridMinDistance(t, cyl, n)
                        - gridResolution(t, n) > 1e-4) {
                    expect(ti.test(t2, cyl2).intersect).toBe(false);
                }
            }, 60);
    });

    it('is invariant under permuting the triangle vertices', () => {
        // The query sorts the vertices by their axial coordinate, so the
        // answer must not depend on the input order.
        check(arbPair, ([t, cyl]) => {
            const expected = ti.test(t, cyl).intersect;
            const orders = [[0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1],
                [2, 1, 0]];
            for (const o of orders) {
                const p = Triangle.fromVertices(t.v[o[0]], t.v[o[1]],
                    t.v[o[2]]);
                expect(ti.test(p, cyl).intersect).toBe(expected);
            }
        });
    });

    it('rejects infinite cylinders', () => {
        const infinite = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(v3(0, 0, 0), v3(0, 0, 1)), 1, -1);
        expect(() => ti.test(tri(v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)),
            infinite)).toThrow('Infinite cylinders are not yet supported.');
    });

    it('does not report a degenerate projected polygon as containing the'
        + ' axis (upstream-bug regression)', () => {
        const cyl = unitCylinder();   // axis (0,0,1), radius 1, |z| <= 1

        // A triangle degenerated to a segment parallel to the axis. Its
        // projection is a single point 10 units from the axis, so every
        // DotPerp of DiskOverlapsPolygon is zero; upstream reports that the
        // polygon contains the origin.
        const degenerate = tri(v3(10, 0, -0.5), v3(10, 0, 0), v3(10, 0, 0.5));
        expect(ti.test(degenerate, cyl).intersect).toBe(false);
        expect(distanceToSolidCylinder(v3(10, 0, 0), cyl)).toBeGreaterThan(8);

        // A nondegenerate triangle in a plane that contains the cylinder
        // axis. Its projection is a segment on a line through the origin, so
        // again every DotPerp is zero. The triangle is 3 units from the axis.
        const coplanar = tri(v3(0, 3, 0), v3(0, 6, 0), v3(0, 4, 0.5));
        expect(ti.test(coplanar, cyl).intersect).toBe(false);
        expect(gridMinDistance(coplanar, cyl, 16)).toBeGreaterThan(1.5);

        // The same configuration moved onto the cylinder still intersects.
        const touching = tri(v3(0, 0.5, 0), v3(0, 6, 0), v3(0, 4, 0.5));
        expect(ti.test(touching, cyl).intersect).toBe(true);

        // A point-degenerate projection inside the disk also still
        // intersects.
        const inside = tri(v3(0.5, 0, -0.5), v3(0.5, 0, 0), v3(0.5, 0, 0.5));
        expect(ti.test(inside, cyl).intersect).toBe(true);
    });

    it('handles the exact slab boundary cases', () => {
        const cyl = unitCylinder();   // |z| <= 1, x^2 + y^2 <= 1
        // A triangle entirely in the plane z = 1 (case 1b/2b of the PDF) that
        // overlaps the disk.
        const onTop = tri(v3(0, 0, 1), v3(2, 0, 1), v3(0, 2, 1));
        expect(ti.test(onTop, cyl).intersect).toBe(true);
        // The same triangle translated so it only touches the rim.
        const rim = tri(v3(1, 0, 1), v3(3, 0, 1), v3(1, 2, 1));
        expect(ti.test(rim, cyl).intersect).toBe(true);
        // Strictly above the slab.
        const above = tri(v3(0, 0, 1.5), v3(2, 0, 1.5), v3(0, 2, 1.5));
        expect(ti.test(above, cyl).intersect).toBe(false);
        // Strictly below the slab.
        const below = tri(v3(0, 0, -1.5), v3(2, 0, -1.5), v3(0, 2, -1.5));
        expect(ti.test(below, cyl).intersect).toBe(false);
        // In the plane z = 1 but outside the disk.
        const offTop = tri(v3(2, 0, 1), v3(4, 0, 1), v3(2, 2, 1));
        expect(ti.test(offTop, cyl).intersect).toBe(false);
    });
});
