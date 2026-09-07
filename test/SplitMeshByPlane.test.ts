import { describe, expect, it } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import { SplitMeshByPlane } from '../src/SplitMeshByPlane.js';
import type { SplitMeshByPlaneResult } from '../src/SplitMeshByPlane.js';
import { Vector, dot, length, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, unitVector
} from './helpers/arbitraries.js';

// A small deterministic pseudorandom generator (mulberry32) so the randomized
// cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function coords(v: Vector): number[] {
    return [v.get(0), v.get(1), v.get(2)];
}

function signedDistance(plane: Hyperplane, p: Vector): number {
    return dot(plane.normal, p) - plane.constant;
}

function triangleArea(p0: Vector, p1: Vector, p2: Vector): number {
    return 0.5 * length(cross(sub(p1, p0), sub(p2, p0)));
}

function meshArea(vertices: readonly Vector[],
    indices: readonly number[]): number {
    let area = 0;
    for (let i = 0; i < indices.length; i += 3) {
        area += triangleArea(vertices[indices[i] as number] as Vector,
            vertices[indices[i + 1] as number] as Vector,
            vertices[indices[i + 2] as number] as Vector);
    }
    return area;
}

// Verify that no triangle of 'indices' has a vertex on the far side of the
// plane, where 'sign' is +1 for the positive output and -1 for the negative
// output. Also verify that no output triangle lies entirely in the plane.
function verifySide(result: SplitMeshByPlaneResult, plane: Hyperplane,
    indices: readonly number[], sign: number, epsilon: number): void {
    expect(indices.length % 3).toBe(0);
    for (let i = 0; i < indices.length; i += 3) {
        let numOffPlane = 0;
        for (let j = 0; j < 3; ++j) {
            const p = result.clipVertices[indices[i + j] as number] as Vector;
            const sDist = sign * signedDistance(plane, p);
            expect(sDist).toBeGreaterThan(-epsilon);
            if (sDist > epsilon) {
                ++numOffPlane;
            }
        }
        expect(numOffPlane).toBeGreaterThan(0);
    }
}

// The directed edges whose two endpoints both lie in the plane.
function onPlaneDirectedEdges(result: SplitMeshByPlaneResult,
    plane: Hyperplane, indices: readonly number[],
    epsilon: number): string[] {
    const onPlane = (index: number): boolean => {
        const p = result.clipVertices[index] as Vector;
        return Math.abs(signedDistance(plane, p)) <= epsilon;
    };
    const directed: string[] = [];
    for (let i = 0; i < indices.length; i += 3) {
        for (let j = 0; j < 3; ++j) {
            const i0 = indices[i + j] as number;
            const i1 = indices[i + (j + 1) % 3] as number;
            if (onPlane(i0) && onPlane(i1)) {
                directed.push(i0 + '->' + i1);
            }
        }
    }
    directed.sort();
    return directed;
}

function reversedEdges(directed: readonly string[]): string[] {
    const reversed = directed.map(e => {
        const parts = e.split('->');
        return parts[1] + '->' + parts[0];
    });
    reversed.sort();
    return reversed;
}

// The unit cube [-1,1]^3 as a triangle mesh with outward-facing triangles.
function makeCube(): { vertices: Vector[], indices: number[] } {
    const vertices = [
        v3(-1, -1, -1), v3(1, -1, -1), v3(1, 1, -1), v3(-1, 1, -1),
        v3(-1, -1, 1), v3(1, -1, 1), v3(1, 1, 1), v3(-1, 1, 1)
    ];
    const indices = [
        0, 3, 2, 0, 2, 1,   // z = -1
        4, 5, 6, 4, 6, 7,   // z = +1
        0, 1, 5, 0, 5, 4,   // y = -1
        1, 2, 6, 1, 6, 5,   // x = +1
        2, 3, 7, 2, 7, 6,   // y = +1
        3, 0, 4, 3, 4, 7    // x = -1
    ];
    return { vertices, indices };
}

describe('SplitMeshByPlane', () => {
    const plusZ = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);

    it('rejects an index count that is not a multiple of 3', () => {
        const split = new SplitMeshByPlane();
        expect(() => split.compute([v3(0, 0, 1), v3(1, 0, 1)], [0, 1], plusZ))
            .toThrowError('The number of indices must be a multiple of 3.');
    });

    it('keeps an unsplit triangle on the correct side', () => {
        const split = new SplitMeshByPlane();

        // All positive.
        let result = split.compute(
            [v3(0, 0, 1), v3(1, 0, 2), v3(0, 1, 3)], [0, 1, 2], plusZ);
        expect(result.posIndices).toEqual([0, 1, 2]);
        expect(result.negIndices).toEqual([]);
        expect(result.clipVertices.length).toBe(3);

        // All negative.
        result = split.compute(
            [v3(0, 0, -1), v3(1, 0, -2), v3(0, 1, -3)], [0, 1, 2], plusZ);
        expect(result.negIndices).toEqual([0, 1, 2]);
        expect(result.posIndices).toEqual([]);

        // Two vertices in the plane, one positive ("+00" and its rotations).
        for (const tri of [[0, 1, 2], [1, 2, 0], [2, 0, 1]]) {
            result = split.compute(
                [v3(0, 0, 1), v3(1, 0, 0), v3(0, 1, 0)], tri, plusZ);
            expect(result.posIndices).toEqual(tri);
            expect(result.negIndices).toEqual([]);
        }

        // Two vertices in the plane, one negative.
        result = split.compute(
            [v3(0, 0, -1), v3(1, 0, 0), v3(0, 1, 0)], [0, 1, 2], plusZ);
        expect(result.negIndices).toEqual([0, 1, 2]);
        expect(result.posIndices).toEqual([]);

        // One vertex in the plane, the others positive ("+0+").
        result = split.compute(
            [v3(0, 0, 1), v3(1, 0, 0), v3(0, 1, 2)], [0, 1, 2], plusZ);
        expect(result.posIndices).toEqual([0, 1, 2]);
        expect(result.negIndices).toEqual([]);
    });

    it('rejects a triangle that lies in the plane', () => {
        const split = new SplitMeshByPlane();
        const result = split.compute(
            [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)], [0, 1, 2], plusZ);
        expect(result.posIndices).toEqual([]);
        expect(result.negIndices).toEqual([]);
        expect(result.clipVertices.length).toBe(3);
    });

    it('splits a "++-" triangle into two positive and one negative', () => {
        const split = new SplitMeshByPlane();
        const result = split.compute(
            [v3(0, 0, 1), v3(1, 0, 1), v3(0, 0, -1)], [0, 1, 2], plusZ);

        // The edge (1,2) is processed before the edge (2,0), so the new
        // vertices are indices 3 and 4 in that order.
        expect(result.clipVertices.length).toBe(5);
        expect(coords(result.clipVertices[3] as Vector)).toEqual([0.5, 0, 0]);
        expect(coords(result.clipVertices[4] as Vector)).toEqual([0, 0, 0]);

        expect(result.posIndices).toEqual([0, 1, 3, 0, 3, 4]);
        expect(result.negIndices).toEqual([2, 4, 3]);

        // The three output triangles tile the input triangle.
        expect(meshArea(result.clipVertices, result.posIndices)
            + meshArea(result.clipVertices, result.negIndices))
            .toBeCloseTo(triangleArea(v3(0, 0, 1), v3(1, 0, 1), v3(0, 0, -1)),
                12);
    });

    it('splits a "--+" triangle into two negative and one positive', () => {
        const split = new SplitMeshByPlane();
        const result = split.compute(
            [v3(0, 0, -1), v3(1, 0, -1), v3(0, 0, 1)], [0, 1, 2], plusZ);
        expect(result.clipVertices.length).toBe(5);
        expect(coords(result.clipVertices[3] as Vector)).toEqual([0.5, 0, 0]);
        expect(coords(result.clipVertices[4] as Vector)).toEqual([0, 0, 0]);
        expect(result.negIndices).toEqual([0, 1, 3, 0, 3, 4]);
        expect(result.posIndices).toEqual([2, 4, 3]);
    });

    it('splits a "+-0" triangle across the vertex in the plane', () => {
        const split = new SplitMeshByPlane();
        // v0 above, v1 below, v2 in the plane.
        const result = split.compute(
            [v3(0, 0, 1), v3(0, 0, -3), v3(2, 0, 0)], [0, 1, 2], plusZ);

        // The only crossing edge is (0,1); t = 1 / (1 + 3) = 1/4.
        expect(result.clipVertices.length).toBe(4);
        expect(coords(result.clipVertices[3] as Vector)).toEqual([0, 0, 0]);
        expect(result.posIndices).toEqual([2, 0, 3]);
        expect(result.negIndices).toEqual([2, 3, 1]);
    });

    it('splits a cube at mid-height into two halves', () => {
        const { vertices, indices } = makeCube();
        const split = new SplitMeshByPlane();
        const result = split.compute(vertices, indices, plusZ);

        // Four vertical cube edges and four face diagonals cross the plane,
        // so eight vertices are appended to the original eight.
        expect(result.clipVertices.length).toBe(16);
        for (let i = 8; i < 16; ++i) {
            expect(signedDistance(plusZ, result.clipVertices[i] as Vector))
                .toBeCloseTo(0, 15);
        }

        // The two z = -1 triangles stay negative, the two z = +1 triangles
        // stay positive, and each of the eight side triangles splits into
        // three, three per side face on each side.
        expect(result.posIndices.length).toBe(3 * (2 + 12));
        expect(result.negIndices.length).toBe(3 * (2 + 12));

        verifySide(result, plusZ, result.posIndices, +1, 1e-12);
        verifySide(result, plusZ, result.negIndices, -1, 1e-12);

        // No surface area is created or destroyed by the split.
        expect(meshArea(result.clipVertices, result.posIndices)
            + meshArea(result.clipVertices, result.negIndices))
            .toBeCloseTo(meshArea(vertices, indices), 12);
        expect(meshArea(vertices, indices)).toBeCloseTo(24, 12);

        // The cut boundary is consistent: each directed on-plane edge of the
        // positive half appears reversed in the negative half.
        const posBoundary = onPlaneDirectedEdges(result, plusZ,
            result.posIndices, 1e-12);
        const negBoundary = onPlaneDirectedEdges(result, plusZ,
            result.negIndices, 1e-12);
        expect(posBoundary.length).toBe(8);
        expect(posBoundary).toEqual(reversedEdges(negBoundary));

        // The input vertices are copied, not aliased.
        expect(result.clipVertices[0]).not.toBe(vertices[0]);
        expect(coords(result.clipVertices[0] as Vector))
            .toEqual(coords(vertices[0] as Vector));
    });

    it('sends the whole cube to one side when the plane misses it', () => {
        const { vertices, indices } = makeCube();
        const split = new SplitMeshByPlane();

        let result = split.compute(vertices, indices,
            Hyperplane.fromNormalConstant(v3(0, 0, 1), -5));
        expect(result.posIndices).toEqual(indices);
        expect(result.negIndices).toEqual([]);
        expect(result.clipVertices.length).toBe(8);

        result = split.compute(vertices, indices,
            Hyperplane.fromNormalConstant(v3(0, 0, 1), 5));
        expect(result.negIndices).toEqual(indices);
        expect(result.posIndices).toEqual([]);
    });

    it('splits the cube by an oblique plane through a corner', () => {
        const { vertices, indices } = makeCube();
        const normal = v3(1, 1, 1);
        normalize(normal);
        // The plane through (1,1,1) with this normal touches only the corner
        // vertex 6, so the cube is entirely on the negative side.
        const plane = Hyperplane.fromNormalOrigin(normal, v3(1, 1, 1));
        const split = new SplitMeshByPlane();
        const result = split.compute(vertices, indices, plane);
        expect(result.clipVertices.length).toBe(8);
        expect(result.posIndices).toEqual([]);
        expect(result.negIndices.length).toBe(indices.length);
        verifySide(result, plane, result.negIndices, -1, 1e-12);
    });

    it('preserves area and sidedness for random meshes and planes', () => {
        const random = makeRandom(0xBADCAFE);
        const { vertices: cube, indices: cubeIndices } = makeCube();

        for (let trial = 0; trial < 60; ++trial) {
            // A random plane that meets the cube.
            const normal = v3(2 * random() - 1, 2 * random() - 1,
                2 * random() - 1);
            if (length(normal) < 1e-3) {
                continue;
            }
            normalize(normal);
            const plane = Hyperplane.fromNormalConstant(normal,
                1.6 * random() - 0.8);

            const split = new SplitMeshByPlane();
            const result = split.compute(cube, cubeIndices, plane);

            // Every generated vertex lies in the plane and is a convex
            // combination of two original vertices.
            for (let i = cube.length; i < result.clipVertices.length; ++i) {
                expect(Math.abs(signedDistance(plane,
                    result.clipVertices[i] as Vector))).toBeLessThan(1e-12);
            }

            verifySide(result, plane, result.posIndices, +1, 1e-9);
            verifySide(result, plane, result.negIndices, -1, 1e-9);

            expect(meshArea(result.clipVertices, result.posIndices)
                + meshArea(result.clipVertices, result.negIndices))
                .toBeCloseTo(meshArea(cube, cubeIndices), 9);

            const posBoundary = onPlaneDirectedEdges(result, plane,
                result.posIndices, 1e-9);
            const negBoundary = onPlaneDirectedEdges(result, plane,
                result.negIndices, 1e-9);
            expect(posBoundary).toEqual(reversedEdges(negBoundary));
        }
    });

    it('preserves area and sidedness for random triangle soups', () => {
        const random = makeRandom(0x1234ABCD);
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);

        for (let trial = 0; trial < 40; ++trial) {
            const vertices: Vector[] = [];
            const indices: number[] = [];
            const numTriangles = 1 + (trial % 8);
            for (let t = 0; t < numTriangles; ++t) {
                for (let j = 0; j < 3; ++j) {
                    vertices.push(v3(4 * random() - 2, 4 * random() - 2,
                        4 * random() - 2));
                    indices.push(3 * t + j);
                }
            }

            const split = new SplitMeshByPlane();
            const result = split.compute(vertices, indices, plane);

            verifySide(result, plane, result.posIndices, +1, 1e-9);
            verifySide(result, plane, result.negIndices, -1, 1e-9);

            expect(meshArea(result.clipVertices, result.posIndices)
                + meshArea(result.clipVertices, result.negIndices))
                .toBeCloseTo(meshArea(vertices, indices), 9);

            // The caller's vertices are untouched.
            for (let i = 0; i < vertices.length; ++i) {
                expect(coords(result.clipVertices[i] as Vector))
                    .toEqual(coords(vertices[i] as Vector));
                expect(result.clipVertices[i]).not.toBe(vertices[i]);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V37): independent re-check against SplitMeshByPlane.h.
// ---------------------------------------------------------------------------

// Six times the signed volume of the cone from the origin over the triangle.
// Summed over a closed mesh this is six times the enclosed volume; summed
// over any mesh it is invariant under a subdivision of the triangles that
// preserves their orientation, which is exactly what the split does.
function sixSignedVolume(vertices: readonly Vector[],
    indices: readonly number[]): number {
    let total = 0;
    for (let i = 0; i < indices.length; i += 3) {
        const p0 = vertices[indices[i] as number] as Vector;
        const p1 = vertices[indices[i + 1] as number] as Vector;
        const p2 = vertices[indices[i + 2] as number] as Vector;
        total += dot(p0, cross(p1, p2));
    }
    return total;
}

// The multiset of directed edges of a triangle list, as counts keyed by
// "a->b".
function directedEdgeCounts(indices: readonly number[]): Map<string, number> {
    const counts = new Map<string, number>();
    for (let i = 0; i < indices.length; i += 3) {
        for (let j = 0; j < 3; ++j) {
            const key = (indices[i + j] as number) + '->'
                + (indices[i + (j + 1) % 3] as number);
            counts.set(key, (counts.get(key) ?? 0) + 1);
        }
    }
    return counts;
}

// The regular octahedron, outward-facing.
function makeOctahedron(): { vertices: Vector[], indices: number[] } {
    const vertices = [
        v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
        v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)
    ];
    const indices = [
        0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
        2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5
    ];
    return { vertices, indices };
}

// A closed mesh whose triangles are consistently outward-facing, scaled and
// translated so that the split plane meets it in a variety of ways.
const closedMeshArb = fc.tuple(fc.boolean(), fc.integer({ min: 1, max: 4 }),
    fc.integer({ min: -2, max: 2 }), fc.integer({ min: -2, max: 2 }),
    fc.integer({ min: -2, max: 2 }))
    .map(([useCube, scale, tx, ty, tz]) => {
        const base = useCube ? makeCube() : makeOctahedron();
        const vertices = base.vertices.map(p => v3(
            scale * p.get(0) + tx, scale * p.get(1) + ty,
            scale * p.get(2) + tz));
        return { vertices, indices: base.indices };
    });

// A plane with a well-scaled unit normal (a subnormal normal would make the
// signed distances meaningless) whose constant is an eighth-integer offset by
// a sixteenth. The mesh coordinates are integers, so for an axis-aligned
// normal every signed distance is at least 1/16 away from zero; that keeps
// fast-check's shrinker from sliding the plane onto a vertex, where the
// output legitimately contains slivers that no fixed epsilon can classify.
const planeArb = fc.tuple(unitVector(3), fc.integer({ min: -32, max: 32 }))
    .map(([normal, k]) => Hyperplane.fromNormalConstant(normal,
        (4 * k + 2) / 16));

describe('SplitMeshByPlane verification', () => {
    it('conserves signed volume and closedness when it splits a closed mesh', () => {
        // A plane that passes within rounding distance of a vertex produces
        // legitimate slivers whose three vertices are all within any fixed
        // epsilon of the plane, which no fixed-epsilon side test can
        // classify. Those draws are skipped; the counters below keep the
        // property from becoming vacuous.
        let trials = 0;
        let skipped = 0;
        let actuallySplit = 0;
        check(fc.tuple(closedMeshArb, planeArb), ([mesh, plane]) => {
            ++trials;
            let scale = 1;
            let closest = Number.MAX_VALUE;
            for (const p of mesh.vertices) {
                scale = Math.max(scale, length(p));
                closest = Math.min(closest,
                    Math.abs(signedDistance(plane, p)));
            }
            if (closest < 1e-6 * scale) {
                ++skipped;
                return;
            }

            const split = new SplitMeshByPlane();
            const result = split.compute(mesh.vertices, mesh.indices, plane);
            if (result.clipVertices.length > mesh.vertices.length) {
                ++actuallySplit;
            }

            // Every generated vertex lies in the plane.
            for (let i = mesh.vertices.length; i < result.clipVertices.length;
                ++i) {
                const p = result.clipVertices[i] as Vector;
                expect(Math.abs(signedDistance(plane, p)))
                    .toBeLessThan(1e-9 * (1 + length(p)));
            }
            // The caller's vertices are copied, not aliased or mutated.
            for (let i = 0; i < mesh.vertices.length; ++i) {
                expect(coords(result.clipVertices[i] as Vector))
                    .toEqual(coords(mesh.vertices[i] as Vector));
                expect(result.clipVertices[i]).not.toBe(mesh.vertices[i]);
            }

            // Sidedness of every output triangle.
            verifySide(result, plane, result.posIndices, +1, 1e-9);
            verifySide(result, plane, result.negIndices, -1, 1e-9);

            // The two halves tile the original surface, so the cone volume
            // and the surface area are both conserved. The only triangles
            // that may be dropped are those lying entirely in the plane,
            // which contribute zero area and, being coplanar with a plane
            // through them, are still counted here because the input mesh is
            // closed and a random plane does not contain a whole triangle.
            const combined = result.posIndices.concat(result.negIndices);
            expectClose(sixSignedVolume(result.clipVertices, combined),
                sixSignedVolume(mesh.vertices, mesh.indices), 1e-8, 1e-9);
            expectClose(meshArea(result.clipVertices, combined),
                meshArea(mesh.vertices, mesh.indices), 1e-8, 1e-9);

            // The union of the two halves is still a closed surface: every
            // directed edge is matched by its reverse. (Each half on its own
            // is open along the cut.)
            const counts = directedEdgeCounts(combined);
            for (const [key, count] of counts) {
                const parts = key.split('->');
                const reverse = parts[1] + '->' + parts[0];
                expect(counts.get(reverse) ?? 0).toBe(count);
            }

            // The cut boundaries of the two halves are reverses of each
            // other, which is what would let a cap be attached.
            expect(onPlaneDirectedEdges(result, plane, result.posIndices, 1e-9))
                .toEqual(reversedEdges(onPlaneDirectedEdges(result, plane,
                    result.negIndices, 1e-9)));

            // Every strictly positive original vertex survives in the
            // positive output and every strictly negative one in the
            // negative output.
            const posSet = new Set<number>(result.posIndices);
            const negSet = new Set<number>(result.negIndices);
            for (let i = 0; i < mesh.vertices.length; ++i) {
                const sDist = signedDistance(plane,
                    mesh.vertices[i] as Vector);
                if (sDist > 0) {
                    expect(posSet.has(i)).toBe(true);
                }
                else if (sDist < 0) {
                    expect(negSet.has(i)).toBe(true);
                }
            }

            // Reusing the query object gives the same answer (mEMap and the
            // signed distances are reset by compute()).
            const again = split.compute(mesh.vertices, mesh.indices, plane);
            expect(again.posIndices).toEqual(result.posIndices);
            expect(again.negIndices).toEqual(result.negIndices);
            expect(again.clipVertices.map(coords)).toEqual(
                result.clipVertices.map(coords));
        }, 150);
        expect(trials - skipped).toBeGreaterThan(140);
        expect(actuallySplit).toBeGreaterThan(40);
    }, 30000);

    it('sends a mesh that misses the plane entirely to one side', () => {
        check(fc.tuple(closedMeshArb, unitVector(3)), ([mesh, normal]) => {
            // Push the plane past the whole mesh in both directions.
            let extreme = 0;
            for (const p of mesh.vertices) {
                extreme = Math.max(extreme, Math.abs(dot(normal, p)));
            }
            const split = new SplitMeshByPlane();

            const below = split.compute(mesh.vertices, mesh.indices,
                Hyperplane.fromNormalConstant(normal, -extreme - 1));
            expect(below.negIndices).toEqual([]);
            expect(below.posIndices).toEqual(Array.from(mesh.indices));
            expect(below.clipVertices.length).toBe(mesh.vertices.length);

            const above = split.compute(mesh.vertices, mesh.indices,
                Hyperplane.fromNormalConstant(normal, extreme + 1));
            expect(above.posIndices).toEqual([]);
            expect(above.negIndices).toEqual(Array.from(mesh.indices));
            expect(above.clipVertices.length).toBe(mesh.vertices.length);
        }, 100);
    });

    it('handles the vertex-in-plane cases without generating a vertex there', () => {
        // A triangle with one vertex exactly in the plane and the other two
        // strictly on opposite sides is split by exactly one new vertex, on
        // the edge joining the two off-plane vertices. This is the "+-0" /
        // "-+0" family, which upstream routes through SplitTrianglePMZ and
        // SplitTriangleMPZ.
        check(fc.tuple(fc.integer({ min: 1, max: 5 }),
            fc.integer({ min: 1, max: 5 }), fc.integer({ min: 0, max: 2 }),
            fc.boolean()), ([above, below, rotate, flip]) => {
            const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);
            const triangle = [
                v3(0, 0, 0),
                v3(1, 0, flip ? above : -below),
                v3(0, 1, flip ? -below : above)
            ];
            const indices = [rotate % 3, (rotate + 1) % 3, (rotate + 2) % 3];

            const split = new SplitMeshByPlane();
            const result = split.compute(triangle, indices, plane);

            // Exactly one new vertex, on the edge between vertices 1 and 2
            // at the parameter t = sDist1 / (sDist1 - sDist2) that upstream
            // computes for that edge.
            expect(result.clipVertices.length).toBe(4);
            const generated = result.clipVertices[3] as Vector;
            const z1 = triangle[1]!.get(2);
            const z2 = triangle[2]!.get(2);
            const t = z1 / (z1 - z2);
            expectVectorClose(generated, v3(1 - t, t, 0), 1e-12, 1e-12);
            expectClose(signedDistance(plane, generated), 0, 1e-12, 1e-12);

            // One triangle on each side, each using the on-plane vertex 0 and
            // the generated vertex.
            expect(result.posIndices.length).toBe(3);
            expect(result.negIndices.length).toBe(3);
            expect(result.posIndices).toContain(0);
            expect(result.posIndices).toContain(3);
            expect(result.negIndices).toContain(0);
            expect(result.negIndices).toContain(3);

            verifySide(result, plane, result.posIndices, +1, 1e-12);
            verifySide(result, plane, result.negIndices, -1, 1e-12);
            expectClose(
                meshArea(result.clipVertices, result.posIndices)
                + meshArea(result.clipVertices, result.negIndices),
                meshArea(triangle, indices), 1e-12, 1e-12);
        }, 120);
    });

    it('drops only the triangles that lie entirely in the plane', () => {
        // A closed mesh whose equator is a ring of vertices in the plane:
        // the two "000" triangles of the flat band are rejected and the rest
        // are classified without any new vertex.
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);
        check(fc.integer({ min: 1, max: 4 }), height => {
            const vertices = [
                v3(1, 0, 0), v3(0, 1, 0), v3(-1, 0, 0), v3(0, -1, 0),
                v3(0, 0, height), v3(0, 0, -height)
            ];
            // The eight side triangles plus two coplanar triangles that
            // tile the equatorial square.
            const indices = [
                0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4,
                1, 0, 5, 2, 1, 5, 3, 2, 5, 0, 3, 5,
                0, 1, 2, 0, 2, 3
            ];

            const split = new SplitMeshByPlane();
            const result = split.compute(vertices, indices, plane);

            // No edge crosses the plane strictly, so no vertex is generated.
            expect(result.clipVertices.length).toBe(vertices.length);
            // Four triangles above, four below, two rejected.
            expect(result.posIndices.length).toBe(12);
            expect(result.negIndices.length).toBe(12);
            expect(result.posIndices).toEqual(indices.slice(0, 12));
            expect(result.negIndices).toEqual(indices.slice(12, 24));
        }, 4);
    });
});
