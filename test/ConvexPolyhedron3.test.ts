import { describe, it, expect } from 'vitest';
import { ConvexPolyhedron3 } from '../src/ConvexPolyhedron3';
import { Vector, dot } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The unit cube [0,1]^3 with counterclockwise-when-viewed-from-outside faces.
function cubeVertices(): Vector[] {
    return [
        v3(0, 0, 0), v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0),
        v3(0, 0, 1), v3(1, 0, 1), v3(1, 1, 1), v3(0, 1, 1)
    ];
}

function cubeIndices(): number[] {
    return [
        // -z face (normal (0,0,-1))
        0, 3, 2, 0, 2, 1,
        // +z face
        4, 5, 6, 4, 6, 7,
        // -y face
        0, 1, 5, 0, 5, 4,
        // +y face
        3, 7, 6, 3, 6, 2,
        // -x face
        0, 4, 7, 0, 7, 3,
        // +x face
        1, 2, 6, 1, 6, 5
    ];
}

// A regular-ish tetrahedron with outward counterclockwise faces.
function tetraVertices(): Vector[] {
    return [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
}

function tetraIndices(): number[] {
    return [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
}

// Evaluate the plane (n0,n1,n2,d) at the point p: Dot(N,p) + d.
function evalPlane(plane: Vector, p: Vector): number {
    return plane.values[0] * p.values[0] + plane.values[1] * p.values[1]
        + plane.values[2] * p.values[2] + plane.values[3];
}

// Compare a plane against expected coefficients, tolerating the signed zeros
// that the cross product and the negated dot product can produce.
function expectPlane(plane: Vector, expected: [number, number, number, number]): void {
    for (let i = 0; i < 4; ++i) {
        expect(plane.values[i]).toBeCloseTo(expected[i], 12);
    }
}

describe('ConvexPolyhedron3', () => {
    it('default-constructs an empty polyhedron with the default aligned box', () => {
        const poly = new ConvexPolyhedron3();
        expect(poly.vertices).toEqual([]);
        expect(poly.indices).toEqual([]);
        expect(poly.planes).toEqual([]);
        expect(poly.alignedBox.dimension).toBe(3);
        expect(poly.alignedBox.min.values).toEqual([-1, -1, -1]);
        expect(poly.alignedBox.max.values).toEqual([1, 1, 1]);
    });

    it('constructs a cube and stores the moved arrays', () => {
        const vertices = cubeVertices();
        const indices = cubeIndices();
        const poly = new ConvexPolyhedron3(vertices, indices, false, false);
        expect(poly.vertices.length).toBe(8);
        expect(poly.indices.length).toBe(36);
        // The arrays are moved (referenced), not copied.
        expect(poly.vertices).toBe(vertices);
        expect(poly.indices).toBe(indices);
        // Neither planes nor the box were requested.
        expect(poly.planes).toEqual([]);
        expect(poly.alignedBox.min.values).toEqual([-1, -1, -1]);
    });

    it('fails construction with fewer than 4 vertices, leaving arrays empty', () => {
        const poly = new ConvexPolyhedron3(
            [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)],
            tetraIndices(), true, true);
        expect(poly.vertices).toEqual([]);
        expect(poly.indices).toEqual([]);
        expect(poly.planes).toEqual([]);
        expect(poly.alignedBox.min.values).toEqual([-1, -1, -1]);
    });

    it('fails construction with fewer than 12 indices', () => {
        const poly = new ConvexPolyhedron3(tetraVertices(), [0, 1, 2, 0, 1, 3],
            true, true);
        expect(poly.vertices).toEqual([]);
        expect(poly.indices).toEqual([]);
    });

    it('generates one plane per triangle for the cube', () => {
        const poly = new ConvexPolyhedron3(cubeVertices(), cubeIndices(), true, false);
        expect(poly.planes.length).toBe(12);
        for (const plane of poly.planes) {
            expect(plane.size).toBe(4);
        }
    });

    it('generates outer-pointing cube planes with the expected values', () => {
        const poly = new ConvexPolyhedron3(cubeVertices(), cubeIndices(), true, false);

        // The two triangles of the -z face share the plane
        // Dot((0,0,-1), X) = 0, stored with an unnormalized normal of
        // length 1 for this axis-aligned cube.
        expectPlane(poly.planes[0], [0, 0, -1, 0]);
        expectPlane(poly.planes[1], [0, 0, -1, 0]);

        // The +z face: Dot((0,0,1), X) - 1 = 0.
        expectPlane(poly.planes[2], [0, 0, 1, -1]);
        expectPlane(poly.planes[3], [0, 0, 1, -1]);

        // The +x face: Dot((1,0,0), X) - 1 = 0.
        expectPlane(poly.planes[10], [1, 0, 0, -1]);
        expectPlane(poly.planes[11], [1, 0, 0, -1]);
    });

    it('has all cube planes negative at the interior and zero on the face', () => {
        const poly = new ConvexPolyhedron3(cubeVertices(), cubeIndices(), true, false);
        const center = v3(0.5, 0.5, 0.5);
        for (let t = 0; t < poly.planes.length; ++t) {
            expect(evalPlane(poly.planes[t], center)).toBeLessThan(0);
            // Every vertex of the triangle lies on its plane.
            for (let k = 0; k < 3; ++k) {
                const V = poly.vertices[poly.indices[3 * t + k]];
                expect(evalPlane(poly.planes[t], V)).toBeCloseTo(0, 12);
            }
        }
        // Every vertex of the cube is inside or on every face plane.
        for (const V of poly.vertices) {
            for (const plane of poly.planes) {
                expect(evalPlane(plane, V)).toBeLessThanOrEqual(1e-12);
            }
        }
    });

    it('generates consistent tetrahedron planes', () => {
        const poly = new ConvexPolyhedron3(tetraVertices(), tetraIndices(), true, false);
        expect(poly.planes.length).toBe(4);

        const centroid = v3(0.25, 0.25, 0.25);
        for (let t = 0; t < 4; ++t) {
            expect(evalPlane(poly.planes[t], centroid)).toBeLessThan(0);
            for (let k = 0; k < 3; ++k) {
                const V = poly.vertices[poly.indices[3 * t + k]];
                expect(evalPlane(poly.planes[t], V)).toBeCloseTo(0, 12);
            }
        }

        // The slanted face x+y+z = 1 has an unnormalized outer normal of
        // (1,1,1) up to the positive scale produced by the cross product.
        const slant = poly.planes[3];
        const n = v3(slant.values[0], slant.values[1], slant.values[2]);
        expect(n.values[0]).toBeGreaterThan(0);
        expect(n.values[1]).toBeCloseTo(n.values[0], 12);
        expect(n.values[2]).toBeCloseTo(n.values[0], 12);
        expect(slant.values[3]).toBeCloseTo(-n.values[0], 12);
    });

    it('the plane normal is the (unnormalized) cross product of the edges', () => {
        const poly = new ConvexPolyhedron3(tetraVertices(), tetraIndices(), true, false);
        for (let t = 0; t < 4; ++t) {
            const V0 = poly.vertices[poly.indices[3 * t + 0]];
            const V1 = poly.vertices[poly.indices[3 * t + 1]];
            const V2 = poly.vertices[poly.indices[3 * t + 2]];
            const e1 = [V1.values[0] - V0.values[0], V1.values[1] - V0.values[1],
                V1.values[2] - V0.values[2]];
            const e2 = [V2.values[0] - V0.values[0], V2.values[1] - V0.values[1],
                V2.values[2] - V0.values[2]];
            const n = v3(
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0]);
            expect(poly.planes[t].values[0]).toBeCloseTo(n.values[0], 12);
            expect(poly.planes[t].values[1]).toBeCloseTo(n.values[1], 12);
            expect(poly.planes[t].values[2]).toBeCloseTo(n.values[2], 12);
            expect(poly.planes[t].values[3]).toBeCloseTo(-dot(n, V0), 12);
        }
    });

    it('generates the aligned box on request', () => {
        const poly = new ConvexPolyhedron3(cubeVertices(), cubeIndices(), false, true);
        expect(poly.alignedBox.min.values).toEqual([0, 0, 0]);
        expect(poly.alignedBox.max.values).toEqual([1, 1, 1]);
    });

    it('generates both planes and the aligned box for a shifted box', () => {
        const vertices = cubeVertices().map(v =>
            v3(2 * v.values[0] - 3, 4 * v.values[1] + 1, v.values[2] - 5));
        const poly = new ConvexPolyhedron3(vertices, cubeIndices(), true, true);
        expect(poly.planes.length).toBe(12);
        expect(poly.alignedBox.min.values).toEqual([-3, 1, -5]);
        expect(poly.alignedBox.max.values).toEqual([-1, 5, -4]);

        const center = v3(-2, 3, -4.5);
        for (const plane of poly.planes) {
            expect(evalPlane(plane, center)).toBeLessThan(0);
        }
    });

    it('recomputes planes and box after the vertices are modified', () => {
        const vertices = cubeVertices();
        const poly = new ConvexPolyhedron3(vertices, cubeIndices(), true, true);
        expect(poly.alignedBox.max.values).toEqual([1, 1, 1]);

        // Scale the cube by 3 in place, then regenerate.
        for (const v of poly.vertices) {
            for (let d = 0; d < 3; ++d) {
                v.values[d] *= 3;
            }
        }
        poly.generatePlanes();
        poly.generateAlignedBox();

        expect(poly.alignedBox.min.values).toEqual([0, 0, 0]);
        expect(poly.alignedBox.max.values).toEqual([3, 3, 3]);
        // The +x plane is now Dot((3,0,0), X) - 9 = 0 (normals scale with the
        // cross product, so they are not unit length).
        expectPlane(poly.planes[10], [9, 0, 0, -27]);
    });

    it('leaves planes and box untouched when the arrays are empty', () => {
        const poly = new ConvexPolyhedron3();
        poly.generatePlanes();
        poly.generateAlignedBox();
        expect(poly.planes).toEqual([]);
        expect(poly.alignedBox.min.values).toEqual([-1, -1, -1]);
    });
});
