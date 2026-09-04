import { describe, it, expect } from 'vitest';
import { ConvexPolyhedron3 } from '../src/ConvexPolyhedron3.js';
import { Vector, dot } from '../src/Vector.js';
import { check, fc, rotationFrame, scaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md).
//
// The class is a container: GeneratePlanes builds one outer-pointing,
// *unnormalized* plane per triangle and GenerateAlignedBox is ComputeExtremes
// over the vertex array. The properties below check both against independent
// computations on random affine images of fixed convex solids (a rotation with
// positive scales preserves convexity and face orientation).
// ---------------------------------------------------------------------------

function octaVertices(): Vector[] {
    return [v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
        v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)];
}

function octaIndices(): number[] {
    return [
        0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
        2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5
    ];
}

interface Solid { vertices: Vector[]; indices: number[]; interior: Vector }

const baseSolids: Solid[] = [
    { vertices: cubeVertices(), indices: cubeIndices(), interior: v3(0.5, 0.5, 0.5) },
    { vertices: tetraVertices(), indices: tetraIndices(), interior: v3(0.25, 0.25, 0.25) },
    { vertices: octaVertices(), indices: octaIndices(), interior: v3(0, 0, 0) }
];

// A random affine image of one of the base solids: X -> R * S * X + T with R a
// rotation and S positive diagonal, so the image is convex, the faces keep
// their counterclockwise-from-outside order and no face degenerates.
const affineSolid = fc.record({
    which: fc.integer({ min: 0, max: baseSolids.length - 1 }),
    frame: rotationFrame(3),
    scale: fc.tuple(scaled(0.5, 3, 16), scaled(0.5, 3, 16), scaled(0.5, 3, 16)),
    translate: fc.tuple(scaled(-5, 5, 32), scaled(-5, 5, 32), scaled(-5, 5, 32))
}).map(({ which, frame, scale, translate }) => {
    const base = baseSolids[which];
    const map = (p: Vector): Vector => {
        const s = [scale[0] * p.get(0), scale[1] * p.get(1), scale[2] * p.get(2)];
        return v3(
            frame[0].get(0) * s[0] + frame[1].get(0) * s[1] + frame[2].get(0) * s[2] + translate[0],
            frame[0].get(1) * s[0] + frame[1].get(1) * s[1] + frame[2].get(1) * s[2] + translate[1],
            frame[0].get(2) * s[0] + frame[1].get(2) * s[1] + frame[2].get(2) * s[2] + translate[2]);
    };
    return {
        vertices: base.vertices.map(map),
        indices: base.indices.slice(),
        interior: map(base.interior)
    };
});

describe('ConvexPolyhedron3 verification', () => {
    it('planes are the unnormalized outer-pointing face planes', () => {
        check(affineSolid, ({ vertices, indices }) => {
            const poly = new ConvexPolyhedron3(vertices.map(p => p.clone()),
                indices.slice(), true, false);
            expect(poly.planes.length).toBe(indices.length / 3);
            for (let t = 0; t < poly.planes.length; ++t) {
                const V0 = poly.vertices[poly.indices[3 * t + 0]];
                const V1 = poly.vertices[poly.indices[3 * t + 1]];
                const V2 = poly.vertices[poly.indices[3 * t + 2]];
                // Independent recomputation of N = (V1-V0) x (V2-V0).
                const e1 = [V1.get(0) - V0.get(0), V1.get(1) - V0.get(1), V1.get(2) - V0.get(2)];
                const e2 = [V2.get(0) - V0.get(0), V2.get(1) - V0.get(1), V2.get(2) - V0.get(2)];
                const n = [e1[1] * e2[2] - e1[2] * e2[1],
                    e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0]];
                for (let i = 0; i < 3; ++i) { expect(poly.planes[t].get(i)).toBe(n[i]); }
                expect(poly.planes[t].get(3))
                    .toBe(-(n[0] * V0.get(0) + n[1] * V0.get(1) + n[2] * V0.get(2)));
                // The plane vanishes on its own three vertices. The evaluation
                // subtracts products of coordinates of size |N|*|V|, so the
                // tolerance is relative to that magnitude.
                const magnitude = Math.hypot(n[0], n[1], n[2])
                    * Math.max(1, Math.hypot(V0.get(0), V0.get(1), V0.get(2)));
                for (const V of [V0, V1, V2]) {
                    expect(Math.abs(evalPlane(poly.planes[t], V)))
                        .toBeLessThanOrEqual(1e-12 * magnitude);
                }
            }
        });
    });

    it('the polyhedron is the intersection of its face halfspaces', () => {
        check(fc.tuple(affineSolid, fc.array(fc.tuple(
            scaled(0, 1, 16), scaled(0, 1, 16), scaled(0.05, 1, 16)),
        { minLength: 4, maxLength: 4 })), ([{ vertices, indices, interior }, blends]) => {
            const poly = new ConvexPolyhedron3(vertices.map(p => p.clone()),
                indices.slice(), true, true);
            const scale = Math.max(1, ...poly.vertices.map(p =>
                Math.hypot(p.get(0), p.get(1), p.get(2))));
            for (let t = 0; t < poly.planes.length; ++t) {
                const plane = poly.planes[t];
                const normalLength = Math.hypot(plane.get(0), plane.get(1), plane.get(2));
                const tol = 1e-12 * normalLength * scale;
                // Every vertex is on or inside every face plane (convexity and
                // outer-pointing normals).
                for (const V of poly.vertices) {
                    expect(evalPlane(plane, V)).toBeLessThanOrEqual(tol);
                }
                // A point strictly inside is strictly inside every halfspace.
                expect(evalPlane(plane, interior)).toBeLessThan(0);
                // Blends of the interior point with a vertex stay inside.
                for (const [u] of blends) {
                    const V = poly.vertices[Math.min(poly.vertices.length - 1,
                        Math.floor(u * poly.vertices.length))];
                    const mid = v3(0.5 * (V.get(0) + interior.get(0)),
                        0.5 * (V.get(1) + interior.get(1)),
                        0.5 * (V.get(2) + interior.get(2)));
                    expect(evalPlane(plane, mid)).toBeLessThanOrEqual(tol);
                }
                // Stepping outward from a face vertex along its own normal
                // leaves that halfspace.
                const V0 = poly.vertices[poly.indices[3 * t + 0]];
                const outward = v3(V0.get(0) + plane.get(0) / normalLength,
                    V0.get(1) + plane.get(1) / normalLength,
                    V0.get(2) + plane.get(2) / normalLength);
                expect(evalPlane(plane, outward)).toBeGreaterThan(0);
            }
        }, 100);
    });

    it('the aligned box is the componentwise extreme of the vertices', () => {
        check(affineSolid, ({ vertices, indices }) => {
            const poly = new ConvexPolyhedron3(vertices.map(p => p.clone()),
                indices.slice(), false, true);
            for (let i = 0; i < 3; ++i) {
                const values = poly.vertices.map(p => p.get(i));
                expect(poly.alignedBox.min.get(i)).toBe(Math.min(...values));
                expect(poly.alignedBox.max.get(i)).toBe(Math.max(...values));
            }
        });
    });

    it('construction succeeds only for >= 4 vertices and >= 12 indices', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 8 }), fc.integer({ min: 0, max: 15 })),
            ([numVertices, numIndices]) => {
                const vertices: Vector[] = [];
                for (let i = 0; i < numVertices; ++i) { vertices.push(v3(i, i * i, 1)); }
                const indices: number[] = [];
                for (let i = 0; i < numIndices; ++i) {
                    indices.push(numVertices > 0 ? i % numVertices : 0);
                }
                const poly = new ConvexPolyhedron3(vertices, indices, true, true);
                const ok = numVertices >= 4 && numIndices >= 12;
                expect(poly.vertices.length).toBe(ok ? numVertices : 0);
                expect(poly.indices.length).toBe(ok ? numIndices : 0);
                // Upstream (issue #175): the constructor never checks
                // indices.size() % 3 == 0, so a trailing partial triangle is
                // silently dropped. Preserved here.
                expect(poly.planes.length).toBe(ok ? Math.floor(numIndices / 3) : 0);
                if (!ok) {
                    // Failed construction leaves the default box, [-1,1]^3.
                    expect(poly.alignedBox.min.values).toEqual([-1, -1, -1]);
                    expect(poly.alignedBox.max.values).toEqual([1, 1, 1]);
                }
            });
    });

    it('the constructor moves (aliases) the input arrays, as documented', () => {
        check(affineSolid, ({ vertices, indices }) => {
            const poly = new ConvexPolyhedron3(vertices, indices, false, false);
            expect(poly.vertices).toBe(vertices);
            expect(poly.indices).toBe(indices);
            // Regenerating after an in-place edit picks the new data up.
            poly.generateAlignedBox();
            const before = poly.alignedBox.max.get(0);
            vertices[0].set(0, before + 10);
            poly.generateAlignedBox();
            expect(poly.alignedBox.max.get(0)).toBe(before + 10);
        });
    });
});
