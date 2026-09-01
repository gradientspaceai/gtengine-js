import { describe, it, expect } from 'vitest';
import { MeshCurvature } from '../src/MeshCurvature';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

function V3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

interface TriMesh {
    vertices: Vector[];
    indices: number[];
}

// An icosphere: an icosahedron subdivided 'level' times and projected onto
// the sphere of radius r. Every vertex is unique, so the mesh is manifold
// and closed.
function icosphere(r: number, level: number): TriMesh {
    const t = (1 + Math.sqrt(5)) / 2;
    let vertices: Vector[] = [
        V3(-1, t, 0), V3(1, t, 0), V3(-1, -t, 0), V3(1, -t, 0),
        V3(0, -1, t), V3(0, 1, t), V3(0, -1, -t), V3(0, 1, -t),
        V3(t, 0, -1), V3(t, 0, 1), V3(-t, 0, -1), V3(-t, 0, 1)
    ];
    let indices: number[] = [
        0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11,
        1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
        3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9,
        4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1
    ];

    for (let s = 0; s < level; ++s) {
        const midpoints = new Map<string, number>();
        const nextIndices: number[] = [];
        const midpoint = (a: number, b: number): number => {
            const key = (a < b ? a + ',' + b : b + ',' + a);
            const existing = midpoints.get(key);
            if (existing !== undefined) {
                return existing;
            }
            const m = mul(add(vertices[a], vertices[b]), 0.5);
            vertices.push(m);
            const index = vertices.length - 1;
            midpoints.set(key, index);
            return index;
        };
        for (let i = 0; i < indices.length; i += 3) {
            const v0 = indices[i], v1 = indices[i + 1], v2 = indices[i + 2];
            const a = midpoint(v0, v1);
            const b = midpoint(v1, v2);
            const c = midpoint(v2, v0);
            nextIndices.push(v0, a, c, v1, b, a, v2, c, b, a, b, c);
        }
        indices = nextIndices;
    }

    vertices = vertices.map(v => {
        const u = v.clone();
        normalize(u);
        return mul(u, r);
    });
    return { vertices, indices };
}

// An open cylinder of radius r and height h with numAxial+1 rings of
// numRadial vertices. The vertices are unique; the mesh wraps around in the
// angular direction.
function cylinder(r: number, h: number, numRadial: number,
    numAxial: number): TriMesh {
    const vertices: Vector[] = [];
    for (let j = 0; j <= numAxial; ++j) {
        const z = h * (j / numAxial - 0.5);
        for (let i = 0; i < numRadial; ++i) {
            const angle = 2 * Math.PI * i / numRadial;
            vertices.push(V3(r * Math.cos(angle), r * Math.sin(angle), z));
        }
    }

    const indices: number[] = [];
    for (let j = 0; j < numAxial; ++j) {
        for (let i = 0; i < numRadial; ++i) {
            const i1 = (i + 1) % numRadial;
            const a = j * numRadial + i;
            const b = j * numRadial + i1;
            const c = (j + 1) * numRadial + i;
            const d = (j + 1) * numRadial + i1;
            // Counterclockwise seen from outside.
            indices.push(a, b, c, b, d, c);
        }
    }
    return { vertices, indices };
}

// A planar grid in the plane z = 0 with normals along +z.
function plane(size: number, n: number): TriMesh {
    const vertices: Vector[] = [];
    for (let j = 0; j <= n; ++j) {
        for (let i = 0; i <= n; ++i) {
            vertices.push(V3(size * (i / n - 0.5), size * (j / n - 0.5), 0));
        }
    }
    const indices: number[] = [];
    for (let j = 0; j < n; ++j) {
        for (let i = 0; i < n; ++i) {
            const a = j * (n + 1) + i;
            const b = a + 1;
            const c = a + (n + 1);
            const d = c + 1;
            indices.push(a, b, c, b, d, c);
        }
    }
    return { vertices, indices };
}

describe('MeshCurvature', () => {
    it('computes area-weighted unit normals', () => {
        const r = 2;
        const { vertices, indices } = icosphere(r, 2);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 0);

        const normals = mc.getNormals();
        expect(normals.length).toBe(vertices.length);
        for (let i = 0; i < vertices.length; ++i) {
            expect(length(normals[i])).toBeCloseTo(1, 12);
            // Outward on the sphere.
            expect(dot(normals[i], vertices[i]) / r).toBeGreaterThan(0.999);
        }
    });

    it('estimates the curvatures of a sphere of radius r as 1/r', () => {
        for (const r of [0.5, 1, 3]) {
            const { vertices, indices } = icosphere(r, 3);
            const mc = new MeshCurvature();
            mc.compute(vertices, indices, 0);

            const kMin = mc.getMinCurvatures();
            const kMax = mc.getMaxCurvatures();
            expect(kMin.length).toBe(vertices.length);
            for (let i = 0; i < vertices.length; ++i) {
                expect(kMin[i]).toBeGreaterThan(0);
                expect(kMin[i]).toBeLessThanOrEqual(kMax[i]);
                expect(Math.abs(kMin[i] * r - 1)).toBeLessThan(0.05);
                expect(Math.abs(kMax[i] * r - 1)).toBeLessThan(0.05);
            }
        }
    });

    it('produces principal directions orthogonal to the normal', () => {
        const { vertices, indices } = icosphere(1.5, 2);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 0);

        const normals = mc.getNormals();
        const dMin = mc.getMinDirections();
        const dMax = mc.getMaxDirections();
        for (let i = 0; i < vertices.length; ++i) {
            expect(length(dMin[i])).toBeCloseTo(1, 10);
            expect(length(dMax[i])).toBeCloseTo(1, 10);
            expect(dot(dMin[i], normals[i])).toBeCloseTo(0, 10);
            expect(dot(dMax[i], normals[i])).toBeCloseTo(0, 10);
        }
    });

    it('estimates the curvatures of a cylinder as (0, 1/r)', () => {
        const r = 2, h = 4, numRadial = 48, numAxial = 12;
        const { vertices, indices } = cylinder(r, h, numRadial, numAxial);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 0);

        const normals = mc.getNormals();
        const kMin = mc.getMinCurvatures();
        const kMax = mc.getMaxCurvatures();
        const dMin = mc.getMinDirections();
        const dMax = mc.getMaxDirections();

        // Interior rings only; the boundary rings have one-sided fans.
        for (let j = 1; j < numAxial; ++j) {
            for (let i = 0; i < numRadial; ++i) {
                const v = j * numRadial + i;

                // The normal is radial and outward.
                const radial = V3(vertices[v].values[0], vertices[v].values[1], 0);
                expect(dot(normals[v], radial) / r).toBeGreaterThan(0.999);

                expect(Math.abs(kMin[v])).toBeLessThan(0.01);
                expect(Math.abs(kMax[v] * r - 1)).toBeLessThan(0.02);

                // The direction of zero curvature is the cylinder axis and
                // the direction of curvature 1/r is circumferential.
                expect(Math.abs(dMin[v].values[2])).toBeGreaterThan(0.99);
                expect(Math.abs(dMax[v].values[2])).toBeLessThan(0.06);
                expect(Math.abs(dot(dMin[v], dMax[v]))).toBeLessThan(1e-8);
            }
        }
    });

    it('estimates zero curvature on a plane', () => {
        const n = 8;
        const { vertices, indices } = plane(4, n);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 0);

        const normals = mc.getNormals();
        const kMin = mc.getMinCurvatures();
        const kMax = mc.getMaxCurvatures();
        for (let i = 0; i < vertices.length; ++i) {
            expect(normals[i].values[0]).toBeCloseTo(0, 12);
            expect(normals[i].values[1]).toBeCloseTo(0, 12);
            expect(Math.abs(normals[i].values[2])).toBeCloseTo(1, 12);
            expect(kMin[i]).toBeCloseTo(0, 10);
            expect(kMax[i]).toBeCloseTo(0, 10);
        }
    });

    it('flags locally planar points via the singularity threshold', () => {
        // On the plane every D*W^T is exactly zero, so any positive
        // threshold flags every vertex; the directions are then the
        // orthogonal complement of the normal.
        const { vertices, indices } = plane(4, 4);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 1e-8);

        const normals = mc.getNormals();
        const dMin = mc.getMinDirections();
        const dMax = mc.getMaxDirections();
        for (let i = 0; i < vertices.length; ++i) {
            expect(mc.getMinCurvatures()[i]).toBe(0);
            expect(mc.getMaxCurvatures()[i]).toBe(0);
            expect(length(dMin[i])).toBeCloseTo(1, 12);
            expect(length(dMax[i])).toBeCloseTo(1, 12);
            expect(dot(dMin[i], normals[i])).toBeCloseTo(0, 12);
            expect(dot(dMax[i], normals[i])).toBeCloseTo(0, 12);
            // {dMin, dMax, N} is right-handed.
            const c = cross(dMin[i], dMax[i]);
            expect(length(sub(c, normals[i]))).toBeCloseTo(0, 12);
        }
    });

    it('flags every vertex of a sphere when the threshold is huge', () => {
        const { vertices, indices } = icosphere(1, 1);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 1e6);
        for (let i = 0; i < vertices.length; ++i) {
            expect(mc.getMinCurvatures()[i]).toBe(0);
            expect(mc.getMaxCurvatures()[i]).toBe(0);
        }
    });

    it('reverses the sign of the curvatures when the mesh is inverted', () => {
        // Flipping the triangle winding flips the estimated normals, which
        // flips the sign of the shape operator.
        const r = 1.5;
        const { vertices, indices } = icosphere(r, 2);
        const flipped: number[] = [];
        for (let i = 0; i < indices.length; i += 3) {
            flipped.push(indices[i], indices[i + 2], indices[i + 1]);
        }

        const outward = new MeshCurvature();
        outward.compute(vertices, indices, 0);
        const inward = new MeshCurvature();
        inward.compute(vertices, flipped, 0);

        for (let i = 0; i < vertices.length; ++i) {
            expect(inward.getMinCurvatures()[i]).toBeCloseTo(
                -outward.getMaxCurvatures()[i], 7);
            expect(inward.getMaxCurvatures()[i]).toBeCloseTo(
                -outward.getMinCurvatures()[i], 7);
            expect(length(add(inward.getNormals()[i],
                outward.getNormals()[i]))).toBeCloseTo(0, 10);
        }
    });

    it('is invariant under rigid motions of the mesh', () => {
        const r = 1.25;
        const { vertices, indices } = icosphere(r, 2);

        // Rotate about the axis (1,1,1)/sqrt(3) by pi/3 and translate.
        const axis = V3(1, 1, 1);
        normalize(axis);
        const angle = Math.PI / 3;
        const cs = Math.cos(angle), sn = Math.sin(angle);
        const translation = V3(3, -2, 7);
        const moved = vertices.map(v => {
            const term0 = mul(v, cs);
            const term1 = mul(cross(axis, v), sn);
            const term2 = mul(axis, dot(axis, v) * (1 - cs));
            return add(add(add(term0, term1), term2), translation);
        });

        const a = new MeshCurvature();
        a.compute(vertices, indices, 0);
        const b = new MeshCurvature();
        b.compute(moved, indices, 0);
        for (let i = 0; i < vertices.length; ++i) {
            expect(b.getMinCurvatures()[i]).toBeCloseTo(a.getMinCurvatures()[i], 6);
            expect(b.getMaxCurvatures()[i]).toBeCloseTo(a.getMaxCurvatures()[i], 6);
        }
    });

    it('handles an empty mesh and a single triangle', () => {
        const mc = new MeshCurvature();
        mc.compute([], [], 0);
        expect(mc.getNormals().length).toBe(0);
        expect(mc.getMinCurvatures().length).toBe(0);

        const single = new MeshCurvature();
        single.compute([V3(0, 0, 0), V3(1, 0, 0), V3(0, 1, 0)], [0, 1, 2], 1e-12);
        // A single flat triangle has zero normal derivatives everywhere.
        for (let i = 0; i < 3; ++i) {
            expect(single.getNormals()[i].values[2]).toBeCloseTo(1, 12);
            expect(single.getMinCurvatures()[i]).toBe(0);
            expect(single.getMaxCurvatures()[i]).toBe(0);
        }
    });
});
