import { describe, it, expect } from 'vitest';
import { ConvexMesh3 } from '../src/ConvexMesh3.js';
import type {
    ConvexMesh3Triangle, ConvexMesh3Vertex
} from '../src/ConvexMesh3.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function vertex(x: number, y: number, z: number): ConvexMesh3Vertex {
    return Vector.fromArray([x, y, z]);
}

describe('ConvexMesh3 construction', () => {
    it('default constructs to the empty configuration', () => {
        const mesh = new ConvexMesh3();
        expect(mesh.configuration).toBe(ConvexMesh3.CFG_EMPTY);
        expect(mesh.vertices).toEqual([]);
        expect(mesh.triangles).toEqual([]);
    });

    it('gives each instance its own containers', () => {
        const mesh0 = new ConvexMesh3();
        const mesh1 = new ConvexMesh3();
        mesh0.vertices.push(vertex(1, 2, 3));
        mesh0.triangles.push([0, 1, 2]);
        expect(mesh1.vertices).toHaveLength(0);
        expect(mesh1.triangles).toHaveLength(0);
    });
});

describe('ConvexMesh3 configuration constants', () => {
    it('are the upstream single-bit flags', () => {
        expect(ConvexMesh3.CFG_EMPTY).toBe(0x00000000);
        expect(ConvexMesh3.CFG_POINT).toBe(0x00000001);
        expect(ConvexMesh3.CFG_SEGMENT).toBe(0x00000002);
        expect(ConvexMesh3.CFG_POLYGON).toBe(0x00000004);
        expect(ConvexMesh3.CFG_POLYHEDRON).toBe(0x00000008);
    });

    it('are distinct and combinable as a bit mask', () => {
        const flags = [
            ConvexMesh3.CFG_POINT,
            ConvexMesh3.CFG_SEGMENT,
            ConvexMesh3.CFG_POLYGON,
            ConvexMesh3.CFG_POLYHEDRON
        ];
        expect(new Set(flags).size).toBe(4);
        for (const flag of flags) {
            // Each nonempty flag has exactly one bit set.
            expect(flag & (flag - 1)).toBe(0);
            expect(flag & ConvexMesh3.CFG_EMPTY).toBe(0);
        }
        const all = flags.reduce((a, b) => a | b, 0);
        expect(all).toBe(0x0000000F);
    });
});

describe('ConvexMesh3 documented degenerate configurations', () => {
    it('represents a point', () => {
        const mesh = new ConvexMesh3();
        mesh.configuration = ConvexMesh3.CFG_POINT;
        mesh.vertices = [vertex(1, 2, 3)];
        expect(mesh.vertices).toHaveLength(1);
        expect(mesh.triangles).toHaveLength(0);
    });

    it('represents a line segment', () => {
        const mesh = new ConvexMesh3();
        mesh.configuration = ConvexMesh3.CFG_SEGMENT;
        mesh.vertices = [vertex(0, 0, 0), vertex(1, 1, 1)];
        expect(mesh.vertices).toHaveLength(2);
        expect(mesh.triangles).toHaveLength(0);
    });

    it('represents a coplanar convex polygon', () => {
        // A square in the plane z = 0, triangulated into two faces.
        const mesh = new ConvexMesh3();
        mesh.configuration = ConvexMesh3.CFG_POLYGON;
        mesh.vertices = [
            vertex(0, 0, 0), vertex(1, 0, 0), vertex(1, 1, 0), vertex(0, 1, 0)
        ];
        mesh.triangles = [[0, 1, 2], [0, 2, 3]];
        expect(mesh.vertices.length).toBeGreaterThanOrEqual(3);
        expect(mesh.triangles.length).toBeGreaterThan(0);

        // The vertices are coplanar.
        const p0 = mesh.vertices[0];
        const normal = cross(sub(mesh.vertices[1], p0), sub(mesh.vertices[2], p0));
        for (const v of mesh.vertices) {
            expect(dot(normal, sub(v, p0))).toBe(0);
        }

        // Every vertex is used by some triangle (contract item 1).
        const used = new Set<number>();
        for (const t of mesh.triangles) {
            t.forEach(i => used.add(i));
        }
        expect(used.size).toBe(mesh.vertices.length);
    });

    it('represents a convex polyhedron with consistently oriented faces', () => {
        // The tetrahedron with vertices at the origin and the three unit
        // axis points, with all faces counterclockwise viewed from outside.
        const mesh = new ConvexMesh3();
        mesh.configuration = ConvexMesh3.CFG_POLYHEDRON;
        mesh.vertices = [
            vertex(0, 0, 0), vertex(1, 0, 0), vertex(0, 1, 0), vertex(0, 0, 1)
        ];
        mesh.triangles = [
            [1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]
        ];

        // The vertices are not coplanar.
        const p0 = mesh.vertices[0];
        const normal = cross(sub(mesh.vertices[1], p0), sub(mesh.vertices[2], p0));
        expect(dot(normal, sub(mesh.vertices[3], p0))).not.toBe(0);

        // Every vertex is used (contract item 1).
        const used = new Set<number>();
        for (const t of mesh.triangles) {
            t.forEach(i => used.add(i));
        }
        expect(used.size).toBe(4);

        // Same chirality: each face normal points away from the mesh
        // centroid, so all faces are counterclockwise viewed from outside
        // (contract item 2).
        const centroid = vertex(0.25, 0.25, 0.25);
        for (const t of mesh.triangles) {
            const a = mesh.vertices[t[0]];
            const b = mesh.vertices[t[1]];
            const c = mesh.vertices[t[2]];
            const faceNormal = cross(sub(b, a), sub(c, a));
            expect(dot(faceNormal, sub(a, centroid))).toBeGreaterThan(0);
        }

        // Each of the 6 edges is shared by exactly 2 faces, traversed in
        // opposite directions - the closed-surface consistency check.
        const directed = new Map<string, number>();
        for (const t of mesh.triangles) {
            for (let i = 0; i < 3; ++i) {
                const key = `${t[i]},${t[(i + 1) % 3]}`;
                directed.set(key, (directed.get(key) ?? 0) + 1);
            }
        }
        expect(directed.size).toBe(12);
        for (const [key, count] of directed) {
            expect(count).toBe(1);
            const [v0, v1] = key.split(',');
            expect(directed.get(`${v1},${v0}`)).toBe(1);
        }
    });
});

describe('ConvexMesh3 vertex and triangle types', () => {
    it('uses 3D vectors for vertices and index triples for triangles', () => {
        const mesh = new ConvexMesh3();
        const v: ConvexMesh3Vertex = vertex(4, 5, 6);
        const t: ConvexMesh3Triangle = [0, 0, 0];
        mesh.vertices.push(v);
        mesh.triangles.push(t);
        expect(mesh.vertices[0].size).toBe(3);
        expect(mesh.vertices[0].values).toEqual([4, 5, 6]);
        expect(mesh.triangles[0]).toHaveLength(3);
    });
});
