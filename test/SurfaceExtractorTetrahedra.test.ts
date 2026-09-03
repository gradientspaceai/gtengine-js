import { describe, expect, it } from 'vitest';
import {
    SurfaceExtractorTetrahedra,
    SurfaceExtractorTetrahedraEdge
} from '../src/SurfaceExtractorTetrahedra.js';
import type { SurfaceExtractorTriangle } from '../src/SurfaceExtractor.js';

// Sample an integer-valued function on a lexicographically ordered grid.
function makeVoxels(xBound: number, yBound: number, zBound: number,
    f: (x: number, y: number, z: number) => number): number[] {
    const voxels: number[] = new Array<number>(xBound * yBound * zBound);
    for (let z = 0; z < zBound; ++z) {
        for (let y = 0; y < yBound; ++y) {
            for (let x = 0; x < xBound; ++x) {
                voxels[x + xBound * (y + yBound * z)] = f(x, y, z);
            }
        }
    }
    return voxels;
}

// The number of times each undirected triangle edge occurs.
function edgeCounts(triangles: readonly SurfaceExtractorTriangle[]):
    Map<string, number> {
    const counts = new Map<string, number>();
    for (const triangle of triangles) {
        for (let i = 0; i < 3; ++i) {
            const a = triangle.v[i];
            const b = triangle.v[(i + 1) % 3];
            const key = (a < b ? `${a}-${b}` : `${b}-${a}`);
            counts.set(key, (counts.get(key) ?? 0) + 1);
        }
    }
    return counts;
}

// Reach the protected getGradient for testing.
function gradientOf(extractor: SurfaceExtractorTetrahedra,
    pos: [number, number, number]): [number, number, number] {
    return (extractor as unknown as {
        getGradient(p: [number, number, number]): [number, number, number];
    }).getGradient(pos);
}

describe('SurfaceExtractorTetrahedraEdge', () => {
    it('stores the vertex indices in increasing order', () => {
        expect(new SurfaceExtractorTetrahedraEdge(3, 1).v).toEqual([1, 3]);
        expect(new SurfaceExtractorTetrahedraEdge(1, 3).v).toEqual([1, 3]);
        expect(new SurfaceExtractorTetrahedraEdge().v).toEqual([0, 0]);
    });

    it('compares as the upstream operator== and operator<', () => {
        const e01 = new SurfaceExtractorTetrahedraEdge(0, 1);
        const e10 = new SurfaceExtractorTetrahedraEdge(1, 0);
        const e02 = new SurfaceExtractorTetrahedraEdge(0, 2);
        const e12 = new SurfaceExtractorTetrahedraEdge(1, 2);

        expect(e01.equals(e10)).toBe(true);
        expect(e01.equals(e02)).toBe(false);
        expect(e01.lessThan(e02)).toBe(true);
        expect(e02.lessThan(e01)).toBe(false);
        expect(e02.lessThan(e12)).toBe(true);
        expect(e01.lessThan(e10)).toBe(false);
        expect(e10.lessThan(e01)).toBe(false);
    });
});

describe('SurfaceExtractorTetrahedra: construction', () => {
    it('rejects bounds smaller than 2', () => {
        expect(() => new SurfaceExtractorTetrahedra(1, 2, 2, [0, 0, 0, 0]))
            .toThrow('Invalid input.');
        expect(() => new SurfaceExtractorTetrahedra(2, 2, 1, [0, 0, 0, 0]))
            .toThrow('Invalid input.');
    });

    it('produces no geometry when the level misses the value range', () => {
        const voxels = makeVoxels(4, 4, 4, () => 7);
        const extractor = new SurfaceExtractorTetrahedra(4, 4, 4, voxels);
        for (const level of [0, 6, 8, 100]) {
            const { vertices, triangles } = extractor.extract(level, true);
            expect(vertices.length).toBe(0);
            expect(triangles.length).toBe(0);
        }
    });
});

describe('SurfaceExtractorTetrahedra: planar level surface', () => {
    // F(x,y,z) = 2x extracted at level 3 gives the plane x = 1.5, which
    // avoids the degenerate case of the level passing through voxel values.
    const xBound = 4, yBound = 3, zBound = 3;
    const voxels = makeVoxels(xBound, yBound, zBound, (x) => 2 * x);
    const extractor = new SurfaceExtractorTetrahedra(xBound, yBound, zBound,
        voxels);

    it('places every vertex exactly on the plane', () => {
        const { vertices, triangles } = extractor.extract(3, true);
        expect(vertices.length).toBeGreaterThan(0);
        expect(triangles.length).toBeGreaterThan(0);
        for (const v of vertices) {
            expect(v[0]).toBe(1.5);
        }
    });

    it('produces the exact rational coordinates', () => {
        const { vertices } = extractor.extractRational(3);
        expect(vertices.length).toBeGreaterThan(0);
        for (const v of vertices) {
            expect(v.xNumer / v.xDenom).toBe(1.5);
            expect(v.xDenom).toBeGreaterThan(0);
            // The tetrahedral decomposition cuts face and box diagonals as
            // well as grid edges, so y and z need not be integers; they stay
            // inside the image, with positive denominators.
            expect(v.yDenom).toBeGreaterThan(0);
            expect(v.zDenom).toBeGreaterThan(0);
            expect(v.yNumer / v.yDenom).toBeGreaterThanOrEqual(0);
            expect(v.yNumer / v.yDenom).toBeLessThanOrEqual(yBound - 1);
            expect(v.zNumer / v.zDenom).toBeGreaterThanOrEqual(0);
            expect(v.zNumer / v.zDenom).toBeLessThanOrEqual(zBound - 1);
        }
    });

    it('the triangulation tiles the plane rectangle without gaps', () => {
        const { vertices, triangles } = extractor.extract(3, true);
        // Each triangle lives in the plane x = 1.5, so the total (unsigned)
        // area equals the area of the rectangle it covers.
        let area = 0;
        for (const t of triangles) {
            const v0 = vertices[t.v[0]];
            const v1 = vertices[t.v[1]];
            const v2 = vertices[t.v[2]];
            area += 0.5 * Math.abs(
                (v1[1] - v0[1]) * (v2[2] - v0[2])
                - (v1[2] - v0[2]) * (v2[1] - v0[1]));
        }
        expect(area).toBeCloseTo((yBound - 1) * (zBound - 1), 12);
    });

    it('vertices are unique and every triangle is non-degenerate', () => {
        const { vertices, triangles } = extractor.extract(3, true);
        const keys = new Set(vertices.map((v) => v.join(',')));
        expect(keys.size).toBe(vertices.length);
        for (const t of triangles) {
            expect(new Set(t.v).size).toBe(3);
        }
    });
});

describe('SurfaceExtractorTetrahedra: degenerate level', () => {
    // F(x,y,z) = x extracted at level 1: the level passes exactly through the
    // voxel values on the plane x = 1, exercising the zero-value cases of
    // processTetrahedron.
    const xBound = 4, yBound = 3, zBound = 3;
    const voxels = makeVoxels(xBound, yBound, zBound, (x) => x);
    const extractor = new SurfaceExtractorTetrahedra(xBound, yBound, zBound,
        voxels);

    it('places every vertex on the grid plane through the level', () => {
        const { vertices, triangles } = extractor.extract(1, true);
        expect(vertices.length).toBeGreaterThan(0);
        expect(triangles.length).toBeGreaterThan(0);
        for (const v of vertices) {
            expect(v[0]).toBe(1);
            // The vertices are grid points of the plane x = 1.
            expect(Number.isInteger(v[1])).toBe(true);
            expect(Number.isInteger(v[2])).toBe(true);
        }
    });

    it('the degenerate triangulation still tiles the plane rectangle', () => {
        const { vertices, triangles } = extractor.extract(1, true);
        let area = 0;
        for (const t of triangles) {
            const v0 = vertices[t.v[0]];
            const v1 = vertices[t.v[1]];
            const v2 = vertices[t.v[2]];
            area += 0.5 * Math.abs(
                (v1[1] - v0[1]) * (v2[2] - v0[2])
                - (v1[2] - v0[2]) * (v2[1] - v0[1]));
        }
        expect(area).toBeCloseTo((yBound - 1) * (zBound - 1), 12);
    });

    it('handles a tetrahedron whose four values are all the level', () => {
        // A constant image extracted at that constant value drives every
        // tetrahedron into the 0000 case, which emits all four faces.
        const constant = makeVoxels(2, 2, 2, () => 5);
        const flat = new SurfaceExtractorTetrahedra(2, 2, 2, constant);
        const { vertices, triangles } = flat.extract(5, true);
        // The only vertices available are the 8 voxel corners.
        expect(vertices.length).toBe(8);
        expect(triangles.length).toBeGreaterThan(0);
        for (const v of vertices) {
            expect(v.every((t) => t === 0 || t === 1)).toBe(true);
        }
        // No triangle appears with both windings.
        const seen = new Set<string>();
        for (const t of triangles) {
            const key = t.v.join(',');
            const reversed = [t.v[0], t.v[2], t.v[1]].join(',');
            expect(seen.has(key)).toBe(false);
            expect(seen.has(reversed)).toBe(false);
            seen.add(key);
        }
    });
});

describe('SurfaceExtractorTetrahedra: spherical level surface', () => {
    const c = 5;
    const bound = 11;
    const r2 = 10;
    // 4*(r2 - |P - C|^2) is a multiple of 4, so the odd level 1 is never a
    // voxel value and no tetrahedron has a zero corner value.
    const voxels = makeVoxels(bound, bound, bound, (x, y, z) =>
        4 * (r2 - ((x - c) * (x - c) + (y - c) * (y - c) + (z - c) * (z - c))));
    const extractor = new SurfaceExtractorTetrahedra(bound, bound, bound,
        voxels);

    it('places vertices within half a voxel of the true sphere', () => {
        const { vertices, triangles } = extractor.extract(1, true);
        expect(vertices.length).toBeGreaterThan(0);
        expect(triangles.length).toBeGreaterThan(0);

        const radius = Math.sqrt(r2 - 0.25);
        for (const v of vertices) {
            const dx = v[0] - c;
            const dy = v[1] - c;
            const dz = v[2] - c;
            const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
            expect(Math.abs(r - radius)).toBeLessThan(0.5);
        }
    });

    it('extracts a closed manifold surface', () => {
        const { triangles } = extractor.extract(1, true);
        for (const count of edgeCounts(triangles).values()) {
            expect(count).toBe(2);
        }
    });

    it('the reported edges are exactly the triangle edges', () => {
        const { triangles } = extractor.extract(1, true);
        const edges = extractor.getEdges();
        // The edge set is built before the duplicate vertices are welded, so
        // it has at least as many edges as the welded triangulation.
        expect(edges.length).toBeGreaterThanOrEqual(
            edgeCounts(triangles).size);

        // getEdges() returns them ordered by the upstream comparison.
        for (let i = 1; i < edges.length; ++i) {
            expect(edges[i - 1].lessThan(edges[i])).toBe(true);
        }
    });

    it('the triangles are returned in the upstream sorted order', () => {
        const { triangles } = extractor.extractRational(1);
        for (let i = 1; i < triangles.length; ++i) {
            expect(triangles[i - 1].lessThan(triangles[i])).toBe(true);
        }
    });

    it('nested levels give nested surfaces', () => {
        const maxRadius = (vertices: [number, number, number][]) => {
            let m = 0;
            for (const v of vertices) {
                const dx = v[0] - c;
                const dy = v[1] - c;
                const dz = v[2] - c;
                m = Math.max(m, Math.sqrt(dx * dx + dy * dy + dz * dz));
            }
            return m;
        };
        const inner = extractor.extract(21, true);
        const outer = extractor.extract(1, true);
        expect(maxRadius(inner.vertices)).toBeLessThan(
            maxRadius(outer.vertices));
    });

    it('orients the triangles consistently with the gradient', () => {
        const { vertices, triangles } = extractor.extract(1, true);
        const signedVolume = () => {
            let sum = 0;
            for (const t of triangles) {
                const v0 = vertices[t.v[0]];
                const v1 = vertices[t.v[1]];
                const v2 = vertices[t.v[2]];
                const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
                const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
                const n = [
                    e1[1] * e2[2] - e1[2] * e2[1],
                    e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0]
                ];
                sum += (v0[0] - c) * n[0] + (v0[1] - c) * n[1]
                    + (v0[2] - c) * n[2];
            }
            return sum;
        };

        // The field decreases outward, so its gradient points inward.
        extractor.orientTriangles(vertices, triangles, true);
        expect(signedVolume()).toBeLessThan(0);
        extractor.orientTriangles(vertices, triangles, false);
        expect(signedVolume()).toBeGreaterThan(0);

        const normals = extractor.computeNormals(vertices, triangles);
        for (let i = 0; i < normals.length; ++i) {
            const n = normals[i];
            expect(Math.sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]))
                .toBeCloseTo(1, 10);
            const v = vertices[i];
            const g = [v[0] - c, v[1] - c, v[2] - c];
            const glen = Math.sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
            expect((g[0] * n[0] + g[1] * n[1] + g[2] * n[2]) / glen)
                .toBeGreaterThan(0.5);
        }
    });
});

describe('SurfaceExtractorTetrahedra: gradient', () => {
    it('is exact for a linear field', () => {
        const voxels = makeVoxels(4, 4, 4,
            (x, y, z) => 2 * x + 3 * y + 5 * z);
        const extractor = new SurfaceExtractorTetrahedra(4, 4, 4, voxels);
        extractor.extractRational(0);
        // A linear field has the same gradient in every tetrahedron of every
        // partition, so all five branches of both parities agree.
        for (const p of [[0.5, 0.5, 0.5], [0.1, 0.1, 0.1], [0.9, 0.9, 0.9],
            [1.5, 0.5, 0.5], [1.9, 0.05, 0.05], [1.9, 0.9, 0.9],
            [2.5, 1.5, 0.5]] as [number, number, number][]) {
            const g = gradientOf(extractor, p);
            expect(g[0]).toBeCloseTo(2, 12);
            expect(g[1]).toBeCloseTo(3, 12);
            expect(g[2]).toBeCloseTo(5, 12);
        }
    });

    it('returns zero outside the voxels of the image', () => {
        const voxels = makeVoxels(4, 4, 4, (x, y, z) => x + y + z);
        const extractor = new SurfaceExtractorTetrahedra(4, 4, 4, voxels);
        extractor.extractRational(0);
        expect(gradientOf(extractor, [-1, 0.5, 0.5])).toEqual([0, 0, 0]);
        expect(gradientOf(extractor, [0.5, -1, 0.5])).toEqual([0, 0, 0]);
        expect(gradientOf(extractor, [0.5, 0.5, -1])).toEqual([0, 0, 0]);
        expect(gradientOf(extractor, [3, 0.5, 0.5])).toEqual([0, 0, 0]);
        expect(gradientOf(extractor, [0.5, 3, 0.5])).toEqual([0, 0, 0]);
        expect(gradientOf(extractor, [0.5, 0.5, 3])).toEqual([0, 0, 0]);
    });

    it('selects the correct tetrahedron in an odd-parity voxel', () => {
        // The voxel with origin (1,0,0) has odd parity. Its corner values are
        // chosen so that the five tetrahedra of the odd-parity partition give
        // five distinguishable gradients.
        //   f000 = v(1,0,0) = 0    f100 = v(2,0,0) = 3
        //   f010 = v(1,1,0) = 1    f110 = v(2,1,0) = 2
        //   f001 = v(1,0,1) = 1    f101 = v(2,0,1) = 0
        //   f011 = v(1,1,1) = 0    f111 = v(2,1,1) = 10
        const xBound = 3, yBound = 2, zBound = 2;
        const value = new Map<string, number>([
            ['1,0,0', 0], ['2,0,0', 3], ['1,1,0', 1], ['2,1,0', 2],
            ['1,0,1', 1], ['2,0,1', 0], ['1,1,1', 0], ['2,1,1', 10],
            ['0,0,0', 0], ['0,1,0', 0], ['0,0,1', 0], ['0,1,1', 0]
        ]);
        const voxels = makeVoxels(xBound, yBound, zBound,
            (x, y, z) => value.get(`${x},${y},${z}`) as number);
        const extractor = new SurfaceExtractorTetrahedra(xBound, yBound,
            zBound, voxels);
        extractor.extractRational(0);

        // Corner tetrahedron 1205 at the (1,0,0) corner of the voxel:
        //   (f100 - f000, f110 - f100, f101 - f100) = (3, -1, -3).
        expect(gradientOf(extractor, [1.9, 0.05, 0.05]))
            .toEqual([3, -1, -3]);

        // Corner tetrahedron 6572 at the (1,1,1) corner of the voxel:
        //   (f111 - f011, f111 - f101, f111 - f110) = (10, 10, 8).
        expect(gradientOf(extractor, [1.9, 0.9, 0.9])).toEqual([10, 10, 8]);

        // The voxel center is in the central tetrahedron 0752, whose gradient
        // is the average of the four differences. Upstream tests
        // 'dx + dy + dz >= 0' for the 6572 branch, which is always true, so
        // this case reached 6572 and returned (10, 10, 8) instead.
        expect(gradientOf(extractor, [1.5, 0.5, 0.5])).toEqual([1, 1, -1]);
    });

    it('selects the correct tetrahedron in an even-parity voxel', () => {
        // The voxel with origin (0,0,0) has even parity. Its central
        // tetrahedron 6314 covers the voxel center.
        //   f000 = 0, f100 = 4, f010 = 2, f110 = 0
        //   f001 = 6, f101 = 0, f011 = 0, f111 = 8
        const value = new Map<string, number>([
            ['0,0,0', 0], ['1,0,0', 4], ['0,1,0', 2], ['1,1,0', 0],
            ['0,0,1', 6], ['1,0,1', 0], ['0,1,1', 0], ['1,1,1', 8]
        ]);
        const voxels = makeVoxels(2, 2, 2,
            (x, y, z) => value.get(`${x},${y},${z}`) as number);
        const extractor = new SurfaceExtractorTetrahedra(2, 2, 2, voxels);
        extractor.extractRational(0);

        // Corner tetrahedron 0134: (f100 - f000, f010 - f000, f001 - f000).
        expect(gradientOf(extractor, [0.05, 0.05, 0.05])).toEqual([4, 2, 6]);

        // Central tetrahedron 6314 at the voxel center:
        //   0.5 * (f111 - f010 + f100 - f001) = 0.5 * (8 - 2 + 4 - 6) = 2
        //   0.5 * (f111 + f010 - f100 - f001) = 0.5 * (8 + 2 - 4 - 6) = 0
        //   0.5 * (f111 - f010 - f100 + f001) = 0.5 * (8 - 2 - 4 + 6) = 4
        expect(gradientOf(extractor, [0.5, 0.5, 0.5])).toEqual([2, 0, 4]);
    });
});
