import { describe, it, expect } from 'vitest';
import { TetrahedraRasterizer } from '../src/TetrahedraRasterizer.js';

// The canonical orientation: {v0,v1,v2,v3} with v0=(0,0,0), v1=(1,0,0),
// v2=(0,1,0), v3=(0,0,1) (positive determinant of edge vectors).

describe('TetrahedraRasterizer', () => {
    it('validates constructor and rasterize arguments', () => {
        const vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]];
        const tetrahedra = [[0, 1, 2, 3]];
        expect(() => new TetrahedraRasterizer([], tetrahedra)).toThrow('Invalid argument.');
        expect(() => new TetrahedraRasterizer(vertices, [])).toThrow('Invalid argument.');
        const rasterizer = new TetrahedraRasterizer(vertices, tetrahedra);
        expect(() => rasterizer.rasterize([0, 0, 0], [1, 1, 1], [1, 2, 2]))
            .toThrow('Invalid argument.');
        expect(() => rasterizer.rasterize([0, 0, 0], [1, 1, 1], [2, 2, 1]))
            .toThrow('Invalid argument.');
    });

    it('rasterizes a single known tetrahedron; boundary points are inside', () => {
        // Scaled canonical tetrahedron with vertices at lattice points. On
        // the grid [0,4]^3 with bound 5, grid coordinates equal vertex
        // coordinates, and (x,y,z) is contained iff x+y+z <= 4 (the boundary
        // planes are inclusive).
        const vertices = [[0, 0, 0], [4, 0, 0], [0, 4, 0], [0, 0, 4]];
        const tetrahedra = [[0, 1, 2, 3]];
        const rasterizer = new TetrahedraRasterizer(vertices, tetrahedra);
        const bound = [5, 5, 5];
        const grid = rasterizer.rasterize([0, 0, 0], [4, 4, 4], bound);

        expect(grid.length).toBe(125);
        for (let z = 0; z < 5; ++z) {
            for (let y = 0; y < 5; ++y) {
                for (let x = 0; x < 5; ++x) {
                    const i = x + bound[0] * (y + bound[1] * z);
                    const expected = (x + y + z <= 4) ? 0 : -1;
                    expect(grid[i]).toBe(expected);
                }
            }
        }

        // Spot-check boundary semantics explicitly: vertices, an edge
        // midpoint and a point interior to the slanted face are all marked.
        expect(grid[0 + 5 * (0 + 5 * 0)]).toBe(0);   // vertex (0,0,0)
        expect(grid[4 + 5 * (0 + 5 * 0)]).toBe(0);   // vertex (4,0,0)
        expect(grid[2 + 5 * (2 + 5 * 0)]).toBe(0);   // edge midpoint (2,2,0)
        expect(grid[1 + 5 * (1 + 5 * 2)]).toBe(0);   // slant face (1,1,2)
        expect(grid[2 + 5 * (2 + 5 * 2)]).toBe(-1);  // just outside (2,2,2)
        expect(grid[4 + 5 * (4 + 5 * 4)]).toBe(-1);  // far corner
    });

    it('applies the region-to-grid transformation', () => {
        // Unit canonical tetrahedron on region [0,2]^3 with bound 5: the
        // multiplier is (5-1)/2 = 2, so grid point (x,y,z) maps to world
        // (x/2,y/2,z/2) and is contained iff x+y+z <= 2.
        const vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]];
        const tetrahedra = [[0, 1, 2, 3]];
        const rasterizer = new TetrahedraRasterizer(vertices, tetrahedra);
        const bound = [5, 5, 5];
        const grid = rasterizer.rasterize([0, 0, 0], [2, 2, 2], bound);

        for (let z = 0; z < 5; ++z) {
            for (let y = 0; y < 5; ++y) {
                for (let x = 0; x < 5; ++x) {
                    const i = x + bound[0] * (y + bound[1] * z);
                    const expected = (x + y + z <= 2) ? 0 : -1;
                    expect(grid[i]).toBe(expected);
                }
            }
        }
    });

    it('later tetrahedra overwrite earlier ones in overlap regions', () => {
        // Two overlapping tetrahedra: the second is the first translated by
        // (1,0,0). Their intersection is nonempty; a contained point gets
        // the larger index because tetrahedra are rasterized in order.
        const vertices = [
            [0, 0, 0], [4, 0, 0], [0, 4, 0], [0, 0, 4],
            [1, 0, 0], [5, 0, 0], [1, 4, 0], [1, 0, 4]
        ];
        const tetrahedra = [[0, 1, 2, 3], [4, 5, 6, 7]];
        const rasterizer = new TetrahedraRasterizer(vertices, tetrahedra);
        const bound = [6, 6, 6];
        const grid = rasterizer.rasterize([0, 0, 0], [5, 5, 5], bound);

        const at = (x: number, y: number, z: number) =>
            grid[x + bound[0] * (y + bound[1] * z)];
        // (0,y,z) can only be in tetrahedron 0.
        expect(at(0, 0, 0)).toBe(0);
        expect(at(0, 2, 1)).toBe(0);
        // (2,1,1): in both (2+1+1 <= 4 and translated 1+1+1 <= 4): index 1.
        expect(at(2, 1, 1)).toBe(1);
        // (5,0,0): only in tetrahedron 1.
        expect(at(5, 0, 0)).toBe(1);
        // Outside both.
        expect(at(5, 5, 5)).toBe(-1);
        expect(at(0, 4, 4)).toBe(-1);
    });

    it('culls tetrahedra whose bounding boxes miss the region', () => {
        const vertices = [[10, 10, 10], [11, 10, 10], [10, 11, 10], [10, 10, 11]];
        const tetrahedra = [[0, 1, 2, 3]];
        const rasterizer = new TetrahedraRasterizer(vertices, tetrahedra);
        const grid = rasterizer.rasterize([0, 0, 0], [4, 4, 4], [5, 5, 5]);
        for (let i = 0; i < grid.length; ++i) {
            expect(grid[i]).toBe(-1);
        }
    });

    it('matches an independent barycentric containment test', () => {
        // A generic (positively oriented) tetrahedron with non-lattice
        // vertices so no grid point lies exactly on a face plane.
        const v0 = [0.3, 0.2, 0.1];
        const v1 = [6.7, 0.4, 0.35];
        const v2 = [0.55, 7.1, 0.6];
        const v3 = [0.8, 0.9, 6.3];
        const rasterizer = new TetrahedraRasterizer([v0, v1, v2, v3], [[0, 1, 2, 3]]);
        const bound = [9, 9, 9];
        // Region [0,8]^3 with bound 9: grid coordinates = world coordinates.
        const grid = rasterizer.rasterize([0, 0, 0], [8, 8, 8], bound);

        // Independent check: solve for barycentric coordinates by Cramer's
        // rule; contained iff all four coordinates are nonnegative.
        const det3 = (a: number[], b: number[], c: number[]) =>
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
        const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        const e3 = [v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]];
        const det = det3(e1, e2, e3);
        expect(det).toBeGreaterThan(0);  // positively oriented

        let insideCount = 0;
        for (let z = 0; z < 9; ++z) {
            for (let y = 0; y < 9; ++y) {
                for (let x = 0; x < 9; ++x) {
                    const p = [x - v0[0], y - v0[1], z - v0[2]];
                    const b1 = det3(p, e2, e3) / det;
                    const b2 = det3(e1, p, e3) / det;
                    const b3 = det3(e1, e2, p) / det;
                    const b0 = 1 - b1 - b2 - b3;
                    const inside = b0 >= 0 && b1 >= 0 && b2 >= 0 && b3 >= 0;
                    if (inside) {
                        ++insideCount;
                    }
                    const i = x + bound[0] * (y + bound[1] * z);
                    expect(grid[i]).toBe(inside ? 0 : -1);
                }
            }
        }
        // Sanity: the comparison covered both classes.
        expect(insideCount).toBeGreaterThan(10);
        expect(insideCount).toBeLessThan(grid.length);
    });

    it('marks nothing for a negatively oriented tetrahedron', () => {
        // The containment test assumes the canonical (positive) orientation;
        // swapping two vertices inverts all face-plane signs, so no grid
        // point is classified as contained. This documents the upstream
        // orientation requirement.
        const vertices = [[0, 0, 0], [4, 0, 0], [0, 4, 0], [0, 0, 4]];
        const rasterizer = new TetrahedraRasterizer(vertices, [[1, 0, 2, 3]]);
        const grid = rasterizer.rasterize([0, 0, 0], [4, 4, 4], [5, 5, 5]);
        let marked = 0;
        for (let i = 0; i < grid.length; ++i) {
            if (grid[i] !== -1) {
                ++marked;
            }
        }
        expect(marked).toBe(0);
    });
});
