import { describe, it, expect } from 'vitest';
import { TetrahedraRasterizer } from '../src/TetrahedraRasterizer.js';
import { check, fc } from './helpers/arbitraries.js';

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

describe('TetrahedraRasterizer verification', () => {
    const BOUND = 9;
    const bound = [BOUND, BOUND, BOUND];
    const regionMin = [0, 0, 0];
    const regionMax = [BOUND - 1, BOUND - 1, BOUND - 1];
    // With this region the grid multiplier is exactly 1, so grid coordinates
    // equal world coordinates and integer vertices keep every determinant
    // exact in binary64. The comparisons below are therefore exact.

    const det3 = (a: readonly number[], b: readonly number[],
        c: readonly number[]): number =>
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);

    const sub = (u: readonly number[], v: readonly number[]): number[] =>
        [u[0] - v[0], u[1] - v[1], u[2] - v[2]];

    // The upstream containment predicate, written directly as four scalar
    // triple products: the point is contained when it is on the negative side
    // of (or exactly on) all four face planes.
    function contains(P: readonly number[], V: readonly number[][]): boolean {
        const [V0, V1, V2, V3] = V;
        const PmV0 = sub(P, V0);
        const V1mV0 = sub(V1, V0);
        const V2mV0 = sub(V2, V0);
        const V3mV0 = sub(V3, V0);
        if (det3(PmV0, V2mV0, V1mV0) > 0) { return false; }
        if (det3(PmV0, V1mV0, V3mV0) > 0) { return false; }
        if (det3(PmV0, V3mV0, V2mV0) > 0) { return false; }
        const PmV1 = sub(P, V1);
        return det3(PmV1, sub(V2, V1), sub(V3, V1)) <= 0;
    }

    // The scan box upstream derives for tetrahedron t: the vertex bounding
    // box clipped to the region, then transformed to grid coordinates with
    // ceil on the low corner and floor on the high corner.
    function scanBox(V: readonly number[][]): { min: number[], max: number[], valid: boolean } {
        const lo = [0, 1, 2].map(i => Math.min(...V.map(v => v[i])));
        const hi = [0, 1, 2].map(i => Math.max(...V.map(v => v[i])));
        const min = [0, 1, 2].map(i => Math.ceil(Math.max(lo[i], regionMin[i]) - regionMin[i]));
        const max = [0, 1, 2].map(i => Math.floor(Math.min(hi[i], regionMax[i]) - regionMin[i]));
        const valid = [0, 1, 2].every(i =>
            Math.max(lo[i], regionMin[i]) <= Math.min(hi[i], regionMax[i]));
        return { min, max, valid };
    }

    // The grid a faithful reference rasterizer produces: for each grid point,
    // the LAST tetrahedron index whose scan box contains the point and whose
    // containment predicate accepts it, or -1.
    function referenceGrid(vertices: number[][], tetrahedra: number[][]): number[] {
        const grid = new Array<number>(BOUND * BOUND * BOUND).fill(-1);
        tetrahedra.forEach((tet, t) => {
            const V = tet.map(i => vertices[i]);
            const box = scanBox(V);
            if (!box.valid) {
                return;
            }
            for (let z = box.min[2]; z <= box.max[2]; ++z) {
                for (let y = box.min[1]; y <= box.max[1]; ++y) {
                    for (let x = box.min[0]; x <= box.max[0]; ++x) {
                        if (contains([x, y, z], V)) {
                            grid[x + BOUND * (y + BOUND * z)] = t;
                        }
                    }
                }
            }
        });
        return grid;
    }

    const lattice = fc.tuple(fc.integer({ min: 0, max: BOUND - 1 }),
        fc.integer({ min: 0, max: BOUND - 1 }),
        fc.integer({ min: 0, max: BOUND - 1 }));

    it('rasterizes exactly the lattice points the predicate accepts', () => {
        // Arbitrary tetrahedra, including negatively oriented and degenerate
        // ones, cross-checked against a brute-force scan of every grid point.
        const mesh = fc.tuple(
            fc.array(lattice, { minLength: 4, maxLength: 10 }),
            fc.array(fc.tuple(fc.nat({ max: 50 }), fc.nat({ max: 50 }),
                fc.nat({ max: 50 }), fc.nat({ max: 50 })),
                { minLength: 1, maxLength: 3 }))
            .map(([vs, raw]) => ({
                vertices: vs.map(v => [v[0], v[1], v[2]]),
                tetrahedra: raw.map(r => r.map(k => k % vs.length))
            }));
        check(mesh, ({ vertices, tetrahedra }) => {
            const rasterizer = new TetrahedraRasterizer(vertices, tetrahedra);
            const grid = rasterizer.rasterize(regionMin, regionMax, bound);
            const expected = referenceGrid(vertices, tetrahedra);
            expect(Array.from(grid)).toEqual(expected);
        }, 100);
    });

    it('agrees with barycentric containment for oriented tetrahedra', () => {
        // For a positively oriented, nondegenerate tetrahedron the upstream
        // four-plane predicate is exactly "all four barycentric coordinates
        // are nonnegative"; the barycentric solve is an independent
        // computation (Cramer's rule on the edge matrix).
        const oriented = fc.tuple(lattice, lattice, lattice, lattice)
            .map(([a, b, c, d]) => [
                [a[0], a[1], a[2]], [b[0], b[1], b[2]],
                [c[0], c[1], c[2]], [d[0], d[1], d[2]]])
            .filter(V => det3(sub(V[1], V[0]), sub(V[2], V[0]), sub(V[3], V[0])) !== 0)
            .map(V => {
                const det = det3(sub(V[1], V[0]), sub(V[2], V[0]), sub(V[3], V[0]));
                // Swap two vertices to make the orientation positive.
                return det > 0 ? V : [V[0], V[2], V[1], V[3]];
            });
        check(oriented, (V) => {
            const rasterizer = new TetrahedraRasterizer(V, [[0, 1, 2, 3]]);
            const grid = rasterizer.rasterize(regionMin, regionMax, bound);

            const e1 = sub(V[1], V[0]), e2 = sub(V[2], V[0]), e3 = sub(V[3], V[0]);
            const det = det3(e1, e2, e3);
            expect(det).toBeGreaterThan(0);
            for (let z = 0; z < BOUND; ++z) {
                for (let y = 0; y < BOUND; ++y) {
                    for (let x = 0; x < BOUND; ++x) {
                        const p = sub([x, y, z], V[0]);
                        const b1 = det3(p, e2, e3);
                        const b2 = det3(e1, p, e3);
                        const b3 = det3(e1, e2, p);
                        const b0 = det - b1 - b2 - b3;
                        const inside = b0 >= 0 && b1 >= 0 && b2 >= 0 && b3 >= 0;
                        expect(grid[x + BOUND * (y + BOUND * z)])
                            .toBe(inside ? 0 : -1);
                    }
                }
            }
        }, 60);
    });

    it('lets the last tetrahedron win on a shared face', () => {
        // Points on a face satisfy the predicate for both tetrahedra sharing
        // it (the tests are "> 0 rejects"), so upstream marks such a point
        // twice and the later index survives. Two tetrahedra of the unit-cube
        // decomposition share the face (0,0,0),(4,0,0),(0,4,0).
        const vertices = [[0, 0, 0], [4, 0, 0], [0, 4, 0], [0, 0, 4], [0, 0, -4]];
        const shared = [[0, 1, 2, 3], [0, 2, 1, 4]];
        const rasterizer = new TetrahedraRasterizer(vertices, shared);
        const grid = rasterizer.rasterize([0, -4, -4], [8, 4, 4], [9, 9, 9]);
        // The grid point (1,1,0) in world coordinates lies on the shared
        // face; its grid coordinates are (1, 5, 4).
        const at = (x: number, y: number, z: number) => grid[x + 9 * (y + 9 * z)];
        expect(at(1, 5, 4)).toBe(1);
        // A point strictly inside the first tetrahedron keeps index 0.
        expect(at(1, 5, 5)).toBe(0);
        // A point strictly inside the second keeps index 1.
        expect(at(1, 5, 3)).toBe(1);
    });

    it('marks exactly one grid point for a fully degenerate tetrahedron', () => {
        // Four coincident vertices make every scalar triple product zero, so
        // the predicate accepts every point; the scan box is the single
        // vertex, so exactly that grid point is marked.
        check(lattice, ([x, y, z]) => {
            const rasterizer = new TetrahedraRasterizer([[x, y, z]],
                [[0, 0, 0, 0]]);
            const grid = rasterizer.rasterize(regionMin, regionMax, bound);
            let marked = 0;
            for (let i = 0; i < grid.length; ++i) {
                if (grid[i] !== -1) {
                    ++marked;
                    expect(i).toBe(x + BOUND * (y + BOUND * z));
                }
            }
            expect(marked).toBe(1);
        });
    });

    it('is invariant under translating the vertices with the region', () => {
        check(fc.tuple(fc.array(lattice, { minLength: 4, maxLength: 6 }),
            fc.integer({ min: -20, max: 20 })), ([vs, shift]) => {
                const vertices = vs.map(v => [v[0], v[1], v[2]]);
                const tetrahedra = [[0, 1, 2, 3]];
                const a = new TetrahedraRasterizer(vertices, tetrahedra)
                    .rasterize(regionMin, regionMax, bound);
                const shifted = vertices.map(v =>
                    [v[0] + shift, v[1] + shift, v[2] + shift]);
                const b = new TetrahedraRasterizer(shifted, tetrahedra)
                    .rasterize(regionMin.map(v => v + shift),
                        regionMax.map(v => v + shift), bound);
                expect(Array.from(b)).toEqual(Array.from(a));
            });
    });

    it('upstream quirk: rasterize clips the stored bounding boxes in place', () => {
        // ClipCullAABBs writes back into mTetraMin/mTetraMax, which are
        // computed once in the constructor. A second call with a larger
        // region therefore scans the box left over from the first call
        // instead of the tetrahedron's own box. The quirk is preserved; this
        // test pins it so a future change is deliberate.
        const vertices = [[0, 0, 0], [8, 0, 0], [0, 8, 0], [0, 0, 8]];
        const rasterizer = new TetrahedraRasterizer(vertices, [[0, 1, 2, 3]]);

        // First call over a quarter of the tetrahedron.
        const small = rasterizer.rasterize([0, 0, 0], [4, 4, 4], [5, 5, 5]);
        expect(small.some(v => v === 0)).toBe(true);

        // Second call over the whole tetrahedron. A fresh rasterizer marks
        // strictly more points than the reused one, whose bounding box is
        // still clipped to [0,4]^3.
        const reused = rasterizer.rasterize([0, 0, 0], [8, 8, 8], [9, 9, 9]);
        const fresh = new TetrahedraRasterizer(vertices, [[0, 1, 2, 3]])
            .rasterize([0, 0, 0], [8, 8, 8], [9, 9, 9]);
        const count = (g: Int32Array) => g.reduce((n, v) => n + (v === 0 ? 1 : 0), 0);
        expect(count(reused)).toBeLessThan(count(fresh));
        // The reused rasterizer sees only the sub-box the first call left.
        for (let z = 0; z < 9; ++z) {
            for (let y = 0; y < 9; ++y) {
                for (let x = 0; x < 9; ++x) {
                    if (x > 4 || y > 4 || z > 4) {
                        expect(reused[x + 9 * (y + 9 * z)]).toBe(-1);
                    }
                }
            }
        }
    });
});
