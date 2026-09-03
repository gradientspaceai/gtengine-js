import { describe, expect, it } from 'vitest';
import { Image3 } from '../src/Image3.js';
import { MarchingCubes } from '../src/MarchingCubes.js';
import { SurfaceExtractorMC } from '../src/SurfaceExtractorMC.js';
import type { Vector } from '../src/Vector.js';

// Build an Image3<number> of the given dimensions from an analytic function
// evaluated at the voxel coordinates.
function makeImage(d0: number, d1: number, d2: number,
    f: (x: number, y: number, z: number) => number): Image3<number> {
    const image = new Image3<number>(d0, d1, d2, () => 0);
    for (let z = 0; z < d2; ++z) {
        for (let y = 0; y < d1; ++y) {
            for (let x = 0; x < d0; ++x) {
                image.set(x, y, z, f(x, y, z));
            }
        }
    }
    return image;
}

// The 8 voxel corners in local coordinates, k -> (k&1, (k&2)>>1, (k&4)>>2).
function corner(k: number): [number, number, number] {
    return [k & 1, (k & 2) >> 1, (k & 4) >> 2];
}

// The number of bits set in a byte.
function popCount(n: number): number {
    let count = 0;
    for (let i = 0; i < 8; ++i) {
        count += (n >> i) & 1;
    }
    return count;
}

describe('SurfaceExtractorMC: table structure', () => {
    const extractor = new SurfaceExtractorMC(makeImage(2, 2, 2, () => 1));

    it('the derived class exposes the inherited 256-entry table', () => {
        expect(extractor).toBeInstanceOf(MarchingCubes);
        for (let entry = 0; entry < 256; ++entry) {
            expect(extractor.getTable(entry)).toBe(
                extractor.getTable(entry));
        }
    });

    it('every entry has a well-formed topology', () => {
        for (let entry = 0; entry < 256; ++entry) {
            const topology = extractor.getTable(entry);
            expect(topology.numVertices).toBeGreaterThanOrEqual(0);
            expect(topology.numVertices).toBeLessThanOrEqual(
                MarchingCubes.maxVertices);
            expect(topology.numTriangles).toBeGreaterThanOrEqual(0);
            expect(topology.numTriangles).toBeLessThanOrEqual(
                MarchingCubes.maxTriangles);

            // Every vertex is on a voxel edge: the two corner indices differ
            // in exactly one bit, and j0 < j1 (relied on by the extraction so
            // that voxels sharing an edge generate the same vertex).
            for (let i = 0; i < topology.numVertices; ++i) {
                const j0 = topology.vpair[i][0];
                const j1 = topology.vpair[i][1];
                expect(j0).toBeGreaterThanOrEqual(0);
                expect(j1).toBeLessThan(8);
                expect(j0).toBeLessThan(j1);
                expect(popCount(j0 ^ j1)).toBe(1);
            }

            // Every triangle references distinct existing vertices.
            for (let t = 0; t < topology.numTriangles; ++t) {
                const triple = topology.itriple[t];
                for (const index of triple) {
                    expect(index).toBeGreaterThanOrEqual(0);
                    expect(index).toBeLessThan(topology.numVertices);
                }
                expect(new Set(triple).size).toBe(3);
            }
        }
    });

    it('the vertex/triangle counts match the documented cases', () => {
        // The all-positive and all-negative voxels produce no geometry.
        expect(extractor.getTable(0).numVertices).toBe(0);
        expect(extractor.getTable(0).numTriangles).toBe(0);
        expect(extractor.getTable(255).numVertices).toBe(0);
        expect(extractor.getTable(255).numTriangles).toBe(0);

        // A single sign bit (or a single missing sign bit) cuts one corner
        // off with one triangle over three edges.
        for (let bit = 0; bit < 8; ++bit) {
            const one = 1 << bit;
            expect(extractor.getTable(one).numVertices).toBe(3);
            expect(extractor.getTable(one).numTriangles).toBe(1);
            expect(extractor.getTable(255 - one).numVertices).toBe(3);
            expect(extractor.getTable(255 - one).numTriangles).toBe(1);
        }

        // Each triangle uses 3 vertices, and every vertex is used at least
        // once whenever there is any geometry at all.
        for (let entry = 0; entry < 256; ++entry) {
            const topology = extractor.getTable(entry);
            const used = new Set<number>();
            for (let t = 0; t < topology.numTriangles; ++t) {
                for (const index of topology.itriple[t]) {
                    used.add(index);
                }
            }
            if (topology.numTriangles > 0) {
                expect(used.size).toBe(topology.numVertices);
            }
        }
    });

    it('complementary entries generate the same set of cut edges', () => {
        // Negating the image negates every corner value, so the level surface
        // is unchanged as a point set: entry and its complement must cut the
        // same voxel edges (only the triangle winding may differ).
        const key = (topology: { numVertices: number, vpair: number[][] }) => {
            const pairs: string[] = [];
            for (let i = 0; i < topology.numVertices; ++i) {
                pairs.push(`${topology.vpair[i][0]}-${topology.vpair[i][1]}`);
            }
            pairs.sort();
            return pairs.join(',');
        };
        for (let entry = 0; entry < 256; ++entry) {
            expect(key(extractor.getTable(entry)))
                .toBe(key(extractor.getTable(255 - entry)));
        }
    });

    it('the cut edges are exactly the sign-changing voxel edges', () => {
        for (let entry = 0; entry < 256; ++entry) {
            // Corner j is negative iff bit j of entry is set.
            const expected: string[] = [];
            for (let j0 = 0; j0 < 8; ++j0) {
                for (let j1 = j0 + 1; j1 < 8; ++j1) {
                    if (popCount(j0 ^ j1) !== 1) {
                        continue;
                    }
                    const n0 = (entry >> j0) & 1;
                    const n1 = (entry >> j1) & 1;
                    if (n0 !== n1) {
                        expected.push(`${j0}-${j1}`);
                    }
                }
            }
            expected.sort();

            const topology = extractor.getTable(entry);
            const actual: string[] = [];
            for (let i = 0; i < topology.numVertices; ++i) {
                actual.push(`${topology.vpair[i][0]}-${topology.vpair[i][1]}`);
            }
            // The vertices are distinct edges.
            expect(new Set(actual).size).toBe(actual.length);
            actual.sort();
            expect(actual).toEqual(expected);
        }
    });
});

describe('SurfaceExtractorMC: single-voxel extraction', () => {
    const extractor = new SurfaceExtractorMC(makeImage(2, 2, 2, () => 1));

    it('cuts a single negative corner with one triangle', () => {
        // Corner 0 is below the level, the other 7 are above.
        const F = [-1, 1, 1, 1, 1, 1, 1, 1];
        const { valid, mesh } = extractor.extractVoxel(0, 0, F);
        expect(valid).toBe(true);
        expect(mesh.topology.numVertices).toBe(3);
        expect(mesh.topology.numTriangles).toBe(1);

        // The three vertices are the midpoints of the edges from corner 0.
        const positions = [];
        for (let i = 0; i < 3; ++i) {
            positions.push(mesh.vertices[i].values.join(','));
        }
        positions.sort();
        expect(positions).toEqual(['0,0,0.5', '0,0.5,0', '0.5,0,0']);
    });

    it('places vertices by linear interpolation along the edge', () => {
        // F(x,y,z) = 4x - 1, so the surface is x = 0.25 and only the four
        // edges parallel to x are cut.
        const F: number[] = new Array<number>(8).fill(0);
        for (let k = 0; k < 8; ++k) {
            F[k] = 4 * corner(k)[0] - 1;
        }
        const { valid, mesh } = extractor.extractVoxel(0, 0, F);
        expect(valid).toBe(true);
        expect(mesh.topology.numVertices).toBe(4);
        expect(mesh.topology.numTriangles).toBe(2);
        for (let i = 0; i < 4; ++i) {
            expect(mesh.vertices[i].values[0]).toBeCloseTo(0.25, 12);
        }
    });

    it('reports no geometry when a corner equals the level and perturb is 0', () => {
        const F = [0, 1, 1, 1, 1, 1, 1, 1];
        const { valid, mesh } = extractor.extractVoxel(0, 0, F);
        expect(valid).toBe(false);
        expect(mesh.topology.numVertices).toBe(0);
        expect(mesh.topology.numTriangles).toBe(0);
    });

    it('produces geometry for the same voxel when perturb is nonzero', () => {
        const F = [0, 1, 1, 1, 1, 1, 1, 1];
        const { valid, mesh } = extractor.extractVoxel(0, -0.5, F);
        expect(valid).toBe(true);
        // The perturbation makes corner 0 negative, so one corner is cut.
        expect(mesh.topology.numVertices).toBe(3);
        expect(mesh.topology.numTriangles).toBe(1);
    });

    it('a positive perturbation leaves the corner positive', () => {
        const F = [0, 1, 1, 1, 1, 1, 1, 1];
        const { valid, mesh } = extractor.extractVoxel(0, 0.5, F);
        expect(valid).toBe(true);
        expect(mesh.topology.numVertices).toBe(0);
        expect(mesh.topology.numTriangles).toBe(0);
    });

    it('agrees with the table for every sign configuration', () => {
        for (let entry = 0; entry < 256; ++entry) {
            const F: number[] = new Array<number>(8).fill(0);
            for (let k = 0; k < 8; ++k) {
                F[k] = ((entry >> k) & 1) ? -1 : 1;
            }
            const { valid, mesh } = extractor.extractVoxel(0, 0, F);
            expect(valid).toBe(true);
            expect(mesh.topology).toBe(extractor.getTable(entry));

            // With +/-1 corner values the vertices are edge midpoints.
            for (let i = 0; i < mesh.topology.numVertices; ++i) {
                const j0 = mesh.topology.vpair[i][0];
                const j1 = mesh.topology.vpair[i][1];
                const c0 = corner(j0);
                const c1 = corner(j1);
                for (let d = 0; d < 3; ++d) {
                    expect(mesh.vertices[i].values[d]).toBeCloseTo(
                        0.5 * (c0[d] + c1[d]), 12);
                }
            }
        }
    });
});

describe('SurfaceExtractorMC: image extraction', () => {
    it('extracts a flat surface from a linear field', () => {
        // F = 2z - 3 vanishes at z = 1.5, which is not a voxel value.
        const image = makeImage(4, 4, 4, (_x, _y, z) => 2 * z - 3);
        const extractor = new SurfaceExtractorMC(image);
        const { vertices, indices } = extractor.extract(0, 0);
        expect(vertices.length).toBeGreaterThan(0);
        expect(indices.length % 3).toBe(0);
        for (const v of vertices) {
            expect(v.values[2]).toBeCloseTo(1.5, 12);
        }

        // The surface covers the whole 3x3 slab of voxels, each of which
        // produces 2 triangles.
        expect(indices.length / 3).toBe(2 * 3 * 3);
    });

    it('extracts a sphere with vertices near the true radius', () => {
        const c = 5;
        const radius = 3.25;
        const image = makeImage(11, 11, 11, (x, y, z) =>
            (x - c) * (x - c) + (y - c) * (y - c) + (z - c) * (z - c)
            - radius * radius);
        const extractor = new SurfaceExtractorMC(image);
        const { vertices, indices } = extractor.extract(0, 0);
        expect(vertices.length).toBeGreaterThan(0);
        expect(indices.length).toBeGreaterThan(0);

        let minR = Number.MAX_VALUE;
        let maxR = 0;
        for (const v of vertices) {
            const dx = v.values[0] - c;
            const dy = v.values[1] - c;
            const dz = v.values[2] - c;
            const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
            minR = Math.min(minR, r);
            maxR = Math.max(maxR, r);
        }
        // Linear interpolation of a quadratic along voxel edges is exact at
        // the endpoints only, so the radii bracket the true radius within
        // half a voxel.
        expect(minR).toBeGreaterThan(radius - 0.5);
        expect(maxR).toBeLessThan(radius + 0.5);
    });

    it('makeUnique welds the shared edge vertices into a closed manifold', () => {
        const c = 5;
        const radius = 3.25;
        const image = makeImage(11, 11, 11, (x, y, z) =>
            (x - c) * (x - c) + (y - c) * (y - c) + (z - c) * (z - c)
            - radius * radius);
        const extractor = new SurfaceExtractorMC(image);
        const raw = extractor.extract(0, 0);
        // The raw extraction stores one vertex per cut edge per voxel, so
        // vertices on edges shared by several voxels are duplicated.
        expect(raw.vertices.length).toBeGreaterThan(0);

        const unique = extractor.makeUnique(raw.vertices, raw.indices);
        expect(unique.indices.length).toBe(raw.indices.length);
        expect(unique.vertices.length).toBeLessThan(raw.vertices.length);

        // Each vertex position occurs once.
        const keys = new Set(unique.vertices.map(
            (v: Vector) => v.values.join(',')));
        expect(keys.size).toBe(unique.vertices.length);

        // Every undirected edge is shared by exactly two triangles: the
        // extracted sphere is a closed surface.
        const edgeCount = new Map<string, number>();
        for (let t = 0; t < unique.indices.length; t += 3) {
            const tri = [unique.indices[t], unique.indices[t + 1],
                unique.indices[t + 2]];
            for (let i = 0; i < 3; ++i) {
                const a = tri[i];
                const b = tri[(i + 1) % 3];
                const key = (a < b ? `${a}-${b}` : `${b}-${a}`);
                edgeCount.set(key, (edgeCount.get(key) ?? 0) + 1);
            }
        }
        for (const count of edgeCount.values()) {
            expect(count).toBe(2);
        }
    });

    it('skips voxels whose corners equal the level when perturb is 0', () => {
        // The level 0 is exactly a voxel value on the plane z = 1, so every
        // voxel touching that plane is skipped and no geometry survives.
        const image = makeImage(4, 4, 4, (_x, _y, z) => z - 1);
        const extractor = new SurfaceExtractorMC(image);
        const zero = extractor.extract(0, 0);
        expect(zero.vertices.length).toBe(0);
        expect(zero.indices.length).toBe(0);

        // A perturbation smaller than the minimum nonzero |F - level| = 1
        // pushes the ambiguous corners off the level.
        const perturbed = extractor.extract(0, -0.25);
        expect(perturbed.indices.length / 3).toBe(2 * 3 * 3);
        for (const v of perturbed.vertices) {
            expect(v.values[2]).toBeCloseTo(1, 12);
        }
    });
});

describe('SurfaceExtractorMC: normals and orientation', () => {
    const c = 5;
    const radius = 3.25;
    const image = makeImage(11, 11, 11, (x, y, z) =>
        (x - c) * (x - c) + (y - c) * (y - c) + (z - c) * (z - c)
        - radius * radius);
    const extractor = new SurfaceExtractorMC(image);

    it('orients triangles with (or against) the image gradient', () => {
        const raw = extractor.extract(0, 0);
        const { vertices, indices } = extractor.makeUnique(
            raw.vertices, raw.indices);

        const triangleDot = (idx: readonly number[], t: number) => {
            const v0 = vertices[idx[t]].values;
            const v1 = vertices[idx[t + 1]].values;
            const v2 = vertices[idx[t + 2]].values;
            const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            const n = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0]
            ];
            // The gradient of the sampled field points away from the center.
            const g = [v0[0] - c, v0[1] - c, v0[2] - c];
            return g[0] * n[0] + g[1] * n[1] + g[2] * n[2];
        };

        const sameDir = indices.slice();
        extractor.orientTriangles(vertices, sameDir, true);
        for (let t = 0; t < sameDir.length; t += 3) {
            expect(triangleDot(sameDir, t)).toBeGreaterThanOrEqual(0);
        }

        const oppDir = indices.slice();
        extractor.orientTriangles(vertices, oppDir, false);
        for (let t = 0; t < oppDir.length; t += 3) {
            expect(triangleDot(oppDir, t)).toBeLessThanOrEqual(0);
        }

        // Reorienting is an involution on the triangle as a point set.
        for (let t = 0; t < indices.length; t += 3) {
            expect(new Set([sameDir[t], sameDir[t + 1], sameDir[t + 2]]))
                .toEqual(new Set([indices[t], indices[t + 1], indices[t + 2]]));
        }
    });

    it('computes unit-length vertex normals pointing radially', () => {
        const raw = extractor.extract(0, 0);
        const { vertices, indices } = extractor.makeUnique(
            raw.vertices, raw.indices);
        const oriented = indices.slice();
        extractor.orientTriangles(vertices, oriented, true);

        const normals = extractor.computeNormals(vertices, oriented);
        expect(normals.length).toBe(vertices.length);
        for (let i = 0; i < normals.length; ++i) {
            const n = normals[i].values;
            const len = Math.sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            expect(len).toBeCloseTo(1, 10);

            // The normal roughly agrees with the outward radial direction.
            const v = vertices[i].values;
            const g = [v[0] - c, v[1] - c, v[2] - c];
            const glen = Math.sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
            const dot = (g[0] * n[0] + g[1] * n[1] + g[2] * n[2]) / glen;
            expect(dot).toBeGreaterThan(0.5);
        }
    });

    it('computeNormals returns zero normals for unreferenced vertices', () => {
        const normals = extractor.computeNormals(
            [], []);
        expect(normals.length).toBe(0);
    });
});

describe('SurfaceExtractorMC: gradient', () => {
    it('matches the analytic gradient of a linear field', () => {
        // F = 2x + 3y + 5z is reproduced exactly by trilinear interpolation.
        const image = makeImage(4, 4, 4, (x, y, z) => 2 * x + 3 * y + 5 * z);
        const extractor = new SurfaceExtractorMC(image);
        // getGradient is protected upstream; reach it for the test.
        const getGradient = (p: number[]) =>
            (extractor as unknown as {
                getGradient(v: { values: number[] }): { values: number[] };
            }).getGradient({ values: p }).values;

        expect(getGradient([0.5, 0.5, 0.5])).toEqual([2, 3, 5]);
        expect(getGradient([1.25, 2.75, 0.125])).toEqual([2, 3, 5]);
    });

    it('returns zero outside the voxels of the image', () => {
        const image = makeImage(4, 4, 4, (x, y, z) => x + y + z);
        const extractor = new SurfaceExtractorMC(image);
        const getGradient = (p: number[]) =>
            (extractor as unknown as {
                getGradient(v: { values: number[] }): { values: number[] };
            }).getGradient({ values: p }).values;

        expect(getGradient([-0.5, 0.5, 0.5])).toEqual([0, 0, 0]);
        expect(getGradient([0.5, -1, 0.5])).toEqual([0, 0, 0]);
        expect(getGradient([0.5, 0.5, -0.001])).toEqual([0, 0, 0]);
        // The last voxel of the image starts at index 2; index 3 is outside.
        expect(getGradient([3, 0.5, 0.5])).toEqual([0, 0, 0]);
        expect(getGradient([0.5, 3, 0.5])).toEqual([0, 0, 0]);
        expect(getGradient([0.5, 0.5, 3])).toEqual([0, 0, 0]);
    });
});
