import { describe, expect, it } from 'vitest';
import { SurfaceExtractorCubes } from '../src/SurfaceExtractorCubes';
import { SurfaceExtractorTetrahedra } from '../src/SurfaceExtractorTetrahedra';
import type { SurfaceExtractorTriangle } from '../src/SurfaceExtractor';

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

describe('SurfaceExtractorCubes: construction', () => {
    it('rejects bounds smaller than 2', () => {
        expect(() => new SurfaceExtractorCubes(1, 2, 2, [0, 0, 0, 0]))
            .toThrow('Invalid input.');
        expect(() => new SurfaceExtractorCubes(2, 1, 2, [0, 0, 0, 0]))
            .toThrow('Invalid input.');
        expect(() => new SurfaceExtractorCubes(2, 2, 1, [0, 0, 0, 0]))
            .toThrow('Invalid input.');
    });

    it('produces no geometry for a constant image', () => {
        const voxels = makeVoxels(4, 4, 4, () => 7);
        const extractor = new SurfaceExtractorCubes(4, 4, 4, voxels);
        for (const level of [0, 6, 7, 8, 100]) {
            const { vertices, triangles } = extractor.extract(level, true);
            expect(vertices.length).toBe(0);
            expect(triangles.length).toBe(0);
        }
    });
});

describe('SurfaceExtractorCubes: planar level surface', () => {
    // The extractor doubles the voxel values and uses the odd level
    // 2*level+1, so the extracted surface of the image F(x,y,z) = x at
    // level = 1 is the plane x = 1.5.
    const xBound = 4, yBound = 3, zBound = 3;
    const voxels = makeVoxels(xBound, yBound, zBound, (x) => x);
    const extractor = new SurfaceExtractorCubes(xBound, yBound, zBound, voxels);

    it('places every vertex exactly on the plane', () => {
        const { vertices, triangles } = extractor.extract(1, true);
        expect(vertices.length).toBeGreaterThan(0);
        expect(triangles.length).toBe(2 * (yBound - 1) * (zBound - 1));
        for (const v of vertices) {
            expect(v[0]).toBe(1.5);
            expect(v[1]).toBeGreaterThanOrEqual(0);
            expect(v[1]).toBeLessThanOrEqual(yBound - 1);
            expect(v[2]).toBeGreaterThanOrEqual(0);
            expect(v[2]).toBeLessThanOrEqual(zBound - 1);
        }
    });

    it('produces the exact rational coordinates', () => {
        const { vertices } = extractor.extractRational(1);
        expect(vertices.length).toBeGreaterThan(0);
        for (const v of vertices) {
            // x = 3/2 after the doubling (numerator 3, denominator 2).
            expect(v.xNumer / v.xDenom).toBe(1.5);
            expect(v.xDenom).toBeGreaterThan(0);
            expect(v.yDenom).toBe(1);
            expect(v.zDenom).toBe(1);
        }
    });

    it('welds duplicated vertices without changing the triangle count', () => {
        const raw = extractor.extract(1, false);
        const unique = extractor.extract(1, true);
        expect(raw.vertices.length).toBe(3 * raw.triangles.length);
        expect(unique.triangles.length).toBe(raw.triangles.length);
        expect(unique.vertices.length).toBeLessThan(raw.vertices.length);

        // The welded vertices are pairwise distinct.
        const keys = new Set(unique.vertices.map((v) => v.join(',')));
        expect(keys.size).toBe(unique.vertices.length);
    });

    it('levels outside the value range produce nothing', () => {
        expect(extractor.extract(-1, true).triangles.length).toBe(0);
        expect(extractor.extract(xBound, true).triangles.length).toBe(0);
    });
});

describe('SurfaceExtractorCubes: spherical level surface', () => {
    const c = 5;
    const bound = 11;
    const r2 = 10;
    // F = r2 - |P - C|^2, extracted at level 0, so the surface is
    // |P - C|^2 = r2 - 0.5.
    const voxels = makeVoxels(bound, bound, bound, (x, y, z) =>
        r2 - ((x - c) * (x - c) + (y - c) * (y - c) + (z - c) * (z - c)));
    const extractor = new SurfaceExtractorCubes(bound, bound, bound, voxels);

    it('places vertices within half a voxel of the true sphere', () => {
        const { vertices, triangles } = extractor.extract(0, true);
        expect(vertices.length).toBeGreaterThan(0);
        expect(triangles.length).toBeGreaterThan(0);

        const radius = Math.sqrt(r2 - 0.5);
        let minR = Number.MAX_VALUE;
        let maxR = 0;
        for (const v of vertices) {
            const dx = v[0] - c;
            const dy = v[1] - c;
            const dz = v[2] - c;
            const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
            minR = Math.min(minR, r);
            maxR = Math.max(maxR, r);
        }
        expect(minR).toBeGreaterThan(radius - 0.5);
        expect(maxR).toBeLessThan(radius + 0.5);
    });

    it('extracts a closed surface: every edge is shared by two triangles', () => {
        const { triangles } = extractor.extract(0, true);
        for (const count of edgeCounts(triangles).values()) {
            expect(count).toBe(2);
        }
    });

    it('is symmetric under the symmetries of the sampled field', () => {
        const { vertices } = extractor.extract(0, true);
        // The field is invariant under reflection about each coordinate
        // plane through the center, so the vertex set must be too. Compare
        // the reflected multisets of coordinates rather than exact keys,
        // because the rational-to-double conversion of a mirrored vertex need
        // not reproduce the double of its mirror image bit for bit.
        for (let d = 0; d < 3; ++d) {
            const forward = vertices.map((v) => v[d]).sort((a, b) => a - b);
            const mirrored = vertices.map((v) => 2 * c - v[d])
                .sort((a, b) => a - b);
            expect(forward.length).toBe(mirrored.length);
            for (let i = 0; i < forward.length; ++i) {
                expect(forward[i]).toBeCloseTo(mirrored[i], 12);
            }
        }
    });

    it('nested levels give nested surfaces', () => {
        const inner = extractor.extract(5, true);
        const outer = extractor.extract(0, true);
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
        expect(maxRadius(inner.vertices)).toBeLessThan(
            maxRadius(outer.vertices));
    });
});

describe('SurfaceExtractorCubes: face ambiguity resolution', () => {
    // A single voxel whose xmin face has all four edges cut. The resolution
    // depends on the sign of det = f00*f11 - f01*f10 for the face values.
    // The extractor uses F = 2*voxel - (2*level + 1), so with level 0 a voxel
    // value v maps to 2v - 1.
    const singleVoxel = (v000: number, v100: number, v010: number,
        v110: number, v001: number, v101: number, v011: number,
        v111: number) => {
        const voxels = [v000, v100, v010, v110, v001, v101, v011, v111];
        return new SurfaceExtractorCubes(2, 2, 2, voxels);
    };

    it('det > 0 pairs the face vertices into a single six-sided loop', () => {
        // f00 = +1, f10 = -1, f11 = +5, f01 = -1: det = 5 - 1 = 4 > 0.
        const extractor = singleVoxel(1, 1, 0, 1, 0, 1, 3, 1);
        const { vertices, triangles } = extractor.extract(0, true);
        expect(vertices.length).toBe(6);
        expect(triangles.length).toBe(4);
        // Every vertex lies on the boundary of the unit voxel.
        for (const v of vertices) {
            expect(Math.min(v[0], v[1], v[2])).toBeGreaterThanOrEqual(0);
            expect(Math.max(v[0], v[1], v[2])).toBeLessThanOrEqual(1);
        }
    });

    it('det < 0 pairs the face vertices into two triangles', () => {
        // f00 = +1, f10 = -3, f11 = +1, f01 = -3: det = 1 - 9 = -8 < 0.
        const extractor = singleVoxel(1, 1, -1, 1, -1, 1, 1, 1);
        const { vertices, triangles } = extractor.extract(0, true);
        expect(vertices.length).toBe(6);
        expect(triangles.length).toBe(2);
    });

    it('det == 0 inserts a branch point on the face', () => {
        // f00 = +1, f10 = -1, f11 = +1, f01 = -1: det = 1 - 1 = 0.
        const extractor = singleVoxel(1, 1, 0, 1, 0, 1, 1, 1);
        const { vertices, triangles } = extractor.extract(0, true);
        // The plus-sign configuration adds the face branch point, so there is
        // one more vertex than in the two hyperbolic resolutions.
        expect(vertices.length).toBe(7);
        expect(triangles.length).toBe(4);

        // The branch point is interior to the xmin face (x = 0).
        const branch = vertices.filter((v) =>
            v[0] === 0 && v[1] > 0 && v[1] < 1 && v[2] > 0 && v[2] < 1);
        expect(branch.length).toBe(1);
    });

    it('all three resolutions cut the same set of voxel edges', () => {
        const runs = [
            singleVoxel(1, 1, 0, 1, 0, 1, 3, 1),
            singleVoxel(1, 1, -1, 1, -1, 1, 1, 1),
            singleVoxel(1, 1, 0, 1, 0, 1, 1, 1)
        ].map((extractor) => extractor.extract(0, true).vertices);

        // The six vertices that lie on voxel edges (exactly one coordinate
        // strictly between 0 and 1) sit on the same six voxel edges in all
        // three runs. The positions along those edges differ because the
        // three runs use different voxel values, so compare the edges, which
        // are named by the varying axis and the two integer coordinates.
        const onEdge = (vertices: [number, number, number][]) => {
            const keys: string[] = [];
            for (const v of vertices) {
                const interior = [0, 1, 2].filter(
                    (d) => v[d] > 0 && v[d] < 1);
                if (interior.length === 1) {
                    const axis = interior[0];
                    const rest = [0, 1, 2].filter((d) => d !== axis)
                        .map((d) => v[d]);
                    keys.push(`${axis}:${rest.join(',')}`);
                }
            }
            keys.sort();
            return keys;
        };
        expect(onEdge(runs[1])).toEqual(onEdge(runs[0]));
        expect(onEdge(runs[2])).toEqual(onEdge(runs[0]));
        expect(onEdge(runs[0]).length).toBe(6);
    });
});

describe('SurfaceExtractorCubes: orientation, normals and gradient', () => {
    const c = 5;
    const bound = 11;
    const r2 = 10;
    const voxels = makeVoxels(bound, bound, bound, (x, y, z) =>
        r2 - ((x - c) * (x - c) + (y - c) * (y - c) + (z - c) * (z - c)));
    const extractor = new SurfaceExtractorCubes(bound, bound, bound, voxels);

    it('orients the triangles consistently with the gradient', () => {
        const { vertices, triangles } = extractor.extract(0, true);
        const getGradient = (p: [number, number, number]) =>
            (extractor as unknown as {
                getGradient(pos: [number, number, number]):
                    [number, number, number];
            }).getGradient(p);

        const signedVolume = (tris: SurfaceExtractorTriangle[]) => {
            let sum = 0;
            for (const t of tris) {
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

        // The sampled field decreases outward, so the gradient points inward
        // and the sameDir orientation gives inward-facing normals.
        extractor.orientTriangles(vertices, triangles, true);
        expect(signedVolume(triangles)).toBeLessThan(0);
        for (const t of triangles) {
            const v0 = vertices[t.v[0]];
            const g = getGradient(v0);
            // The gradient of r2 - |P - C|^2 points toward the center.
            const radial = (v0[0] - c) * g[0] + (v0[1] - c) * g[1]
                + (v0[2] - c) * g[2];
            expect(radial).toBeLessThanOrEqual(0);
        }

        extractor.orientTriangles(vertices, triangles, false);
        expect(signedVolume(triangles)).toBeGreaterThan(0);
    });

    it('computes unit-length vertex normals', () => {
        const { vertices, triangles } = extractor.extract(0, true);
        extractor.orientTriangles(vertices, triangles, false);
        const normals = extractor.computeNormals(vertices, triangles);
        expect(normals.length).toBe(vertices.length);
        for (let i = 0; i < normals.length; ++i) {
            const n = normals[i];
            const len = Math.sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            expect(len).toBeCloseTo(1, 10);

            // The outward orientation makes the normal roughly radial.
            const v = vertices[i];
            const g = [v[0] - c, v[1] - c, v[2] - c];
            const glen = Math.sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
            expect((g[0] * n[0] + g[1] * n[1] + g[2] * n[2]) / glen)
                .toBeGreaterThan(0.5);
        }
    });

    it('the gradient is the exact gradient of a linear field', () => {
        // The internal voxels are 2*(2x + 3y + 5z) - (2*level + 1), and the
        // constant offset does not change the gradient.
        const linear = makeVoxels(4, 4, 4, (x, y, z) => 2 * x + 3 * y + 5 * z);
        const lin = new SurfaceExtractorCubes(4, 4, 4, linear);
        lin.extractRational(0);
        const getGradient = (p: [number, number, number]) =>
            (lin as unknown as {
                getGradient(pos: [number, number, number]):
                    [number, number, number];
            }).getGradient(p);

        expect(getGradient([0.5, 0.5, 0.5])).toEqual([4, 6, 10]);
        expect(getGradient([2.25, 1.75, 0.125])).toEqual([4, 6, 10]);
        // Outside the voxels of the image.
        expect(getGradient([-1, 0.5, 0.5])).toEqual([0, 0, 0]);
        expect(getGradient([0.5, -1, 0.5])).toEqual([0, 0, 0]);
        expect(getGradient([0.5, 0.5, -1])).toEqual([0, 0, 0]);
        expect(getGradient([3, 0.5, 0.5])).toEqual([0, 0, 0]);
        expect(getGradient([0.5, 3, 0.5])).toEqual([0, 0, 0]);
        expect(getGradient([0.5, 0.5, 3])).toEqual([0, 0, 0]);
    });
});

describe('SurfaceExtractorCubes vs SurfaceExtractorTetrahedra', () => {
    // The cubes extractor solves 2*V - (2*L + 1) = 0 and the tetrahedra
    // extractor solves V' - L' = 0, so doubling the voxels and using the odd
    // level 2*L+1 makes the two extract the same level surface.
    it('agree exactly on a planar level surface', () => {
        const xBound = 4, yBound = 3, zBound = 3;
        const voxels = makeVoxels(xBound, yBound, zBound, (x) => x);
        const doubled = voxels.map((v) => 2 * v);

        const cubes = new SurfaceExtractorCubes(xBound, yBound, zBound, voxels);
        const tetra = new SurfaceExtractorTetrahedra(xBound, yBound, zBound,
            doubled);

        const fromCubes = cubes.extract(1, true);
        const fromTetra = tetra.extract(3, true);

        expect(fromCubes.vertices.length).toBeGreaterThan(0);
        expect(fromTetra.vertices.length).toBeGreaterThan(0);
        for (const v of fromCubes.vertices) {
            expect(v[0]).toBe(1.5);
        }
        for (const v of fromTetra.vertices) {
            expect(v[0]).toBe(1.5);
        }

        // Both meshes cover the same rectangle of the plane.
        const bounds = (vertices: [number, number, number][], d: number) => {
            let lo = Number.MAX_VALUE;
            let hi = -Number.MAX_VALUE;
            for (const v of vertices) {
                lo = Math.min(lo, v[d]);
                hi = Math.max(hi, v[d]);
            }
            return [lo, hi];
        };
        expect(bounds(fromTetra.vertices, 1)).toEqual(
            bounds(fromCubes.vertices, 1));
        expect(bounds(fromTetra.vertices, 2)).toEqual(
            bounds(fromCubes.vertices, 2));
    });

    it('agree on a spherical level surface to within half a voxel', () => {
        const c = 5;
        const bound = 11;
        const r2 = 10;
        const voxels = makeVoxels(bound, bound, bound, (x, y, z) =>
            r2 - ((x - c) * (x - c) + (y - c) * (y - c) + (z - c) * (z - c)));
        const doubled = voxels.map((v) => 2 * v);

        const cubes = new SurfaceExtractorCubes(bound, bound, bound, voxels);
        const tetra = new SurfaceExtractorTetrahedra(bound, bound, bound,
            doubled);
        const fromCubes = cubes.extract(0, true);
        const fromTetra = tetra.extract(1, true);

        const radius = Math.sqrt(r2 - 0.5);
        for (const vertices of [fromCubes.vertices, fromTetra.vertices]) {
            expect(vertices.length).toBeGreaterThan(0);
            for (const v of vertices) {
                const dx = v[0] - c;
                const dy = v[1] - c;
                const dz = v[2] - c;
                const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
                expect(Math.abs(r - radius)).toBeLessThan(0.5);
            }
        }

        // The two extractors cut the same voxel edges, so every vertex that
        // one of them places in the interior of a voxel edge is placed by the
        // other as well.
        const onEdge = (vertices: [number, number, number][]) => {
            const keys = new Set<string>();
            for (const v of vertices) {
                if (v.filter((t) => t !== Math.floor(t)).length === 1) {
                    keys.add(v.join(','));
                }
            }
            return keys;
        };
        const cubeKeys = onEdge(fromCubes.vertices);
        const tetraKeys = onEdge(fromTetra.vertices);
        expect(cubeKeys.size).toBeGreaterThan(0);
        expect(tetraKeys.size).toBe(cubeKeys.size);
        for (const key of cubeKeys) {
            expect(tetraKeys.has(key)).toBe(true);
        }
    });
});
