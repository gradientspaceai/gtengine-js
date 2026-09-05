import { describe, it, expect } from 'vitest';
import {
    AdaptiveSkeletonClimbing3,
    type AdaptiveSkeletonClimbing3Vertex
} from '../src/AdaptiveSkeletonClimbing3.js';
import { TriangleKey } from '../src/TriangleKey.js';
import { SurfaceExtractorCubes } from '../src/SurfaceExtractorCubes.js';
import { check, fc } from './helpers/arbitraries.js';

type Vertex = AdaptiveSkeletonClimbing3Vertex;

// Build a (2^N+1)^3 integer image from a function f(x, y, z), stored in
// lexicographic order for (x, y, z).
function makeImage(N: number, f: (x: number, y: number, z: number) => number): number[] {
    const size = (1 << N) + 1;
    const voxels = new Array<number>(size * size * size);
    for (let z = 0; z < size; ++z) {
        for (let y = 0; y < size; ++y) {
            for (let x = 0; x < size; ++x) {
                voxels[x + size * (y + size * z)] = f(x, y, z);
            }
        }
    }
    return voxels;
}

// A sphere level set: f(p) = round(scale * (R^2 - |p - center|^2)), which is
// positive inside the sphere of radius R and negative outside. Extracting
// level 0.5 gives (approximately) that sphere. The image gradient points
// inward.
function makeSphereImage(N: number, center: number, R: number, scale: number): number[] {
    return makeImage(N, (x, y, z) => Math.round(scale
        * (R * R - ((x - center) ** 2 + (y - center) ** 2 + (z - center) ** 2))));
}

// A deterministic linear congruential generator, so the randomized checks
// are reproducible.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (1664525 * s + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

const isInteger = (t: number): boolean => Math.abs(t - Math.round(t)) < 1.0e-12;

// Every vertex produced by the box-edge interpolation lies on a grid edge
// (exactly two integer coordinates) and must interpolate the image to the
// extraction level exactly. Vertices with fewer integer coordinates are the
// plus-sign branch points and the merged-box centroids, which are not
// required to lie on the level surface; they are counted and returned.
function checkVerticesOnLevelSet(voxels: number[], size: number, level: number,
    vertices: Vertex[]): number {
    let numComposite = 0;
    let numBad = 0;
    for (const v of vertices) {
        if (!Number.isFinite(v[0]) || !Number.isFinite(v[1]) || !Number.isFinite(v[2])) {
            ++numBad;
            continue;
        }
        const isInt = [isInteger(v[0]), isInteger(v[1]), isInteger(v[2])];
        const numInt = isInt.filter(Boolean).length;
        if (numInt === 3) {
            // All three integer would mean the level equals an image
            // sample, which the non-integer level forbids.
            ++numBad;
            continue;
        }
        if (numInt !== 2) {
            ++numComposite;
            continue;
        }

        const k = isInt.indexOf(false);
        const p = [Math.round(v[0]), Math.round(v[1]), Math.round(v[2])];
        p[k] = Math.floor(v[k]);
        const t = v[k] - p[k];
        const base = p[0] + size * (p[1] + size * p[2]);
        const step = k === 0 ? 1 : (k === 1 ? size : size * size);
        const f0 = voxels[base];
        const f1 = voxels[base + step];
        if (Math.abs(f0 + t * (f1 - f0) - level) > 1.0e-9) {
            ++numBad;
        }
    }
    expect(numBad).toBe(0);
    return numComposite;
}

// The vertices that are not on a grid edge: the plus-sign branch points and
// the merged-box fan centroids.
function countCompositeVertices(vertices: Vertex[]): number {
    let count = 0;
    for (const v of vertices) {
        const numInt = [isInteger(v[0]), isInteger(v[1]), isInteger(v[2])].filter(Boolean).length;
        if (numInt !== 2) {
            ++count;
        }
    }
    return count;
}

// Every triangle index must be a valid vertex index and the three indices
// must be distinct.
function checkIndices(vertices: Vertex[], triangles: TriangleKey[]): void {
    let numBad = 0;
    for (const t of triangles) {
        for (let i = 0; i < 3; ++i) {
            if (!Number.isInteger(t.V[i]) || t.V[i] < 0 || t.V[i] >= vertices.length) {
                ++numBad;
            }
        }
        if (t.V[0] === t.V[1] || t.V[1] === t.V[2] || t.V[0] === t.V[2]) {
            ++numBad;
        }
    }
    expect(numBad).toBe(0);
}

// The number of triangles sharing each undirected mesh edge.
function edgeUseCounts(triangles: TriangleKey[]): Map<string, number> {
    const counts = new Map<string, number>();
    for (const t of triangles) {
        const pairs: [number, number][] = [[t.V[0], t.V[1]], [t.V[1], t.V[2]], [t.V[2], t.V[0]]];
        for (const [a, b] of pairs) {
            const key = Math.min(a, b) + '-' + Math.max(a, b);
            counts.set(key, (counts.get(key) ?? 0) + 1);
        }
    }
    return counts;
}

// Six times the signed volume enclosed by the (closed) triangle mesh.
function signedVolume(vertices: Vertex[], triangles: TriangleKey[]): number {
    let volume = 0;
    for (const t of triangles) {
        const a = vertices[t.V[0]];
        const b = vertices[t.V[1]];
        const c = vertices[t.V[2]];
        volume += (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
    }
    return volume;
}

describe('AdaptiveSkeletonClimbing3', () => {
    it('throws on invalid construction input', () => {
        expect(() => new AdaptiveSkeletonClimbing3(0, new Array<number>(8).fill(0)))
            .toThrow('Invalid input.');
        expect(() => new AdaptiveSkeletonClimbing3(-1, new Array<number>(8).fill(0)))
            .toThrow('Invalid input.');
    });

    it('extracts a plane from a linear x-ramp at full resolution', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const voxels = makeImage(N, (x) => x);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const level = 3.5;

        const result = asc.extract(level, N);
        // One box per voxel; the level plane crosses the 8x8 columns of
        // voxels at x = 3, each contributing a quad split into 2 triangles.
        expect(asc.getNumBoxes()).toBe(512);
        asc.makeUnique(result.vertices, result.triangles);
        expect(result.triangles.length).toBe(128);
        // The 9x9 lattice of crossing points in the (y, z) plane.
        expect(result.vertices.length).toBe(81);
        for (const v of result.vertices) {
            expect(v[0]).toBeCloseTo(3.5, 12);
        }
        checkIndices(result.vertices, result.triangles);
        expect(checkVerticesOnLevelSet(voxels, size, level, result.vertices)).toBe(0);
    });

    it('merges the x-ramp into a coarse plane at depth 1', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const voxels = makeImage(N, (x) => x);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const level = 3.5;

        const result = asc.extract(level, 1);
        // The root's eight children each merge into a single 4x4x4 box.
        expect(asc.getNumBoxes()).toBe(8);
        asc.makeUnique(result.vertices, result.triangles);
        expect(result.triangles.length).toBe(8);
        // The 3x3 lattice at x = 3.5 with y, z in {0, 4, 8}.
        expect(result.vertices.length).toBe(9);
        const ys = new Set<number>();
        const zs = new Set<number>();
        for (const v of result.vertices) {
            expect(v[0]).toBeCloseTo(3.5, 12);
            ys.add(v[1]);
            zs.add(v[2]);
        }
        expect([...ys].sort((a, b) => a - b)).toEqual([0, 4, 8]);
        expect([...zs].sort((a, b) => a - b)).toEqual([0, 4, 8]);
        expect(checkVerticesOnLevelSet(voxels, size, level, result.vertices)).toBe(0);
    });

    it('extracts the same plane from ramps along y and along z', () => {
        const N = 3;
        for (const axis of [1, 2]) {
            const voxels = makeImage(N, (x, y, z) => (axis === 1 ? y : z));
            const asc = new AdaptiveSkeletonClimbing3(N, voxels);
            const result = asc.extract(3.5, 1);
            asc.makeUnique(result.vertices, result.triangles);
            expect(result.triangles.length).toBe(8);
            expect(result.vertices.length).toBe(9);
            for (const v of result.vertices) {
                expect(v[axis]).toBeCloseTo(3.5, 12);
            }
        }
    });

    it('extracts a closed manifold surface for a sphere', () => {
        for (const N of [3, 4]) {
            const size = (1 << N) + 1;
            const center = (size - 1) / 2;
            const R = (size - 1) / 4;
            const voxels = makeSphereImage(N, center, R, 16);
            const asc = new AdaptiveSkeletonClimbing3(N, voxels);
            const result = asc.extract(0.5, N);
            asc.makeUnique(result.vertices, result.triangles);

            expect(result.triangles.length).toBeGreaterThan(0);
            checkIndices(result.vertices, result.triangles);
            // The sphere lies strictly inside the image, so the surface is
            // closed: every edge is shared by exactly two triangles.
            for (const count of edgeUseCounts(result.triangles).values()) {
                expect(count).toBe(2);
            }
            // Euler characteristic of a sphere: V - E + F = 2.
            const numEdges = edgeUseCounts(result.triangles).size;
            expect(result.vertices.length - numEdges + result.triangles.length).toBe(2);
            expect(checkVerticesOnLevelSet(voxels, size, 0.5, result.vertices)).toBe(0);
        }
    });

    it('encloses approximately the correct volume for a sphere', () => {
        const N = 4;
        const size = (1 << N) + 1;
        const center = (size - 1) / 2;
        const R = 5;
        const voxels = makeSphereImage(N, center, R, 16);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const result = asc.extract(0.5, N);
        asc.makeUnique(result.vertices, result.triangles);
        // The image gradient points inward for this sphere, so 'sameDir'
        // false gives outward normals and a positive signed volume.
        asc.orientTriangles(result.vertices, result.triangles, false);
        const volume = signedVolume(result.vertices, result.triangles);
        const exact = (4 / 3) * Math.PI * R * R * R;
        expect(volume).toBeGreaterThan(0);
        // The polygonal approximation inscribes the sphere, so it slightly
        // underestimates the volume.
        expect(volume).toBeLessThan(exact);
        expect(volume / exact).toBeGreaterThan(0.95);
    });

    it('orients triangles with respect to the image gradient', () => {
        const N = 4;
        const size = (1 << N) + 1;
        const center = (size - 1) / 2;
        const voxels = makeSphereImage(N, center, 5, 16);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);

        const withGradient = asc.extract(0.5, N);
        asc.makeUnique(withGradient.vertices, withGradient.triangles);
        asc.orientTriangles(withGradient.vertices, withGradient.triangles, true);

        const againstGradient = asc.extract(0.5, N);
        asc.makeUnique(againstGradient.vertices, againstGradient.triangles);
        asc.orientTriangles(againstGradient.vertices, againstGradient.triangles, false);

        // The gradient of R^2 - r^2 points inward, so 'sameDir' true gives
        // inward normals (negative volume) and false gives outward normals.
        const volumeIn = signedVolume(withGradient.vertices, withGradient.triangles);
        const volumeOut = signedVolume(againstGradient.vertices, againstGradient.triangles);
        expect(volumeIn).toBeLessThan(0);
        expect(volumeOut).toBeGreaterThan(0);
        expect(volumeIn).toBeCloseTo(-volumeOut, 9);

        // Re-orienting an already oriented mesh is a no-op.
        const before = withGradient.triangles.map((t) => t.V.join(','));
        asc.orientTriangles(withGradient.vertices, withGradient.triangles, true);
        expect(withGradient.triangles.map((t) => t.V.join(','))).toEqual(before);
    });

    it('computes unit-length normals that follow the sphere radius', () => {
        const N = 4;
        const size = (1 << N) + 1;
        const center = (size - 1) / 2;
        const voxels = makeSphereImage(N, center, 5, 16);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const result = asc.extract(0.5, N);
        asc.makeUnique(result.vertices, result.triangles);
        asc.orientTriangles(result.vertices, result.triangles, false);

        const normals = asc.computeNormals(result.vertices, result.triangles);
        expect(normals.length).toBe(result.vertices.length);
        for (let i = 0; i < normals.length; ++i) {
            const n = normals[i];
            expect(Math.hypot(n[0], n[1], n[2])).toBeCloseTo(1, 12);
            // Outward orientation: the normal has a positive component
            // along the radial direction.
            const v = result.vertices[i];
            const radial: Vertex = [v[0] - center, v[1] - center, v[2] - center];
            const length = Math.hypot(radial[0], radial[1], radial[2]);
            const dot = (n[0] * radial[0] + n[1] * radial[1] + n[2] * radial[2]) / length;
            expect(dot).toBeGreaterThan(0.5);
        }
    });

    it('computes zero normals for vertices with no triangles', () => {
        const N = 3;
        const voxels = makeSphereImage(N, 4, 2, 16);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const vertices: Vertex[] = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [5, 5, 5]];
        const normals = asc.computeNormals(vertices, [new TriangleKey(true, 0, 1, 2)]);
        expect(normals.length).toBe(4);
        expect(normals[3]).toEqual([0, 0, 0]);
        // The triangle lies in the z = 0 plane, so its vertex normals are
        // (0, 0, 1) with the (0,1,2) winding.
        for (let i = 0; i < 3; ++i) {
            expect(normals[i][0]).toBeCloseTo(0, 12);
            expect(normals[i][1]).toBeCloseTo(0, 12);
            expect(normals[i][2]).toBeCloseTo(1, 12);
        }
    });

    it('produces coarser meshes as the merge depth decreases', () => {
        const N = 4;
        const size = (1 << N) + 1;
        const center = (size - 1) / 2;
        const voxels = makeSphereImage(N, center, 5, 16);

        let previous = Number.POSITIVE_INFINITY;
        let previousBoxes = Number.POSITIVE_INFINITY;
        for (let depth = N; depth >= -1; --depth) {
            const asc = new AdaptiveSkeletonClimbing3(N, voxels);
            const result = asc.extract(0.5, depth);
            asc.makeUnique(result.vertices, result.triangles);
            expect(result.triangles.length).toBeLessThanOrEqual(previous);
            expect(asc.getNumBoxes()).toBeLessThanOrEqual(previousBoxes);
            previous = result.triangles.length;
            previousBoxes = asc.getNumBoxes();
        }
        // The coarsest extraction is far cheaper than the finest.
        const fine = new AdaptiveSkeletonClimbing3(N, voxels);
        const fineResult = fine.extract(0.5, N);
        fine.makeUnique(fineResult.vertices, fineResult.triangles);
        expect(previous).toBeLessThan(fineResult.triangles.length / 2);
    });

    it('forces full resolution when the boundary is fixed', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const center = (size - 1) / 2;
        const voxels = makeSphereImage(N, center, 2, 16);

        const fixed = new AdaptiveSkeletonClimbing3(N, voxels, true);
        const fixedResult = fixed.extract(0.5, -1);
        fixed.makeUnique(fixedResult.vertices, fixedResult.triangles);
        // Every leaf becomes a 1x1x1 box, so no merging happens at all.
        expect(fixed.getNumBoxes()).toBe(512);

        const full = new AdaptiveSkeletonClimbing3(N, voxels);
        const fullResult = full.extract(0.5, N);
        full.makeUnique(fullResult.vertices, fullResult.triangles);
        expect(fixedResult.vertices.length).toBe(fullResult.vertices.length);
        expect(fixedResult.triangles.length).toBe(fullResult.triangles.length);

        // Without the fixed boundary the same image merges.
        const adaptive = new AdaptiveSkeletonClimbing3(N, voxels);
        adaptive.extract(0.5, -1);
        expect(adaptive.getNumBoxes()).toBeLessThan(512);
    });

    it('returns nothing for a constant image or an out-of-range level', () => {
        const N = 3;
        const constant = makeImage(N, () => 5);
        const asc = new AdaptiveSkeletonClimbing3(N, constant);
        for (const level of [2.5, 5.5, -1.5]) {
            const result = asc.extract(level, N);
            expect(result.vertices.length).toBe(0);
            expect(result.triangles.length).toBe(0);
        }

        const size = (1 << N) + 1;
        const sphere = makeSphereImage(N, (size - 1) / 2, 2, 16);
        const asc2 = new AdaptiveSkeletonClimbing3(N, sphere);
        for (const level of [-1000.5, 1000.5]) {
            const result = asc2.extract(level, N);
            expect(result.vertices.length).toBe(0);
            expect(result.triangles.length).toBe(0);
        }
    });

    it('makeUnique removes duplicate vertices and keeps the triangles valid', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const voxels = makeSphereImage(N, (size - 1) / 2, 3, 16);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const result = asc.extract(0.5, N);

        // The raw extraction emits three positions per triangle.
        expect(result.vertices.length).toBe(3 * result.triangles.length);
        const rawTriangles = result.triangles.length;

        asc.makeUnique(result.vertices, result.triangles);
        expect(result.vertices.length).toBeLessThan(3 * rawTriangles);
        expect(result.triangles.length).toBe(rawTriangles);
        checkIndices(result.vertices, result.triangles);
        // The packed vertices are pairwise distinct.
        const keys = new Set(result.vertices.map((v) => v.join(',')));
        expect(keys.size).toBe(result.vertices.length);

        // makeUnique is idempotent.
        const vertexKeys = result.vertices.map((v) => v.join(','));
        const triangleKeys = result.triangles.map((t) => t.V.join(','));
        asc.makeUnique(result.vertices, result.triangles);
        expect(result.vertices.map((v) => v.join(','))).toEqual(vertexKeys);
        expect(result.triangles.map((t) => t.V.join(','))).toEqual(triangleKeys);
    });

    it('makeUnique does nothing for empty input', () => {
        const N = 3;
        const asc = new AdaptiveSkeletonClimbing3(N, makeImage(N, () => 0));
        const vertices: Vertex[] = [];
        const triangles: TriangleKey[] = [];
        asc.makeUnique(vertices, triangles);
        expect(vertices.length).toBe(0);
        expect(triangles.length).toBe(0);

        const someVertices: Vertex[] = [[0, 0, 0], [1, 0, 0]];
        asc.makeUnique(someVertices, []);
        expect(someVertices.length).toBe(2);
    });

    it('accepts typed-array images', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const center = (size - 1) / 2;
        const plain = makeSphereImage(N, center, 3, 16);
        const typed = Int32Array.from(plain);

        const a = new AdaptiveSkeletonClimbing3(N, plain);
        const ra = a.extract(0.5, N);
        a.makeUnique(ra.vertices, ra.triangles);

        const b = new AdaptiveSkeletonClimbing3(N, typed);
        const rb = b.extract(0.5, N);
        b.makeUnique(rb.vertices, rb.triangles);

        expect(rb.vertices).toEqual(ra.vertices);
        expect(rb.triangles.map((t) => t.V.join(','))).toEqual(ra.triangles.map((t) => t.V.join(',')));
    });

    it('extracts the two components of a two-sphere image', () => {
        const N = 4;
        const size = (1 << N) + 1;
        const voxels = makeImage(N, (x, y, z) => {
            const d0 = 9 - ((x - 4) ** 2 + (y - 4) ** 2 + (z - 4) ** 2);
            const d1 = 9 - ((x - 12) ** 2 + (y - 12) ** 2 + (z - 12) ** 2);
            return Math.round(16 * Math.max(d0, d1));
        });
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const result = asc.extract(0.5, N);
        asc.makeUnique(result.vertices, result.triangles);

        checkIndices(result.vertices, result.triangles);
        for (const count of edgeUseCounts(result.triangles).values()) {
            expect(count).toBe(2);
        }
        // Two closed components: V - E + F = 4.
        const numEdges = edgeUseCounts(result.triangles).size;
        expect(result.vertices.length - numEdges + result.triangles.length).toBe(4);
        expect(checkVerticesOnLevelSet(voxels, size, 0.5, result.vertices)).toBe(0);

        // The vertices split cleanly into the two spheres.
        let near0 = 0;
        let near1 = 0;
        for (const v of result.vertices) {
            const r0 = Math.hypot(v[0] - 4, v[1] - 4, v[2] - 4);
            const r1 = Math.hypot(v[0] - 12, v[1] - 12, v[2] - 12);
            if (r0 < 4) {
                ++near0;
            } else if (r1 < 4) {
                ++near1;
            }
        }
        expect(near0 + near1).toBe(result.vertices.length);
        expect(near0).toBe(near1);
    });

    it('disambiguates the four-intersection face configuration by the determinant', () => {
        // A 1x1x1 box whose zmin face has a sign change on all four edges.
        // The sign of f00 * f11 - f01 * f10 selects which pair of face
        // intersections is joined.
        const N = 1;
        const size = 3;
        const build = (f: (x: number, y: number, z: number) => number) => makeImage(N, f);

        // The corner values on the z = 0 face form a saddle. The two
        // determinant signs must produce different edge pairings, hence
        // different triangulations of the same twelve-crossing box.
        const positiveDet = build((x, y, z) => {
            if (z > 0) {
                return -1;
            }
            return (x === y) ? 4 : -1;
        });
        const negativeDet = build((x, y, z) => {
            if (z > 0) {
                return -1;
            }
            return (x !== y) ? 4 : -1;
        });

        for (const voxels of [positiveDet, negativeDet]) {
            const asc = new AdaptiveSkeletonClimbing3(N, voxels);
            const result = asc.extract(0.5, N);
            asc.makeUnique(result.vertices, result.triangles);
            checkIndices(result.vertices, result.triangles);
            checkVerticesOnLevelSet(voxels, size, 0.5, result.vertices);
            expect(result.triangles.length).toBeGreaterThan(0);
        }

        // The two configurations differ, which is the point of the
        // determinant test.
        const a = new AdaptiveSkeletonClimbing3(N, positiveDet);
        const ra = a.extract(0.5, N);
        a.makeUnique(ra.vertices, ra.triangles);
        const b = new AdaptiveSkeletonClimbing3(N, negativeDet);
        const rb = b.extract(0.5, N);
        b.makeUnique(rb.vertices, rb.triangles);
        expect(ra.vertices.map((v) => v.join(','))).not.toEqual(rb.vertices.map((v) => v.join(',')));
    });

    it('adds a face branch point for the plus-sign configuration (det = 0)', () => {
        // The z = 0 face of the box at the origin has the corner values
        // f00 = f11 = 2 and f01 = f10 = -2, so all four face edges are
        // crossed and f00 * f11 - f01 * f10 = 0. Upstream then inserts a
        // branch point at the center of the plus sign rather than choosing
        // one of the two hyperbolic pairings.
        const N = 1;
        const size = 3;
        const voxels = makeImage(N, (x, y, z) => {
            if (z > 0) {
                return -2;
            }
            return (x + y) % 2 === 0 ? 2 : -2;
        });

        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const result = asc.extract(0.5, N);
        asc.makeUnique(result.vertices, result.triangles);
        checkIndices(result.vertices, result.triangles);

        // The z = 0 layer is a 3x3 checkerboard, so all four voxels of that
        // layer are plus-sign configurations. Each contributes one branch
        // point, the only vertices that are not on a grid edge: they have
        // exactly one integer coordinate (z = 0) and take their x from the
        // ymin-zmin crossing and their y from the xmin-zmin crossing.
        const branch = result.vertices.filter((v) =>
            [isInteger(v[0]), isInteger(v[1]), isInteger(v[2])].filter(Boolean).length === 1);
        expect(branch.length).toBe(4);
        for (const b of branch) {
            expect(b[2]).toBe(0);
        }
        // The voxel at the origin crosses x at 0 + (0.5 - 2) / (-2 - 2) and
        // y at the same fraction.
        const keys = branch.map((b) => b.join(','));
        expect(keys).toContain([0.375, 0.375, 0].join(','));
        expect(checkVerticesOnLevelSet(voxels, size, 0.5, result.vertices)).toBe(4);
    });

    it('produces well-formed meshes for randomized blob images', () => {
        const N = 4;
        const size = (1 << N) + 1;
        let numComposite = 0;
        let numOutOfDomain = 0;
        let maxEdgeUse = 0;
        for (let seed = 1; seed <= 8; ++seed) {
            const random = makeRandom(seed);
            const blobs: number[][] = [];
            for (let b = 0; b < 3; ++b) {
                blobs.push([random() * (size - 1), random() * (size - 1),
                    random() * (size - 1), 2 + random() * 3]);
            }
            const voxels = makeImage(N, (x, y, z) => {
                let value = -1000;
                for (const [bx, by, bz, br] of blobs) {
                    value = Math.max(value,
                        br * br - ((x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2));
                }
                return Math.round(4 * value);
            });

            for (const depth of [-1, 2, N]) {
                const asc = new AdaptiveSkeletonClimbing3(N, voxels);
                const result = asc.extract(0.5, depth);
                asc.makeUnique(result.vertices, result.triangles);
                checkIndices(result.vertices, result.triangles);
                for (const v of result.vertices) {
                    // Every vertex is finite and lies inside the image
                    // domain.
                    for (let i = 0; i < 3; ++i) {
                        if (!Number.isFinite(v[i]) || v[i] < 0 || v[i] > size - 1) {
                            ++numOutOfDomain;
                        }
                    }
                }
                for (const count of edgeUseCounts(result.triangles).values()) {
                    maxEdgeUse = Math.max(maxEdgeUse, count);
                }
                numComposite += countCompositeVertices(result.vertices);
            }
        }
        expect(numOutOfDomain).toBe(0);
        // Every mesh edge is used at most twice.
        expect(maxEdgeUse).toBe(2);
        // These images exercise the merged-box fan tessellation, whose
        // centroid vertices are the only ones not on a grid edge.
        expect(numComposite).toBeGreaterThan(0);
    });

    it('reuses the extractor across calls without accumulating state', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const voxels = makeSphereImage(N, (size - 1) / 2, 3, 16);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);

        const first = asc.extract(0.5, N);
        const firstBoxes = asc.getNumBoxes();
        asc.makeUnique(first.vertices, first.triangles);

        // A different level in between must not disturb the repeat.
        asc.extract(100.5, N);

        const second = asc.extract(0.5, N);
        expect(asc.getNumBoxes()).toBe(firstBoxes);
        asc.makeUnique(second.vertices, second.triangles);
        expect(second.vertices).toEqual(first.vertices);
        expect(second.triangles.map((t) => t.V.join(','))).toEqual(first.triangles.map((t) => t.V.join(',')));
    });

    it('drops the level surface when the whole image merges into one box (upstream bug)', () => {
        // Upstream AdaptiveSkeletonClimbing3::Merge adds a monobox only from
        // the parent of the octree node that produced it. When the root node
        // itself merges into a single box the recursion returns true to a
        // caller that ignores the value, so the box is never added and the
        // level surface is lost. A linear ramp is strongly monotone in every
        // direction, so the whole image merges. The port preserves this
        // behavior; see the PR discussion.
        const N = 3;
        const voxels = makeImage(N, (x) => x);
        for (const depth of [0, -1, -5]) {
            const asc = new AdaptiveSkeletonClimbing3(N, voxels);
            const result = asc.extract(3.5, depth);
            expect(asc.getNumBoxes()).toBe(0);
            expect(result.vertices.length).toBe(0);
            expect(result.triangles.length).toBe(0);
        }

        // At depth 1 the root's children are the largest boxes that can be
        // added, and the surface reappears.
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const result = asc.extract(3.5, 1);
        expect(asc.getNumBoxes()).toBe(8);
        expect(result.triangles.length).toBe(8);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V27): property-based cross-checks against the upstream
// header AdaptiveSkeletonClimbing3.h.
// ---------------------------------------------------------------------------

describe('AdaptiveSkeletonClimbing3 verification', () => {
    const N = 3;
    const size = (1 << N) + 1;   // 9

    // A field of two spheres, the configuration the extractor is designed
    // for: large monotone regions that merge plus a nontrivial surface.
    const twoSphereArb = fc.tuple(
        fc.integer({ min: 2, max: 4 }), fc.integer({ min: 2, max: 4 }),
        fc.integer({ min: 2, max: 4 }), fc.integer({ min: 4, max: 6 }),
        fc.integer({ min: 4, max: 6 }), fc.integer({ min: 4, max: 6 }),
        fc.integer({ min: 2, max: 5 }), fc.integer({ min: 2, max: 5 }))
        .map(([ax, ay, az, bx, by, bz, ra2, rb2]) =>
            makeImage(N, (x, y, z) => Math.max(
                ra2 - ((x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2),
                rb2 - ((x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2))));

    // The grid edge a vertex lies on, as "x,y,z|d" with (x,y,z) the lattice
    // endpoint of smaller coordinate along the direction d; null when the
    // vertex is not in the interior of a grid edge.
    const edgeKey = (v: readonly number[]): string | null => {
        let d = -1;
        for (let i = 0; i < 3; ++i) {
            if (!isInteger(v[i])) {
                if (d >= 0) {
                    return null;
                }
                d = i;
            }
        }
        if (d < 0) {
            return null;
        }
        return Math.round(v[0] - (d === 0 ? v[0] - Math.floor(v[0]) : 0))
            + ',' + Math.round(v[1] - (d === 1 ? v[1] - Math.floor(v[1]) : 0))
            + ',' + Math.round(v[2] - (d === 2 ? v[2] - Math.floor(v[2]) : 0))
            + '|' + d;
    };

    it('cuts the same grid edges as SurfaceExtractorCubes at full resolution',
        () => {
            // AdaptiveSkeletonClimbing3 extracts the surface F = level for a
            // real level; SurfaceExtractorCubes doubles the voxel values and
            // subtracts the odd integer 2*L+1, so it extracts F = L + 1/2.
            // With level = 1/2 and L = 0 the two solve the same equation by
            // completely different algorithms.
            check(twoSphereArb, (voxels) => {
                const asc = new AdaptiveSkeletonClimbing3(N, voxels);
                const fine = asc.extract(0.5, N);
                asc.makeUnique(fine.vertices, fine.triangles);

                const cubes = new SurfaceExtractorCubes(size, size, size,
                    voxels).extract(0, true);

                const collect = (vertices: readonly (readonly number[])[]) => {
                    const map = new Map<string, number>();
                    for (const v of vertices) {
                        const key = edgeKey(v);
                        if (key !== null) {
                            map.set(key, v[Number(key.split('|')[1])]);
                        }
                    }
                    return map;
                };
                const a = collect(fine.vertices);
                const b = collect(cubes.vertices);
                expect([...a.keys()].sort()).toEqual([...b.keys()].sort());
                for (const [key, value] of a) {
                    // The two evaluate the same crossing by different
                    // sequences of divisions.
                    expect(value).toBeCloseTo(b.get(key) as number, 9);
                }
            }, 25);
        });

    it('every vertex is on the level set or is a branch point or centroid',
        () => {
            check(fc.tuple(twoSphereArb, fc.integer({ min: 1, max: N })),
                ([voxels, depth]) => {
                    const asc = new AdaptiveSkeletonClimbing3(N, voxels);
                    const result = asc.extract(0.5, depth);
                    asc.makeUnique(result.vertices, result.triangles);
                    checkIndices(result.vertices, result.triangles);
                    checkVerticesOnLevelSet(voxels, size, 0.5,
                        result.vertices);
                    for (const v of result.vertices) {
                        for (let i = 0; i < 3; ++i) {
                            expect(v[i]).toBeGreaterThanOrEqual(0);
                            expect(v[i]).toBeLessThanOrEqual(size - 1);
                        }
                    }
                }, 25);
        });

    it('merging reduces the box count and the mesh size', () => {
        check(twoSphereArb, (voxels) => {
            const fine = new AdaptiveSkeletonClimbing3(N, voxels);
            const fineResult = fine.extract(0.5, N);
            fine.makeUnique(fineResult.vertices, fineResult.triangles);
            const fineBoxes = fine.getNumBoxes();

            let previousBoxes = fineBoxes;
            for (let depth = N - 1; depth >= 1; --depth) {
                const asc = new AdaptiveSkeletonClimbing3(N, voxels);
                const result = asc.extract(0.5, depth);
                asc.makeUnique(result.vertices, result.triangles);
                // Allowing more merging never adds boxes.
                expect(asc.getNumBoxes()).toBeLessThanOrEqual(previousBoxes);
                expect(result.vertices.length)
                    .toBeLessThanOrEqual(fineResult.vertices.length);
                previousBoxes = asc.getNumBoxes();
            }
            // The two-sphere field has large monotone regions, so the
            // coarsest allowed depth really does merge.
            expect(previousBoxes).toBeLessThan(fineBoxes);
        }, 20);
    });

    it('extracts closed manifold spheres at every depth', () => {
        check(fc.tuple(fc.integer({ min: 3, max: 5 }),
            fc.integer({ min: 1, max: N })), ([r2, depth]) => {
            const center = (size - 1) / 2;
            const voxels = makeImage(N, (x, y, z) =>
                r2 - ((x - center) ** 2 + (y - center) ** 2
                    + (z - center) ** 2));
            const asc = new AdaptiveSkeletonClimbing3(N, voxels);
            const result = asc.extract(0.5, depth);
            asc.makeUnique(result.vertices, result.triangles);
            expect(result.triangles.length).toBeGreaterThan(0);

            const counts = edgeUseCounts(result.triangles);
            for (const count of counts.values()) {
                expect(count).toBe(2);
            }
            // A single closed sphere: V - E + F = 2.
            expect(result.vertices.length - counts.size
                + result.triangles.length).toBe(2);
            // After a consistent orientation the mesh encloses a positive
            // volume inside the image. The merged boxes make the coarse
            // meshes much smaller than the true sphere, so only the sign and
            // the containment are asserted here; the deterministic test
            // above pins the volume at full resolution.
            asc.orientTriangles(result.vertices, result.triangles, false);
            const volume = Math.abs(signedVolume(result.vertices,
                result.triangles));
            expect(volume).toBeGreaterThan(0);
            expect(volume).toBeLessThan((size - 1) ** 3);
        }, 25);
    });

    it('reuse leaves no stale merge-tree or box state', () => {
        check(fc.tuple(twoSphereArb, fc.integer({ min: 1, max: N })),
            ([voxels, depth]) => {
                const shared = new AdaptiveSkeletonClimbing3(N, voxels);
                const first = shared.extract(0.5, depth);
                const boxes = shared.getNumBoxes();
                shared.makeUnique(first.vertices, first.triangles);

                // Interleave a different extraction.
                shared.extract(-0.5, 1);
                shared.extract(1000.5, N);

                const second = shared.extract(0.5, depth);
                expect(shared.getNumBoxes()).toBe(boxes);
                shared.makeUnique(second.vertices, second.triangles);
                expect(second.vertices).toEqual(first.vertices);
                expect(second.triangles.map((t) => t.V.join(',')))
                    .toEqual(first.triangles.map((t) => t.V.join(',')));

                // A fresh extractor must agree with the reused one.
                const fresh = new AdaptiveSkeletonClimbing3(N, voxels);
                const result = fresh.extract(0.5, depth);
                fresh.makeUnique(result.vertices, result.triangles);
                expect(result.vertices).toEqual(first.vertices);
            }, 20);
    });

    it('makeUnique packs first-encounter order and is idempotent', () => {
        check(fc.tuple(twoSphereArb, fc.integer({ min: 1, max: N })),
            ([voxels, depth]) => {
                const asc = new AdaptiveSkeletonClimbing3(N, voxels);
                const result = asc.extract(0.5, depth);
                const raw = result.vertices.map((v) => v.join(','));
                const rawTriangles = result.triangles.length;

                asc.makeUnique(result.vertices, result.triangles);

                // The packed vertices are the distinct raw positions in
                // first-encounter order, exactly as the upstream std::map
                // insertion indices give (this mirrors the 2D sibling).
                const expected: string[] = [];
                const seen = new Set<string>();
                for (const key of raw) {
                    if (!seen.has(key)) {
                        seen.add(key);
                        expected.push(key);
                    }
                }
                expect(result.vertices.map((v) => v.join(','))).toEqual(
                    expected);
                expect(result.triangles.length)
                    .toBeLessThanOrEqual(rawTriangles);
                checkIndices(result.vertices, result.triangles);

                // A second call changes nothing.
                const vertices = result.vertices.map((v) =>
                    [...v] as Vertex);
                const triangles = result.triangles.map((t) => t.V.join(','));
                asc.makeUnique(result.vertices, result.triangles);
                expect(result.vertices).toEqual(vertices);
                expect(result.triangles.map((t) => t.V.join(',')))
                    .toEqual(triangles);
            }, 20);
    });

    it('orientTriangles aligns every normal with the image gradient', () => {
        check(twoSphereArb, (voxels) => {
            const asc = new AdaptiveSkeletonClimbing3(N, voxels);
            const result = asc.extract(0.5, N);
            asc.makeUnique(result.vertices, result.triangles);
            const before = result.triangles.map((t) =>
                [...t.V].sort((a, b) => a - b).join(','));
            const gradient = (p: Vertex) => (asc as unknown as {
                getGradient(q: Vertex): Vertex;
            }).getGradient(p);

            // The dot product of the triangle normal with the average image
            // gradient at its vertices, the quantity orientTriangles tests.
            const dotWithGradient = (t: TriangleKey) => {
                const v = [result.vertices[t.V[0]], result.vertices[t.V[1]],
                    result.vertices[t.V[2]]];
                const e1 = [0, 1, 2].map((k) => v[1][k] - v[0][k]);
                const e2 = [0, 1, 2].map((k) => v[2][k] - v[0][k]);
                const n = [e1[1] * e2[2] - e1[2] * e2[1],
                    e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0]];
                const g = v.map(gradient);
                const avr = [0, 1, 2].map((k) =>
                    (g[0][k] + g[1][k] + g[2][k]) / 3);
                return avr[0] * n[0] + avr[1] * n[1] + avr[2] * n[2];
            };

            asc.orientTriangles(result.vertices, result.triangles, true);
            const same = result.triangles.map((t) => t.V.join(','));
            for (const t of result.triangles) {
                expect(dotWithGradient(t)).toBeGreaterThanOrEqual(0);
            }
            // Reorienting an already oriented mesh is a no-op.
            asc.orientTriangles(result.vertices, result.triangles, true);
            expect(result.triangles.map((t) => t.V.join(','))).toEqual(same);

            asc.orientTriangles(result.vertices, result.triangles, false);
            for (const t of result.triangles) {
                expect(dotWithGradient(t)).toBeLessThanOrEqual(0);
            }
            const opposite = result.triangles.map((t) => t.V.join(','));
            asc.orientTriangles(result.vertices, result.triangles, false);
            expect(result.triangles.map((t) => t.V.join(','))).toEqual(
                opposite);

            for (let i = 0; i < same.length; ++i) {
                const a = same[i].split(',');
                // Each triangle either keeps its winding, when the dot
                // product is exactly zero and neither branch swaps, or has
                // its last two indices exchanged.
                expect([same[i], [a[0], a[2], a[1]].join(',')])
                    .toContain(opposite[i]);
            }
            // The winding changes but the vertex set of each triangle does
            // not.
            expect(result.triangles.map((t) =>
                [...t.V].sort((a, b) => a - b).join(','))).toEqual(before);
        }, 20);
    });

    it('computeNormals returns unit normals, or zero for unused vertices',
        () => {
            check(twoSphereArb, (voxels) => {
                const asc = new AdaptiveSkeletonClimbing3(N, voxels);
                const result = asc.extract(0.5, N);
                asc.makeUnique(result.vertices, result.triangles);
                asc.orientTriangles(result.vertices, result.triangles, true);
                const normals = asc.computeNormals(result.vertices,
                    result.triangles);
                expect(normals.length).toBe(result.vertices.length);
                for (const n of normals) {
                    const length = Math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2);
                    expect(Number.isFinite(length)).toBe(true);
                    expect(length === 0 || Math.abs(length - 1) < 1e-12)
                        .toBe(true);
                }
            }, 20);
        });

    it('the face determinant is exact for large voxel values', () => {
        // The upstream disambiguation computes the int64_t determinant
        // f00 * f11 - f01 * f10 of the four face values. The port uses
        // bigint: for voxel magnitudes near 2^30 the two products exceed
        // 2^53 and a double computation loses the difference completely.
        const a = 2 ** 30 + 1;
        const b = 2 ** 30 + 2;
        const c = 2 ** 30;
        // The exact determinant is a*a - c*b = 1, but in doubles it is 0.
        expect(a * a - c * b).toBe(0);

        // The zmin face of the box at the origin has the four corner values
        // (f00, f10, f11, f01) = (a, -b, a, -c), so all four of its edges
        // are crossed and the determinant is the tiny positive number above.
        // The exact sign selects the disjoint hyperbolic pairing; a
        // determinant of zero would instead insert a face branch point.
        const voxels = makeImage(1, (x, y, z) => {
            if (z > 0 || x > 1 || y > 1) {
                return -1;
            }
            if (x === 0 && y === 0) {
                return a;
            }
            if (x === 1 && y === 0) {
                return -b;
            }
            if (x === 1 && y === 1) {
                return a;
            }
            return -c;   // (0,1,0)
        });

        const asc = new AdaptiveSkeletonClimbing3(1, voxels);
        const result = asc.extract(0.5, 1);
        asc.makeUnique(result.vertices, result.triangles);
        checkIndices(result.vertices, result.triangles);
        // No branch point: every vertex is on a grid edge.
        expect(countCompositeVertices(result.vertices)).toBe(0);
        expect(result.triangles.length).toBeGreaterThan(0);
    });

    it('getGradient is zero outside the voxels and trilinear inside', () => {
        const voxels = makeImage(N, (x, y, z) => 2 * x + 3 * y + 5 * z);
        const asc = new AdaptiveSkeletonClimbing3(N, voxels);
        const gradient = (p: Vertex) => (asc as unknown as {
            getGradient(q: Vertex): Vertex;
        }).getGradient(p);

        // A linear field is reproduced exactly by trilinear interpolation.
        expect(gradient([0.5, 0.5, 0.5])).toEqual([2, 3, 5]);
        expect(gradient([1.25, 2.75, 7.5])).toEqual([2, 3, 5]);

        // Outside the voxels of the image the gradient is zero. Upstream
        // truncates toward zero and rejects only strictly negative integer
        // parts, so a position in (-1, 0) still extrapolates from the first
        // voxel; the port keeps that.
        expect(gradient([-1.5, 0.5, 0.5])).toEqual([0, 0, 0]);
        expect(gradient([0.5, -1.5, 0.5])).toEqual([0, 0, 0]);
        expect(gradient([0.5, 0.5, -1.5])).toEqual([0, 0, 0]);
        expect(gradient([-0.5, 0.5, 0.5])).toEqual([2, 3, 5]);
        expect(gradient([size - 1, 0.5, 0.5])).toEqual([0, 0, 0]);
        expect(gradient([0.5, size - 1, 0.5])).toEqual([0, 0, 0]);
        expect(gradient([0.5, 0.5, size - 1])).toEqual([0, 0, 0]);
    });

    it('accepts a typed-array image and gives the same mesh as an array',
        () => {
            check(twoSphereArb, (voxels) => {
                const typed = Int32Array.from(voxels);
                const fromArray = new AdaptiveSkeletonClimbing3(N, voxels)
                    .extract(0.5, N);
                const fromTyped = new AdaptiveSkeletonClimbing3(N, typed)
                    .extract(0.5, N);
                expect(fromTyped.vertices).toEqual(fromArray.vertices);
                expect(fromTyped.triangles.map((t) => t.V.join(',')))
                    .toEqual(fromArray.triangles.map((t) => t.V.join(',')));
            }, 20);
        });
});
