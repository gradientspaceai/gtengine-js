import { describe, it, expect } from 'vitest';
import { MarchingCubes } from '../src/MarchingCubes.js';
import { check, fc } from './helpers/arbitraries.js';

describe('MarchingCubes', () => {
    const mc = new MarchingCubes();

    it('has constants matching upstream Topology limits', () => {
        expect(MarchingCubes.maxVertices).toBe(12);
        expect(MarchingCubes.maxTriangles).toBe(5);
    });

    it('generates a table that matches the upstream prebuilt table exactly', () => {
        const prebuilt = MarchingCubes.getPrebuiltTable();
        expect(prebuilt.length).toBe(256);
        for (let entry = 0; entry < 256; ++entry) {
            const topology = mc.getTable(entry);
            const row = prebuilt[entry];
            expect(row.length).toBe(41);
            const flat: number[] = [topology.numVertices, topology.numTriangles];
            for (const pair of topology.vpair) {
                flat.push(pair[0], pair[1]);
            }
            for (const triple of topology.itriple) {
                flat.push(triple[0], triple[1], triple[2]);
            }
            expect(flat, `entry ${entry}`).toEqual([...row]);
        }
    });

    it('getFlatTable returns the 256 x 41 table flattened row-major', () => {
        const flat = mc.getFlatTable();
        expect(flat.length).toBe(256 * 41);
        const prebuilt = MarchingCubes.getPrebuiltTable();
        for (let entry = 0; entry < 256; ++entry) {
            expect(flat.slice(41 * entry, 41 * (entry + 1))).toEqual([...prebuilt[entry]]);
        }
    });

    it('satisfies the row invariants for all 256 configurations', () => {
        for (let entry = 0; entry < 256; ++entry) {
            const t = mc.getTable(entry);

            // Counts in valid ranges.
            expect(t.numVertices).toBeGreaterThanOrEqual(0);
            expect(t.numVertices).toBeLessThanOrEqual(12);
            expect(t.numTriangles).toBeGreaterThanOrEqual(0);
            expect(t.numTriangles).toBeLessThanOrEqual(5);
            if (entry === 0 || entry === 255) {
                expect(t.numVertices).toBe(0);
                expect(t.numTriangles).toBe(0);
            } else {
                expect(t.numVertices).toBeGreaterThanOrEqual(3);
                expect(t.numTriangles).toBeGreaterThanOrEqual(1);
            }

            expect(t.vpair.length).toBe(12);
            for (let i = 0; i < 12; ++i) {
                const [a, b] = t.vpair[i];
                if (i < t.numVertices) {
                    // Corner indices in bounds, minimum first.
                    expect(a).toBeGreaterThanOrEqual(0);
                    expect(b).toBeLessThanOrEqual(7);
                    expect(a).toBeLessThan(b);
                    // The corners must be joined by a voxel edge (they
                    // differ in exactly one coordinate bit).
                    expect([1, 2, 4]).toContain(a ^ b);
                    // A surface vertex lies on an edge whose endpoints have
                    // opposite signs in this configuration.
                    expect((entry >> a) & 1, `entry ${entry} pair ${i}`)
                        .not.toBe((entry >> b) & 1);
                } else {
                    expect([a, b]).toEqual([0, 0]);
                }
            }

            expect(t.itriple.length).toBe(5);
            for (let i = 0; i < 5; ++i) {
                const [v0, v1, v2] = t.itriple[i];
                if (i < t.numTriangles) {
                    // Vertex indices in bounds and distinct.
                    for (const v of [v0, v1, v2]) {
                        expect(v).toBeGreaterThanOrEqual(0);
                        expect(v).toBeLessThan(t.numVertices);
                    }
                    expect(v0).not.toBe(v1);
                    expect(v1).not.toBe(v2);
                    expect(v0).not.toBe(v2);
                } else {
                    expect([v0, v1, v2]).toEqual([0, 0, 0]);
                }
            }
        }
    });

    it('matches known configurations (spot checks)', () => {
        // Entry 1: only corner 0 is negative; three vertices on the edges
        // from corner 0, one triangle.
        const t1 = mc.getTable(1);
        expect(t1.numVertices).toBe(3);
        expect(t1.numTriangles).toBe(1);
        expect(t1.vpair.slice(0, 3)).toEqual([[0, 1], [0, 4], [0, 2]]);
        expect(t1.itriple[0]).toEqual([0, 1, 2]);

        // Entry 254 is the complement of entry 1: same vertices, triangle
        // with opposite orientation.
        const t254 = mc.getTable(254);
        expect(t254.numVertices).toBe(3);
        expect(t254.numTriangles).toBe(1);
        expect(t254.vpair.slice(0, 3)).toEqual([[0, 1], [0, 4], [0, 2]]);
        expect(t254.itriple[0]).toEqual([0, 2, 1]);

        // Entry 3 (corners 0,1 negative): a quad split into two triangles.
        const t3 = mc.getTable(3);
        expect(t3.numVertices).toBe(4);
        expect(t3.numTriangles).toBe(2);
        expect(t3.vpair.slice(0, 4)).toEqual([[0, 4], [0, 2], [1, 3], [1, 5]]);
        expect(t3.itriple.slice(0, 2)).toEqual([[0, 1, 2], [0, 2, 3]]);

        // Entry 105 (0b01101001, corners 0,3,5,6): four separated corners,
        // 12 vertices and 4 triangles.
        const t105 = mc.getTable(105);
        expect(t105.numVertices).toBe(12);
        expect(t105.numTriangles).toBe(4);

        // Entry 15 (bottom face negative): four edge-crossings, two
        // triangles.
        const t15 = mc.getTable(15);
        expect(t15.numVertices).toBe(4);
        expect(t15.numTriangles).toBe(2);
        expect(t15.vpair.slice(0, 4)).toEqual([[0, 4], [2, 6], [3, 7], [1, 5]]);
    });

    it('reports configuration types by name', () => {
        expect(MarchingCubes.getConfigurationType(0)).toBe('Bits0');
        expect(MarchingCubes.getConfigurationType(255)).toBe('Bits0');
        expect(MarchingCubes.getConfigurationType(1)).toBe('Bits1');
        expect(MarchingCubes.getConfigurationType(254)).toBe('Bits7');
        expect(MarchingCubes.getConfigurationType(3)).toBe('Bits2Edge');
        expect(MarchingCubes.getConfigurationType(0b00011000)).toBe('Bits2BoxDiag');
        expect(MarchingCubes.getConfigurationType(0b00000110)).toBe('Bits2FaceDiag');
        expect(MarchingCubes.getConfigurationType(15)).toBe('Bits4SameFace');
        expect(MarchingCubes.getConfigurationType(105)).toBe('Bits4EdgeEdgePerp');
        expect(MarchingCubes.getConfigurationType(0b01100110)).toBe('Bits4EdgeEdgePara');
        expect(MarchingCubes.getConfigurationType(256)).toBe('');
        expect(MarchingCubes.getConfigurationType(-1)).toBe('');
    });

    it('complementary configurations have the same vertex sets', () => {
        // Entry and 255-entry have the same sign-change edges, so the same
        // vertex pair multisets (triangulation may differ due to ambiguous
        // face handling).
        const key = (pairs: number[][], n: number): string =>
            pairs.slice(0, n).map((p) => `${p[0]}-${p[1]}`).sort().join(',');
        for (let entry = 0; entry < 256; ++entry) {
            const t = mc.getTable(entry);
            const tc = mc.getTable(255 - entry);
            expect(t.numVertices).toBe(tc.numVertices);
            expect(key(t.vpair, t.numVertices))
                .toBe(key(tc.vpair, tc.numVertices));
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream MarchingCubes.h lookup table.
// ---------------------------------------------------------------------------

// Extract a level surface from a sampled scalar field with the port's table,
// welding vertices by the grid edge they lie on. The bit assignment is the
// upstream one: bit i of the entry is set when the corner (x + (i & 1),
// y + ((i >> 1) & 1), z + ((i >> 2) & 1)) has a negative value.
function extractSurface(mc: MarchingCubes, n: number, f: Float64Array) {
    const index = (x: number, y: number, z: number) => x + n * (y + n * z);
    const vmap = new Map<string, number>();
    const vertices: [number, number, number][] = [];
    const triangles: [number, number, number][] = [];

    for (let z = 0; z + 1 < n; ++z) {
        for (let y = 0; y + 1 < n; ++y) {
            for (let x = 0; x + 1 < n; ++x) {
                const corner = (i: number): [number, number, number] =>
                    [x + (i & 1), y + ((i >> 1) & 1), z + ((i >> 2) & 1)];
                let entry = 0;
                for (let i = 0; i < 8; ++i) {
                    const c = corner(i);
                    if (f[index(c[0], c[1], c[2])] < 0) {
                        entry |= (1 << i);
                    }
                }
                const topology = mc.getTable(entry);
                const local: number[] = [];
                for (let i = 0; i < topology.numVertices; ++i) {
                    const [a, b] = topology.vpair[i];
                    const ca = corner(a), cb = corner(b);
                    const ia = index(ca[0], ca[1], ca[2]);
                    const ib = index(cb[0], cb[1], cb[2]);
                    const key = ia < ib ? `${ia}_${ib}` : `${ib}_${ia}`;
                    let vi = vmap.get(key);
                    if (vi === undefined) {
                        const t = f[ia] / (f[ia] - f[ib]);
                        vi = vertices.length;
                        vertices.push([
                            ca[0] + t * (cb[0] - ca[0]),
                            ca[1] + t * (cb[1] - ca[1]),
                            ca[2] + t * (cb[2] - ca[2])]);
                        vmap.set(key, vi);
                    }
                    local.push(vi);
                }
                for (let i = 0; i < topology.numTriangles; ++i) {
                    const [p, q, r] = topology.itriple[i];
                    triangles.push([local[p], local[q], local[r]]);
                }
            }
        }
    }
    return { vertices, triangles };
}

// Edge bookkeeping for a triangle soup: a closed oriented manifold uses each
// directed edge exactly once and each undirected edge exactly twice.
function edgeCounts(triangles: readonly [number, number, number][]) {
    const directed = new Map<string, number>();
    const undirected = new Map<string, number>();
    for (const [a, b, c] of triangles) {
        for (const [u, v] of [[a, b], [b, c], [c, a]] as [number, number][]) {
            directed.set(`${u}_${v}`, (directed.get(`${u}_${v}`) ?? 0) + 1);
            const key = u < v ? `${u}_${v}` : `${v}_${u}`;
            undirected.set(key, (undirected.get(key) ?? 0) + 1);
        }
    }
    return {
        directed, undirected,
        repeatedDirected: [...directed.values()].filter(v => v !== 1).length,
        nonManifold: [...undirected.values()].filter(v => v !== 2).length
    };
}

function signedVolume(vertices: readonly [number, number, number][],
    triangles: readonly [number, number, number][]): number {
    let volume = 0;
    for (const [a, b, c] of triangles) {
        const p = vertices[a], q = vertices[b], r = vertices[c];
        volume += (p[0] * (q[1] * r[2] - q[2] * r[1])
            - p[1] * (q[0] * r[2] - q[2] * r[0])
            + p[2] * (q[0] * r[1] - q[1] * r[0])) / 6;
    }
    return volume;
}

describe('MarchingCubes verification', () => {
    const mc = new MarchingCubes();

    it('lists exactly the sign-changing voxel edges, once each', () => {
        // The vertex pairs of an entry must be the complete set of the 12
        // voxel edges whose two corner signs differ -- no missing edge and
        // no duplicate. This pins every row of the configuration table and
        // every Bits* generator against the sign pattern it belongs to.
        for (let entry = 0; entry < 256; ++entry) {
            const expected: string[] = [];
            for (let a = 0; a < 8; ++a) {
                for (const bit of [1, 2, 4]) {
                    const b = a ^ bit;
                    if (b > a && ((entry >> a) & 1) !== ((entry >> b) & 1)) {
                        expected.push(`${a},${b}`);
                    }
                }
            }
            const topology = mc.getTable(entry);
            const got: string[] = [];
            for (let i = 0; i < topology.numVertices; ++i) {
                got.push(`${topology.vpair[i][0]},${topology.vpair[i][1]}`);
            }
            expect(got.sort(), `entry ${entry}`).toEqual(expected.sort());
        }
    });

    it('gives every configuration a manifold triangle patch', () => {
        // Within one voxel no directed edge may repeat (consistent winding)
        // and no undirected edge may be shared by more than two triangles.
        for (let entry = 0; entry < 256; ++entry) {
            const topology = mc.getTable(entry);
            const triangles = topology.itriple.slice(0, topology.numTriangles)
                .map(t => [t[0], t[1], t[2]] as [number, number, number]);
            const counts = edgeCounts(triangles);
            expect(counts.repeatedDirected, `entry ${entry}`).toBe(0);
            expect([...counts.undirected.values()].every(v => v <= 2),
                `entry ${entry}`).toBe(true);
            // Every vertex of the entry is used by at least one triangle.
            const used = new Set<number>();
            for (const [a, b, c] of triangles) { used.add(a); used.add(b); used.add(c); }
            expect(used.size, `entry ${entry}`).toBe(topology.numVertices);
        }
    });

    it('closes each patch on the voxel faces', () => {
        // A boundary edge of the patch (one used by a single triangle) must
        // join two surface vertices that lie on voxel edges of a common
        // voxel face. That is what lets the patches of adjacent voxels meet,
        // and it fails immediately if a triangle triple indexes the wrong
        // vertex pair.
        const inFace = (edge: readonly number[], axis: number, side: number) =>
            edge.every(c => ((c >> axis) & 1) === side);
        for (let entry = 0; entry < 256; ++entry) {
            const topology = mc.getTable(entry);
            const triangles = topology.itriple.slice(0, topology.numTriangles)
                .map(t => [t[0], t[1], t[2]] as [number, number, number]);
            const counts = edgeCounts(triangles);
            for (const [key, count] of counts.undirected) {
                if (count !== 1) { continue; }
                const [u, v] = key.split('_').map(Number);
                const e0 = topology.vpair[u], e1 = topology.vpair[v];
                let onFace = false;
                for (let axis = 0; axis < 3 && !onFace; ++axis) {
                    for (const side of [0, 1]) {
                        if (inFace(e0, axis, side) && inFace(e1, axis, side)) {
                            onFace = true;
                            break;
                        }
                    }
                }
                expect(onFace, `entry ${entry} boundary edge ${key}`).toBe(true);
            }
        }
    });

    it('extracts a closed, inward-oriented sphere with Euler characteristic 2', () => {
        const n = 12;
        const center = (n - 1) / 2;
        const radius = 4.3;
        const f = new Float64Array(n * n * n);
        for (let z = 0; z < n; ++z) {
            for (let y = 0; y < n; ++y) {
                for (let x = 0; x < n; ++x) {
                    f[x + n * (y + n * z)] = (x - center) ** 2 + (y - center) ** 2
                        + (z - center) ** 2 - radius * radius;
                }
            }
        }
        // The table assumes no sample is zero.
        for (const value of f) { expect(value).not.toBe(0); }

        const { vertices, triangles } = extractSurface(mc, n, f);
        const counts = edgeCounts(triangles);
        expect(counts.repeatedDirected).toBe(0);
        expect(counts.nonManifold).toBe(0);
        expect(vertices.length - counts.undirected.size + triangles.length).toBe(2);

        // Upstream orients triangles counterclockwise as seen from the
        // negative side, that is with normals pointing into the sphere, so
        // the divergence-theorem volume is negative.
        const volume = signedVolume(vertices, triangles);
        expect(volume).toBeLessThan(0);
        const exact = (4 / 3) * Math.PI * radius ** 3;
        // Marching cubes on a quadratic field under-resolves the sphere; the
        // relative error at this resolution is a few percent.
        expect(Math.abs(-volume - exact) / exact).toBeLessThan(0.06);
    });

    it('extracts closed oriented surfaces from random two-sphere fields', () => {
        // Two overlapping or disjoint balls produce saddle configurations
        // that exercise the ambiguous-face entries of the table.
        const coordinate = fc.integer({ min: 30, max: 70 }).map(v => v / 10);
        check(fc.tuple(coordinate, coordinate, coordinate,
            coordinate, coordinate, coordinate,
            fc.integer({ min: 15, max: 30 }), fc.integer({ min: 15, max: 30 })),
            ([ax, ay, az, bx, by, bz, ra, rb]) => {
                const n = 12;
                const radiusA = ra / 10, radiusB = rb / 10;
                const f = new Float64Array(n * n * n);
                for (let z = 0; z < n; ++z) {
                    for (let y = 0; y < n; ++y) {
                        for (let x = 0; x < n; ++x) {
                            const da = (x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2
                                - radiusA * radiusA;
                            const db = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
                                - radiusB * radiusB;
                            // The offset keeps the level surface off the
                            // lattice; a sample that still lands on zero
                            // violates the table's precondition and the
                            // draw is skipped below.
                            f[x + n * (y + n * z)] = Math.min(da, db) + 1 / 3;
                        }
                    }
                }
                if (f.some(value => value === 0)) { return; }

                const { vertices, triangles } = extractSurface(mc, n, f);
                if (triangles.length === 0) { return; }
                const counts = edgeCounts(triangles);
                // Closed (no boundary edge) and consistently oriented.
                expect(counts.repeatedDirected).toBe(0);
                expect(counts.nonManifold).toBe(0);
                const euler = vertices.length - counts.undirected.size + triangles.length;
                expect(euler % 2).toBe(0);
                expect(euler).toBeGreaterThan(0);
                expect(signedVolume(vertices, triangles)).toBeLessThan(0);
            }, 40);
    });


});
