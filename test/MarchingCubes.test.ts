import { describe, it, expect } from 'vitest';
import { MarchingCubes } from '../src/MarchingCubes';

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
