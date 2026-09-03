import { describe, it, expect } from 'vitest';
import {
    SurfaceExtractor,
    SurfaceExtractorTriangle,
    SurfaceExtractorVertex
} from '../src/SurfaceExtractor.js';
import { MarchingCubes } from '../src/MarchingCubes.js';

// Concrete extractor built on the ported MarchingCubes tables, mirroring the
// structure of the upstream derived classes. Voxel corner k of the cube with
// origin (x,y,z) is at (x + (k&1), y + ((k>>1)&1), z + ((k>>2)&1)), matching
// the MarchingCubes bit assignment. A level-crossing on a lattice edge with
// integer endpoint values F0, F1 (level already subtracted) is the exact
// rational point where the linear interpolant vanishes.
class MCExtractor extends SurfaceExtractor {
    private readonly mMC = new MarchingCubes();

    constructor(xBound: number, yBound: number, zBound: number,
        inputVoxels: ArrayLike<number>) {
        super(xBound, yBound, zBound, inputVoxels);
    }

    extractRational(level: number): {
        vertices: SurfaceExtractorVertex[];
        triangles: SurfaceExtractorTriangle[];
    } {
        const vertices: SurfaceExtractorVertex[] = [];
        const triangles: SurfaceExtractorTriangle[] = [];

        for (let i = 0; i < this.mVoxels.length; ++i) {
            this.mVoxels[i] = this.mInputVoxels[i] - level;
        }

        const F = new Array<number>(8);
        for (let z = 0; z + 1 < this.mZBound; ++z) {
            for (let y = 0; y + 1 < this.mYBound; ++y) {
                nextVoxel:
                for (let x = 0; x + 1 < this.mXBound; ++x) {
                    let entry = 0;
                    for (let k = 0; k < 8; ++k) {
                        const cx = x + (k & 1);
                        const cy = y + ((k >> 1) & 1);
                        const cz = z + ((k >> 2) & 1);
                        F[k] = this.mVoxels[cx + this.mXBound * cy + this.mXYBound * cz];
                        if (F[k] === 0) {
                            // The tables assume nonzero corner values.
                            continue nextVoxel;
                        }
                        if (F[k] < 0) {
                            entry |= (1 << k);
                        }
                    }

                    const topology = this.mMC.getTable(entry);
                    const vbase = vertices.length;
                    for (let i = 0; i < topology.numVertices; ++i) {
                        const j0 = topology.vpair[i][0];
                        const j1 = topology.vpair[i][1];
                        const d = F[j0] - F[j1];
                        const c0 = [x + (j0 & 1), y + ((j0 >> 1) & 1), z + ((j0 >> 2) & 1)];
                        const c1 = [x + (j1 & 1), y + ((j1 >> 1) & 1), z + ((j1 >> 2) & 1)];
                        const numer = [0, 0, 0];
                        const denom = [1, 1, 1];
                        for (let axis = 0; axis < 3; ++axis) {
                            if (c1[axis] === c0[axis]) {
                                numer[axis] = c0[axis];
                                denom[axis] = 1;
                            } else if (c1[axis] === c0[axis] + 1) {
                                // c0 + F0/(F0 - F1)
                                numer[axis] = c0[axis] * d + F[j0];
                                denom[axis] = d;
                            } else {
                                // c0 - F0/(F0 - F1)
                                numer[axis] = c0[axis] * d - F[j0];
                                denom[axis] = d;
                            }
                        }
                        vertices.push(new SurfaceExtractorVertex(
                            numer[0], denom[0], numer[1], denom[1], numer[2], denom[2]));
                    }

                    for (let i = 0; i < topology.numTriangles; ++i) {
                        triangles.push(new SurfaceExtractorTriangle(
                            vbase + topology.itriple[i][0],
                            vbase + topology.itriple[i][1],
                            vbase + topology.itriple[i][2]));
                    }
                }
            }
        }

        return { vertices, triangles };
    }

    // Trilinear-interpolation gradient of the input image, as in the
    // upstream SurfaceExtractorMC::GetGradient.
    protected getGradient(pos: [number, number, number]): [number, number, number] {
        const x = Math.floor(pos[0]);
        const y = Math.floor(pos[1]);
        const z = Math.floor(pos[2]);
        if (pos[0] < 0 || x + 1 >= this.mXBound
            || pos[1] < 0 || y + 1 >= this.mYBound
            || pos[2] < 0 || z + 1 >= this.mZBound) {
            return [0, 0, 0];
        }

        const dx = pos[0] - x, dy = pos[1] - y, dz = pos[2] - z;
        const omx = 1 - dx, omy = 1 - dy, omz = 1 - dz;
        const at = (ix: number, iy: number, iz: number) =>
            this.mInputVoxels[ix + this.mXBound * iy + this.mXYBound * iz];
        const f000 = at(x, y, z), f100 = at(x + 1, y, z);
        const f010 = at(x, y + 1, z), f110 = at(x + 1, y + 1, z);
        const f001 = at(x, y, z + 1), f101 = at(x + 1, y, z + 1);
        const f011 = at(x, y + 1, z + 1), f111 = at(x + 1, y + 1, z + 1);

        let tmp0 = omy * (f100 - f000) + dy * (f110 - f010);
        let tmp1 = omy * (f101 - f001) + dy * (f111 - f011);
        const gx = omz * tmp0 + dz * tmp1;

        tmp0 = omx * (f010 - f000) + dx * (f110 - f100);
        tmp1 = omx * (f011 - f001) + dx * (f111 - f101);
        const gy = omz * tmp0 + dz * tmp1;

        tmp0 = omx * (f001 - f000) + dx * (f101 - f100);
        tmp1 = omx * (f011 - f010) + dx * (f111 - f110);
        const gz = omy * tmp0 + dy * tmp1;

        return [gx, gy, gz];
    }
}

// A trivial concrete class for testing base-class utilities directly with a
// prescribed gradient field.
class StubExtractor extends MCExtractor {
    gradient: [number, number, number] = [0, 0, 1];

    protected override getGradient(): [number, number, number] {
        return this.gradient;
    }
}

// Integer sphere-like field f = 4*((x-c)^2 + (y-c)^2 + (z-c)^2) - 25, whose
// zero set is the sphere of radius 2.5 about (c,c,c). All samples are odd
// integers, so no lattice value equals 0.
function sphereField(bound: number, c: number): number[] {
    const data = new Array<number>(bound * bound * bound);
    let i = 0;
    for (let z = 0; z < bound; ++z) {
        for (let y = 0; y < bound; ++y) {
            for (let x = 0; x < bound; ++x, ++i) {
                data[i] = 4 * ((x - c) * (x - c) + (y - c) * (y - c) + (z - c) * (z - c)) - 25;
            }
        }
    }
    return data;
}

function triangleNormal(v0: number[], v1: number[], v2: number[]): number[] {
    const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    return [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0]
    ];
}

describe('SurfaceExtractorVertex', () => {
    it('normalizes signs so denominators are positive', () => {
        const v = new SurfaceExtractorVertex(-3, -4, 5, 2, -7, -8);
        expect(v.xNumer).toBe(3);
        expect(v.xDenom).toBe(4);
        expect(v.yNumer).toBe(5);
        expect(v.yDenom).toBe(2);
        expect(v.zNumer).toBe(7);
        expect(v.zDenom).toBe(8);
    });

    it('compares rationals with different denominators', () => {
        const a = new SurfaceExtractorVertex(1, 2, 3, 4, 5, 6);
        const b = new SurfaceExtractorVertex(2, 4, 6, 8, 10, 12);
        const c = new SurfaceExtractorVertex(1, 2, 3, 4, 6, 6);
        expect(a.equals(b)).toBe(true);
        expect(a.equals(c)).toBe(false);
        expect(a.lessThan(b)).toBe(false);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);   // 5/6 < 6/6
        expect(c.lessThan(a)).toBe(false);
        // Lexicographic: x dominates.
        const d = new SurfaceExtractorVertex(1, 3, 100, 1, 100, 1);
        expect(d.lessThan(a)).toBe(true);   // 1/3 < 1/2
    });
});

describe('SurfaceExtractorTriangle', () => {
    it('stores the min-first cyclic permutation, preserving winding', () => {
        expect(new SurfaceExtractorTriangle(1, 2, 3).v).toEqual([1, 2, 3]);
        expect(new SurfaceExtractorTriangle(3, 1, 2).v).toEqual([1, 2, 3]);
        expect(new SurfaceExtractorTriangle(2, 3, 1).v).toEqual([1, 2, 3]);
        // Opposite winding is preserved, not sorted.
        expect(new SurfaceExtractorTriangle(2, 1, 3).v).toEqual([1, 3, 2]);
        expect(new SurfaceExtractorTriangle(3, 2, 1).v).toEqual([1, 3, 2]);
    });

    it('compares lexicographically', () => {
        const a = new SurfaceExtractorTriangle(1, 2, 3);
        const b = new SurfaceExtractorTriangle(2, 3, 1);
        const c = new SurfaceExtractorTriangle(1, 3, 2);
        expect(a.equals(b)).toBe(true);
        expect(a.equals(c)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
    });
});

describe('SurfaceExtractor', () => {
    const voxels222 = [1, 1, 1, 1, 1, 1, 1, -1];

    it('validates constructor input', () => {
        expect(() => new MCExtractor(1, 2, 2, voxels222)).toThrow('Invalid input.');
        expect(() => new MCExtractor(2, 2, 1, voxels222)).toThrow('Invalid input.');
        expect(() => new MCExtractor(2, 2, 2, voxels222)).not.toThrow();
    });

    it('makeUnique removes duplicate vertices and triangles and remaps indices', () => {
        const ex = new MCExtractor(2, 2, 2, voxels222);
        const A = new SurfaceExtractorVertex(1, 2, 0, 1, 0, 1);
        const B = new SurfaceExtractorVertex(1, 1, 1, 2, 0, 1);
        const A2 = new SurfaceExtractorVertex(2, 4, 0, 2, 0, 3);  // equals A
        const C = new SurfaceExtractorVertex(0, 1, 1, 1, 1, 2);
        const vertices = [A, B, A2, C];
        const triangles = [
            new SurfaceExtractorTriangle(0, 1, 3),
            new SurfaceExtractorTriangle(2, 1, 3),  // same as (0,1,3) after remap
            new SurfaceExtractorTriangle(3, 1, 0)   // opposite winding: distinct
        ];
        ex.makeUnique(vertices, triangles);

        expect(vertices.length).toBe(3);
        expect(vertices[0]).toBe(A);
        expect(vertices[1]).toBe(B);
        expect(vertices[2]).toBe(C);
        expect(triangles.length).toBe(2);
        expect(triangles[0].v).toEqual([0, 1, 2]);
        expect(triangles[1].v).toEqual([0, 2, 1]);
    });

    it('convert produces floating-point triples', () => {
        const ex = new MCExtractor(2, 2, 2, voxels222);
        const out = ex.convert([
            new SurfaceExtractorVertex(1, 2, 3, 4, 9, 3),
            new SurfaceExtractorVertex(0, 1, 7, 2, 1, 8)
        ]);
        expect(out).toEqual([[0.5, 0.75, 3], [0, 3.5, 0.125]]);
    });

    it('computeNormals produces unit normals; zero for unused vertices', () => {
        const ex = new MCExtractor(2, 2, 2, voxels222);
        const vertices: [number, number, number][] = [
            [0, 0, 0], [2, 0, 0], [0, 2, 0], [5, 5, 5]
        ];
        const triangles = [new SurfaceExtractorTriangle(0, 1, 2)];
        const normals = ex.computeNormals(vertices, triangles);
        expect(normals.length).toBe(4);
        // Counterclockwise in the xy-plane: +z normal at each used vertex.
        for (let i = 0; i < 3; ++i) {
            expect(normals[i][0]).toBeCloseTo(0, 14);
            expect(normals[i][1]).toBeCloseTo(0, 14);
            expect(normals[i][2]).toBeCloseTo(1, 14);
        }
        expect(normals[3]).toEqual([0, 0, 0]);
    });

    it('orientTriangles flips windings against the gradient field', () => {
        const ex = new StubExtractor(2, 2, 2, voxels222);
        const vertices: [number, number, number][] = [[0, 0, 0], [1, 0, 0], [0, 1, 0]];
        // Normal of (0,1,2) is +z; the stub gradient is +z.
        const tSame = [new SurfaceExtractorTriangle(0, 1, 2)];
        ex.orientTriangles(vertices, tSame, true);
        expect(tSame[0].v).toEqual([0, 1, 2]);  // already aligned
        ex.orientTriangles(vertices, tSame, false);
        expect(tSame[0].v).toEqual([0, 2, 1]);  // flipped to oppose
        ex.orientTriangles(vertices, tSame, true);
        expect(tSame[0].v).toEqual([0, 1, 2]);  // flipped back
    });

    describe('sphere level-surface extraction', () => {
        const bound = 9;
        const center = 4;
        const radius = 2.5;
        const field = sphereField(bound, center);
        const ex = new MCExtractor(bound, bound, bound, field);
        const { vertices, triangles } = ex.extract(0, true);

        it('produces a nonempty mesh with valid indices', () => {
            expect(vertices.length).toBeGreaterThan(0);
            expect(triangles.length).toBeGreaterThan(0);
            for (const tri of triangles) {
                expect(tri.v[0]).not.toBe(tri.v[1]);
                expect(tri.v[1]).not.toBe(tri.v[2]);
                expect(tri.v[2]).not.toBe(tri.v[0]);
                for (const idx of tri.v) {
                    expect(idx).toBeGreaterThanOrEqual(0);
                    expect(idx).toBeLessThan(vertices.length);
                }
            }
        });

        it('vertices lie near the true sphere', () => {
            for (const v of vertices) {
                const r = Math.hypot(v[0] - center, v[1] - center, v[2] - center);
                expect(Math.abs(r - radius)).toBeLessThan(0.25);
            }
        });

        it('removes duplicate vertices shared by adjacent voxels', () => {
            const raw = ex.extract(0, false);
            expect(raw.triangles.length).toBe(triangles.length);
            expect(raw.vertices.length).toBeGreaterThan(vertices.length);
        });

        it('is a closed manifold with sphere topology (V - E + F = 2)', () => {
            const edges = new Map<string, number>();
            for (const tri of triangles) {
                for (let i = 0; i < 3; ++i) {
                    const a = tri.v[i], b = tri.v[(i + 1) % 3];
                    const key = a < b ? `${a},${b}` : `${b},${a}`;
                    edges.set(key, (edges.get(key) ?? 0) + 1);
                }
            }
            for (const count of edges.values()) {
                expect(count).toBe(2);
            }
            expect(vertices.length - edges.size + triangles.length).toBe(2);
        });

        it('has consistent winding: all normals point inward (negative side)', () => {
            // MarchingCubes triangles are counterclockwise for an observer on
            // the negative side of the level surface; the field is negative
            // inside the sphere, so every normal points toward the center.
            for (const tri of triangles) {
                const v0 = vertices[tri.v[0]];
                const v1 = vertices[tri.v[1]];
                const v2 = vertices[tri.v[2]];
                const n = triangleNormal(v0, v1, v2);
                const cx = (v0[0] + v1[0] + v2[0]) / 3 - center;
                const cy = (v0[1] + v1[1] + v2[1]) / 3 - center;
                const cz = (v0[2] + v1[2] + v2[2]) / 3 - center;
                const dot = n[0] * cx + n[1] * cy + n[2] * cz;
                expect(dot).toBeLessThan(0);
            }
        });

        it('orientTriangles(sameDir=true) aligns windings with the gradient', () => {
            const oriented = ex.extract(0, true);
            ex.orientTriangles(oriented.vertices, oriented.triangles, true);
            for (const tri of oriented.triangles) {
                const v0 = oriented.vertices[tri.v[0]];
                const v1 = oriented.vertices[tri.v[1]];
                const v2 = oriented.vertices[tri.v[2]];
                const n = triangleNormal(v0, v1, v2);
                // The gradient of the field is radially outward.
                const cx = (v0[0] + v1[0] + v2[0]) / 3 - center;
                const cy = (v0[1] + v1[1] + v2[1]) / 3 - center;
                const cz = (v0[2] + v1[2] + v2[2]) / 3 - center;
                const dot = n[0] * cx + n[1] * cy + n[2] * cz;
                expect(dot).toBeGreaterThan(0);
            }
        });

        it('computeNormals on the oriented mesh points outward', () => {
            const oriented = ex.extract(0, true);
            ex.orientTriangles(oriented.vertices, oriented.triangles, true);
            const normals = ex.computeNormals(oriented.vertices, oriented.triangles);
            for (let i = 0; i < normals.length; ++i) {
                const v = oriented.vertices[i];
                const rx = v[0] - center, ry = v[1] - center, rz = v[2] - center;
                const len = Math.hypot(rx, ry, rz);
                const dot = (normals[i][0] * rx + normals[i][1] * ry + normals[i][2] * rz) / len;
                // Unit normal nearly radial for a well-resolved sphere.
                expect(dot).toBeGreaterThan(0.8);
                expect(Math.hypot(...normals[i])).toBeCloseTo(1, 12);
            }
        });

        it('extracts a different level of the same field', () => {
            // level 10: 4 rho^2 - 35 = 0, radius sqrt(35)/2 ~ 2.958.
            const shifted = ex.extract(10, true);
            const target = Math.sqrt(35) / 2;
            expect(shifted.vertices.length).toBeGreaterThan(0);
            for (const v of shifted.vertices) {
                const r = Math.hypot(v[0] - center, v[1] - center, v[2] - center);
                expect(Math.abs(r - target)).toBeLessThan(0.25);
            }
        });

        it('rational vertices are exact on lattice edges', () => {
            const { vertices: rational } = ex.extractRational(0);
            for (const v of rational) {
                expect(v.xDenom).toBeGreaterThan(0);
                expect(v.yDenom).toBeGreaterThan(0);
                expect(v.zDenom).toBeGreaterThan(0);
                // Exactly two of the three coordinates are integers (the
                // vertex lies in the interior of a lattice edge; corner
                // values are never zero for this field).
                let integers = 0;
                if (v.xNumer % v.xDenom === 0) { ++integers; }
                if (v.yNumer % v.yDenom === 0) { ++integers; }
                if (v.zNumer % v.zDenom === 0) { ++integers; }
                expect(integers).toBe(2);
            }
        });
    });
});
