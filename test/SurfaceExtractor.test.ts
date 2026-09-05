import { describe, it, expect } from 'vitest';
import {
    SurfaceExtractor,
    SurfaceExtractorTriangle,
    SurfaceExtractorVertex
} from '../src/SurfaceExtractor.js';
import { MarchingCubes } from '../src/MarchingCubes.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

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
            new SurfaceExtractorTriangle(0, 1, 3),  // stored [0,1,3] -> [0,1,2]
            new SurfaceExtractorTriangle(2, 1, 3),  // stored [1,3,2] -> [1,2,0]
            new SurfaceExtractorTriangle(3, 1, 0)   // stored [0,3,1] -> [0,2,1]
        ];
        ex.makeUnique(vertices, triangles);

        expect(vertices.length).toBe(3);
        expect(vertices[0]).toBe(A);
        expect(vertices[1]).toBe(B);
        expect(vertices[2]).toBe(C);
        // Upstream dedups on the remapped triple as it stands and does not
        // rebuild the min-first rotation, so the second triangle -- a
        // rotation of the first, hence the same oriented triangle -- is kept
        // as distinct. The third has the opposite winding and is distinct
        // under any convention.
        expect(triangles.length).toBe(3);
        expect(triangles[0].v).toEqual([0, 1, 2]);
        expect(triangles[1].v).toEqual([1, 2, 0]);
        expect(triangles[2].v).toEqual([0, 2, 1]);
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

describe('SurfaceExtractor verification', () => {
    // Small nonnegative rationals, so duplicates of the same value written
    // with different numerators and denominators occur often.
    const rational = fc.tuple(fc.integer({ min: 0, max: 3 }),
        fc.integer({ min: 1, max: 3 }));
    const vertexArb = fc.tuple(rational, rational, rational).map(
        ([[xn, xd], [yn, yd], [zn, zd]]) =>
            new SurfaceExtractorVertex(xn, xd, yn, yd, zn, zd));

    const valueOf = (v: SurfaceExtractorVertex): [number, number, number] =>
        [v.xNumer / v.xDenom, v.yNumer / v.yDenom, v.zNumer / v.zDenom];

    it('Vertex comparisons agree with the rational values they encode', () => {
        check(fc.tuple(vertexArb, vertexArb), ([a, b]) => {
            const va = valueOf(a), vb = valueOf(b);
            const lexLess = (va[0] !== vb[0] ? va[0] < vb[0]
                : va[1] !== vb[1] ? va[1] < vb[1]
                    : va[2] < vb[2]);
            expect(a.lessThan(b)).toBe(lexLess);
            expect(a.equals(b)).toBe(va[0] === vb[0] && va[1] === vb[1]
                && va[2] === vb[2]);
            // Trichotomy: exactly one of <, ==, > holds.
            const count = Number(a.lessThan(b)) + Number(a.equals(b))
                + Number(b.lessThan(a));
            expect(count).toBe(1);
        });
    });

    it('Vertex lessThan is a strict weak ordering', () => {
        check(fc.tuple(vertexArb, vertexArb, vertexArb), ([a, b, c]) => {
            expect(a.lessThan(a)).toBe(false);
            if (a.lessThan(b) && b.lessThan(c)) {
                expect(a.lessThan(c)).toBe(true);
            }
            // Incomparability is transitive because it coincides with equals.
            const incomparable = (p: SurfaceExtractorVertex,
                q: SurfaceExtractorVertex) => !p.lessThan(q) && !q.lessThan(p);
            if (incomparable(a, b) && incomparable(b, c)) {
                expect(incomparable(a, c)).toBe(true);
            }
        });
    });

    it('Triangle stores the min-first rotation and preserves the winding', () => {
        check(fc.tuple(fc.integer({ min: -5, max: 5 }),
            fc.integer({ min: -5, max: 5 }),
            fc.integer({ min: -5, max: 5 })), ([v0, v1, v2]) => {
                const t = new SurfaceExtractorTriangle(v0, v1, v2);
                const rotations = [[v0, v1, v2], [v1, v2, v0], [v2, v0, v1]];
                // The stored triple is one of the three cyclic rotations, so
                // the winding is unchanged.
                expect(rotations.some(r => r[0] === t.v[0] && r[1] === t.v[1]
                    && r[2] === t.v[2])).toBe(true);
                // With distinct indices the first entry is the minimum.
                if (v0 !== v1 && v1 !== v2 && v0 !== v2) {
                    expect(t.v[0]).toBe(Math.min(v0, v1, v2));
                }
            });
    });

    // An independent reimplementation of upstream MakeUnique: unique vertices
    // in first-occurrence order, triangle indices remapped, and duplicate
    // remapped triples dropped WITHOUT re-canonicalizing the rotation (which
    // is what upstream's std::map<Triangle,int32_t> does, since
    // Triangle::operator< is lexicographic on the stored triple).
    function referenceMakeUnique(vertices: SurfaceExtractorVertex[],
        triangles: SurfaceExtractorTriangle[]): {
            vertices: SurfaceExtractorVertex[];
            triangles: number[][];
        } {
        if (vertices.length === 0 || triangles.length === 0) {
            return {
                vertices: vertices.slice(),
                triangles: triangles.map(t => [t.v[0], t.v[1], t.v[2]])
            };
        }
        const unique: SurfaceExtractorVertex[] = [];
        const remap: number[] = [];
        for (const v of vertices) {
            let found = unique.findIndex(u => u.equals(v));
            if (found < 0) {
                found = unique.length;
                unique.push(v);
            }
            remap.push(found);
        }
        const seen = new Set<string>();
        const outTriangles: number[][] = [];
        for (const t of triangles) {
            const v = [remap[t.v[0]], remap[t.v[1]], remap[t.v[2]]];
            const key = v.join(',');
            if (!seen.has(key)) {
                seen.add(key);
                outTriangles.push(v);
            }
        }
        return { vertices: unique, triangles: outTriangles };
    }

    it('makeUnique matches an independent implementation', () => {
        const mesh = fc.tuple(
            fc.array(vertexArb, { minLength: 1, maxLength: 8 }),
            fc.array(fc.tuple(fc.nat({ max: 100 }), fc.nat({ max: 100 }),
                fc.nat({ max: 100 })), { minLength: 1, maxLength: 6 }))
            .map(([vs, raw]) => ({
                vertices: vs,
                triples: raw.map(([a, b, c]) =>
                    [a % vs.length, b % vs.length, c % vs.length] as
                    [number, number, number])
            }));
        check(mesh, ({ vertices, triples }) => {
            const makeTriangles = () => triples.map(([a, b, c]) =>
                new SurfaceExtractorTriangle(a, b, c));
            const actualVertices = vertices.slice();
            const actualTriangles = makeTriangles();
            const extractor = new StubExtractor(2, 2, 2, new Array<number>(8).fill(0));
            extractor.makeUnique(actualVertices, actualTriangles);

            const expected = referenceMakeUnique(vertices.slice(), makeTriangles());
            expect(actualVertices.length).toBe(expected.vertices.length);
            for (let i = 0; i < expected.vertices.length; ++i) {
                expect(actualVertices[i].equals(expected.vertices[i])).toBe(true);
            }
            expect(actualTriangles.map(t => [t.v[0], t.v[1], t.v[2]]))
                .toEqual(expected.triangles);

            // Post-conditions: the surviving vertices are pairwise distinct
            // and every triangle index is in range.
            for (let i = 0; i < actualVertices.length; ++i) {
                for (let j = i + 1; j < actualVertices.length; ++j) {
                    expect(actualVertices[i].equals(actualVertices[j])).toBe(false);
                }
            }
            for (const t of actualTriangles) {
                for (const v of t.v) {
                    expect(v >= 0 && v < actualVertices.length).toBe(true);
                }
            }
        });
    });

    it('makeUnique keeps the remapped rotation upstream produces', () => {
        // Upstream remaps triangle.v[] in place and then dedups on the
        // remapped triple as it stands; it does not rebuild the min-first
        // rotation. Two triangles whose remapped triples are rotations of one
        // another therefore both survive, and each keeps its own rotation.
        const A = new SurfaceExtractorVertex(0, 1, 0, 1, 0, 1);
        const B = new SurfaceExtractorVertex(1, 2, 0, 1, 0, 1);
        const C = new SurfaceExtractorVertex(0, 1, 1, 2, 0, 1);
        // The same three points written with different numerators and
        // denominators, which the rational comparison must see as equal.
        const A2 = new SurfaceExtractorVertex(0, 3, 0, 5, 0, 7);
        const B2 = new SurfaceExtractorVertex(2, 4, 0, 1, 0, 1);
        const C2 = new SurfaceExtractorVertex(0, 1, 3, 6, 0, 1);
        const vertices = [A, B, C, A2, B2, C2];
        // (0,1,2) is stored as [0,1,2]; (3,1,5) is stored as [1,5,3] and
        // remaps to [1,2,0], the rotation of [0,1,2] starting at 1.
        const triangles = [
            new SurfaceExtractorTriangle(0, 1, 2),
            new SurfaceExtractorTriangle(3, 1, 5)
        ];
        const extractor = new StubExtractor(2, 2, 2, new Array<number>(8).fill(0));
        extractor.makeUnique(vertices, triangles);

        expect(vertices.length).toBe(3);
        expect(triangles.length).toBe(2);
        expect(triangles[0].v).toEqual([0, 1, 2]);
        expect(triangles[1].v).toEqual([1, 2, 0]);
    });

    it('convert divides each rational component exactly', () => {
        check(fc.array(vertexArb, { minLength: 1, maxLength: 6 }), (vs) => {
            const extractor = new StubExtractor(2, 2, 2, new Array<number>(8).fill(0));
            const out = extractor.convert(vs);
            expect(out.length).toBe(vs.length);
            for (let i = 0; i < vs.length; ++i) {
                expect(out[i]).toEqual(valueOf(vs[i]));
            }
        });
    });

    it('extracts a closed orientable sphere for several radii and centers', () => {
        // A level surface of a sphere is a closed manifold, so every edge is
        // shared by exactly two triangles and V - E + F = 2. The extraction
        // itself uses no topology, so this checks the vertex identification
        // in makeUnique together with the marching-cubes triangulation.
        const bound = 9;
        for (const c of [4, 4.5]) {
            for (const r2x4 of [17, 25, 41, 61]) {
                const data = new Array<number>(bound * bound * bound);
                let i = 0;
                for (let z = 0; z < bound; ++z) {
                    for (let y = 0; y < bound; ++y) {
                        for (let x = 0; x < bound; ++x, ++i) {
                            data[i] = 4 * ((x - c) * (x - c) + (y - c) * (y - c)
                                + (z - c) * (z - c)) - r2x4;
                        }
                    }
                }
                // The tables exclude zero samples, and a sphere that pokes
                // out of the grid is clipped into an open surface.
                if (data.some(v => v === 0)
                    || c + Math.sqrt(r2x4 / 4) > bound - 1.05) {
                    continue;
                }
                const extractor = new MCExtractor(bound, bound, bound, data);
                const { vertices, triangles } = extractor.extract(0, true);
                expect(triangles.length).toBeGreaterThan(0);

                const edges = new Map<string, number>();
                for (const t of triangles) {
                    for (let k = 0; k < 3; ++k) {
                        const a = t.v[k], b = t.v[(k + 1) % 3];
                        const key = (a < b ? a + ':' + b : b + ':' + a);
                        edges.set(key, (edges.get(key) ?? 0) + 1);
                    }
                }
                for (const count of edges.values()) {
                    expect(count).toBe(2);
                }
                expect(vertices.length - edges.size + triangles.length).toBe(2);

                // Every vertex lies on the extracted level surface: the
                // trilinear field along the lattice edge vanishes there, so
                // the point is within half a voxel of the true sphere.
                const radius = Math.sqrt(r2x4 / 4);
                for (const v of vertices) {
                    const d = Math.hypot(v[0] - c, v[1] - c, v[2] - c);
                    expect(Math.abs(d - radius)).toBeLessThan(0.5);
                }

                // orientTriangles(sameDir = true) makes every triangle normal
                // agree with the average image gradient, which points outward
                // for this field.
                extractor.orientTriangles(vertices, triangles, true);
                for (const t of triangles) {
                    const n = triangleNormal(vertices[t.v[0]], vertices[t.v[1]],
                        vertices[t.v[2]]);
                    const centroid = [0, 1, 2].map(k =>
                        (vertices[t.v[0]][k] + vertices[t.v[1]][k]
                            + vertices[t.v[2]][k]) / 3);
                    const outward = [centroid[0] - c, centroid[1] - c, centroid[2] - c];
                    expect(n[0] * outward[0] + n[1] * outward[1] + n[2] * outward[2])
                        .toBeGreaterThanOrEqual(0);
                }

                const normals = extractor.computeNormals(vertices, triangles);
                for (const n of normals) {
                    expectClose(Math.hypot(n[0], n[1], n[2]), 1, 1e-12, 1e-12);
                }
            }
        }
    }, 30000);
});
