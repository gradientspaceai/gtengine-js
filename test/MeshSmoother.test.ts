import { describe, it, expect } from 'vitest';
import { MeshSmoother } from '../src/MeshSmoother.js';
import { Vector, add, dot, length, mul as mulVector, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function copyVerts(v: Vector[]): number[][] {
    return v.map(p => p.values.slice());
}

// A triangulated rectangular grid in the plane z = 0 with the given number
// of samples per direction and spacing 1.
function makeGrid(n: number, height: (i: number, j: number) => number):
    { vertices: Vector[], indices: number[], interior: number[], boundary: number[] } {
    const vertices: Vector[] = [];
    const interior: number[] = [];
    const boundary: number[] = [];
    for (let j = 0; j < n; ++j) {
        for (let i = 0; i < n; ++i) {
            vertices.push(v3(i, j, height(i, j)));
            const index = i + n * j;
            if (i > 0 && i < n - 1 && j > 0 && j < n - 1) {
                interior.push(index);
            } else {
                boundary.push(index);
            }
        }
    }
    const indices: number[] = [];
    for (let j = 0; j + 1 < n; ++j) {
        for (let i = 0; i + 1 < n; ++i) {
            const i00 = i + n * j;
            const i10 = i00 + 1;
            const i01 = i00 + n;
            const i11 = i10 + n;
            // Counterclockwise when viewed from +z.
            indices.push(i00, i10, i11);
            indices.push(i00, i11, i01);
        }
    }
    return { vertices, indices, interior, boundary };
}

// An independent reimplementation of one smoothing step, written directly
// from the mathematics rather than from the port.
function referenceUpdate(vertices: number[][], indices: number[],
    tanWeight: number, norWeight: number): void {
    const n = vertices.length;
    const normals: number[][] = [];
    const means: number[][] = [];
    const counts: number[] = new Array(n).fill(0);
    for (let i = 0; i < n; ++i) {
        normals.push([0, 0, 0]);
        means.push([0, 0, 0]);
    }
    for (let t = 0; t * 3 < indices.length; ++t) {
        counts[indices[3 * t]] += 2;
        counts[indices[3 * t + 1]] += 2;
        counts[indices[3 * t + 2]] += 2;
    }
    for (let t = 0; t * 3 < indices.length; ++t) {
        const a = indices[3 * t];
        const b = indices[3 * t + 1];
        const c = indices[3 * t + 2];
        const A = vertices[a], B = vertices[b], C = vertices[c];
        const e1 = [B[0] - A[0], B[1] - A[1], B[2] - A[2]];
        const e2 = [C[0] - A[0], C[1] - A[1], C[2] - A[2]];
        const nrm = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0]
        ];
        for (let k = 0; k < 3; ++k) {
            normals[a][k] += nrm[k];
            normals[b][k] += nrm[k];
            normals[c][k] += nrm[k];
            means[a][k] += B[k] + C[k];
            means[b][k] += C[k] + A[k];
            means[c][k] += A[k] + B[k];
        }
    }
    for (let i = 0; i < n; ++i) {
        const len = Math.hypot(normals[i][0], normals[i][1], normals[i][2]);
        if (len > 0) {
            for (let k = 0; k < 3; ++k) {
                normals[i][k] /= len;
            }
        }
        for (let k = 0; k < 3; ++k) {
            means[i][k] /= counts[i];
        }
    }
    for (let i = 0; i < n; ++i) {
        const diff = [
            means[i][0] - vertices[i][0],
            means[i][1] - vertices[i][1],
            means[i][2] - vertices[i][2]
        ];
        const dotDN = diff[0] * normals[i][0] + diff[1] * normals[i][1]
            + diff[2] * normals[i][2];
        for (let k = 0; k < 3; ++k) {
            const tangent = diff[k] - dotDN * normals[i][k];
            vertices[i][k] += tanWeight * tangent + norWeight * normals[i][k];
        }
    }
}

describe('MeshSmoother construction', () => {
    it('starts empty', () => {
        const smoother = new MeshSmoother();
        expect(smoother.getNumVertices()).toBe(0);
        expect(smoother.getNumTriangles()).toBe(0);
        expect(smoother.getVertices()).toEqual([]);
        expect(smoother.getNormals()).toEqual([]);
        expect(smoother.getMeans()).toEqual([]);
        expect(smoother.getNeighborCounts()).toEqual([]);
    });

    it('rejects meshes with too few vertices or no triangles', () => {
        const smoother = new MeshSmoother();
        expect(() => smoother.initialize([v3(0, 0, 0), v3(1, 0, 0)], [0, 1, 0]))
            .toThrow('Invalid input.');
        expect(() => smoother.initialize([v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)], []))
            .toThrow('Invalid input.');
    });

    it('records the mesh and counts vertex neighbors as twice the incidence', () => {
        // Two triangles sharing the edge (1,2).
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0)];
        const indices = [0, 1, 2, 1, 3, 2];
        const smoother = new MeshSmoother();
        smoother.initialize(vertices, indices);

        expect(smoother.getNumVertices()).toBe(4);
        expect(smoother.getNumTriangles()).toBe(2);
        expect(smoother.getVertices()).toBe(vertices);
        expect(smoother.getIndices()).toBe(indices);
        // Vertex 0 and vertex 3 are in one triangle; 1 and 2 are in both.
        expect(smoother.getNeighborCounts()).toEqual([2, 4, 4, 2]);
    });
});

describe('MeshSmoother.update known values', () => {
    it('matches a hand computation for a single triangle', () => {
        // Triangle in the plane z = 0, counterclockwise from +z.
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const indices = [0, 1, 2];
        const smoother = new MeshSmoother();
        smoother.initialize(vertices, indices);
        smoother.update();

        // Every vertex has neighbor count 2, so mean[i] is the midpoint of
        // the other two vertices. The unit normal is (0,0,1) everywhere, and
        // diff lies in the plane, so tangent = diff and the vertex moves
        // halfway to the midpoint of the opposite edge (the centroid-ward
        // half-step of a Laplacian smoother).
        expect(smoother.getNeighborCounts()).toEqual([2, 2, 2]);
        expect(smoother.getNormals()[0].values).toEqual([0, 0, 1]);
        expect(smoother.getMeans()[0].values).toEqual([0.5, 0.5, 0]);
        expect(smoother.getMeans()[1].values).toEqual([0, 0.5, 0]);
        expect(smoother.getMeans()[2].values).toEqual([0.5, 0, 0]);

        expect(vertices[0].values).toEqual([0.25, 0.25, 0]);
        expect(vertices[1].values).toEqual([0.5, 0.25, 0]);
        expect(vertices[2].values).toEqual([0.25, 0.5, 0]);
    });

    it('shrinks a triangle toward its centroid while preserving the centroid', () => {
        const vertices = [v3(0, 0, 0), v3(4, 0, 0), v3(0, 4, 0)];
        const indices = [0, 1, 2];
        const smoother = new MeshSmoother();
        smoother.initialize(vertices, indices);

        const centroidBefore = [4 / 3, 4 / 3, 0];
        for (let iter = 0; iter < 20; ++iter) {
            smoother.update();
        }
        const centroidAfter = [0, 1, 2].map(k =>
            (vertices[0].values[k] + vertices[1].values[k] + vertices[2].values[k]) / 3);
        for (let k = 0; k < 3; ++k) {
            expect(centroidAfter[k]).toBeCloseTo(centroidBefore[k], 12);
        }
        // The vertices collapse toward the centroid: v -> v + 0.5*(mid - v)
        // and mid - v = -1.5*(v - centroid), so v - centroid shrinks by 1/4
        // per iteration.
        for (const v of vertices) {
            const d = Math.hypot(v.values[0] - centroidBefore[0],
                v.values[1] - centroidBefore[1], v.values[2] - centroidBefore[2]);
            expect(d).toBeLessThan(1e-10);
        }
    });

    it('applies the normal weight along the unit normal', () => {
        class Inflater extends MeshSmoother {
            protected override getTangentWeight(): number { return 0; }
            protected override getNormalWeight(): number { return 0.25; }
        }
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const indices = [0, 1, 2];
        const smoother = new Inflater();
        smoother.initialize(vertices, indices);
        smoother.update();
        // The unit normal is (0,0,1), so every vertex rises by 0.25 and the
        // x and y components are unchanged.
        expect(vertices[0].values).toEqual([0, 0, 0.25]);
        expect(vertices[1].values).toEqual([1, 0, 0.25]);
        expect(vertices[2].values).toEqual([0, 1, 0.25]);
    });
});

describe('MeshSmoother.update cross-check', () => {
    it('agrees with an independent implementation over many iterations', () => {
        let seed = 24680;
        const next = () => {
            seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF;
            return seed / 0x7FFFFFFF;
        };
        const n = 5;
        const grid = makeGrid(n, () => 0);
        // Jitter all three coordinates so the normals are genuinely varying.
        for (const v of grid.vertices) {
            v.values[0] += 0.3 * (next() - 0.5);
            v.values[1] += 0.3 * (next() - 0.5);
            v.values[2] += 0.3 * (next() - 0.5);
        }
        const reference = copyVerts(grid.vertices);

        const smoother = new MeshSmoother();
        smoother.initialize(grid.vertices, grid.indices);
        for (let iter = 0; iter < 10; ++iter) {
            smoother.update();
            referenceUpdate(reference, grid.indices, 0.5, 0.0);
            for (let i = 0; i < grid.vertices.length; ++i) {
                for (let k = 0; k < 3; ++k) {
                    expect(grid.vertices[i].values[k]).toBeCloseTo(reference[i][k], 10);
                }
            }
        }
    });

    it('agrees with the reference for nonzero tangent and normal weights', () => {
        class Weighted extends MeshSmoother {
            protected override getTangentWeight(): number { return 0.3; }
            protected override getNormalWeight(): number { return 0.05; }
        }
        const grid = makeGrid(4, (i, j) => 0.2 * Math.sin(i + j));
        const reference = copyVerts(grid.vertices);
        const smoother = new Weighted();
        smoother.initialize(grid.vertices, grid.indices);
        for (let iter = 0; iter < 5; ++iter) {
            smoother.update();
            referenceUpdate(reference, grid.indices, 0.3, 0.05);
            for (let i = 0; i < grid.vertices.length; ++i) {
                for (let k = 0; k < 3; ++k) {
                    expect(grid.vertices[i].values[k]).toBeCloseTo(reference[i][k], 10);
                }
            }
        }
    });
});

describe('MeshSmoother planarity and convergence', () => {
    it('keeps an exactly planar mesh planar', () => {
        const grid = makeGrid(6, () => 0);
        let seed = 4242;
        const next = () => {
            seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF;
            return seed / 0x7FFFFFFF;
        };
        for (const v of grid.vertices) {
            v.values[0] += 0.25 * (next() - 0.5);
            v.values[1] += 0.25 * (next() - 0.5);
        }
        const smoother = new MeshSmoother();
        smoother.initialize(grid.vertices, grid.indices);
        for (let iter = 0; iter < 50; ++iter) {
            smoother.update();
        }
        // The normals are all (0,0,1) and the displacements are tangential,
        // so z stays exactly 0.
        for (const v of grid.vertices) {
            expect(v.values[2]).toBe(0);
        }
    });

    it('reduces the in-plane roughness of a jittered planar grid', () => {
        const grid = makeGrid(7, () => 0);
        let seed = 777;
        const next = () => {
            seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF;
            return seed / 0x7FFFFFFF;
        };
        for (const index of grid.interior) {
            grid.vertices[index].values[0] += 0.3 * (next() - 0.5);
            grid.vertices[index].values[1] += 0.3 * (next() - 0.5);
        }
        // Pin the boundary so the target lattice is fixed; the interior
        // update is then exactly the Laplacian v += (1/2)(mean - v), which
        // is a contraction toward the harmonic (regular) configuration.
        const boundarySet = new Set(grid.boundary);
        class Pinned extends MeshSmoother {
            protected override vertexInfluenced(i: number): boolean {
                return !boundarySet.has(i);
            }
        }
        const smoother = new Pinned();
        smoother.initialize(grid.vertices, grid.indices);

        const roughness = () => {
            let sum = 0;
            const means = smoother.getMeans();
            for (const index of grid.interior) {
                const v = grid.vertices[index].values;
                const m = means[index].values;
                sum += (m[0] - v[0]) ** 2 + (m[1] - v[1]) ** 2 + (m[2] - v[2]) ** 2;
            }
            return sum;
        };

        smoother.update();
        let previous = roughness();
        for (let iter = 0; iter < 30; ++iter) {
            smoother.update();
            const current = roughness();
            expect(current).toBeLessThanOrEqual(previous + 1e-12);
            previous = current;
        }
        // The interior settles onto the regular lattice positions.
        expect(previous).toBeLessThan(1e-6);
    });

    it('is tangential: it converges to a fixed point that keeps the height field', () => {
        const grid = makeGrid(9, () => 0);
        let seed = 31337;
        const next = () => {
            seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF;
            return seed / 0x7FFFFFFF;
        };
        for (const index of grid.interior) {
            grid.vertices[index].values[2] = 0.3 * (next() - 0.5);
        }
        const initialMax = Math.max(...grid.interior.map(i =>
            Math.abs(grid.vertices[i].values[2])));

        // Pin the boundary so the smoothing is driven by the interior noise
        // alone; this is the standard way to keep a patch's border fixed.
        const boundarySet = new Set(grid.boundary);
        class Pinned extends MeshSmoother {
            protected override vertexInfluenced(i: number): boolean {
                return !boundarySet.has(i);
            }
        }
        const smoother = new Pinned();
        smoother.initialize(grid.vertices, grid.indices);
        for (let iter = 0; iter < 400; ++iter) {
            smoother.update();
        }
        // The iteration has reached a fixed point: one more step moves
        // nothing measurable.
        const before = copyVerts(grid.vertices);
        smoother.update();
        let maxStep = 0;
        for (let i = 0; i < grid.vertices.length; ++i) {
            maxStep = Math.max(maxStep, Math.hypot(
                grid.vertices[i].values[0] - before[i][0],
                grid.vertices[i].values[1] - before[i][1],
                grid.vertices[i].values[2] - before[i][2]));
        }
        expect(maxStep).toBeLessThan(1e-9);

        // The displacement is the component of (mean - v) tangent to the
        // estimated surface normal, so the smoother relaxes the in-surface
        // vertex distribution and leaves the height field essentially
        // unchanged. This is the upstream behavior with the default weights
        // (tangent 1/2, normal 0); to move the surface itself, override
        // getNormalWeight.
        const finalMax = Math.max(...grid.interior.map(i =>
            Math.abs(grid.vertices[i].values[2])));
        expect(finalMax).toBeLessThanOrEqual(initialMax);
        expect(finalMax).toBeGreaterThan(0.9 * initialMax);

        // The boundary is untouched.
        for (const index of grid.boundary) {
            expect(grid.vertices[index].values[2]).toBe(0);
        }
    });
});

describe('MeshSmoother vertex influence', () => {
    it('leaves uninfluenced vertices exactly where they were', () => {
        class PinFirst extends MeshSmoother {
            protected override vertexInfluenced(i: number, t: number): boolean {
                // The time argument is forwarded by update().
                return i !== 0 || t > 5;
            }
        }
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const indices = [0, 1, 2];
        const smoother = new PinFirst();
        smoother.initialize(vertices, indices);
        smoother.update(0);
        expect(vertices[0].values).toEqual([0, 0, 0]);
        expect(vertices[1].values).not.toEqual([1, 0, 0]);

        // With t > 5 the first vertex moves too.
        smoother.update(10);
        expect(vertices[0].values).not.toEqual([0, 0, 0]);
    });
});

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// A triangulated n-by-n grid with the vertices jittered in all three
// coordinates. wellScaledVector keeps the jitter away from subnormals, which
// matters because update() normalizes a cross product of vertex differences.
const jitteredGrid = fc.tuple(fc.integer({ min: 3, max: 5 }),
    fc.array(wellScaledVector(3, -0.4, 0.4), { minLength: 25, maxLength: 25 }))
    .map(([n, jitter]) => {
        const grid = makeGrid(n, () => 0);
        for (let i = 0; i < grid.vertices.length; ++i) {
            const j = jitter[i % jitter.length];
            for (let k = 0; k < 3; ++k) {
                grid.vertices[i].values[k] += j.values[k];
            }
        }
        return grid;
    });

// An independent accumulation of the un-normalized means, in the order in
// which update() visits the triangles.
function referenceMeanSums(vertices: Vector[], indices: number[]): number[][] {
    const sums: number[][] = vertices.map(() => [0, 0, 0]);
    for (let t = 0; 3 * t < indices.length; ++t) {
        const a = indices[3 * t], b = indices[3 * t + 1], c = indices[3 * t + 2];
        const A = vertices[a].values, B = vertices[b].values, C = vertices[c].values;
        for (let k = 0; k < 3; ++k) {
            sums[a][k] += B[k] + C[k];
            sums[b][k] += C[k] + A[k];
            sums[c][k] += A[k] + B[k];
        }
    }
    return sums;
}

describe('MeshSmoother verification', () => {
    it('counts each vertex twice per incident triangle', () => {
        check(jitteredGrid, grid => {
            const smoother = new MeshSmoother();
            smoother.initialize(grid.vertices, grid.indices);
            const expected = new Array<number>(grid.vertices.length).fill(0);
            for (const index of grid.indices) { expected[index] += 2; }
            expect(smoother.getNeighborCounts()).toEqual(expected);
        });
    });

    it('divides the mean by the reciprocal of the neighbor count, as Vector operator/= does', () => {
        // Upstream 'mMeans[i] /= T(count)' multiplies by 1/count; a plain
        // componentwise division differs in the last ulp. The reference
        // accumulates in the same order, so the comparison is exact.
        check(jitteredGrid, grid => {
            const before = grid.vertices.map(v => v.clone());
            const sums = referenceMeanSums(before, grid.indices);
            const counts = new Array<number>(before.length).fill(0);
            for (const index of grid.indices) { counts[index] += 2; }

            const smoother = new MeshSmoother();
            smoother.initialize(grid.vertices, grid.indices);
            smoother.update();
            const means = smoother.getMeans();
            for (let i = 0; i < before.length; ++i) {
                const inv = 1 / counts[i];
                for (let k = 0; k < 3; ++k) {
                    expect(means[i].values[k] + 0).toBe(sums[i][k] * inv + 0);
                }
            }
        });
    });

    it('gives a vertex with no incident triangle a zero mean, not NaN', () => {
        // Regression for the port defect fixed in V16. Upstream's
        // 'mMeans[i] /= T(0)' is Vector operator/= with a zero divisor, which
        // sets the vector to ZERO. A componentwise 0/0 gives NaN, which then
        // propagates into the vertex position.
        check(fc.tuple(wellScaledVector(3, -5, 5), fc.integer({ min: 1, max: 4 })),
            ([orphanPosition, numOrphans]) => {
                const vertices = [
                    Vector.fromArray([0, 0, 0]),
                    Vector.fromArray([1, 0, 0]),
                    Vector.fromArray([0, 1, 0])
                ];
                for (let i = 0; i < numOrphans; ++i) {
                    vertices.push(orphanPosition.clone());
                }
                const smoother = new MeshSmoother();
                smoother.initialize(vertices, [0, 1, 2]);
                smoother.update();

                for (let i = 3; i < vertices.length; ++i) {
                    expect(smoother.getNeighborCounts()[i]).toBe(0);
                    expect(smoother.getMeans()[i].values).toEqual([0, 0, 0]);
                    expect(smoother.getNormals()[i].values).toEqual([0, 0, 0]);
                    // mean - v = -v, the normal is zero so the whole
                    // displacement is tangential: v -> v + 0.5*(-v).
                    for (let k = 0; k < 3; ++k) {
                        expect(Number.isNaN(vertices[i].values[k])).toBe(false);
                        expectClose(vertices[i].values[k],
                            0.5 * orphanPosition.values[k], 1e-15, 1e-15);
                    }
                }
            });
    });

    it('mutates the caller vertices in place and leaves the topology alone', () => {
        check(jitteredGrid, grid => {
            const objects = grid.vertices.slice();
            const indices = grid.indices.slice();
            const smoother = new MeshSmoother();
            smoother.initialize(grid.vertices, grid.indices);
            smoother.update();
            smoother.update();
            expect(smoother.getVertices()).toBe(grid.vertices);
            for (let i = 0; i < objects.length; ++i) {
                expect(grid.vertices[i]).toBe(objects[i]);
                expect(grid.vertices[i].size).toBe(3);
            }
            expect(grid.indices).toEqual(indices);
            expect(smoother.getNumTriangles()).toBe(indices.length / 3);
            expect(smoother.getNumVertices()).toBe(objects.length);
        });
    });

    it('moves every vertex tangentially when the normal weight is zero', () => {
        // The default weights are 1/2 for the tangent and 0 for the normal, so
        // the displacement is orthogonal to the (unit) vertex normal and is
        // exactly half the tangential part of mean - position.
        check(jitteredGrid, grid => {
            const before = grid.vertices.map(v => v.clone());
            const smoother = new MeshSmoother();
            smoother.initialize(grid.vertices, grid.indices);
            smoother.update();
            const normals = smoother.getNormals();
            const means = smoother.getMeans();
            for (let i = 0; i < before.length; ++i) {
                const delta = sub(grid.vertices[i], before[i]);
                expectClose(dot(delta, normals[i]), 0, 1e-12, 1e-12);
                const diff = sub(means[i], before[i]);
                const tangent = sub(diff,
                    mulVector(normals[i], dot(diff, normals[i])));
                expectVectorClose(delta, mulVector(tangent, 0.5), 1e-12, 1e-12);
            }
        });
    });

    it('does not move anything when both weights are zero', () => {
        class Frozen extends MeshSmoother {
            protected override getTangentWeight(): number { return 0; }
            protected override getNormalWeight(): number { return 0; }
        }
        check(jitteredGrid, grid => {
            const before = grid.vertices.map(v => v.values.slice());
            const smoother = new Frozen();
            smoother.initialize(grid.vertices, grid.indices);
            smoother.update();
            for (let i = 0; i < before.length; ++i) {
                expect(grid.vertices[i].values).toEqual(before[i]);
            }
        });
    });

    it('commutes with a rigid motion of the mesh', () => {
        check(fc.tuple(jitteredGrid, rotationFrame(3),
            wellScaledVector(3, -3, 3)), ([grid, R, t]) => {
            const rot = (x: Vector): Vector => {
                const r = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    r.values[i] = R[0].values[i] * x.values[0]
                        + R[1].values[i] * x.values[1]
                        + R[2].values[i] * x.values[2];
                }
                return r;
            };
            const moved = grid.vertices.map(v => add(rot(v), t));

            const a = new MeshSmoother();
            a.initialize(grid.vertices, grid.indices);
            a.update();

            const b = new MeshSmoother();
            b.initialize(moved, grid.indices);
            b.update();

            for (let i = 0; i < moved.length; ++i) {
                // The rotation of a smoothed vertex and the smoothing of a
                // rotated vertex differ only by the round-off of the rotation
                // itself, which is a few ulps of the coordinate magnitudes.
                expectVectorClose(moved[i], add(rot(grid.vertices[i]), t),
                    1e-11, 1e-11);
            }
        });
    });

    it('commutes with a positive uniform scale when the normal weight is zero', () => {
        check(fc.tuple(jitteredGrid, positive(8, 0.05)), ([grid, s]) => {
            const scaled = grid.vertices.map(v => mulVector(v, s));

            const a = new MeshSmoother();
            a.initialize(grid.vertices, grid.indices);
            a.update();

            const b = new MeshSmoother();
            b.initialize(scaled, grid.indices);
            b.update();

            for (let i = 0; i < scaled.length; ++i) {
                expectVectorClose(scaled[i], mulVector(grid.vertices[i], s),
                    1e-11, 1e-11);
            }
        });
    });

    it('normalizes the accumulated triangle normals, or zeroes them', () => {
        check(jitteredGrid, grid => {
            const smoother = new MeshSmoother();
            smoother.initialize(grid.vertices, grid.indices);
            smoother.update();
            for (const n of smoother.getNormals()) {
                const len = length(n);
                expect(len === 0 || Math.abs(len - 1) <= 1e-12).toBe(true);
            }
        });
    });

    it('rejects meshes with fewer than three vertices or no triangle', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 2 }),
            fc.integer({ min: 0, max: 6 })), ([numVertices, numIndices]) => {
            const vertices: Vector[] = [];
            for (let i = 0; i < numVertices; ++i) { vertices.push(new Vector(3)); }
            const indices = new Array<number>(numIndices).fill(0);
            expect(() => new MeshSmoother().initialize(vertices, indices))
                .toThrow('Invalid input.');
        });
        check(fc.integer({ min: 3, max: 10 }), numVertices => {
            const vertices: Vector[] = [];
            for (let i = 0; i < numVertices; ++i) { vertices.push(new Vector(3)); }
            for (const numIndices of [0, 1, 2]) {
                expect(() => new MeshSmoother()
                    .initialize(vertices, new Array<number>(numIndices).fill(0)))
                    .toThrow('Invalid input.');
            }
        });
    });
});
