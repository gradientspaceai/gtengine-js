import { describe, it, expect } from 'vitest';
import { MeshCurvature } from '../src/MeshCurvature.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3, cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    seededRandom, wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function V3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

interface TriMesh {
    vertices: Vector[];
    indices: number[];
}

// An icosphere: an icosahedron subdivided 'level' times and projected onto
// the sphere of radius r. Every vertex is unique, so the mesh is manifold
// and closed.
function icosphere(r: number, level: number): TriMesh {
    const t = (1 + Math.sqrt(5)) / 2;
    let vertices: Vector[] = [
        V3(-1, t, 0), V3(1, t, 0), V3(-1, -t, 0), V3(1, -t, 0),
        V3(0, -1, t), V3(0, 1, t), V3(0, -1, -t), V3(0, 1, -t),
        V3(t, 0, -1), V3(t, 0, 1), V3(-t, 0, -1), V3(-t, 0, 1)
    ];
    let indices: number[] = [
        0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11,
        1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
        3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9,
        4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1
    ];

    for (let s = 0; s < level; ++s) {
        const midpoints = new Map<string, number>();
        const nextIndices: number[] = [];
        const midpoint = (a: number, b: number): number => {
            const key = (a < b ? a + ',' + b : b + ',' + a);
            const existing = midpoints.get(key);
            if (existing !== undefined) {
                return existing;
            }
            const m = mul(add(vertices[a], vertices[b]), 0.5);
            vertices.push(m);
            const index = vertices.length - 1;
            midpoints.set(key, index);
            return index;
        };
        for (let i = 0; i < indices.length; i += 3) {
            const v0 = indices[i], v1 = indices[i + 1], v2 = indices[i + 2];
            const a = midpoint(v0, v1);
            const b = midpoint(v1, v2);
            const c = midpoint(v2, v0);
            nextIndices.push(v0, a, c, v1, b, a, v2, c, b, a, b, c);
        }
        indices = nextIndices;
    }

    vertices = vertices.map(v => {
        const u = v.clone();
        normalize(u);
        return mul(u, r);
    });
    return { vertices, indices };
}

// An open cylinder of radius r and height h with numAxial+1 rings of
// numRadial vertices. The vertices are unique; the mesh wraps around in the
// angular direction.
function cylinder(r: number, h: number, numRadial: number,
    numAxial: number): TriMesh {
    const vertices: Vector[] = [];
    for (let j = 0; j <= numAxial; ++j) {
        const z = h * (j / numAxial - 0.5);
        for (let i = 0; i < numRadial; ++i) {
            const angle = 2 * Math.PI * i / numRadial;
            vertices.push(V3(r * Math.cos(angle), r * Math.sin(angle), z));
        }
    }

    const indices: number[] = [];
    for (let j = 0; j < numAxial; ++j) {
        for (let i = 0; i < numRadial; ++i) {
            const i1 = (i + 1) % numRadial;
            const a = j * numRadial + i;
            const b = j * numRadial + i1;
            const c = (j + 1) * numRadial + i;
            const d = (j + 1) * numRadial + i1;
            // Counterclockwise seen from outside.
            indices.push(a, b, c, b, d, c);
        }
    }
    return { vertices, indices };
}

// A planar grid in the plane z = 0 with normals along +z.
function plane(size: number, n: number): TriMesh {
    const vertices: Vector[] = [];
    for (let j = 0; j <= n; ++j) {
        for (let i = 0; i <= n; ++i) {
            vertices.push(V3(size * (i / n - 0.5), size * (j / n - 0.5), 0));
        }
    }
    const indices: number[] = [];
    for (let j = 0; j < n; ++j) {
        for (let i = 0; i < n; ++i) {
            const a = j * (n + 1) + i;
            const b = a + 1;
            const c = a + (n + 1);
            const d = c + 1;
            indices.push(a, b, c, b, d, c);
        }
    }
    return { vertices, indices };
}

describe('MeshCurvature', () => {
    it('computes area-weighted unit normals', () => {
        const r = 2;
        const { vertices, indices } = icosphere(r, 2);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 0);

        const normals = mc.getNormals();
        expect(normals.length).toBe(vertices.length);
        for (let i = 0; i < vertices.length; ++i) {
            expect(length(normals[i])).toBeCloseTo(1, 12);
            // Outward on the sphere.
            expect(dot(normals[i], vertices[i]) / r).toBeGreaterThan(0.999);
        }
    });

    it('estimates the curvatures of a sphere of radius r as 1/r', () => {
        for (const r of [0.5, 1, 3]) {
            const { vertices, indices } = icosphere(r, 3);
            const mc = new MeshCurvature();
            mc.compute(vertices, indices, 0);

            const kMin = mc.getMinCurvatures();
            const kMax = mc.getMaxCurvatures();
            expect(kMin.length).toBe(vertices.length);
            for (let i = 0; i < vertices.length; ++i) {
                expect(kMin[i]).toBeGreaterThan(0);
                expect(kMin[i]).toBeLessThanOrEqual(kMax[i]);
                expect(Math.abs(kMin[i] * r - 1)).toBeLessThan(0.05);
                expect(Math.abs(kMax[i] * r - 1)).toBeLessThan(0.05);
            }
        }
    });

    it('produces principal directions orthogonal to the normal', () => {
        const { vertices, indices } = icosphere(1.5, 2);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 0);

        const normals = mc.getNormals();
        const dMin = mc.getMinDirections();
        const dMax = mc.getMaxDirections();
        for (let i = 0; i < vertices.length; ++i) {
            expect(length(dMin[i])).toBeCloseTo(1, 10);
            expect(length(dMax[i])).toBeCloseTo(1, 10);
            expect(dot(dMin[i], normals[i])).toBeCloseTo(0, 10);
            expect(dot(dMax[i], normals[i])).toBeCloseTo(0, 10);
        }
    });

    it('estimates the curvatures of a cylinder as (0, 1/r)', () => {
        const r = 2, h = 4, numRadial = 48, numAxial = 12;
        const { vertices, indices } = cylinder(r, h, numRadial, numAxial);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 0);

        const normals = mc.getNormals();
        const kMin = mc.getMinCurvatures();
        const kMax = mc.getMaxCurvatures();
        const dMin = mc.getMinDirections();
        const dMax = mc.getMaxDirections();

        // Interior rings only; the boundary rings have one-sided fans.
        for (let j = 1; j < numAxial; ++j) {
            for (let i = 0; i < numRadial; ++i) {
                const v = j * numRadial + i;

                // The normal is radial and outward.
                const radial = V3(vertices[v].values[0], vertices[v].values[1], 0);
                expect(dot(normals[v], radial) / r).toBeGreaterThan(0.999);

                expect(Math.abs(kMin[v])).toBeLessThan(0.01);
                expect(Math.abs(kMax[v] * r - 1)).toBeLessThan(0.02);

                // The direction of zero curvature is the cylinder axis and
                // the direction of curvature 1/r is circumferential.
                expect(Math.abs(dMin[v].values[2])).toBeGreaterThan(0.99);
                expect(Math.abs(dMax[v].values[2])).toBeLessThan(0.06);
                expect(Math.abs(dot(dMin[v], dMax[v]))).toBeLessThan(1e-8);
            }
        }
    });

    it('estimates zero curvature on a plane', () => {
        const n = 8;
        const { vertices, indices } = plane(4, n);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 0);

        const normals = mc.getNormals();
        const kMin = mc.getMinCurvatures();
        const kMax = mc.getMaxCurvatures();
        for (let i = 0; i < vertices.length; ++i) {
            expect(normals[i].values[0]).toBeCloseTo(0, 12);
            expect(normals[i].values[1]).toBeCloseTo(0, 12);
            expect(Math.abs(normals[i].values[2])).toBeCloseTo(1, 12);
            expect(kMin[i]).toBeCloseTo(0, 10);
            expect(kMax[i]).toBeCloseTo(0, 10);
        }
    });

    it('flags locally planar points via the singularity threshold', () => {
        // On the plane every D*W^T is exactly zero, so any positive
        // threshold flags every vertex; the directions are then the
        // orthogonal complement of the normal.
        const { vertices, indices } = plane(4, 4);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 1e-8);

        const normals = mc.getNormals();
        const dMin = mc.getMinDirections();
        const dMax = mc.getMaxDirections();
        for (let i = 0; i < vertices.length; ++i) {
            expect(mc.getMinCurvatures()[i]).toBe(0);
            expect(mc.getMaxCurvatures()[i]).toBe(0);
            expect(length(dMin[i])).toBeCloseTo(1, 12);
            expect(length(dMax[i])).toBeCloseTo(1, 12);
            expect(dot(dMin[i], normals[i])).toBeCloseTo(0, 12);
            expect(dot(dMax[i], normals[i])).toBeCloseTo(0, 12);
            // {dMin, dMax, N} is right-handed.
            const c = cross(dMin[i], dMax[i]);
            expect(length(sub(c, normals[i]))).toBeCloseTo(0, 12);
        }
    });

    it('flags every vertex of a sphere when the threshold is huge', () => {
        const { vertices, indices } = icosphere(1, 1);
        const mc = new MeshCurvature();
        mc.compute(vertices, indices, 1e6);
        for (let i = 0; i < vertices.length; ++i) {
            expect(mc.getMinCurvatures()[i]).toBe(0);
            expect(mc.getMaxCurvatures()[i]).toBe(0);
        }
    });

    it('reverses the sign of the curvatures when the mesh is inverted', () => {
        // Flipping the triangle winding flips the estimated normals, which
        // flips the sign of the shape operator.
        const r = 1.5;
        const { vertices, indices } = icosphere(r, 2);
        const flipped: number[] = [];
        for (let i = 0; i < indices.length; i += 3) {
            flipped.push(indices[i], indices[i + 2], indices[i + 1]);
        }

        const outward = new MeshCurvature();
        outward.compute(vertices, indices, 0);
        const inward = new MeshCurvature();
        inward.compute(vertices, flipped, 0);

        for (let i = 0; i < vertices.length; ++i) {
            expect(inward.getMinCurvatures()[i]).toBeCloseTo(
                -outward.getMaxCurvatures()[i], 7);
            expect(inward.getMaxCurvatures()[i]).toBeCloseTo(
                -outward.getMinCurvatures()[i], 7);
            expect(length(add(inward.getNormals()[i],
                outward.getNormals()[i]))).toBeCloseTo(0, 10);
        }
    });

    it('is invariant under rigid motions of the mesh', () => {
        const r = 1.25;
        const { vertices, indices } = icosphere(r, 2);

        // Rotate about the axis (1,1,1)/sqrt(3) by pi/3 and translate.
        const axis = V3(1, 1, 1);
        normalize(axis);
        const angle = Math.PI / 3;
        const cs = Math.cos(angle), sn = Math.sin(angle);
        const translation = V3(3, -2, 7);
        const moved = vertices.map(v => {
            const term0 = mul(v, cs);
            const term1 = mul(cross(axis, v), sn);
            const term2 = mul(axis, dot(axis, v) * (1 - cs));
            return add(add(add(term0, term1), term2), translation);
        });

        const a = new MeshCurvature();
        a.compute(vertices, indices, 0);
        const b = new MeshCurvature();
        b.compute(moved, indices, 0);
        for (let i = 0; i < vertices.length; ++i) {
            expect(b.getMinCurvatures()[i]).toBeCloseTo(a.getMinCurvatures()[i], 6);
            expect(b.getMaxCurvatures()[i]).toBeCloseTo(a.getMaxCurvatures()[i], 6);
        }
    });

    it('handles an empty mesh and a single triangle', () => {
        const mc = new MeshCurvature();
        mc.compute([], [], 0);
        expect(mc.getNormals().length).toBe(0);
        expect(mc.getMinCurvatures().length).toBe(0);

        const single = new MeshCurvature();
        single.compute([V3(0, 0, 0), V3(1, 0, 0), V3(0, 1, 0)], [0, 1, 2], 1e-12);
        // A single flat triangle has zero normal derivatives everywhere.
        for (let i = 0; i < 3; ++i) {
            expect(single.getNormals()[i].values[2]).toBeCloseTo(1, 12);
            expect(single.getMinCurvatures()[i]).toBe(0);
            expect(single.getMaxCurvatures()[i]).toBe(0);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// A closed icosphere with a random radius, optionally jittered so the mesh is
// a generic (non-symmetric) surface. wellScaled keeps the jitter away from
// subnormals, which matters because the algorithm squares vertex differences.
const roundMesh = fc.tuple(positive(4, 0.25),
    fc.array(wellScaled(-0.05, 0.05), { minLength: 3 * 42, maxLength: 3 * 42 }))
    .map(([r, jitter]) => {
        const mesh = icosphere(r, 1);
        for (let i = 0; i < mesh.vertices.length; ++i) {
            for (let k = 0; k < 3; ++k) {
                mesh.vertices[i].values[k] += r * jitter[3 * i + k];
            }
        }
        return mesh;
    });

// Rotate by the frame whose columns are the frame vectors.
function rotateBy(R: Vector[], x: Vector): Vector {
    const r = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        r.values[i] = R[0].values[i] * x.values[0]
            + R[1].values[i] * x.values[1]
            + R[2].values[i] * x.values[2];
    }
    return r;
}

function computeOn(mesh: TriMesh, threshold: number): MeshCurvature {
    const query = new MeshCurvature();
    query.compute(mesh.vertices, mesh.indices, threshold);
    return query;
}

describe('MeshCurvature verification', () => {
    it('produces unit normals, ordered curvatures and tangential directions', () => {
        // The directions are unit length except at an exactly umbilic vertex,
        // where upstream builds both candidate eigenvectors from S(0,1) and
        // curvature - S(0,0), both of which are zero; see the dedicated test
        // below. Whatever their length, they lie in the tangent plane.
        check(roundMesh, mesh => {
            const query = computeOn(mesh, 0);
            const normals = query.getNormals();
            const kmin = query.getMinCurvatures();
            const kmax = query.getMaxCurvatures();
            const dmin = query.getMinDirections();
            const dmax = query.getMaxDirections();
            expect(normals.length).toBe(mesh.vertices.length);
            for (let i = 0; i < mesh.vertices.length; ++i) {
                expectClose(length(normals[i]), 1, 1e-12, 1e-12);
                expect(kmin[i]).toBeLessThanOrEqual(kmax[i]);
                for (const d of [dmin[i], dmax[i]]) {
                    const len = length(d);
                    expect(len === 0 || Math.abs(len - 1) <= 1e-9).toBe(true);
                    expectClose(dot(d, normals[i]), 0, 1e-9, 1e-9);
                }
            }
        });
    });

    it('returns zero principal directions at an exactly umbilic vertex (upstream quirk)', () => {
        // On the unit icosahedron one vertex has a shape matrix that is
        // exactly a multiple of the identity. Upstream then forms
        //   W0 = (S(0,1), k - S(0,0)) and W1 = (k - S(1,1), S(1,0)),
        // both of which are the zero vector, takes the first branch because
        // Dot(W0,W0) >= Dot(W1,W1) holds for 0 >= 0, and Normalize leaves the
        // zero vector alone. The direction returned is therefore (0,0,0)
        // rather than an arbitrary unit tangent. The port preserves this.
        const mesh = icosphere(1, 0);
        const query = computeOn(mesh, 0);
        const degenerate: number[] = [];
        for (let i = 0; i < mesh.vertices.length; ++i) {
            if (length(query.getMinDirections()[i]) === 0) { degenerate.push(i); }
        }
        expect(degenerate.length).toBeGreaterThan(0);
        for (const i of degenerate) {
            expect(query.getMinDirections()[i].values).toEqual([0, 0, 0]);
            expect(query.getMaxDirections()[i].values).toEqual([0, 0, 0]);
            // The curvatures themselves are still the correct 1/r = 1.
            expectClose(query.getMinCurvatures()[i], 1, 0, 0.1);
            expectClose(query.getMaxCurvatures()[i], 1, 0, 0.1);
            expectClose(query.getMinCurvatures()[i],
                query.getMaxCurvatures()[i], 0, 0);
        }
    });

    it('computes the normals as the normalized area-weighted triangle sum', () => {
        check(roundMesh, mesh => {
            const expected = mesh.vertices.map(() => new Vector(3));
            for (let t = 0; 3 * t < mesh.indices.length; ++t) {
                const a = mesh.indices[3 * t];
                const b = mesh.indices[3 * t + 1];
                const c = mesh.indices[3 * t + 2];
                const n = cross(sub(mesh.vertices[b], mesh.vertices[a]),
                    sub(mesh.vertices[c], mesh.vertices[a]));
                for (const v of [a, b, c]) {
                    for (let k = 0; k < 3; ++k) {
                        expected[v].values[k] += n.values[k];
                    }
                }
            }
            const query = computeOn(mesh, 0);
            const normals = query.getNormals();
            for (let i = 0; i < expected.length; ++i) {
                const len = length(expected[i]);
                expect(len).toBeGreaterThan(0);
                for (let k = 0; k < 3; ++k) {
                    expectClose(normals[i].values[k],
                        expected[i].values[k] / len, 1e-12, 1e-12);
                }
            }
        });
    });

    it('leaves the caller vertex array untouched and is repeatable', () => {
        check(roundMesh, mesh => {
            const before = mesh.vertices.map(v => v.values.slice());
            const first = computeOn(mesh, 0);
            const firstMin = first.getMinCurvatures().slice();
            for (let i = 0; i < before.length; ++i) {
                expect(mesh.vertices[i].values).toEqual(before[i]);
            }
            const second = computeOn(mesh, 0);
            expect(second.getMinCurvatures()).toEqual(firstMin);
        });
    });

    it('estimates 1/r at every vertex of a sphere of radius r', () => {
        check(positive(6, 0.25), r => {
            const mesh = icosphere(r, 2);
            const query = computeOn(mesh, 0);
            for (let i = 0; i < mesh.vertices.length; ++i) {
                // The icosphere is only an approximation of the sphere, so
                // the discrete estimate carries several percent of error.
                expectClose(query.getMinCurvatures()[i], 1 / r, 0, 0.08);
                expectClose(query.getMaxCurvatures()[i], 1 / r, 0, 0.08);
            }
        });
    }, 30000);

    it('estimates (0, 1/r) away from the rim of a cylinder of radius r', () => {
        check(positive(5, 0.5), r => {
            const numRadial = 24;
            const numAxial = 8;
            const mesh = cylinder(r, 4 * r, numRadial, numAxial);
            const query = computeOn(mesh, 0);
            // The first and last rings are on the boundary, where the
            // one-sided neighborhood biases the estimate.
            for (let j = 1; j < numAxial; ++j) {
                for (let i = 0; i < numRadial; ++i) {
                    const v = j * numRadial + i;
                    // The zero curvature is only approached: the polygonal
                    // cross-section makes the axial direction slightly
                    // non-straight in the discrete estimate, so the tolerance
                    // is relative to the other principal curvature 1/r.
                    expectClose(query.getMinCurvatures()[v], 0, 0.01 / r, 0);
                    expectClose(query.getMaxCurvatures()[v], 1 / r, 0, 0.02);
                }
            }
        });
    }, 30000);

    it('is invariant under a rigid motion of the mesh', () => {
        check(fc.tuple(roundMesh, rotationFrame(3), wellScaledVector(3, -3, 3)),
            ([mesh, R, t]) => {
                const a = computeOn(mesh, 0);
                const moved: TriMesh = {
                    vertices: mesh.vertices.map(v => add(rotateBy(R, v), t)),
                    indices: mesh.indices.slice()
                };
                const b = computeOn(moved, 0);
                for (let i = 0; i < mesh.vertices.length; ++i) {
                    // The rotation re-rounds every coordinate, so the
                    // accumulated W*W^T and D*W^T sums and the 3x3 inverse
                    // differ by a few ulps of the curvature scale.
                    expectClose(a.getMinCurvatures()[i], b.getMinCurvatures()[i],
                        1e-9, 1e-6);
                    expectClose(a.getMaxCurvatures()[i], b.getMaxCurvatures()[i],
                        1e-9, 1e-6);
                    // The normals rotate with the mesh.
                        expectVectorClose(b.getNormals()[i],
                        rotateBy(R, a.getNormals()[i]), 1e-12, 1e-9);
                }
            });
    });

    it('scales the curvatures by 1/s when the mesh is scaled by s', () => {
        check(fc.tuple(roundMesh, positive(6, 0.25)), ([mesh, s]) => {
            const a = computeOn(mesh, 0);
            const scaled: TriMesh = {
                vertices: mesh.vertices.map(v => mul(v, s)),
                indices: mesh.indices.slice()
            };
            const b = computeOn(scaled, 0);
            for (let i = 0; i < mesh.vertices.length; ++i) {
                expectClose(b.getMinCurvatures()[i],
                    a.getMinCurvatures()[i] / s, 1e-8, 1e-7);
                expectClose(b.getMaxCurvatures()[i],
                    a.getMaxCurvatures()[i] / s, 1e-8, 1e-7);
                expectVectorClose(b.getNormals()[i], a.getNormals()[i],
                    1e-9, 1e-9);
            }
        });
    });

    it('permutes the outputs when the vertices are relabelled', () => {
        // Relabelling keeps the triangle visit order, so every per-vertex
        // accumulation happens in the same order and the results must be
        // bit-identical after the permutation.
        check(fc.tuple(roundMesh, fc.integer({ min: 1, max: 1000 })),
            ([mesh, seed]) => {
                const n = mesh.vertices.length;
                const perm = Array.from({ length: n }, (_, i) => i);
                const rand = seededRandom(seed);
                for (let i = n - 1; i > 0; --i) {
                    const j = Math.floor(rand() * (i + 1));
                    [perm[i], perm[j]] = [perm[j], perm[i]];
                }
                const relabelled: TriMesh = {
                    vertices: new Array<Vector>(n),
                    indices: mesh.indices.map(v => perm[v])
                };
                for (let i = 0; i < n; ++i) {
                    relabelled.vertices[perm[i]] = mesh.vertices[i].clone();
                }

                const a = computeOn(mesh, 0);
                const b = computeOn(relabelled, 0);
                for (let i = 0; i < n; ++i) {
                    expect(b.getMinCurvatures()[perm[i]] + 0)
                        .toBe(a.getMinCurvatures()[i] + 0);
                    expect(b.getMaxCurvatures()[perm[i]] + 0)
                        .toBe(a.getMaxCurvatures()[i] + 0);
                    for (let k = 0; k < 3; ++k) {
                        expect(b.getNormals()[perm[i]].values[k] + 0)
                            .toBe(a.getNormals()[i].values[k] + 0);
                        expect(b.getMinDirections()[perm[i]].values[k] + 0)
                            .toBe(a.getMinDirections()[i].values[k] + 0);
                        expect(b.getMaxDirections()[perm[i]].values[k] + 0)
                            .toBe(a.getMaxDirections()[i].values[k] + 0);
                    }
                }
            });
    });

    it('returns the tangent frame itself at points flagged as locally planar', () => {
        // A threshold above every |D*W^T| entry flags every vertex, and the
        // upstream planar branch returns (0, 0) with directions U and V, the
        // orthogonal complement of the vertex normal.
        check(roundMesh, mesh => {
            const query = computeOn(mesh, Number.MAX_VALUE);
            for (let i = 0; i < mesh.vertices.length; ++i) {
                expect(query.getMinCurvatures()[i]).toBe(0);
                expect(query.getMaxCurvatures()[i]).toBe(0);
                const basis = [query.getNormals()[i].clone(),
                    new Vector(3), new Vector(3)];
                computeOrthogonalComplement3(1, basis);
                expectVectorClose(query.getMinDirections()[i], basis[1], 0, 0);
                expectVectorClose(query.getMaxDirections()[i], basis[2], 0, 0);
            }
        });
    });

    it('flips the sign of the curvatures when the triangles are reversed', () => {
        // Reversing every triangle reverses the normals, which reverses the
        // shape operator; the curvatures negate and swap roles.
        check(roundMesh, mesh => {
            const a = computeOn(mesh, 0);
            const flipped: TriMesh = {
                vertices: mesh.vertices.map(v => v.clone()),
                indices: []
            };
            for (let t = 0; 3 * t < mesh.indices.length; ++t) {
                flipped.indices.push(mesh.indices[3 * t],
                    mesh.indices[3 * t + 2], mesh.indices[3 * t + 1]);
            }
            const b = computeOn(flipped, 0);
            for (let i = 0; i < mesh.vertices.length; ++i) {
                expectVectorClose(b.getNormals()[i],
                    mul(a.getNormals()[i], -1), 1e-12, 1e-12);
                // The per-vertex W*W^T and D*W^T sums are accumulated in a
                // different order once the triangles are reversed, so the
                // curvatures agree only to within the conditioning of the
                // 3x3 inverse.
                expectClose(b.getMinCurvatures()[i],
                    -a.getMaxCurvatures()[i], 1e-9, 1e-6);
                expectClose(b.getMaxCurvatures()[i],
                    -a.getMinCurvatures()[i], 1e-9, 1e-6);
            }
        });
    });
});
