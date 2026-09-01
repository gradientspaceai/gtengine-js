// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SurfaceExtractorMC.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Level-surface extraction using the Marching Cubes lookup table. The class
// derives from MarchingCubes, which owns the 256-entry table of voxel sign
// configurations.
//
// Port notes: T and IndexType are both number. The nested struct Mesh is
// exported as SurfaceExtractorMCMesh. The two upstream Extract overloads are
// extractVoxel(level, perturb, F) (returning the per-voxel mesh and whether
// geometry was produced) and extract(level, perturb) (returning the vertices
// and indices for the whole image). The upstream output reference parameters
// become returned object literals. Vector3<T> is the runtime-dimension
// Vector of size 3.

import { MarchingCubes } from './MarchingCubes';
import type { MarchingCubesTopology } from './MarchingCubes';
import { Image3 } from './Image3';
import { UniqueVerticesSimplices } from './UniqueVerticesSimplices';
import { Vector } from './Vector';
import { cross } from './Vector3';

// The triangle mesh generated for a single voxel. The topology is the
// Marching Cubes table entry and the vertices are in the local coordinates
// of the voxel, where the corners are (0,0,0) through (1,1,1).
export interface SurfaceExtractorMCMesh {
    topology: MarchingCubesTopology;
    // MarchingCubes.maxVertices entries; only the first
    // topology.numVertices are meaningful.
    vertices: Vector[];
}

// Create a mesh with a zeroed topology and zero vertices.
function makeMesh(): SurfaceExtractorMCMesh {
    const vertices: Vector[] = [];
    for (let i = 0; i < MarchingCubes.maxVertices; ++i) {
        vertices.push(new Vector(3));
    }
    return {
        topology: {
            numVertices: 0,
            numTriangles: 0,
            vpair: Array.from({ length: MarchingCubes.maxVertices }, () => [0, 0]),
            itriple: Array.from({ length: MarchingCubes.maxTriangles }, () => [0, 0, 0])
        },
        vertices
    };
}

export class SurfaceExtractorMC extends MarchingCubes {
    protected mImage: Image3<number>;

    constructor(image: Image3<number>) {
        super();
        this.mImage = image;
    }

    // Extract the triangle mesh approximating F = level for a single voxel
    // whose origin corner is (x,y,z). The image has dimensions d0, d1 and
    // d2, and the origin corner satisfies 0 <= x < d0-1, 0 <= y < d1-1 and
    // 0 <= z < d2-1. The input function values must be stored as
    //  F[0] = image(x  ,y  ,z), F[4] = image(x  ,y  ,z+1),
    //  F[1] = image(x+1,y  ,z), F[5] = image(x+1,y  ,z+1),
    //  F[2] = image(x  ,y+1,z), F[6] = image(x  ,y+1,z+1),
    //  F[3] = image(x+1,y+1,z), F[7] = image(x+1,y+1,z+1)
    // In local coordinates where the corners are (0,0,0), (1,0,0), (0,1,0),
    // (1,1,0), (0,0,1), (1,0,1), (0,1,1) and (1,1,1), we have
    // F[k] = imageLocal(k & 1, (k & 2) >> 1, (k & 4) >> 2). The caller must
    // add the (x,y,z) origin corner to mesh.vertices to obtain the global
    // coordinates of the vertices.
    //
    // The returned 'valid' is true iff the F[] values are all not equal to
    // 'level'. If at least one of F[] is 'level', the returned mesh has no
    // vertices and no triangles. If you want this behavior, set 'perturb' to
    // zero.
    //
    // If you want to avoid the case when F[i] = level for some i, set
    // 'perturb' to a small nonzero number whose absolute value is smaller
    // than the minimum absolute value of the differences between voxel
    // values and 'level'.
    extractVoxel(level: number, perturb: number, F: readonly number[]):
        { valid: boolean, mesh: SurfaceExtractorMCMesh } {
        const mesh = makeMesh();
        const localF: number[] = new Array<number>(8).fill(0);

        let entry = 0;
        for (let i = 0, mask = 1; i < 8; ++i, mask <<= 1) {
            localF[i] = F[i] - level;
            if (localF[i] === 0) {
                localF[i] += perturb;
            }

            if (localF[i] < 0) {
                entry |= mask;
            } else if (localF[i] === 0) {
                // If 'perturb' is zero, the function reports that no geometry
                // is generated for this voxel. If 'perturb' is not zero, the
                // comparison to zero still needs to be made after the
                // perturbation in case floating-point rounding errors cause
                // localF[i] still to be zero.
                return { valid: false, mesh };
            }
        }

        mesh.topology = this.getTable(entry);

        for (let i = 0; i < mesh.topology.numVertices; ++i) {
            const j0 = mesh.topology.vpair[i][0];
            const j1 = mesh.topology.vpair[i][1];

            // The vertex can be computed with 3D-only computations as
            //   V = (F[j0] * k1 - F[j1] * k0) / (F[j1] - F[j0])
            // but floating-point rounding errors can cause at least one of
            // the integer-valued components of V not to be 0 or 1 as the case
            // may be. Such errors in turn can lead to 2 nearly identical
            // vertices in the mesh, and makeUnique will not be able to
            // characterize them as the same vertex. The componentwise
            // computations avoid these floating-point rounding errors. It is
            // guaranteed that j0 < j1, so multiple voxels sharing the same
            // edge will generate the same vertex.
            const vertex = mesh.vertices[i];
            const k0 = [j0 & 1, (j0 & 2) >> 1, (j0 & 4) >> 2];
            const k1 = [j1 & 1, (j1 & 2) >> 1, (j1 & 4) >> 2];
            for (let index = 0; index < 3; ++index) {
                if (k0[index] === 0) {
                    if (k1[index] === 0) {
                        vertex.values[index] = 0;
                    } else { // k1[index] = 1
                        vertex.values[index] = F[j0] / (F[j0] - F[j1]);
                    }
                } else { // k0[index] = 1
                    if (k1[index] === 0) {
                        vertex.values[index] = F[j1] / (F[j1] - F[j0]);
                    } else { // k1[index] = 1
                        vertex.values[index] = 1;
                    }
                }
            }
        }
        return { valid: true, mesh };
    }

    // Extract the triangle mesh approximating F = level for all the voxels in
    // the 3D image. The output 'indices' consists of indices.length/3
    // triangles, each a triple of indices into 'vertices'.
    //
    // The triangle table lookups depend on voxel values never being exactly
    // equal to 'level'. Set 'perturb' to zero so that any voxel cube with at
    // least one corner value equal to 'level' is ignored in the final mesh;
    // that is, such a voxel does not generate any triangles in the final
    // mesh. If you want triangles from such voxels, set 'perturb' to a small
    // nonzero number whose absolute value is smaller than the minimum
    // absolute value of the differences between voxel values and 'level'.
    extract(level: number, perturb: number):
        { vertices: Vector[], indices: number[] } {
        const vertices: Vector[] = [];
        const indices: number[] = [];

        for (let z0 = 0, z1 = 1; z1 < this.mImage.getDimension(2); z0 = z1++) {
            for (let y0 = 0, y1 = 1; y1 < this.mImage.getDimension(1); y0 = y1++) {
                for (let x0 = 0, x1 = 1; x1 < this.mImage.getDimension(0); x0 = x1++) {
                    const F: number[] = [
                        this.mImage.get(x0, y0, z0),
                        this.mImage.get(x1, y0, z0),
                        this.mImage.get(x0, y1, z0),
                        this.mImage.get(x1, y1, z0),
                        this.mImage.get(x0, y0, z1),
                        this.mImage.get(x1, y0, z1),
                        this.mImage.get(x0, y1, z1),
                        this.mImage.get(x1, y1, z1)
                    ];

                    const { valid, mesh } = this.extractVoxel(level, perturb, F);
                    if (valid) {
                        const vbase = vertices.length;
                        for (let i = 0; i < mesh.topology.numVertices; ++i) {
                            const position = mesh.vertices[i];
                            position.values[0] += x0;
                            position.values[1] += y0;
                            position.values[2] += z0;
                            vertices.push(position);
                        }

                        for (let i = 0; i < mesh.topology.numTriangles; ++i) {
                            for (let j = 0; j < 3; ++j) {
                                indices.push(vbase + mesh.topology.itriple[i][j]);
                            }
                        }
                    }
                }
            }
        }

        return { vertices, indices };
    }

    // The extraction has duplicate vertices on edges shared by voxels. This
    // function eliminates the duplication.
    makeUnique(vertices: readonly Vector[], indices: readonly number[]):
        { vertices: Vector[], indices: number[] } {
        const uvt = new UniqueVerticesSimplices<Vector>(3,
            (v) => `${v.values[0]},${v.values[1]},${v.values[2]}`);
        return uvt.removeDuplicateVertices(vertices, indices);
    }

    // The extraction does not use any topological information about the
    // level surface. The triangles can be a mixture of clockwise-ordered and
    // counterclockwise-ordered. This function is an attempt to give the
    // triangles a consistent ordering by selecting a normal in approximately
    // the same direction as the average gradient at the vertices (when
    // sameDir is true), or in the opposite direction (when sameDir is
    // false). This might not always produce a consistent order, but is fast.
    // A consistent order can be computed if you build a table of vertex, edge
    // and face adjacencies, but the resulting data structure is somewhat
    // expensive to process to reorient triangles.
    orientTriangles(vertices: readonly Vector[], indices: number[],
        sameDir: boolean): void {
        const numTriangles = Math.floor(indices.length / 3);
        for (let t = 0; t < numTriangles; ++t) {
            const base = 3 * t;

            // Get triangle vertices.
            const v0 = vertices[indices[base]];
            const v1 = vertices[indices[base + 1]];
            const v2 = vertices[indices[base + 2]];

            // Construct triangle normal based on current orientation.
            const edge1 = Vector.fromArray([
                v1.values[0] - v0.values[0],
                v1.values[1] - v0.values[1],
                v1.values[2] - v0.values[2]
            ]);
            const edge2 = Vector.fromArray([
                v2.values[0] - v0.values[0],
                v2.values[1] - v0.values[1],
                v2.values[2] - v0.values[2]
            ]);
            const normal = cross(edge1, edge2);

            // Get the image gradient at the vertices.
            const gradient0 = this.getGradient(v0);
            const gradient1 = this.getGradient(v1);
            const gradient2 = this.getGradient(v2);

            // Compute the average gradient.
            const gradientAvr = [
                (gradient0.values[0] + gradient1.values[0] + gradient2.values[0]) / 3,
                (gradient0.values[1] + gradient1.values[1] + gradient2.values[1]) / 3,
                (gradient0.values[2] + gradient1.values[2] + gradient2.values[2]) / 3
            ];

            // Compute the dot product of normal and average gradient.
            const dot = gradientAvr[0] * normal.values[0]
                + gradientAvr[1] * normal.values[1]
                + gradientAvr[2] * normal.values[2];

            // Choose triangle orientation based on gradient direction.
            const wrong = (sameDir ? dot < 0 : dot > 0);
            if (wrong) {
                // Wrong orientation, reorder it.
                const save = indices[base + 1];
                indices[base + 1] = indices[base + 2];
                indices[base + 2] = save;
            }
        }
    }

    // Compute vertex normals for the mesh. A vertex normal is the normalized
    // sum of the area-weighted normals to the triangles that share the
    // vertex.
    computeNormals(vertices: readonly Vector[],
        indices: readonly number[]): Vector[] {
        // Maintain a running sum of triangle normals at each vertex.
        const normals: Vector[] = [];
        for (let i = 0; i < vertices.length; ++i) {
            normals.push(new Vector(3));
        }

        const numTriangles = Math.floor(indices.length / 3);
        for (let t = 0; t < numTriangles; ++t) {
            const base = 3 * t;
            const i0 = indices[base];
            const i1 = indices[base + 1];
            const i2 = indices[base + 2];

            // Get triangle vertices.
            const v0 = vertices[i0];
            const v1 = vertices[i1];
            const v2 = vertices[i2];

            // Construct triangle normal.
            const edge1 = Vector.fromArray([
                v1.values[0] - v0.values[0],
                v1.values[1] - v0.values[1],
                v1.values[2] - v0.values[2]
            ]);
            const edge2 = Vector.fromArray([
                v2.values[0] - v0.values[0],
                v2.values[1] - v0.values[1],
                v2.values[2] - v0.values[2]
            ]);
            const normal = cross(edge1, edge2);

            // Maintain the sum of normals at each vertex.
            for (const i of [i0, i1, i2]) {
                for (let j = 0; j < 3; ++j) {
                    normals[i].values[j] += normal.values[j];
                }
            }
        }

        // The normal vector storage was used to accumulate the sum of
        // triangle normals. Now these vectors must be rescaled to be unit
        // length.
        for (const normal of normals) {
            const sqrLength = normal.values[0] * normal.values[0]
                + normal.values[1] * normal.values[1]
                + normal.values[2] * normal.values[2];
            const len = Math.sqrt(sqrLength);
            if (len > 0) {
                for (let i = 0; i < 3; ++i) {
                    normal.values[i] /= len;
                }
            } else {
                for (let i = 0; i < 3; ++i) {
                    normal.values[i] = 0;
                }
            }
        }

        return normals;
    }

    // The gradient of the trilinearly interpolated image at 'position',
    // which is in image (voxel) coordinates. The zero vector is returned when
    // the position is not in a voxel of the image.
    protected getGradient(inPosition: Vector): Vector {
        const position = [
            inPosition.values[0], inPosition.values[1], inPosition.values[2]
        ];

        let x = 0;
        if (position[0] >= 0) {
            x = Math.floor(position[0]);
            if (x + 1 >= this.mImage.getDimension(0)) {
                return new Vector(3);
            }
        } else {
            return new Vector(3);
        }

        let y = 0;
        if (position[1] >= 0) {
            y = Math.floor(position[1]);
            if (y + 1 >= this.mImage.getDimension(1)) {
                return new Vector(3);
            }
        } else {
            return new Vector(3);
        }

        let z = 0;
        if (position[2] >= 0) {
            z = Math.floor(position[2]);
            if (z + 1 >= this.mImage.getDimension(2)) {
                return new Vector(3);
            }
        } else {
            return new Vector(3);
        }

        position[0] -= x;
        position[1] -= y;
        position[2] -= z;
        const oneMX = 1 - position[0];
        const oneMY = 1 - position[1];
        const oneMZ = 1 - position[2];

        // Get image values at corners of voxel.
        const corners = this.mImage.getCorners8(x, y, z);
        const f000 = this.mImage.get(corners[0]);
        const f100 = this.mImage.get(corners[1]);
        const f010 = this.mImage.get(corners[2]);
        const f110 = this.mImage.get(corners[3]);
        const f001 = this.mImage.get(corners[4]);
        const f101 = this.mImage.get(corners[5]);
        const f011 = this.mImage.get(corners[6]);
        const f111 = this.mImage.get(corners[7]);

        const gradient = new Vector(3);

        let tmp0 = oneMY * (f100 - f000) + position[1] * (f110 - f010);
        let tmp1 = oneMY * (f101 - f001) + position[1] * (f111 - f011);
        gradient.values[0] = oneMZ * tmp0 + position[2] * tmp1;

        tmp0 = oneMX * (f010 - f000) + position[0] * (f110 - f100);
        tmp1 = oneMX * (f011 - f001) + position[0] * (f111 - f101);
        gradient.values[1] = oneMZ * tmp0 + position[2] * tmp1;

        tmp0 = oneMX * (f001 - f000) + position[0] * (f101 - f100);
        tmp1 = oneMX * (f011 - f010) + position[0] * (f111 - f110);
        gradient.values[2] = oneMY * tmp0 + position[1] * tmp1;

        return gradient;
    }
}
