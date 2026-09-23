// Replays oracle/cpp/cases/v11-compgeom.cpp (ConstrainedDelaunay2.h,
// ConvexHull3.h, Delaunay2Mesh.h, Delaunay3Mesh.h, MinimumAreaBox2.h,
// MinimumWidthPoints2.h, SeparatePoints2.h). Keep the two files in the same
// order.
//
// Canonicalization: upstream enumerates mesh features in std::unordered_map
// order and inserted edges in std::unordered_set order, while the port
// enumerates them in sorted-key order, so every such list is sorted on both
// sides and adjacency indices are remapped through the same permutation. See
// the C++ file's header comment.
import { describe } from 'vitest';
import { ConstrainedDelaunay2 } from '../../src/ConstrainedDelaunay2.js';
import { ConvexHull3 } from '../../src/ConvexHull3.js';
import { Delaunay2 } from '../../src/Delaunay2.js';
import { Delaunay2Mesh } from '../../src/Delaunay2Mesh.js';
import { Delaunay3 } from '../../src/Delaunay3.js';
import { Delaunay3Mesh } from '../../src/Delaunay3Mesh.js';
import { MinimumAreaBox2 } from '../../src/MinimumAreaBox2.js';
import { MinimumWidthPoints2 } from '../../src/MinimumWidthPoints2.js';
import { SeparatePoints2 } from '../../src/SeparatePoints2.js';
import type { OrientedBox2 } from '../../src/OrientedBox.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

function points(io: OracleIO, n: number, dimension: number): Vector[] {
    const pts: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        pts.push(io.vec(dimension));
    }
    return pts;
}

// order[r] is the index of the feature whose stored vertex tuple has rank r.
function canonicalOrder(indices: readonly number[], count: number,
    width: number): number[] {
    const order: number[] = [];
    for (let i = 0; i < count; ++i) {
        order.push(i);
    }
    order.sort((a, b) => {
        for (let j = 0; j < width; ++j) {
            if (indices[width * a + j] !== indices[width * b + j]) {
                return indices[width * a + j] - indices[width * b + j];
            }
        }
        return 0;
    });
    return order;
}

function rankOf(order: readonly number[]): number[] {
    const rank = new Array<number>(order.length).fill(0);
    for (let r = 0; r < order.length; ++r) {
        rank[order[r]] = r;
    }
    return rank;
}

function emitFeatures(io: OracleIO, indices: readonly number[],
    adjacencies: readonly number[], order: readonly number[],
    rank: readonly number[], width: number): void {
    io.outInt(order.length);
    for (let r = 0; r < order.length; ++r) {
        const t = order[r];
        for (let j = 0; j < width; ++j) {
            io.outInt(indices[width * t + j]);
        }
        for (let j = 0; j < width; ++j) {
            const a = adjacencies[width * t + j];
            io.outInt(a < 0 ? -1 : rank[a]);
        }
    }
}

function emitSortedTuples(io: OracleIO, flat: readonly number[],
    width: number): void {
    const count = flat.length / width;
    const tuples: number[][] = [];
    for (let i = 0; i < count; ++i) {
        tuples.push(flat.slice(width * i, width * i + width));
    }
    tuples.sort((a, b) => {
        for (let j = 0; j < width; ++j) {
            if (a[j] !== b[j]) {
                return a[j] - b[j];
            }
        }
        return 0;
    });
    io.outInt(count);
    for (const tuple of tuples) {
        for (let j = 0; j < width; ++j) {
            io.outInt(tuple[j]);
        }
    }
}

// The oriented box, field by field, in upstream's member order.
function emitBox(io: OracleIO, box: OrientedBox2): void {
    io.outVec(box.center);
    io.outVec(box.axis[0]);
    io.outVec(box.axis[1]);
    io.outVec(box.extent);
}

describe('oracle: v11-compgeom', () => {
    const family = new OracleFamily('v11-compgeom');

    // ---- ConvexHull3 -----------------------------------------------------

    family.case('ConvexHull3.compute', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const hull = new ConvexHull3();
        hull.compute(pts);
        const dimension = hull.getDimension();
        io.outInt(dimension);

        const vertices = [...hull.getVertices()];
        const indices = hull.getHull();
        if (dimension === 3) {
            vertices.sort((a, b) => a - b);
            io.outInt(vertices.length);
            for (const v of vertices) {
                io.outInt(v);
            }
            emitSortedTuples(io, indices, 3);
        }
        else {
            io.outInt(vertices.length);
            for (const v of vertices) {
                io.outInt(v);
            }
            io.outInt(indices.length);
            for (const v of indices) {
                io.outInt(v);
            }
        }

        const mesh = hull.getHullMesh();
        io.outInt(mesh.getNumVertices());
        io.outInt(mesh.getNumEdges());
        io.outInt(mesh.getNumTriangles());
    }, { exact: true });

    // ---- Delaunay2Mesh ---------------------------------------------------

    family.case('Delaunay2Mesh.query', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const del = new Delaunay2();
        io.outBool(del.compute(pts));

        const mesh = new Delaunay2Mesh(del);
        io.outInt(mesh.getNumVertices());
        io.outInt(mesh.getNumTriangles());

        const numTriangles = mesh.getNumTriangles();
        const vertices = mesh.getVertices();
        for (let i = 0; i < n; ++i) {
            io.outVec(vertices[i]);
        }

        const indices = mesh.getIndices();
        const adjacencies = mesh.getAdjacencies();
        const order = canonicalOrder(indices, numTriangles, 3);
        const rank = rankOf(order);
        emitFeatures(io, indices, adjacencies, order, rank, 3);

        io.outInt(mesh.getInvalidIndex());

        const q = io.vec(2);
        const found = mesh.getContainingTriangle(q);
        io.outInt(found === mesh.getInvalidIndex() ? -1 : rank[found]);

        const r = io.integer();
        const t = order[r];

        const triVertices = mesh.getTriangleVertices(t);
        io.outBool(triVertices !== null);
        for (let j = 0; j < 3; ++j) {
            io.outVec(triVertices === null ? Vector.zero(2) : triVertices[j]);
        }

        const triIndices = mesh.getTriangleIndices(t);
        io.outBool(triIndices !== null);
        for (let j = 0; j < 3; ++j) {
            io.outInt(triIndices === null ? 0 : triIndices[j]);
        }

        const triAdjacencies = mesh.getTriangleAdjacencies(t);
        io.outBool(triAdjacencies !== null);
        for (let j = 0; j < 3; ++j) {
            const a = triAdjacencies === null ? 0 : triAdjacencies[j];
            io.outInt(a < 0 ? -1 : rank[a]);
        }

        const bary = mesh.getBarycentrics(t, q);
        io.outBool(bary !== null);
        for (let j = 0; j < 3; ++j) {
            io.outReal(bary === null ? 0 : bary[j]);
        }

        // The out-of-range branch. Upstream zeroes the output arrays of
        // GetVertices and GetBarycentrics; the port returns null.
        const badVertices = mesh.getTriangleVertices(numTriangles);
        io.outBool(badVertices !== null);
        for (let j = 0; j < 3; ++j) {
            io.outVec(badVertices === null ? Vector.zero(2) : badVertices[j]);
        }
        io.outBool(mesh.getTriangleIndices(numTriangles) !== null);
        io.outBool(mesh.getTriangleAdjacencies(numTriangles) !== null);
        const badBary = mesh.getBarycentrics(numTriangles, q);
        io.outBool(badBary !== null);
        for (let j = 0; j < 3; ++j) {
            io.outReal(badBary === null ? 0 : badBary[j]);
        }
    }, { exact: true });

    // Throw parity: the constructor rejects a dimension other than 2.
    family.case('Delaunay2Mesh.constructorThrows', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const del = new Delaunay2();
        const built = del.compute(pts);
        const mesh = new Delaunay2Mesh(del);
        io.outBool(built);
        io.outInt(mesh.getNumTriangles());
    }, { exact: true });

    // ---- Delaunay3Mesh ---------------------------------------------------

    family.case('Delaunay3Mesh.query', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        io.outBool(del.compute(pts));

        const mesh = new Delaunay3Mesh(del);
        io.outInt(mesh.getNumVertices());
        io.outInt(mesh.getNumTetrahedra());

        const numTetrahedra = mesh.getNumTetrahedra();
        const vertices = mesh.getVertices();
        for (let i = 0; i < n; ++i) {
            io.outVec(vertices[i]);
        }

        const indices = mesh.getIndices();
        const adjacencies = mesh.getAdjacencies();
        const order = canonicalOrder(indices, numTetrahedra, 4);
        const rank = rankOf(order);
        emitFeatures(io, indices, adjacencies, order, rank, 4);

        io.outInt(mesh.getInvalidIndex());

        const q = io.vec(3);
        const found = mesh.getContainingTetrahedron(q);
        io.outInt(found === mesh.getInvalidIndex() ? -1 : rank[found]);

        const r = io.integer();
        const t = order[r];

        const tetVertices = mesh.getTetrahedronVertices(t);
        io.outBool(tetVertices !== null);
        for (let j = 0; j < 4; ++j) {
            io.outVec(tetVertices === null ? Vector.zero(3) : tetVertices[j]);
        }

        const tetIndices = mesh.getTetrahedronIndices(t);
        io.outBool(tetIndices !== null);
        for (let j = 0; j < 4; ++j) {
            io.outInt(tetIndices === null ? 0 : tetIndices[j]);
        }

        const tetAdjacencies = mesh.getTetrahedronAdjacencies(t);
        io.outBool(tetAdjacencies !== null);
        for (let j = 0; j < 4; ++j) {
            const a = tetAdjacencies === null ? 0 : tetAdjacencies[j];
            io.outInt(a < 0 ? -1 : rank[a]);
        }

        const bary = mesh.getBarycentrics(t, q);
        io.outBool(bary !== null);
        for (let j = 0; j < 4; ++j) {
            io.outReal(bary === null ? 0 : bary[j]);
        }

        const badVertices = mesh.getTetrahedronVertices(numTetrahedra);
        io.outBool(badVertices !== null);
        for (let j = 0; j < 4; ++j) {
            io.outVec(badVertices === null ? Vector.zero(3) : badVertices[j]);
        }
        io.outBool(mesh.getTetrahedronIndices(numTetrahedra) !== null);
        io.outBool(mesh.getTetrahedronAdjacencies(numTetrahedra) !== null);
        const badBary = mesh.getBarycentrics(numTetrahedra, q);
        io.outBool(badBary !== null);
        for (let j = 0; j < 4; ++j) {
            io.outReal(badBary === null ? 0 : badBary[j]);
        }
    }, { exact: true, timeout: 180000 });

    // Throw parity: the constructor rejects a dimension other than 3.
    family.case('Delaunay3Mesh.constructorThrows', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        const built = del.compute(pts);
        const mesh = new Delaunay3Mesh(del);
        io.outBool(built);
        io.outInt(mesh.getNumTetrahedra());
    }, { exact: true });

    family.finish();
});
