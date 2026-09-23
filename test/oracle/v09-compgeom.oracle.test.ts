// Replays oracle/cpp/cases/v09-compgeom.cpp (ConvexHull2.h, Delaunay2.h,
// Delaunay3.h). Keep the two files in the same order.
//
// Canonicalization: upstream enumerates triangles / tetrahedra in
// std::unordered_map order and the port enumerates them in sorted-key order,
// so every feature list is sorted by the stored vertex tuple on both sides
// and the adjacency indices are remapped through the same permutation. Hull
// index lists, whose element order upstream leaves unspecified, are sorted
// lexicographically. See the C++ file's header comment.
import { describe } from 'vitest';
import { ConvexHull2 } from '../../src/ConvexHull2.js';
import { Delaunay2, Delaunay2SearchInfo } from '../../src/Delaunay2.js';
import { Delaunay3, Delaunay3SearchInfo } from '../../src/Delaunay3.js';
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

describe('oracle: v09-compgeom', () => {
    const family = new OracleFamily('v09-compgeom');

    family.case('ConvexHull2.compute', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const hull = new ConvexHull2();
        io.outBool(hull.compute(pts));
        io.outInt(hull.getDimension());
        io.outInt(hull.getNumPoints());
        io.outInt(hull.getNumUniquePoints());
        const indices = hull.getHull();
        io.outInt(indices.length);
        for (const index of indices) {
            io.outInt(index);
        }
        io.outVec(hull.getLine().origin);
        io.outVec(hull.getLine().direction);
    }, { exact: true });

    // The full 2D result, canonicalized, plus the accessor overloads.
    function emitDelaunay2(io: OracleIO, del: Delaunay2, n: number): void {
        io.outInt(del.getDimension());
        io.outInt(del.getNumVertices());
        io.outInt(del.getNumUniqueVertices());
        io.outInt(del.getNumTriangles());
        const duplicates = del.getDuplicates();
        io.outInt(duplicates.length);
        for (const duplicate of duplicates) {
            io.outInt(duplicate);
        }
        io.outInt(n);

        const numTriangles = del.getNumTriangles();
        const indices = del.getIndices();
        const adjacencies = del.getAdjacencies();
        const order = canonicalOrder(indices, numTriangles, 3);
        const rank = rankOf(order);
        emitFeatures(io, indices, adjacencies, order, rank, 3);

        const hull = del.getHull();
        io.outBool(true);
        emitSortedTuples(io, hull, 2);

        io.outVec(del.getLine().origin);
        io.outVec(del.getLine().direction);
        io.outInt(del.getGraph().getNumTriangles());
        io.outInt(del.getGraph().getNumEdges());

        const tri = del.getTriangleIndices(order[0]);
        io.outBool(tri !== null);
        for (let j = 0; j < 3; ++j) {
            io.outInt(tri === null ? 0 : tri[j]);
        }
        const adj = del.getTriangleAdjacencies(order[0]);
        io.outBool(adj !== null);
        for (let j = 0; j < 3; ++j) {
            const a = adj === null ? 0 : adj[j];
            io.outInt(a < 0 ? -1 : rank[a]);
        }
        io.outBool(del.getTriangleIndices(numTriangles) !== null);
        io.outBool(del.getTriangleAdjacencies(numTriangles) !== null);
    }

    family.case('Delaunay2.compute', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const del = new Delaunay2();
        io.outBool(del.compute(pts));
        emitDelaunay2(io, del, pts.length);
    }, { exact: true });

    family.case('Delaunay2.compute.lowDimension', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const del = new Delaunay2();
        io.outBool(del.compute(pts));
        io.outInt(del.getDimension());
        io.outInt(del.getNumUniqueVertices());
        io.outInt(del.getNumTriangles());
        io.outInt(del.getDuplicates().length);
        io.outVec(del.getLine().origin);
        io.outVec(del.getLine().direction);
    }, { exact: true });

    // The port returns the stored vertex count for every dimension; upstream
    // returns mIRVertices.size(), which is 0 after a dimension-0 or -1 input.
    family.case('Delaunay2.compute.deviation.numVertices', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const del = new Delaunay2();
        io.outBool(del.compute(pts));
        io.outInt(del.getDimension());
        io.outInt(del.getNumVertices());
    }, { exact: true, deviation: '#277 (GetNumVertices after degenerate input)' });

    // Throw parity: getHull() and getContainingTriangle() both reject a
    // dimension other than 2, so every record of this case is a throw record.
    family.case('Delaunay2.degenerateAccessorsThrow', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const del = new Delaunay2();
        const result = del.compute(pts);
        if (io.index % 2 === 0) {
            del.getHull();
            io.outBool(true);
        } else {
            const found = del.getContainingTriangle(Vector.fromArray([0.5, 0.25]),
                new Delaunay2SearchInfo());
            io.outInt(found);
        }
        io.outBool(result);
        io.outInt(del.getDimension());
    }, { exact: true });

    // The port classifies the intrinsic dimension with its own exact ToLine
    // predicate instead of the intrinsics' epsilon = 0 floating-point test.
    family.case('Delaunay2.compute.deviation.epsilon', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const del = new Delaunay2();
        io.outBool(del.compute(pts));
        io.outInt(del.getDimension());
        io.outInt(del.getNumTriangles());
    }, { exact: true, deviation: '#391 (hardcoded epsilon = 0 in Delaunay2)' });

    function emitSearch2(io: OracleIO, del: Delaunay2, q: Vector,
        info: Delaunay2SearchInfo, rank: readonly number[]): void {
        const found = del.getContainingTriangle(q, info);
        io.outInt(found < 0 ? -1 : rank[found]);
        io.outInt(info.numPath);
        for (let i = 0; i < info.numPath; ++i) {
            io.outInt(rank[info.path[i]]);
        }
        io.outInt(rank[info.finalTriangle]);
        for (let j = 0; j < 3; ++j) {
            io.outInt(info.finalV[j]);
        }
        io.outInt(rank[info.initialTriangle]);
    }

    family.case('Delaunay2.getContainingTriangle', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const del = new Delaunay2();
        io.outBool(del.compute(pts));
        io.outInt(del.getNumTriangles());
        const order = canonicalOrder(del.getIndices(), del.getNumTriangles(), 3);
        const rank = rankOf(order);
        const q = io.vec(2);
        const info = new Delaunay2SearchInfo();
        info.initialTriangle = order[0];
        emitSearch2(io, del, q, info, rank);
    }, { exact: true });

    family.case('Delaunay2.getContainingTriangle.defaultStart', (io) => {
        const pts = points(io, 3, 2);
        const del = new Delaunay2();
        io.outBool(del.compute(pts));
        io.outInt(del.getNumTriangles());
        const order = canonicalOrder(del.getIndices(), del.getNumTriangles(), 3);
        const rank = rankOf(order);
        const q = io.vec(2);
        emitSearch2(io, del, q, new Delaunay2SearchInfo(), rank);
    }, { exact: true });

    // The full 3D result, canonicalized, plus the accessor overloads.
    function emitDelaunay3(io: OracleIO, del: Delaunay3, n: number): void {
        io.outInt(del.getDimension());
        io.outInt(del.getNumVertices());
        io.outInt(del.getNumUniqueVertices());
        io.outInt(del.getNumTetrahedra());
        const duplicates = del.getDuplicates();
        io.outInt(duplicates.length);
        for (const duplicate of duplicates) {
            io.outInt(duplicate);
        }
        io.outInt(n);

        const numTetrahedra = del.getNumTetrahedra();
        const indices = del.getIndices();
        const adjacencies = del.getAdjacencies();
        const order = canonicalOrder(indices, numTetrahedra, 4);
        const rank = rankOf(order);
        emitFeatures(io, indices, adjacencies, order, rank, 4);

        const hull = del.getHull();
        io.outBool(true);
        emitSortedTuples(io, hull, 3);

        io.outVec(del.getLine().origin);
        io.outVec(del.getLine().direction);
        io.outVec(del.getPlane().normal);
        io.outReal(del.getPlane().constant);
        io.outInt(del.getGraph().getNumTetrahedra());
        io.outInt(del.getGraph().getNumTriangles());

        const tetra = del.getTetrahedronIndices(order[0]);
        io.outBool(tetra !== null);
        for (let j = 0; j < 4; ++j) {
            io.outInt(tetra === null ? 0 : tetra[j]);
        }
        const adj = del.getTetrahedronAdjacencies(order[0]);
        io.outBool(adj !== null);
        for (let j = 0; j < 4; ++j) {
            const a = adj === null ? 0 : adj[j];
            io.outInt(a < 0 ? -1 : rank[a]);
        }
        io.outBool(del.getTetrahedronIndices(numTetrahedra) !== null);
        io.outBool(del.getTetrahedronAdjacencies(numTetrahedra) !== null);
    }

    family.case('Delaunay3.compute', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        io.outBool(del.compute(pts));
        emitDelaunay3(io, del, pts.length);
    }, { exact: true, timeout: 120000 });

    family.case('Delaunay3.compute.lowDimension', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        io.outBool(del.compute(pts));
        io.outInt(del.getDimension());
        io.outInt(del.getNumUniqueVertices());
        io.outInt(del.getNumTetrahedra());
        io.outInt(del.getDuplicates().length);
        io.outVec(del.getLine().origin);
        io.outVec(del.getLine().direction);
        io.outVec(del.getPlane().normal);
        io.outReal(del.getPlane().constant);
    }, { exact: true });

    function emitSearch3(io: OracleIO, del: Delaunay3, q: Vector,
        info: Delaunay3SearchInfo, rank: readonly number[]): void {
        const found = del.getContainingTetrahedron(q, info);
        io.outInt(found < 0 ? -1 : rank[found]);
        io.outInt(info.numPath);
        for (let i = 0; i < info.numPath; ++i) {
            io.outInt(rank[info.path[i]]);
        }
        io.outInt(rank[info.finalTetrahedron]);
        for (let j = 0; j < 4; ++j) {
            io.outInt(info.finalV[j]);
        }
        io.outInt(rank[info.initialTetrahedron]);
    }

    family.case('Delaunay3.getContainingTetrahedron', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        io.outBool(del.compute(pts));
        io.outInt(del.getNumTetrahedra());
        const order = canonicalOrder(del.getIndices(), del.getNumTetrahedra(), 4);
        const rank = rankOf(order);
        const q = io.vec(3);
        const info = new Delaunay3SearchInfo();
        info.initialTetrahedron = order[0];
        emitSearch3(io, del, q, info, rank);
    }, { exact: true, timeout: 120000 });

    family.case('Delaunay3.getContainingTetrahedron.defaultStart', (io) => {
        const pts = points(io, 4, 3);
        const del = new Delaunay3();
        io.outBool(del.compute(pts));
        io.outInt(del.getNumTetrahedra());
        const order = canonicalOrder(del.getIndices(), del.getNumTetrahedra(), 4);
        const rank = rankOf(order);
        const q = io.vec(3);
        emitSearch3(io, del, q, new Delaunay3SearchInfo(), rank);
    }, { exact: true });

    // Throw parity: getHull() and getContainingTetrahedron() both reject a
    // dimension other than 3.
    family.case('Delaunay3.degenerateAccessorsThrow', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        const result = del.compute(pts);
        if (io.index % 2 === 0) {
            del.getHull();
            io.outBool(true);
        } else {
            const found = del.getContainingTetrahedron(
                Vector.fromArray([0.5, 0.25, 0.125]), new Delaunay3SearchInfo());
            io.outInt(found);
        }
        io.outBool(result);
        io.outInt(del.getDimension());
    }, { exact: true });

    // The port's ProcessedVertex hashes and compares the vertex only, so a
    // repeated input vertex is recognized; upstream's also compares the
    // location and therefore never matches.
    family.case('Delaunay3.compute.deviation.duplicates', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        io.outBool(del.compute(pts));
        io.outInt(del.getDimension());
        io.outInt(del.getNumUniqueVertices());
        io.outInt(del.getDuplicates().length);
        for (const duplicate of del.getDuplicates()) {
            io.outInt(duplicate);
        }
        io.outInt(del.getNumTetrahedra());
    }, { exact: true, deviation: '#283 (ProcessedVertex compares the location)' });

    family.case('Delaunay3.compute.deviation.epsilon', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        io.outBool(del.compute(pts));
        io.outInt(del.getDimension());
        io.outInt(del.getNumTetrahedra());
    }, { exact: true, deviation: '#391 (hardcoded epsilon = 0 in Delaunay3)' });

    family.case('Delaunay3.compute.deviation.numVertices', (io) => {
        const n = io.integer();
        const pts = points(io, n, 3);
        const del = new Delaunay3();
        io.outBool(del.compute(pts));
        io.outInt(del.getDimension());
        io.outInt(del.getNumVertices());
    }, { exact: true, deviation: '#283 (GetNumVertices after degenerate input)' });

    family.finish();
});
