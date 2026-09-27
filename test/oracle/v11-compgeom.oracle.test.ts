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

    // ---- MinimumAreaBox2 -------------------------------------------------

    family.case('MinimumAreaBox2.compute', (io) => {
        const useRotatingCalipers = io.boolean();
        const n = io.integer();
        const pts = points(io, n, 2);
        const mab = new MinimumAreaBox2();
        emitBox(io, mab.compute(pts, useRotatingCalipers));
        io.outReal(mab.getArea());
        const support = mab.getSupportIndices();
        for (let i = 0; i < 4; ++i) {
            io.outInt(support[i]);
        }
        const hull = mab.getHull();
        io.outInt(hull.length);
        for (const h of hull) {
            io.outInt(h);
        }
        io.outInt(mab.getNumPoints());
        io.outBool(mab.getPoints().length > 0);
    }, { exact: true });

    // Everything but mHull agrees on the hull-dimension-1 path; mHull is the
    // subject of MinimumAreaBox2.deviation.dimension1Extremes.
    family.case('MinimumAreaBox2.compute.dimension1', (io) => {
        const useRotatingCalipers = io.boolean();
        const n = io.integer();
        const pts = points(io, n, 2);
        const mab = new MinimumAreaBox2();
        emitBox(io, mab.compute(pts, useRotatingCalipers));
        io.outReal(mab.getArea());
        const support = mab.getSupportIndices();
        for (let i = 0; i < 4; ++i) {
            io.outInt(support[i]);
        }
        io.outInt(mab.getHull().length);
        io.outInt(mab.getNumPoints());
    }, { exact: true });

    // Overload 3: the caller supplies a counterclockwise convex polygon,
    // either directly or as an index subset. getNumPoints() / getPoints() are
    // not emitted here; they are the subject of the #402 deviation case.
    family.case('MinimumAreaBox2.computeConvexPolygon', (io) => {
        const useRotatingCalipers = io.boolean();
        const numPoints = io.integer();
        const pts = points(io, numPoints, 2);
        const numIndices = io.integer();
        const indices: number[] = [];
        for (let i = 0; i < numIndices; ++i) {
            indices.push(io.integer());
        }
        const mab = new MinimumAreaBox2();
        emitBox(io, mab.computeConvexPolygon(pts, indices, useRotatingCalipers));
        io.outReal(mab.getArea());
        const support = mab.getSupportIndices();
        for (let i = 0; i < 4; ++i) {
            io.outInt(support[i]);
        }
        const hull = mab.getHull();
        io.outInt(hull.length);
        for (const h of hull) {
            io.outInt(h);
        }
    }, { exact: true });

    // The port seeds the dimension-1 t-extremes from point 0 instead of from
    // the line origin's index 0, so getHull() names the real extremes.
    family.case('MinimumAreaBox2.deviation.dimension1Extremes', (io) => {
        const useRotatingCalipers = io.boolean();
        const n = io.integer();
        const pts = points(io, n, 2);
        const mab = new MinimumAreaBox2();
        const box = mab.compute(pts, useRotatingCalipers);
        io.outReal(box.extent.get(0));
        const hull = mab.getHull();
        io.outInt(hull.length);
        for (const h of hull) {
            io.outInt(h);
        }
    }, { exact: true, deviation: '#328 (dimension-1 extreme indices)' });

    // The port resets the area and the support indices on every query;
    // upstream's degenerate branches return the previous query's values.
    family.case('MinimumAreaBox2.deviation.staleState', (io) => {
        const useRotatingCalipers = io.boolean();
        const n = io.integer();
        const first = points(io, n, 2);
        const px = io.real();
        const py = io.real();
        const second = [Vector.fromArray([px, py]), Vector.fromArray([px, py])];
        const mab = new MinimumAreaBox2();
        io.outReal(mab.compute(first, useRotatingCalipers).extent.get(0));
        emitBox(io, mab.compute(second, useRotatingCalipers));
        io.outReal(mab.getArea());
        const support = mab.getSupportIndices();
        for (let i = 0; i < 4; ++i) {
            io.outInt(support[i]);
        }
    }, { exact: true, deviation: '#328 (stale area and support indices)' });

    // The port assigns mNumPoints/mPoints on the convex-polygon overload too.
    family.case('MinimumAreaBox2.deviation.polygonPoints', (io) => {
        const useRotatingCalipers = io.boolean();
        const numPoints = io.integer();
        const polygon = points(io, numPoints, 2);
        const mab = new MinimumAreaBox2();
        const box = mab.computeConvexPolygon(polygon, [], useRotatingCalipers);
        io.outReal(box.extent.get(0));
        io.outInt(mab.getNumPoints());
        io.outBool(mab.getPoints().length > 0);
    }, { exact: true, deviation: '#402 (convex-polygon overload never assigns mPoints)' });

    // The port's removeCollinearPoints compares against the most recent
    // nonzero edge, so a duplicated polygon vertex does not take the next
    // genuine corner with it.
    family.case('MinimumAreaBox2.deviation.removeCollinear', (io) => {
        const useRotatingCalipers = io.boolean();
        const numPoints = io.integer();
        const polygon = points(io, numPoints, 2);
        const mab = new MinimumAreaBox2();
        emitBox(io, mab.computeConvexPolygon(polygon, [], useRotatingCalipers));
        io.outReal(mab.getArea());
        const support = mab.getSupportIndices();
        for (let i = 0; i < 4; ++i) {
            io.outInt(support[i]);
        }
    }, { exact: true, deviation: '#286 (duplicate point drops the next corner)' });

    // ---- MinimumWidthPoints2 ---------------------------------------------

    family.case('MinimumWidthPoints2.compute', (io) => {
        const useRotatingCalipers = io.boolean();
        const n = io.integer();
        const pts = points(io, n, 2);
        const mwp = new MinimumWidthPoints2();
        emitBox(io, mwp.compute(pts, useRotatingCalipers));
    }, { exact: true });

    family.case('MinimumWidthPoints2.computeIndexed', (io) => {
        const useRotatingCalipers = io.boolean();
        const n = io.integer();
        const pts = points(io, n, 2);
        const numIndices = io.integer();
        const indices: number[] = [];
        for (let i = 0; i < numIndices; ++i) {
            indices.push(io.integer());
        }
        const mwp = new MinimumWidthPoints2();
        emitBox(io, mwp.computeIndexed(pts, indices, useRotatingCalipers));
    }, { exact: true });

    // Throw parity for the input validation of both overloads.
    family.case('MinimumWidthPoints2.invalidInputThrows', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const numIndices = io.integer();
        const indices: number[] = [];
        for (let i = 0; i < numIndices; ++i) {
            indices.push(io.integer());
        }
        const mwp = new MinimumWidthPoints2();
        emitBox(io, io.index % 2 === 0
            ? mwp.compute(pts, true)
            : mwp.computeIndexed(pts, indices, true));
    }, { exact: true });

    // ---- SeparatePoints2 -------------------------------------------------

    // The separating line is emitted only on a true return: upstream leaves
    // the caller's line holding the last tested candidate on a false return,
    // while the port returns a default-constructed line there.
    family.case('SeparatePoints2.compute', (io) => {
        const n0 = io.integer();
        const n1 = io.integer();
        const pts0 = points(io, n0, 2);
        const pts1 = points(io, n1, 2);
        const query = new SeparatePoints2();
        const result = query.compute(pts0, pts1);
        io.outBool(result.separated);
        if (result.separated) {
            io.outVec(result.separatingLine.origin);
            io.outVec(result.separatingLine.direction);
        }
    }, { exact: true });

    // The port evaluates the side of a point relative to a candidate hull
    // edge with an exact orientation predicate on the unnormalized edge
    // normal, so a rounded plane constant can no longer report overlapping
    // point sets as separated.
    family.case('SeparatePoints2.deviation.roundoff', (io) => {
        const n0 = io.integer();
        const n1 = io.integer();
        const pts0 = points(io, n0, 2);
        const pts1 = points(io, n1, 2);
        const query = new SeparatePoints2();
        io.outBool(query.compute(pts0, pts1).separated);
    }, { exact: true, deviation: '#328 (round-off in the side tests)' });

    // ---- ConstrainedDelaunay2 --------------------------------------------

    function emitPartition(io: OracleIO, partitionedEdge: readonly number[]): void {
        io.outInt(partitionedEdge.length);
        for (const v of partitionedEdge) {
            io.outInt(v);
        }
    }

    function emitInsertedEdges(io: OracleIO, cdt: ConstrainedDelaunay2): void {
        const inserted: number[] = [];
        for (const ekey of cdt.getInsertedEdges()) {
            inserted.push(ekey.V[0], ekey.V[1]);
        }
        emitSortedTuples(io, inserted, 2);
    }

    // Mirrors EmitCDT2: inserted edges, canonical compact arrays after
    // updateIndicesAdjacencies, hull edges, graph sizes.
    function emitCDT2(io: OracleIO, cdt: ConstrainedDelaunay2): void {
        emitInsertedEdges(io, cdt);
        cdt.updateIndicesAdjacencies();
        const numTriangles = cdt.getNumTriangles();
        const indices = cdt.getIndices();
        const adjacencies = cdt.getAdjacencies();
        const order = canonicalOrder(indices, numTriangles, 3);
        const rank = rankOf(order);
        emitFeatures(io, indices, adjacencies, order, rank, 3);
        const hull = cdt.getHull();
        io.outBool(true);
        emitSortedTuples(io, hull, 2);
        io.outInt(cdt.getGraph().getNumTriangles());
        io.outInt(cdt.getGraph().getNumEdges());
    }

    function readEdge(io: OracleIO): [number, number] {
        const e0 = io.integer();
        const e1 = io.integer();
        return [e0, e1];
    }

    // Both upstream operator() overloads forward to Delaunay2::operator();
    // the port has the one compute().
    family.case('ConstrainedDelaunay2.insert', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const numEdges = io.integer();
        const edges: [number, number][] = [];
        for (let k = 0; k < numEdges; ++k) {
            edges.push(readEdge(io));
        }
        const cdt = new ConstrainedDelaunay2();
        io.outBool(cdt.compute(pts));
        for (const edge of edges) {
            emitPartition(io, cdt.insert(edge));
        }
        emitCDT2(io, cdt);
    }, { exact: true });

    family.case('ConstrainedDelaunay2.insert.strip', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const cross = io.integer() !== 0;
        const left = io.integer();
        const right = io.integer();
        const cdt = new ConstrainedDelaunay2();
        io.outBool(cdt.compute(pts));
        emitPartition(io, cdt.insert([0, 1]));
        if (cross) {
            emitPartition(io, cdt.insert([left, right]));
            emitPartition(io, cdt.insert([1, 0]));
        }
        emitCDT2(io, cdt);
    }, { exact: true });

    // Throw parity for 'Invalid edge.' after the duplicate substitution.
    family.case('ConstrainedDelaunay2.insert.invalidEdgeThrows', (io) => {
        const n = io.integer();
        const pts = points(io, n, 2);
        const edge = readEdge(io);
        const cdt = new ConstrainedDelaunay2();
        const built = cdt.compute(pts);
        const partitionedEdge = cdt.insert(edge);
        io.outBool(built);
        emitPartition(io, partitionedEdge);
    }, { exact: true });

    // The port clears the inserted-edge set in compute().
    family.case('ConstrainedDelaunay2.deviation.staleInsertedEdges', (io) => {
        const nA = io.integer();
        const ptsA = points(io, nA, 2);
        const nB = io.integer();
        const ptsB = points(io, nB, 2);
        const numEdgesA = io.integer();
        const edgesA: [number, number][] = [];
        for (let k = 0; k < numEdgesA; ++k) {
            edgesA.push(readEdge(io));
        }
        const edgeB = readEdge(io);
        const cdt = new ConstrainedDelaunay2();
        cdt.compute(ptsA);
        for (const edge of edgesA) {
            cdt.insert(edge);
        }
        io.outBool(cdt.compute(ptsB));
        emitInsertedEdges(io, cdt);
        emitPartition(io, cdt.insert(edgeB));
        emitCDT2(io, cdt);
    }, { exact: true, deviation: '#325 (mInsertedEdges survives operator())' });

    // The port range-checks the raw indices before the duplicate
    // substitution and throws; upstream reads the stale element j of
    // mDuplicates left by the previous data set and inserts <r, s>.
    family.case('ConstrainedDelaunay2.deviation.duplicatesRead', (io) => {
        const nB = io.integer();
        const m = io.integer();
        const ptsA = points(io, nB + m, 2);
        const edge = readEdge(io);
        const ptsB = ptsA.slice(0, nB);
        const cdt = new ConstrainedDelaunay2();
        cdt.compute(ptsA);
        io.outBool(cdt.compute(ptsB));
        emitPartition(io, cdt.insert(edge));
        emitCDT2(io, cdt);
    }, { exact: true, deviation: '#325 (duplicates[] read before the range check)' });

    family.finish();
});
