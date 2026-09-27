// Replays oracle/cpp/cases/v10-compgeom.cpp. Keep the two files in the same
// order. The C++ file documents what is comparable for each header.
import { describe } from 'vitest';
import { ExtremalQuery3BSP } from '../../src/ExtremalQuery3BSP.js';
import { IncrementalDelaunay2, IncrementalDelaunay2SearchInfo }
    from '../../src/IncrementalDelaunay2.js';
import { MinimumAreaCircle2 } from '../../src/MinimumAreaCircle2.js';
import { MinimumVolumeSphere3 } from '../../src/MinimumVolumeSphere3.js';
import { Polyhedron3 } from '../../src/Polyhedron3.js';
import { RotatingCalipers } from '../../src/RotatingCalipers.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// Read a recorded point list: the count first, then the coordinates.
function readPoints(io: OracleIO, dim: number): Vector[] {
    const count = io.integer();
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        points.push(io.vec(dim));
    }
    return points;
}

function emitAntipodes(io: OracleIO, polygon: readonly Vector[]): void {
    const antipodes = RotatingCalipers.computeAntipodes(polygon);
    io.outInt(antipodes.length);
    for (const a of antipodes) {
        io.outInt(a.vertex);
        io.outInt(a.edge[0]);
        io.outInt(a.edge[1]);
    }
}

describe('oracle: v10-compgeom', () => {
    const family = new OracleFamily('v10-compgeom');

    // ---- RotatingCalipers ------------------------------------------------

    family.case('RotatingCalipers.computeAntipodes', (io) => {
        emitAntipodes(io, readPoints(io, 2));
    }, { exact: true });

    family.case('RotatingCalipers.computeAntipodes.collinear', (io) => {
        emitAntipodes(io, readPoints(io, 2));
    }, { exact: true });

    // The port compares each candidate corner against the most recent
    // *nonzero* edge, so a duplicated vertex no longer discards the next real
    // corner (issue #286). Upstream keeps one corner fewer, and throws when
    // that leaves fewer than three.
    family.case('RotatingCalipers.computeAntipodes.deviation.duplicate', (io) => {
        emitAntipodes(io, readPoints(io, 2));
    }, { deviation: '#286 (RotatingCalipers::CreatePolygon drops the corner after a duplicate)' });

    family.case('RotatingCalipers.computeAntipodes.collinearThrows', (io) => {
        emitAntipodes(io, readPoints(io, 2));
    }, { exact: true });

    // ---- MinimumAreaCircle2 / MinimumVolumeSphere3 -----------------------

    // The C++ side makes upstream run the port's own shuffle permutation (a
    // lockstep std::default_random_engine; see PortOrderResult in the C++
    // file), so the result, the support set and the support order are all
    // compared as they are. A fresh query per record: the port's first
    // compute() applies the permutation the C++ side reproduced.
    function emitCircle(io: OracleIO, points: readonly Vector[]): void {
        const query = new MinimumAreaCircle2();
        const { minimal, success } = query.compute(points);
        io.outBool(success);
        io.outReal(minimal.center.get(0));
        io.outReal(minimal.center.get(1));
        io.outReal(minimal.radius);
        io.outInt(query.numSupport);
        for (let i = 0; i < query.numSupport; ++i) {
            io.outInt(query.support[i]);
        }
    }

    function emitSphere(io: OracleIO, points: readonly Vector[]): void {
        const query = new MinimumVolumeSphere3();
        const { minimal, success } = query.compute(points);
        io.outBool(success);
        io.outReal(minimal.center.get(0));
        io.outReal(minimal.center.get(1));
        io.outReal(minimal.center.get(2));
        io.outReal(minimal.radius);
        io.outInt(query.numSupport);
        for (let i = 0; i < query.numSupport; ++i) {
            io.outInt(query.support[i]);
        }
    }

    family.case('MinimumAreaCircle2.compute', (io) => {
        emitCircle(io, readPoints(io, 2));
    }, { exact: true });

    // The port passes the whole input array to getContainerCircle2 where
    // upstream passes the unique-point count with the full array and so
    // bounds only a prefix (issue #286).
    family.case('MinimumAreaCircle2.compute.deviation.trappedFallback', (io) => {
        emitCircle(io, readPoints(io, 2));
    }, { deviation: '#286 (the trapped-failure fallback bounds only a prefix)' });

    family.case('MinimumAreaCircle2.compute.empty', (io) => {
        io.vec(2);
        const r = new MinimumAreaCircle2().compute([]);
        io.outBool(r.success);
    }, { exact: true });

    // Aimed at the three-point support: an acute lattice triangle plus
    // points the exact in-circle predicate places strictly inside its
    // circumcircle, so exactCircle3's output is the emitted result.
    family.case('MinimumAreaCircle2.compute.circumcircle', (io) => {
        io.integer();
        emitCircle(io, readPoints(io, 2));
    }, { exact: true });

    family.case('MinimumVolumeSphere3.compute', (io) => {
        emitSphere(io, readPoints(io, 3));
    }, { exact: true });

    family.case('MinimumVolumeSphere3.compute.deviation.trappedFallback', (io) => {
        emitSphere(io, readPoints(io, 3));
    }, { deviation: '#286 (the trapped-failure fallback bounds only a prefix)' });

    family.case('MinimumVolumeSphere3.compute.empty', (io) => {
        io.vec(3);
        const r = new MinimumVolumeSphere3().compute([]);
        io.outBool(r.success);
    }, { exact: true });

    // ---- ExtremalQuery3BSP -----------------------------------------------

    // The base polytopes of the C++ case, indexed by the recorded selector.
    const BASE_INDICES: readonly (readonly number[])[] = [
        [0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2],
        [0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
            2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5],
        [0, 2, 3, 0, 3, 4, 0, 4, 2, 1, 3, 2, 1, 4, 3, 1, 2, 4]
    ];
    const BASE_NUM_VERTICES: readonly number[] = [4, 6, 5];

    family.case('ExtremalQuery3BSP.getExtremeVertices', (io) => {
        const which = io.integer();
        const vertices: Vector[] = [];
        for (let i = 0; i < BASE_NUM_VERTICES[which]; ++i) {
            vertices.push(io.vec(3));
        }
        const indices = BASE_INDICES[which].slice();
        const polytope = new Polyhedron3(vertices, indices.length, indices, true);
        const query = new ExtremalQuery3BSP(polytope);

        const normals = query.getFaceNormals();
        io.outInt(normals.length);
        for (const n of normals) {
            io.outVec(n);
        }

        for (let k = 0; k < 6; ++k) {
            const direction = io.vec(3);
            const r = query.getExtremeVertices(direction);
            io.outInt(r.positiveDirection);
            io.outInt(r.negativeDirection);
        }
    }, { exact: true });

    // ---- IncrementalDelaunay2 --------------------------------------------

    // The triangle numbering is hash-table order upstream and sorted order
    // in the port, so the features are sorted by their stored vertex tuple
    // and the adjacency indices are remapped through the same permutation.
    function emitTriangles(io: OracleIO, del: IncrementalDelaunay2): void {
        const tris = del.getTriangles();
        const adjs = del.getAdjacencies();
        const order = tris.map((_, i) => i);
        order.sort((a, b) => {
            for (let j = 0; j < 3; ++j) {
                if (tris[a][j] !== tris[b][j]) { return tris[a][j] - tris[b][j]; }
            }
            return 0;
        });
        const rank = new Array<number>(order.length);
        order.forEach((t, r) => { rank[t] = r; });

        io.outInt(tris.length);
        for (const t of order) {
            for (let j = 0; j < 3; ++j) { io.outInt(tris[t][j]); }
            for (let j = 0; j < 3; ++j) {
                const a = adjs[t][j];
                io.outReal(a < 0 ? -1 : rank[a]);
            }
        }
    }

    function compareTuples(a: readonly number[], b: readonly number[]): number {
        for (let j = 0; j < 3; ++j) {
            if (a[j] !== b[j]) { return a[j] - b[j]; }
        }
        return 0;
    }

    function emitHull(io: OracleIO, del: IncrementalDelaunay2): void {
        const hull = del.getHull();
        io.outInt(hull.length);
        for (const v of hull) { io.outInt(v); }
    }

    function makeDelaunay(io: OracleIO): IncrementalDelaunay2 {
        const xMin = io.real();
        const yMin = io.real();
        const xMax = io.real();
        const yMax = io.real();
        return new IncrementalDelaunay2(xMin, yMin, xMax, yMax);
    }

    family.case('IncrementalDelaunay2.insert', (io) => {
        const del = makeDelaunay(io);
        const n = io.integer();
        for (let i = 0; i < n; ++i) {
            io.outReal(del.insert(io.vec(2)));
        }
        io.outInt(del.getNumVertices());
        io.outInt(del.getNumTriangles());
        emitTriangles(io, del);
        emitHull(io, del);
    }, { exact: true });

    family.case('IncrementalDelaunay2.remove', (io) => {
        const del = makeDelaunay(io);
        const n = io.integer();
        const numRemove = io.integer();
        for (let i = 0; i < n; ++i) {
            del.insert(io.vec(2));
        }
        for (let k = 0; k < numRemove; ++k) {
            io.outReal(del.remove(io.vec(2)));
        }
        const numAfter = io.integer();
        for (let k = 0; k < numAfter; ++k) {
            io.outReal(del.insert(io.vec(2)));
        }
        io.outInt(del.getNumVertices());
        io.outInt(del.getNumTriangles());
        emitTriangles(io, del);
        emitHull(io, del);
    }, { exact: true });

    // The walk starts at a recorded canonical rank (triangles sorted by
    // their stored vertex tuple), and every triangle index of the SearchInfo
    // is emitted as a canonical rank; see the C++ case.
    family.case('IncrementalDelaunay2.getContainingTriangle', (io) => {
        const del = makeDelaunay(io);
        const n = io.integer();
        for (let i = 0; i < n; ++i) {
            del.insert(io.vec(2));
        }

        const { vertices, triangles } = del.getTriangulation();
        io.outInt(vertices.length);
        for (const v of vertices) { io.outVec(v); }
        const sorted = triangles.slice().sort(compareTuples);
        io.outInt(sorted.length);
        for (const t of sorted) {
            io.outInt(t[0]);
            io.outInt(t[1]);
            io.outInt(t[2]);
        }

        const tris = del.getTriangles();
        const order = tris.map((_, i) => i);
        order.sort((a, b) => compareTuples(tris[a], tris[b]));
        const rank = new Array<number>(order.length);
        order.forEach((t, r) => { rank[t] = r; });
        const asRank = (t: number): number => (t < 0 ? -1 : rank[t]);

        for (let k = 0; k < 4; ++k) {
            const q = io.vec(2);
            const startRank = io.integer();
            const info = new IncrementalDelaunay2SearchInfo();
            info.initialTriangle = order[startRank];
            const t = del.getContainingTriangle(q, info);
            io.outReal(asRank(t));
            io.outReal(asRank(info.initialTriangle));
            io.outReal(asRank(info.finalTriangle));
            io.outInt(info.finalV[0]);
            io.outInt(info.finalV[1]);
            io.outInt(info.finalV[2]);
            io.outInt(info.numPath);
            for (let i = 0; i < info.numPath; ++i) {
                io.outReal(asRank(info.path[i]));
            }
            if (t >= 0) {
                const triangle = del.getTriangle(t);
                io.outBool(triangle !== null);
                io.outInt((triangle as number[])[0]);
                io.outInt((triangle as number[])[1]);
                io.outInt((triangle as number[])[2]);
                const adjacent = del.getAdjacent(t);
                io.outBool(adjacent !== null);
                for (let j = 0; j < 3; ++j) {
                    io.outReal(asRank((adjacent as number[])[j]));
                }
            }
        }

        io.outBool(del.getTriangle(tris.length) !== null);
        io.outBool(del.getAdjacent(tris.length) !== null);
    }, { exact: true });

    family.case('IncrementalDelaunay2.finalizeTriangulation', (io) => {
        const del = makeDelaunay(io);
        const n = io.integer();
        const points: Vector[] = [];
        for (let i = 0; i < n; ++i) {
            const p = io.vec(2);
            points.push(p);
            del.insert(p);
        }
        io.outBool(del.finalizeTriangulation());
        io.outInt(del.getNumVertices());
        io.outInt(del.getNumTriangles());
        emitTriangles(io, del);
        emitHull(io, del);
        io.outBool(del.finalizeTriangulation());
        io.outReal(del.insert(points[0]));
        io.outReal(del.remove(points[0]));
    }, { exact: true });

    family.case('IncrementalDelaunay2.domainAsserts', (io) => {
        const mode = io.index % 4;
        const xMin = io.real();
        const yMin = io.real();
        let xMax = 0;
        let yMax = 0;
        if (mode === 0) {
            const flatX = io.boolean();
            const lo = io.real();
            xMax = flatX ? xMin + lo : xMin + 4;
            yMax = flatX ? yMin + 4 : yMin + lo;
        } else {
            xMax = xMin + 8;
            yMax = yMin + 8;
        }
        const del = new IncrementalDelaunay2(xMin, yMin, xMax, yMax);
        const first = del.insert(Vector.fromArray([xMin + 3, yMin + 5]));

        const t = io.real();
        const side = io.integer();
        const off = io.real();
        const p = side === 0 ? Vector.fromArray([xMin - off, yMin + t])
            : side === 1 ? Vector.fromArray([xMax + off, yMin + t])
                : side === 2 ? Vector.fromArray([xMin + t, yMin - off])
                    : Vector.fromArray([xMin + t, yMax + off]);
        if (mode === 1) {
            const index = del.insert(p);
            io.outReal(first);
            io.outReal(index);
        } else if (mode === 2) {
            const index = del.remove(p);
            io.outReal(first);
            io.outReal(index);
        } else if (mode === 3) {
            io.outReal(first);
            io.outBool(del.finalizeTriangulation());
            io.outReal(del.remove(p));
        }
    }, { exact: true });

    // The port throws where upstream's hull walk returns to its start early
    // (a hull[] padded with supervertex 0) or runs off the end of hull[]
    // (issue #290). The C++ side calls the real GetHull when that is safe
    // and emits the replica's prefix when it is not; see the C++ case.
    family.case('IncrementalDelaunay2.getHull.deviation.collinear', (io) => {
        const del = makeDelaunay(io);
        const n = io.integer();
        for (let i = 0; i < n; ++i) {
            del.insert(io.vec(2));
        }
        io.outBool(del.finalizeTriangulation());
        io.outInt(del.getNumTriangles());
        const hull = del.getHull();
        io.outBool(false);
        io.outInt(hull.length);
        for (const v of hull) { io.outInt(v); }
    }, { deviation: '#290 (GetHull on a collinear triangulation: early return padded with 0, or an unbounded walk)' });

    family.finish();
});
