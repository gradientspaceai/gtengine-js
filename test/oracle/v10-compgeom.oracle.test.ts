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

    // ---- MinimumAreaCircle2 / MinimumVolumeSphere3 -----------------------

    // Only point sets whose result does not depend on the implementation
    // defined std::shuffle permutation reach this case; see the C++ file.
    function emitCircle(io: OracleIO, points: readonly Vector[]): void {
        const query = new MinimumAreaCircle2();
        const { minimal, success } = query.compute(points);
        io.outBool(success);
        io.outReal(minimal.center.get(0));
        io.outReal(minimal.center.get(1));
        io.outReal(minimal.radius);
        io.outInt(query.numSupport);
        const support = query.support.slice(0, query.numSupport).slice().sort(
            (a, b) => a - b);
        for (const s of support) {
            io.outInt(s);
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
        const support = query.support.slice(0, query.numSupport).slice().sort(
            (a, b) => a - b);
        for (const s of support) {
            io.outInt(s);
        }
    }

    family.case('MinimumAreaCircle2.compute', (io) => {
        const n = io.integer();
        const points: Vector[] = [];
        for (let i = 0; i < n; ++i) {
            points.push(io.vec(2));
        }
        emitCircle(io, points);
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
        const n = io.integer();
        const points: Vector[] = [];
        for (let i = 0; i < n; ++i) {
            points.push(io.vec(3));
        }
        emitSphere(io, points);
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
        io.outInt(del.getNumVertices());
        io.outInt(del.getNumTriangles());
        emitTriangles(io, del);
        emitHull(io, del);
    }, { exact: true });

    family.case('IncrementalDelaunay2.getContainingTriangle', (io) => {
        const del = makeDelaunay(io);
        const n = io.integer();
        for (let i = 0; i < n; ++i) {
            del.insert(io.vec(2));
        }

        const { vertices, triangles } = del.getTriangulation();
        io.outInt(vertices.length);
        for (const v of vertices) { io.outVec(v); }
        const sorted = triangles.slice().sort((a, b) => {
            for (let j = 0; j < 3; ++j) {
                if (a[j] !== b[j]) { return a[j] - b[j]; }
            }
            return 0;
        });
        io.outInt(sorted.length);
        for (const t of sorted) {
            io.outInt(t[0]);
            io.outInt(t[1]);
            io.outInt(t[2]);
        }

        for (let k = 0; k < 4; ++k) {
            const q = io.vec(2);
            const info = new IncrementalDelaunay2SearchInfo();
            const t = del.getContainingTriangle(q, info);
            io.outBool(t >= 0);
            if (t >= 0) {
                const triangle = del.getTriangle(t);
                io.outBool(triangle !== null);
                io.outInt((triangle as number[])[0]);
                io.outInt((triangle as number[])[1]);
                io.outInt((triangle as number[])[2]);
            }
        }

        io.outBool(del.getTriangle(del.getNumTriangles()) !== null);
        io.outBool(del.getAdjacent(del.getNumTriangles()) !== null);
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

    // The port's bounded hull walk throws where upstream's unbounded walk
    // runs off the end of hull[] (issue #290). The C++ side runs a verbatim
    // copy of upstream's walk with a step cap and emits the prefix upstream
    // would write, so every record whose edges do not form a closed cycle is
    // a disagreement.
    family.case('IncrementalDelaunay2.getHull.deviation.collinear', (io) => {
        const del = makeDelaunay(io);
        const n = io.integer();
        for (let i = 0; i < n; ++i) {
            del.insert(io.vec(2));
        }
        io.outBool(del.finalizeTriangulation());
        io.outInt(del.getNumTriangles());
        emitHull(io, del);
    }, { deviation: '#290 (GetHull walks an open edge path out of bounds)' });

    family.finish();
});
