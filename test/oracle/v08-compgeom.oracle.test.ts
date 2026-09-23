// Replays oracle/cpp/cases/v08-compgeom.cpp. Keep the two files in the same
// order; each case here reads the inputs the C++ case recorded, in the same
// order, and emits the outputs the C++ case recorded, in the same order.
//
// See the comment at the top of the C++ file for the comparison policy: the
// only libm-dependent paths (SortPointsOnCircle's atan2 and
// InscribedFixedAspectRectInQuad's quadrant index) are pinned by the
// generators, so every case is declared exact.
import { describe } from 'vitest';
import { ConvexHullSimplePolygon } from '../../src/ConvexHullSimplePolygon.js';
import { ConvexPolyhedron3 } from '../../src/ConvexPolyhedron3.js';
import { DisjointIntervals } from '../../src/DisjointIntervals.js';
import { DisjointRectangles } from '../../src/DisjointRectangles.js';
import { ExtremalQuery3PRJ } from '../../src/ExtremalQuery3PRJ.js';
import { Polyhedron3 } from '../../src/Polyhedron3.js';
import { PrimalQuery2 } from '../../src/PrimalQuery2.js';
import { PrimalQuery3 } from '../../src/PrimalQuery3.js';
import { SortPointsOnCircle } from '../../src/SortPointsOnCircle.js';
import {
    circleThroughPointSpecifiedTangentAndRadius
} from '../../src/CircleThroughPointSpecifiedTangentAndRadius.js';
import {
    circleThroughTwoPointsSpecifiedRadius
} from '../../src/CircleThroughTwoPointsSpecifiedRadius.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Vector, dot, sub } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- DisjointIntervals / DisjointRectangles ----------------------------

function inIntervalSet(s: DisjointIntervals, t: number): boolean {
    const n = s.getNumIntervals();
    for (let i = 0; i < n; ++i) {
        const iv = s.getInterval(i);
        if (iv !== null && iv.xmin <= t && t < iv.xmax) {
            return true;
        }
    }
    return false;
}

function intervalSetIsCanonical(s: DisjointIntervals): boolean {
    const n = s.getNumIntervals();
    let prevMax = -Number.MAX_VALUE;
    for (let i = 0; i < n; ++i) {
        const iv = s.getInterval(i);
        if (iv === null) {
            return false;
        }
        if (!(iv.xmin < iv.xmax) || !(prevMax < iv.xmin)) {
            return false;
        }
        prevMax = iv.xmax;
    }
    return true;
}

function outIntervalSet(io: OracleIO, s: DisjointIntervals): void {
    const n = s.getNumIntervals();
    io.outInt(n);
    for (let i = 0; i < n; ++i) {
        const iv = s.getInterval(i);
        io.outBool(iv !== null);
        io.outReal(iv !== null ? iv.xmin : 0);
        io.outReal(iv !== null ? iv.xmax : 0);
    }

    const high = s.getInterval(n);
    io.outBool(high !== null);
    io.outReal(high !== null ? high.xmin : 0);
    io.outReal(high !== null ? high.xmax : 0);

    const low = s.getInterval(-1);
    io.outBool(low !== null);
    io.outReal(low !== null ? low.xmin : 0);
    io.outReal(low !== null ? low.xmax : 0);

    io.outBool(intervalSetIsCanonical(s));
}

function intervalOpAgreesWithReference(s0: DisjointIntervals, s1: DisjointIntervals,
    result: DisjointIntervals, op: number): boolean {
    for (let k = -40; k <= 40; ++k) {
        const t = 0.25 * k;
        const in0 = inIntervalSet(s0, t);
        const in1 = inIntervalSet(s1, t);
        const expected = (op === 0 ? (in0 || in1)
            : (op === 1 ? (in0 && in1)
                : (op === 2 ? (in0 && !in1) : (in0 !== in1))));
        if (inIntervalSet(result, t) !== expected) {
            return false;
        }
    }
    return true;
}

// Every generator mode records exactly two doubles per interval.
function drawInterval(io: OracleIO): [number, number] {
    const xmin = io.real();
    const xmax = io.real();
    return [xmin, xmax];
}

function makeIntervalSet(io: OracleIO, numOps: number): DisjointIntervals {
    const s = new DisjointIntervals();
    for (let k = 0; k < numOps; ++k) {
        const [xmin, xmax] = drawInterval(io);
        s.insert(xmin, xmax);
    }
    return s;
}

function inRectangleSet(s: DisjointRectangles, x: number, y: number): boolean {
    const n = s.getNumRectangles();
    for (let i = 0; i < n; ++i) {
        const r = s.getRectangle(i);
        if (r !== null && r.xmin <= x && x < r.xmax && r.ymin <= y && y < r.ymax) {
            return true;
        }
    }
    return false;
}

function rectangleSetIsCanonical(s: DisjointRectangles): boolean {
    const n = s.getNumStrips();
    let prevMax = -Number.MAX_VALUE;
    let total = 0;
    for (let i = 0; i < n; ++i) {
        const strip = s.getStrip(i);
        if (strip === null) {
            return false;
        }
        if (!(strip.ymin < strip.ymax) || !(prevMax <= strip.ymin)) {
            return false;
        }
        if (!intervalSetIsCanonical(strip.intervalSet)) {
            return false;
        }
        prevMax = strip.ymax;
        total += strip.intervalSet.getNumIntervals();
    }
    return total === s.getNumRectangles();
}

function outRectangleSet(io: OracleIO, s: DisjointRectangles): void {
    const numRectangles = s.getNumRectangles();
    io.outInt(numRectangles);
    for (let i = 0; i < numRectangles; ++i) {
        const r = s.getRectangle(i);
        io.outBool(r !== null);
        io.outReal(r !== null ? r.xmin : 0);
        io.outReal(r !== null ? r.xmax : 0);
        io.outReal(r !== null ? r.ymin : 0);
        io.outReal(r !== null ? r.ymax : 0);
    }

    const high = s.getRectangle(numRectangles);
    io.outBool(high !== null);
    io.outReal(high !== null ? high.xmin : 0);
    io.outReal(high !== null ? high.xmax : 0);
    io.outReal(high !== null ? high.ymin : 0);
    io.outReal(high !== null ? high.ymax : 0);

    const numStrips = s.getNumStrips();
    io.outInt(numStrips);
    for (let i = 0; i < numStrips; ++i) {
        const strip = s.getStrip(i);
        io.outBool(strip !== null);
        io.outReal(strip !== null ? strip.ymin : 0);
        io.outReal(strip !== null ? strip.ymax : 0);
        outIntervalSet(io, strip !== null ? strip.intervalSet : new DisjointIntervals());
    }

    io.outBool(s.getStrip(numStrips) !== null);
    io.outBool(s.getStrip(-1) !== null);
    io.outBool(rectangleSetIsCanonical(s));
}

function rectangleOpAgreesWithReference(s0: DisjointRectangles, s1: DisjointRectangles,
    result: DisjointRectangles, op: number): boolean {
    for (let kx = -14; kx <= 14; ++kx) {
        const x = 0.5 * kx;
        for (let ky = -14; ky <= 14; ++ky) {
            const y = 0.5 * ky;
            const in0 = inRectangleSet(s0, x, y);
            const in1 = inRectangleSet(s1, x, y);
            const expected = (op === 0 ? (in0 || in1)
                : (op === 1 ? (in0 && in1)
                    : (op === 2 ? (in0 && !in1) : (in0 !== in1))));
            if (inRectangleSet(result, x, y) !== expected) {
                return false;
            }
        }
    }
    return true;
}

function drawRectangle(io: OracleIO): [number, number, number, number] {
    const [xmin, xmax] = drawInterval(io);
    const [ymin, ymax] = drawInterval(io);
    return [xmin, xmax, ymin, ymax];
}

function makeRectangleSet(io: OracleIO, numOps: number): DisjointRectangles {
    const s = new DisjointRectangles();
    for (let k = 0; k < numOps; ++k) {
        const [xmin, xmax, ymin, ymax] = drawRectangle(io);
        s.insert(xmin, xmax, ymin, ymax);
    }
    return s;
}

// ---- PrimalQuery2 / PrimalQuery3 ---------------------------------------
// Every generator mode records exactly two doubles per 2D point and three
// per 3D point, so the replay reads coordinates without reproducing the
// construction.

function readPoints(io: OracleIO, n: number, dimension: number): Vector[] {
    const P: Vector[] = new Array<Vector>(n);
    for (let i = 0; i < n; ++i) {
        P[i] = io.vec(dimension);
    }
    return P;
}

// ---- the two circle constructions --------------------------------------

function outCircles(io: OracleIO, numCircles: number, circle: Hypersphere[]): void {
    io.outInt(numCircles);
    for (let i = 0; i < 2; ++i) {
        io.outVec(circle[i].center);
        io.outReal(circle[i].radius);
    }
}

// ---- ConvexHullSimplePolygon -------------------------------------------

function hullIsValid(polygon: readonly Vector[], hull: readonly number[]): boolean {
    const m = hull.length;
    if (m < 3) {
        return false;
    }
    let area = 0;
    for (let i = 0; i < m; ++i) {
        const a = polygon[hull[i]];
        const b = polygon[hull[(i + 1) % m]];
        area += a.values[0] * b.values[1] - b.values[0] * a.values[1];
        const c = polygon[hull[(i + 2) % m]];
        const turn = (b.values[0] - a.values[0]) * (c.values[1] - b.values[1])
            - (b.values[1] - a.values[1]) * (c.values[0] - b.values[0]);
        if (turn < 0) {
            return false;
        }
    }
    if (!(area > 0)) {
        return false;
    }
    for (const p of polygon) {
        for (let i = 0; i < m; ++i) {
            const a = polygon[hull[i]];
            const b = polygon[hull[(i + 1) % m]];
            const side = (b.values[0] - a.values[0]) * (p.values[1] - a.values[1])
                - (b.values[1] - a.values[1]) * (p.values[0] - a.values[0]);
            if (side < 0) {
                return false;
            }
        }
    }
    return true;
}

// ---- meshes shared by ConvexPolyhedron3 and the extremal queries -------

function readMesh(io: OracleIO): { vertices: Vector[], indices: number[] } {
    const numVertices = io.integer();
    const vertices = readPoints(io, numVertices, 3);
    const numIndices = io.integer();
    const indices: number[] = new Array<number>(numIndices);
    for (let i = 0; i < numIndices; ++i) {
        indices[i] = io.integer();
    }
    return { vertices, indices };
}

describe('oracle: v08-compgeom', () => {
    const family = new OracleFamily('v08-compgeom');

    family.case('DisjointIntervals.insertRemove', (io) => {
        io.integer();
        const numOps = io.integer();
        const s = new DisjointIntervals();
        for (let k = 0; k < numOps; ++k) {
            const op = io.integer();
            const [xmin, xmax] = drawInterval(io);
            const success = (op === 0 ? s.insert(xmin, xmax) : s.remove(xmin, xmax));
            io.outBool(success);
            io.outInt(s.getNumIntervals());
        }
        outIntervalSet(io, s);
    }, { exact: true });

    family.case('DisjointIntervals.operators', (io) => {
        io.integer();
        const numOps0 = io.integer();
        const numOps1 = io.integer();
        const s0 = makeIntervalSet(io, numOps0);
        const s1 = makeIntervalSet(io, numOps1);

        outIntervalSet(io, s0);
        outIntervalSet(io, s1);

        const rUnion = DisjointIntervals.union(s0, s1);
        outIntervalSet(io, rUnion);
        io.outBool(intervalOpAgreesWithReference(s0, s1, rUnion, 0));

        const rIntersect = DisjointIntervals.intersection(s0, s1);
        outIntervalSet(io, rIntersect);
        io.outBool(intervalOpAgreesWithReference(s0, s1, rIntersect, 1));

        const rDifference = DisjointIntervals.difference(s0, s1);
        outIntervalSet(io, rDifference);
        io.outBool(intervalOpAgreesWithReference(s0, s1, rDifference, 2));

        const rXor = DisjointIntervals.exclusiveOr(s0, s1);
        outIntervalSet(io, rXor);
        io.outBool(intervalOpAgreesWithReference(s0, s1, rXor, 3));

        const empty = new DisjointIntervals();
        outIntervalSet(io, DisjointIntervals.union(empty, s0));
        outIntervalSet(io, DisjointIntervals.intersection(s0, empty));
        outIntervalSet(io, DisjointIntervals.difference(empty, s0));
        outIntervalSet(io, DisjointIntervals.exclusiveOr(s0, empty));
    }, { exact: true });

    family.case('DisjointRectangles.insertRemove', (io) => {
        io.integer();
        const numOps = io.integer();
        const s = new DisjointRectangles();
        for (let k = 0; k < numOps; ++k) {
            const op = io.integer();
            const [xmin, xmax, ymin, ymax] = drawRectangle(io);
            const success = (op === 0 ? s.insert(xmin, xmax, ymin, ymax)
                : s.remove(xmin, xmax, ymin, ymax));
            io.outBool(success);
            io.outInt(s.getNumRectangles());
        }
        outRectangleSet(io, s);
    }, { exact: true });

    family.case('DisjointRectangles.operators', (io) => {
        io.integer();
        const numOps0 = io.integer();
        const numOps1 = io.integer();
        const s0 = makeRectangleSet(io, numOps0);
        const s1 = makeRectangleSet(io, numOps1);

        outRectangleSet(io, s0);
        outRectangleSet(io, s1);

        const rUnion = DisjointRectangles.union(s0, s1);
        outRectangleSet(io, rUnion);
        io.outBool(rectangleOpAgreesWithReference(s0, s1, rUnion, 0));

        const rIntersect = DisjointRectangles.intersection(s0, s1);
        outRectangleSet(io, rIntersect);
        io.outBool(rectangleOpAgreesWithReference(s0, s1, rIntersect, 1));

        const rDifference = DisjointRectangles.difference(s0, s1);
        outRectangleSet(io, rDifference);
        io.outBool(rectangleOpAgreesWithReference(s0, s1, rDifference, 2));

        const rXor = DisjointRectangles.exclusiveOr(s0, s1);
        outRectangleSet(io, rXor);
        io.outBool(rectangleOpAgreesWithReference(s0, s1, rXor, 3));

        const empty = new DisjointRectangles();
        outRectangleSet(io, DisjointRectangles.union(empty, s0));
        outRectangleSet(io, DisjointRectangles.intersection(s0, empty));
        outRectangleSet(io, DisjointRectangles.difference(empty, s0));
        outRectangleSet(io, DisjointRectangles.exclusiveOr(s0, empty));
    }, { exact: true });

    family.case('PrimalQuery2.toLine', (io) => {
        io.integer();
        const n = io.integer();
        const P = readPoints(io, n, 2);
        const query = new PrimalQuery2(n, P);
        io.outInt(query.getNumVertices());

        for (let q = 0; q < 4; ++q) {
            const i = io.integer();
            const v0 = io.integer();
            const v1 = io.integer();
            io.outInt(query.toLine(i, v0, v1));

            const test = io.vec(2);
            io.outInt(query.toLine(test, v0, v1));
        }
    }, { exact: true });

    family.case('PrimalQuery2.toLineWithOrder', (io) => {
        io.integer();
        const n = io.integer();
        const P = readPoints(io, n, 2);
        const query = new PrimalQuery2(n, P);

        for (let q = 0; q < 4; ++q) {
            const i = io.integer();
            const v0 = io.integer();
            const v1 = io.integer();
            const ri = query.toLineWithOrder(i, v0, v1);
            io.outInt(ri.sign);
            io.outInt(ri.order);

            const test = io.vec(2);
            const rt = query.toLineWithOrder(test, v0, v1);
            io.outInt(rt.sign);
            io.outInt(rt.order);
        }
    }, { exact: true });

    family.case('PrimalQuery2.toTriangle', (io) => {
        io.integer();
        const n = io.integer();
        const P = readPoints(io, n, 2);
        const query = new PrimalQuery2(n, P);

        for (let q = 0; q < 4; ++q) {
            const i = io.integer();
            const v0 = io.integer();
            const v1 = io.integer();
            const v2 = io.integer();
            io.outInt(query.toTriangle(i, v0, v1, v2));

            const test = io.vec(2);
            io.outInt(query.toTriangle(test, v0, v1, v2));
        }
    }, { exact: true });

    family.case('PrimalQuery2.toCircumcircle', (io) => {
        io.integer();
        const n = io.integer();
        const P = readPoints(io, n, 2);
        const query = new PrimalQuery2(n, P);

        for (let q = 0; q < 4; ++q) {
            const i = io.integer();
            const v0 = io.integer();
            const v1 = io.integer();
            const v2 = io.integer();
            io.outInt(query.toCircumcircle(i, v0, v1, v2));

            const test = io.vec(2);
            io.outInt(query.toCircumcircle(test, v0, v1, v2));
        }
    }, { exact: true });

    family.case('PrimalQuery2.toLineExtended', (io) => {
        io.integer();
        const n = io.integer();
        const P = readPoints(io, n, 2);
        const query = new PrimalQuery2(n, P);

        for (let q = 0; q < 5; ++q) {
            const i = io.integer();
            const v0 = io.integer();
            const v1 = io.integer();
            io.outInt(query.toLineExtended(P[i], P[v0], P[v1]));
        }
    }, { exact: true });

    family.case('PrimalQuery3.toPlane', (io) => {
        io.integer();
        const n = io.integer();
        const P = readPoints(io, n, 3);
        const query = new PrimalQuery3(n, P);
        io.outInt(query.getNumVertices());

        for (let q = 0; q < 4; ++q) {
            const i = io.integer();
            const v0 = io.integer();
            const v1 = io.integer();
            const v2 = io.integer();
            io.outInt(query.toPlane(i, v0, v1, v2));

            const test = io.vec(3);
            io.outInt(query.toPlane(test, v0, v1, v2));
        }
    }, { exact: true });

    family.case('PrimalQuery3.toTetrahedron', (io) => {
        io.integer();
        const n = io.integer();
        const P = readPoints(io, n, 3);
        const query = new PrimalQuery3(n, P);

        for (let q = 0; q < 4; ++q) {
            const i = io.integer();
            const v0 = io.integer();
            const v1 = io.integer();
            const v2 = io.integer();
            const v3 = io.integer();
            io.outInt(query.toTetrahedron(i, v0, v1, v2, v3));

            const test = io.vec(3);
            io.outInt(query.toTetrahedron(test, v0, v1, v2, v3));
        }
    }, { exact: true });

    family.case('PrimalQuery3.toCircumsphere', (io) => {
        io.integer();
        const n = io.integer();
        const P = readPoints(io, n, 3);
        const query = new PrimalQuery3(n, P);

        for (let q = 0; q < 4; ++q) {
            const i = io.integer();
            const v0 = io.integer();
            const v1 = io.integer();
            const v2 = io.integer();
            const v3 = io.integer();
            io.outInt(query.toCircumsphere(i, v0, v1, v2, v3));

            const test = io.vec(3);
            io.outInt(query.toCircumsphere(test, v0, v1, v2, v3));
        }
    }, { exact: true });

    family.case('SortPointsOnCircle.byAngleAndByGeometry', (io) => {
        io.integer();
        const n = io.integer();
        const P: [number, number][] = [];
        for (let i = 0; i < n; ++i) {
            const x = io.real();
            const y = io.real();
            P.push([x, y]);
        }
        const cx = io.real();
        const cy = io.real();
        const dx = io.real();
        const dy = io.real();
        const sortCCW = io.boolean();
        const C: [number, number] = [cx, cy];
        const D: [number, number] = [dx, dy];

        const byAngle = SortPointsOnCircle.byAngle(P, C, D, sortCCW);
        const byGeometry = SortPointsOnCircle.byGeometry(P, C, D, sortCCW);

        io.outInt(byAngle.length);
        for (const i of byAngle) {
            io.outInt(i);
        }
        for (const i of byGeometry) {
            io.outInt(i);
        }
        io.outBool(byAngle.length === byGeometry.length
            && byAngle.every((v, i) => v === byGeometry[i]));
    }, { exact: true });

    family.case('SortPointsOnCircle.byGeometry.ties', (io) => {
        const n = io.integer();
        const P: [number, number][] = [];
        for (let i = 0; i < n; ++i) {
            const x = io.real();
            const y = io.real();
            P.push([x, y]);
        }
        const cx = io.real();
        const cy = io.real();
        const dx = io.real();
        const dy = io.real();
        const sortCCW = io.boolean();

        const byGeometry = SortPointsOnCircle.byGeometry(P, [cx, cy], [dx, dy], sortCCW);
        io.outInt(byGeometry.length);
        for (const i of byGeometry) {
            io.outInt(i);
        }
    }, { exact: true });

    family.case('CircleThroughTwoPointsSpecifiedRadius.compute', (io) => {
        io.integer();
        const P = io.vec(2);
        const Q = io.vec(2);
        const r = io.real();

        const result = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
        outCircles(io, result.numCircles, result.circle);

        let referenceOk = true;
        for (let i = 0; i < result.numCircles; ++i) {
            const dp = sub(P, result.circle[i].center);
            const dq = sub(Q, result.circle[i].center);
            const scale = Math.max(1, r * r);
            referenceOk = referenceOk
                && Math.abs(dot(dp, dp) - r * r) <= 1e-9 * scale
                && Math.abs(dot(dq, dq) - r * r) <= 1e-9 * scale;
        }
        io.outBool(referenceOk);
    }, { exact: true });

    family.case('CircleThroughPointSpecifiedTangentAndRadius.compute', (io) => {
        io.integer();
        io.integer();
        const P = io.vec(2);
        const A = io.vec(2);
        const N = io.vec(2);
        const r = io.real();

        const result = circleThroughPointSpecifiedTangentAndRadius(P, A, N, r);
        outCircles(io, result.numCircles, result.circle);

        let referenceOk = true;
        for (let i = 0; i < result.numCircles; ++i) {
            const dp = sub(P, result.circle[i].center);
            const scale = Math.max(1, r * r);
            referenceOk = referenceOk
                && Math.abs(dot(dp, dp) - r * r) <= 1e-9 * scale
                && Math.abs(Math.abs(dot(N, sub(result.circle[i].center, A))) - r)
                    <= 1e-9 * Math.max(1, r);
        }
        io.outBool(referenceOk);
    }, { exact: true });

    family.case('ConvexHullSimplePolygon.compute', (io) => {
        io.integer();
        const n = io.integer();
        const polygon = readPoints(io, n, 2);

        const hull = new ConvexHullSimplePolygon().compute(polygon);
        io.outInt(hull.length);
        for (const i of hull) {
            io.outInt(i);
        }
        io.outBool(hullIsValid(polygon, hull));
    }, { exact: true });

    family.case('ConvexPolyhedron3.construct', (io) => {
        const mode = io.integer();
        const mesh = readMesh(io);
        const wantPlanes = io.boolean();
        const wantAlignedBox = io.boolean();

        const polyhedron = new ConvexPolyhedron3(mesh.vertices.map((v) => v.clone()),
            mesh.indices.slice(), wantPlanes, wantAlignedBox);

        io.outInt(polyhedron.vertices.length);
        io.outInt(polyhedron.indices.length);
        io.outInt(polyhedron.planes.length);
        for (const plane of polyhedron.planes) {
            io.outVec(plane);
        }
        io.outVec(polyhedron.alignedBox.min);
        io.outVec(polyhedron.alignedBox.max);

        let outward = true;
        if (mode <= 2 && wantPlanes) {
            for (const plane of polyhedron.planes) {
                for (const v of polyhedron.vertices) {
                    const value = plane.values[0] * v.values[0]
                        + plane.values[1] * v.values[1]
                        + plane.values[2] * v.values[2] + plane.values[3];
                    outward = outward && (value <= 0);
                }
            }
        }
        io.outBool(outward);

        polyhedron.generatePlanes();
        polyhedron.generateAlignedBox();
        io.outInt(polyhedron.planes.length);
        for (const plane of polyhedron.planes) {
            io.outVec(plane);
        }
        io.outVec(polyhedron.alignedBox.min);
        io.outVec(polyhedron.alignedBox.max);
    }, { exact: true });

    family.case('ExtremalQuery3PRJ.getExtremeVertices', (io) => {
        io.integer();
        const mesh = readMesh(io);

        const polytope = new Polyhedron3(mesh.vertices, mesh.indices.length,
            mesh.indices, true);
        const query = new ExtremalQuery3PRJ(polytope);

        const normals = query.getFaceNormals();
        io.outInt(normals.length);
        for (const n of normals) {
            io.outVec(n);
        }
        io.outInt(polytope.getUniqueIndices().length);

        for (let q = 0; q < 4; ++q) {
            io.integer();
            const direction = io.vec(3);
            const r = query.getExtremeVertices(direction);
            io.outInt(r.positiveDirection);
            io.outInt(r.negativeDirection);

            let best = -Number.MAX_VALUE;
            let worst = Number.MAX_VALUE;
            for (const i of polytope.getUniqueIndices()) {
                const d = dot(direction, mesh.vertices[i]);
                best = Math.max(best, d);
                worst = Math.min(worst, d);
            }
            const dPos = dot(direction, mesh.vertices[r.positiveDirection]);
            const dNeg = dot(direction, mesh.vertices[r.negativeDirection]);
            io.outBool(dPos === best && dNeg === worst);
        }
    }, { exact: true });

    family.finish();
});
