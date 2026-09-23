// Replays oracle/cpp/cases/v08-compgeom.cpp. Keep the two files in the same
// order; each case here reads the inputs the C++ case recorded, in the same
// order, and emits the outputs the C++ case recorded, in the same order.
//
// See the comment at the top of the C++ file for the comparison policy: the
// only libm-dependent paths (SortPointsOnCircle's atan2 and
// InscribedFixedAspectRectInQuad's quadrant index) are pinned by the
// generators, so every case is declared exact.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { BoxManager } from '../../src/BoxManager.js';
import { ConvexHullSimplePolygon } from '../../src/ConvexHullSimplePolygon.js';
import { ConvexPolyhedron3 } from '../../src/ConvexPolyhedron3.js';
import { DisjointIntervals } from '../../src/DisjointIntervals.js';
import { DisjointRectangles } from '../../src/DisjointRectangles.js';
import { ExtremalQuery3PRJ } from '../../src/ExtremalQuery3PRJ.js';
import {
    InscribedFixedAspectRectInQuad
} from '../../src/InscribedFixedAspectRectInQuad.js';
import {
    MinimumVolumeBox3FloatingPoint
} from '../../src/MinimumVolumeBox3FloatingPoint.js';
import { MinimumVolumeBox3Rational } from '../../src/MinimumVolumeBox3Rational.js';
import { NearestNeighborQuery, PositionSite } from '../../src/NearestNeighborQuery.js';
import type { OrientedBox3 } from '../../src/OrientedBox.js';
import { PolygonTree } from '../../src/PolygonTree.js';
import { Polyhedron3 } from '../../src/Polyhedron3.js';
import { TriangulateEC } from '../../src/TriangulateEC.js';
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

// ---- NearestNeighborQuery ----------------------------------------------
// std::nth_element leaves an implementation-defined permutation, so the
// sorted-point array is never compared. What is compared is the tree shape,
// the split values and, per leaf, the SORTED list of the original site
// indices, all of which follow from the partition postconditions because the
// generator gives every axis pairwise distinct coordinates.

function outTree(io: OracleIO, query: NearestNeighborQuery, n: number): void {
    io.outInt(query.getMaxLeafSize());
    io.outInt(query.getMaxLevel());
    io.outInt(query.getDepth());
    io.outInt(query.getLargestNodeSize());
    io.outInt(query.getNumNodes());

    const sortedPoints = query.getSortedPoints();
    io.outInt(sortedPoints.length);

    for (const node of query.getNodes()) {
        io.outReal(node.split);
        io.outInt(node.axis);
        io.outInt(node.numSites);
        io.outInt(node.siteOffset);
        io.outInt(node.left);
        io.outInt(node.right);
        if (node.siteOffset !== -1) {
            const leaf: number[] = [];
            for (let k = 0; k < node.numSites; ++k) {
                leaf.push(sortedPoints[node.siteOffset + k].index);
            }
            leaf.sort((a, b) => a - b);
            for (const v of leaf) {
                io.outInt(v);
            }
        }
    }
    io.outInt(n);
}

function readSites(io: OracleIO, n: number): { points: Vector[], sites: PositionSite[] } {
    const points = readPoints(io, n, 3);
    return { points, sites: points.map((p) => new PositionSite(p)) };
}

// ---- BoxManager --------------------------------------------------------
// The overlap container is a std::set<EdgeKey<false>> upstream; the port's
// getOverlap() returns its keys sorted lexicographically by (V[0], V[1]),
// which is the std::set iteration order, so the two lists are directly
// comparable.

function readBox(io: OracleIO): AlignedBox {
    const min = new Vector(3);
    const max = new Vector(3);
    for (let d = 0; d < 3; ++d) {
        const a = io.real();
        const b = io.real();
        min.values[d] = Math.min(a, b);
        max.values[d] = Math.max(a, b);
    }
    return AlignedBox.fromMinMax(min, max);
}

function boxesOverlap(b0: AlignedBox, b1: AlignedBox): boolean {
    for (let d = 0; d < 3; ++d) {
        if (b0.max.values[d] < b1.min.values[d] || b0.min.values[d] > b1.max.values[d]) {
            return false;
        }
    }
    return true;
}

function outOverlap(io: OracleIO, manager: BoxManager, boxes: readonly AlignedBox[]): void {
    const overlap = manager.getOverlap();
    io.outInt(overlap.length);
    for (const key of overlap) {
        io.outInt(key.V[0]);
        io.outInt(key.V[1]);
    }

    const expected = new Set<string>();
    for (let i = 0; i < boxes.length; ++i) {
        for (let j = i + 1; j < boxes.length; ++j) {
            if (boxesOverlap(boxes[i], boxes[j])) {
                expected.add(`${i},${j}`);
            }
        }
    }
    const reported = new Set<string>();
    for (const key of overlap) {
        reported.add(`${key.V[0]},${key.V[1]}`);
    }
    let same = reported.size === expected.size;
    if (same) {
        for (const k of reported) {
            if (!expected.has(k)) {
                same = false;
                break;
            }
        }
    }
    io.outBool(same);
}

// ---- TriangulateEC -----------------------------------------------------
// Upstream is TriangulateEC<InputType, ComputeType>, instantiated on the
// C++ side as <double, double> to match this number-only port. Ear clipping
// produces its triangle list in a deterministic order with no container
// whose order C++ leaves unspecified, so the lists are compared triple by
// triple in the order produced.

function readPolygon(io: OracleIO): number[] {
    const n = io.integer();
    const polygon: number[] = new Array<number>(n);
    for (let i = 0; i < n; ++i) {
        polygon[i] = io.integer();
    }
    return polygon;
}

function shoelaceArea(points: readonly Vector[], polygon: readonly number[]): number {
    let twiceArea = 0;
    const m = polygon.length;
    for (let i = 0; i < m; ++i) {
        const a = points[polygon[i]];
        const b = points[polygon[(i + 1) % m]];
        twiceArea += a.values[0] * b.values[1] - b.values[0] * a.values[1];
    }
    return twiceArea;
}

function triangulationIsValid(points: readonly Vector[],
    triangles: readonly [number, number, number][], expectedTwiceArea: number): boolean {
    let twiceArea = 0;
    for (const t of triangles) {
        const a = points[t[0]];
        const b = points[t[1]];
        const c = points[t[2]];
        const d = (b.values[0] - a.values[0]) * (c.values[1] - a.values[1])
            - (b.values[1] - a.values[1]) * (c.values[0] - a.values[0]);
        if (d < 0) {
            return false;
        }
        twiceArea += d;
    }
    return twiceArea === expectedTwiceArea;
}

function outTriangles(io: OracleIO, triangles: readonly [number, number, number][]): void {
    io.outInt(triangles.length);
    for (const t of triangles) {
        io.outInt(t[0]);
        io.outInt(t[1]);
        io.outInt(t[2]);
    }
}

// ---- MinimumVolumeBox3 -------------------------------------------------

function containmentViolation(box: OrientedBox3, points: readonly Vector[]): number {
    let worst = 0;
    for (const p of points) {
        const d = sub(p, box.center);
        for (let i = 0; i < 3; ++i) {
            const t = Math.abs(dot(d, box.axis[i])) - box.extent.values[i];
            worst = Math.max(worst, t);
        }
    }
    return worst;
}

function pointScale(points: readonly Vector[]): number {
    let scale = 1;
    for (const p of points) {
        for (let i = 0; i < 3; ++i) {
            scale = Math.max(scale, Math.abs(p.values[i]));
        }
    }
    return scale;
}

function readCloud(io: OracleIO): Vector[] {
    const n = io.integer();
    return readPoints(io, n, 3);
}

// The box is emitted in a form invariant under the two freedoms the
// algorithm leaves unspecified: the sign of each axis and the order of the
// three axes. Which candidate wins a volume tie depends on the order in
// which ExtractMeshTopology numbers the mesh edges and triangles, and that
// order comes out of a std::unordered_map, so it is not comparable between
// the two builds; the main cases additionally reject inputs on which
// upstream's own answer changes with the mesh order. The invariants below
// are the centre, the sorted extents, the six distinct entries of
// M = sum_i extent[i]^2 * axis[i] * axis[i]^T, and the volume.
function outBox(io: OracleIO, dimension: number, box: OrientedBox3, volume: number): void {
    io.outInt(dimension);
    io.outVec(box.center);

    const extent = [box.extent.values[0], box.extent.values[1], box.extent.values[2]];
    extent.sort((a, b) => a - b);
    for (let i = 0; i < 3; ++i) {
        io.outReal(extent[i]);
    }

    const m = [0, 0, 0, 0, 0, 0];
    for (let i = 0; i < 3; ++i) {
        const w = box.extent.values[i] * box.extent.values[i];
        const a = box.axis[i].values;
        m[0] += w * a[0] * a[0];
        m[1] += w * a[0] * a[1];
        m[2] += w * a[0] * a[2];
        m[3] += w * a[1] * a[1];
        m[4] += w * a[1] * a[2];
        m[5] += w * a[2] * a[2];
    }
    for (let i = 0; i < 6; ++i) {
        io.outReal(m[i]);
    }

    io.outReal(volume);
}

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

    family.case('NearestNeighborQuery.build', (io) => {
        const n = io.integer();
        const maxLeafSize = io.integer();
        const maxLevel = io.integer();
        const { sites } = readSites(io, n);
        const query = new NearestNeighborQuery(sites, maxLeafSize, maxLevel);
        outTree(io, query, n);
    }, { exact: true });

    family.case('NearestNeighborQuery.build.tiedCoordinates', (io) => {
        const n = io.integer();
        const maxLeafSize = io.integer();
        const maxLevel = io.integer();
        const { sites } = readSites(io, n);
        const query = new NearestNeighborQuery(sites, maxLeafSize, maxLevel);
        // Tied split coordinates: std::nth_element's placement of elements
        // equal to the median is unspecified, so the leaf memberships and
        // every split value below the root differ between the two builds.
        // Only what the site counts determine is compared here.
        io.outInt(query.getDepth());
        io.outInt(query.getLargestNodeSize());
        io.outInt(query.getNumNodes());
        io.outReal(query.getNodes()[0].split);

        for (const node of query.getNodes()) {
            io.outInt(node.axis);
            io.outInt(node.numSites);
            io.outInt(node.siteOffset);
            io.outInt(node.left);
            io.outInt(node.right);
        }
    }, { exact: true });

    family.case('NearestNeighborQuery.build.maxLevelAssert', (io) => {
        const n = io.integer();
        const maxLeafSize = io.integer();
        const maxLevel = io.integer();
        const { sites } = readSites(io, n);
        const query = new NearestNeighborQuery(sites, maxLeafSize, maxLevel);
        outTree(io, query, n);
    }, { exact: true });

    family.case('NearestNeighborQuery.findNeighbors.all', (io) => {
        const n = io.integer();
        const maxLeafSize = io.integer();
        const maxLevel = io.integer();
        const { points, sites } = readSites(io, n);
        const query = new NearestNeighborQuery(sites, maxLeafSize, maxLevel);
        outTree(io, query, n);

        for (let q = 0; q < 4; ++q) {
            const point = io.vec(3);
            const radius = io.real();
            const neighbors = query.findNeighbors(point, radius, 12);
            io.outInt(neighbors.length);
            for (const k of neighbors) {
                io.outInt(k);
            }

            let expected = 0;
            for (let i = 0; i < n; ++i) {
                const diff = sub(points[i], point);
                if (dot(diff, diff) <= radius * radius) {
                    ++expected;
                }
            }
            io.outBool(neighbors.length === expected);
        }
    }, { exact: true });

    family.case('NearestNeighborQuery.findNeighbors.limited', (io) => {
        const n = io.integer();
        const maxLeafSize = io.integer();
        const maxLevel = io.integer();
        const maxNeighbors = io.integer();
        const { points, sites } = readSites(io, n);
        const query = new NearestNeighborQuery(sites, maxLeafSize, maxLevel);

        for (let q = 0; q < 3; ++q) {
            const usedPoint = io.vec(3);
            const usedRadius = io.real();
            const neighbors = query.findNeighbors(usedPoint, usedRadius, maxNeighbors);
            io.outInt(neighbors.length);
            for (const k of neighbors) {
                io.outInt(k);
            }

            const all: { d2: number, i: number }[] = [];
            for (let i = 0; i < n; ++i) {
                const diff = sub(points[i], usedPoint);
                const d2 = dot(diff, diff);
                if (d2 <= usedRadius * usedRadius) {
                    all.push({ d2, i });
                }
            }
            all.sort((a, b) => (a.d2 !== b.d2 ? a.d2 - b.d2 : a.i - b.i));
            const expected = Math.min(all.length, maxNeighbors);
            let referenceOk = neighbors.length === expected;
            if (referenceOk) {
                const got = neighbors.slice().sort((a, b) => a - b);
                const want = all.slice(0, expected).map((e) => e.i).sort((a, b) => a - b);
                referenceOk = got.every((v, i) => v === want[i]);
            }
            io.outBool(referenceOk);
        }
    }, { exact: true });

    family.case('BoxManager.initializeAndUpdate', (io) => {
        const n = io.integer();
        const boxes: AlignedBox[] = [];
        for (let i = 0; i < n; ++i) {
            boxes.push(readBox(io));
        }

        const manager = new BoxManager(boxes);
        outOverlap(io, manager, boxes);

        const numRounds = io.integer();
        for (let round = 0; round < numRounds; ++round) {
            const numMoves = io.integer();
            for (let m = 0; m < numMoves; ++m) {
                const index = io.integer();
                const box = readBox(io);
                manager.setBox(index, box);
            }
            manager.update();
            outOverlap(io, manager, boxes);

            const probe = io.integer();
            const got = manager.getBox(probe);
            io.outVec(got.min);
            io.outVec(got.max);
        }

        manager.initialize();
        outOverlap(io, manager, boxes);
    }, { exact: true });

    family.case('InscribedFixedAspectRectInQuad.execute', (io) => {
        // Upstream's std::atan2 call only picks the quadrant index of each
        // inner edge normal; the generator keeps every normal away from the
        // quadrant boundaries, so the index is decided by the geometry and
        // the rest of the solve is arithmetic. Records on which upstream's
        // alpha assertion or its untoleranced interval test fires (finding
        // #395) are kept and compared for throw parity: the port preserves
        // both failures.
        io.integer();
        const quad = readPoints(io, 4, 2);
        const aspectRatio = io.real();

        const r = InscribedFixedAspectRectInQuad.execute(quad, aspectRatio);
        io.outBool(r.isUnique);
        io.outVec(r.rectOrigin);
        io.outReal(r.rectWidth);
        io.outReal(r.rectHeight);

        let referenceOk = r.rectWidth >= 0;
        const rect: Vector[] = [
            r.rectOrigin,
            Vector.fromArray([r.rectOrigin.values[0] + r.rectWidth, r.rectOrigin.values[1]]),
            Vector.fromArray([r.rectOrigin.values[0] + r.rectWidth,
                r.rectOrigin.values[1] + r.rectHeight]),
            Vector.fromArray([r.rectOrigin.values[0], r.rectOrigin.values[1] + r.rectHeight])
        ];
        let scale = 0;
        for (let i = 0; i < 4; ++i) {
            scale = Math.max(scale, Math.abs(quad[i].values[0]));
            scale = Math.max(scale, Math.abs(quad[i].values[1]));
        }
        for (let i = 0; i < 4; ++i) {
            const e = sub(quad[(i + 1) % 4], quad[i]);
            for (let k = 0; k < 4; ++k) {
                const d = sub(rect[k], quad[i]);
                referenceOk = referenceOk
                    && (e.values[0] * d.values[1] - e.values[1] * d.values[0]
                        >= -1e-9 * scale * scale);
            }
        }
        io.outBool(referenceOk);
    }, { exact: true });

    family.case('TriangulateEC.triangulate', (io) => {
        io.integer();
        const n = io.integer();
        const points = readPoints(io, n, 2);

        const triangulator = new TriangulateEC(points);
        triangulator.triangulate();
        outTriangles(io, triangulator.getTriangles());

        const polygon: number[] = [];
        for (let i = 0; i < points.length; ++i) {
            polygon.push(i);
        }
        const referenceOk = triangulator.getTriangles().length === points.length - 2
            && triangulationIsValid(points, triangulator.getTriangles(),
                shoelaceArea(points, polygon));
        io.outBool(referenceOk);
    }, { exact: true });

    family.case('TriangulateEC.triangulatePolygon', (io) => {
        io.integer();
        const n = io.integer();
        const points = readPoints(io, n, 2);
        const polygon = readPolygon(io);

        const triangulator = new TriangulateEC(points);
        triangulator.triangulatePolygon(polygon);
        outTriangles(io, triangulator.getTriangles());

        const referenceOk = triangulator.getTriangles().length === polygon.length - 2
            && triangulationIsValid(points, triangulator.getTriangles(),
                shoelaceArea(points, polygon));
        io.outBool(referenceOk);
    }, { exact: true });

    family.case('TriangulateEC.triangulateWithHole', (io) => {
        const n = io.integer();
        const points = readPoints(io, n, 2);
        const outer = readPolygon(io);
        const inner = readPolygon(io);

        const triangulator = new TriangulateEC(points);
        triangulator.triangulateWithHole(outer, inner);
        outTriangles(io, triangulator.getTriangles());

        const expected = shoelaceArea(points, outer) + shoelaceArea(points, inner);
        io.outBool(triangulationIsValid(points, triangulator.getTriangles(), expected));
    }, { exact: true });

    family.case('TriangulateEC.triangulateWithHoles', (io) => {
        const n = io.integer();
        const points = readPoints(io, n, 2);
        const outer = readPolygon(io);
        const inner0 = readPolygon(io);
        const inner1 = readPolygon(io);

        const triangulator = new TriangulateEC(points);
        triangulator.triangulateWithHoles(outer, [inner0, inner1]);
        outTriangles(io, triangulator.getTriangles());

        const expected = shoelaceArea(points, outer) + shoelaceArea(points, inner0)
            + shoelaceArea(points, inner1);
        io.outBool(triangulationIsValid(points, triangulator.getTriangles(), expected));
    }, { exact: true });

    family.case('TriangulateEC.triangulateTree', (io) => {
        const n = io.integer();
        const points = readPoints(io, n, 2);
        const outer = readPolygon(io);
        const hole = readPolygon(io);
        const innerOuter = readPolygon(io);

        const root = new PolygonTree();
        root.polygon = outer;
        const holeNode = new PolygonTree();
        holeNode.polygon = hole;
        const innerNode = new PolygonTree();
        innerNode.polygon = innerOuter;
        holeNode.child.push(innerNode);
        root.child.push(holeNode);

        const triangulator = new TriangulateEC(points);
        triangulator.triangulateTree(root);
        outTriangles(io, triangulator.getTriangles());

        const expected = shoelaceArea(points, outer) + shoelaceArea(points, hole)
            + shoelaceArea(points, innerOuter);
        io.outBool(triangulationIsValid(points, triangulator.getTriangles(), expected));
    }, { exact: true });

    family.case('MinimumVolumeBox3FloatingPoint.compute', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const points = readCloud(io);

        const query = new MinimumVolumeBox3FloatingPoint(0);
        const r = query.compute(points, lgMaxSample);

        // Only the hull dimension, the minimum volume and the reference
        // checks are comparable; see the C++ case comment.
        io.outInt(r.dimension);
        io.outReal(r.volume);

        const scale = pointScale(points);
        io.outBool(containmentViolation(r.box, points) <= 1e-9 * scale);
        const product = 8 * r.box.extent.values[0] * r.box.extent.values[1]
            * r.box.extent.values[2];
        io.outBool(Math.abs(product - r.volume) <= 1e-9 * Math.max(1, Math.abs(r.volume)));
        // Tolerance 1e-14, measured maximum scaled error 1.8e-16 (one record
        // in twenty at the committed size). Cause: the support vertex picked
        // among vertices whose double projections are equal depends on the
        // path the GetExtreme hill climb takes, and that path follows the
        // vertex adjacency order, which ExtractVertexAdjacencies builds by
        // iterating ETManifoldMesh's std::unordered_map. Two vertices with
        // equal double projections can have different exact projections, so
        // the exact rational GetMinimumVolumeBox then places the centre an
        // ulp apart. Everything else in this case is bit-identical, and the
        // sibling computeHull case, which fixes the mesh, is exact.
    }, { tol: 1e-14 });

    family.case('MinimumVolumeBox3FloatingPoint.compute.lowDimension', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const points = readCloud(io);

        const query = new MinimumVolumeBox3FloatingPoint(0);
        const r = query.compute(points, lgMaxSample);
        outBox(io, r.dimension, r.box, r.volume);
    }, { exact: true });

    family.case('MinimumVolumeBox3FloatingPoint.computeHull', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        const mesh = readMesh(io);

        // 'usable' is recorded by the C++ side: it is true when upstream's
        // answer does not change with the mesh order and its box contains
        // the polytope. Only then is the box comparable.
        const usable = io.integer();

        const query = new MinimumVolumeBox3FloatingPoint(0);
        const r = query.computeHull(mesh.vertices, mesh.indices, lgMaxSample);
        if (usable !== 0) {
            const scale = pointScale(mesh.vertices);
            outBox(io, 3, r.box, r.volume);
            io.outBool(containmentViolation(r.box, mesh.vertices) <= 1e-9 * scale);
        }
    }, { exact: true });

    family.case('MinimumVolumeBox3FloatingPoint.compute.coplanar', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const points = readCloud(io);

        const query = new MinimumVolumeBox3FloatingPoint(0);
        const r = query.compute(points, lgMaxSample);
        outBox(io, r.dimension, r.box, r.volume);
    }, { deviation: '#352 (dimension-2 Newell normal loop drops the wrap-around term)' });

    family.case('MinimumVolumeBox3FloatingPoint.compute.nonContaining', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const points = readCloud(io);

        const query = new MinimumVolumeBox3FloatingPoint(0);
        const r = query.compute(points, lgMaxSample);
        outBox(io, r.dimension, r.box, r.volume);
        io.outBool(containmentViolation(r.box, points) <= 1e-9 * pointScale(points));
    }, { deviation: '#405 (ComputeVolume axis minima) and #426 (GetExtreme plateau)' });

    family.case('MinimumVolumeBox3Rational.compute', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const points = readCloud(io);

        const query = new MinimumVolumeBox3Rational(0);
        const r = query.compute(points, lgMaxSample);

        // Only the hull dimension, the minimum volume and containment are
        // comparable; see the C++ case comment.
        io.outInt(r.dimension);
        io.outReal(r.volume);
        const scale = pointScale(points);
        io.outBool(containmentViolation(r.box, points) <= 1e-9 * scale);
        // Tolerance 1e-14, measured maximum scaled error 6.3e-16 (one record
        // in twenty at the committed size), in the extents and the axis
        // matrix only; the centre and the volume are bit-identical. Cause:
        // the rational pipeline never normalizes a candidate axis, so the
        // same geometric box reached through a different edge pair is
        // represented by a differently scaled axis triple, and the final
        // conversion extent = rScaledExtent / sqrt(rSqrLengthAxis) then
        // rounds differently. Which edge pair wins an exact volume tie
        // depends on the std::unordered_map enumeration order, which is not
        // comparable; the generator already rejects clouds whose answer
        // changes under three alternative input orders.
    }, { tol: 1e-14, timeout: 120000 });

    family.case('MinimumVolumeBox3Rational.compute.lowDimension', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const points = readCloud(io);

        const query = new MinimumVolumeBox3Rational(0);
        const r = query.compute(points, lgMaxSample);
        outBox(io, r.dimension, r.box, r.volume);
    }, { exact: true, timeout: 120000 });

    family.case('MinimumVolumeBox3Rational.computeHull', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        const mesh = readMesh(io);
        const usable = io.integer();

        const query = new MinimumVolumeBox3Rational(0);
        const r = query.computeHull(mesh.vertices, mesh.indices, lgMaxSample);
        if (usable !== 0) {
            const scale = pointScale(mesh.vertices);
            outBox(io, 3, r.box, r.volume);
            io.outBool(containmentViolation(r.box, mesh.vertices) <= 1e-9 * scale);
        }
    }, { exact: true, timeout: 120000 });

    family.case('MinimumVolumeBox3Rational.compute.coplanar', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const points = readCloud(io);

        const query = new MinimumVolumeBox3Rational(0);
        const r = query.compute(points, lgMaxSample);
        outBox(io, r.dimension, r.box, r.volume);
    }, {
        deviation: '#355 (dimension-2 Newell normal loop drops the wrap-around term)',
        timeout: 120000
    });

    family.finish();
});
