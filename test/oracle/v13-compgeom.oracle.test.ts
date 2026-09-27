// Replays oracle/cpp/cases/v13-compgeom.cpp (SeparatePoints3.h,
// TriangulateCDT.h). Keep the two files in the same order.
//
// Canonicalization: PolygonTreeEx.allTriangles and outsideTriangles follow
// std::unordered_map order upstream and sorted-key order in the port, so both
// lists are sorted by their stored vertex tuple on both sides. Everything
// else is emitted in the order the algorithm produces it. See the C++ file's
// header comment.
//
// Each main case ends with an independent reference check that each side
// computes on its own output (exact integer arithmetic on the common dyadic
// scale of the inputs, test/helpers/exact.ts). The checks are also counted
// here, and the last test requires that none of them failed, so a record on
// which both implementations agree on a wrong answer is still reported.
import { describe, expect, it } from 'vitest';
import { ConvexHull3 } from '../../src/ConvexHull3.js';
import { PolygonTree, PolygonTreeEx } from '../../src/PolygonTree.js';
import { SeparatePoints3 } from '../../src/SeparatePoints3.js';
import type { Plane3 } from '../../src/Hyperplane.js';
import { TriangulateCDT } from '../../src/TriangulateCDT.js';
import { Vector } from '../../src/Vector.js';
import { exactDyadic, orient2, orient3 } from '../helpers/exact.js';
import { OracleFamily, type OracleIO } from './harness.js';

// Reference checks that failed on the port's side (case, record).
const referenceFailures: string[] = [];

function points3(io: OracleIO, n: number): Vector[] {
    const pts: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        pts.push(io.vec(3));
    }
    return pts;
}

// ---- exact 3D reference (SeparatePoints3) ---------------------------------

type P3 = readonly bigint[];

function sgn(x: bigint): number {
    return x > 0n ? 1 : (x < 0n ? -1 : 0);
}

function orient2Drop(a: P3, b: P3, c: P3, drop: number): number {
    const i0 = (drop + 1) % 3, i1 = (drop + 2) % 3;
    return sgn((b[i0] - a[i0]) * (c[i1] - a[i1]) - (b[i1] - a[i1]) * (c[i0] - a[i0]));
}

function onSegmentDrop(p: P3, q: P3, x: P3, drop: number): boolean {
    for (let k = 0; k < 3; ++k) {
        if (k === drop) {
            continue;
        }
        const lo = p[k] < q[k] ? p[k] : q[k];
        const hi = p[k] < q[k] ? q[k] : p[k];
        if (x[k] < lo || hi < x[k]) {
            return false;
        }
    }
    return true;
}

function segmentsMeetDrop(p: P3, q: P3, a: P3, b: P3, drop: number): boolean {
    const d0 = orient2Drop(p, q, a, drop);
    const d1 = orient2Drop(p, q, b, drop);
    const d2 = orient2Drop(a, b, p, drop);
    const d3 = orient2Drop(a, b, q, drop);
    if (d0 * d1 < 0 && d2 * d3 < 0) {
        return true;
    }
    return (d0 === 0 && onSegmentDrop(p, q, a, drop))
        || (d1 === 0 && onSegmentDrop(p, q, b, drop))
        || (d2 === 0 && onSegmentDrop(a, b, p, drop))
        || (d3 === 0 && onSegmentDrop(a, b, q, drop));
}

// Closed segment <p,q> against the closed nondegenerate triangle <a,b,c>.
function segmentMeetsTriangle(p: P3, q: P3, a: P3, b: P3, c: P3): boolean {
    const sP = orient3(a, b, c, p);
    const sQ = orient3(a, b, c, q);
    if (sP * sQ > 0) {
        return false;
    }
    if (sP === 0 && sQ === 0) {
        const u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        const v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        const nx = u[1] * v[2] - u[2] * v[1];
        const ny = u[2] * v[0] - u[0] * v[2];
        const drop = nx !== 0n ? 0 : (ny !== 0n ? 1 : 2);
        for (const end of [p, q]) {
            const o0 = orient2Drop(a, b, end, drop);
            const o1 = orient2Drop(b, c, end, drop);
            const o2 = orient2Drop(c, a, end, drop);
            if ((o0 >= 0 && o1 >= 0 && o2 >= 0) || (o0 <= 0 && o1 <= 0 && o2 <= 0)) {
                return true;
            }
        }
        return segmentsMeetDrop(p, q, a, b, drop) || segmentsMeetDrop(p, q, b, c, drop)
            || segmentsMeetDrop(p, q, c, a, drop);
    }
    const v0 = orient3(p, q, a, b);
    const v1 = orient3(p, q, b, c);
    const v2 = orient3(p, q, c, a);
    return (v0 >= 0 && v1 >= 0 && v2 >= 0) || (v0 <= 0 && v1 <= 0 && v2 <= 0);
}

function pointInClosedHull(q: P3, pts: readonly P3[], hull: readonly number[]): boolean {
    for (let t = 0; t + 2 < hull.length; t += 3) {
        if (orient3(pts[hull[t]], pts[hull[t + 1]], pts[hull[t + 2]], q) > 0) {
            return false;
        }
    }
    return true;
}

// The closed convex hulls meet: a point of one set in the other hull, or a
// hull edge of one meeting a hull triangle of the other (complete for convex
// polytopes).
function closedHullsMeet(pts0: readonly P3[], hull0: readonly number[],
    pts1: readonly P3[], hull1: readonly number[]): boolean {
    if (pts0.some(q => pointInClosedHull(q, pts1, hull1))
        || pts1.some(q => pointInClosedHull(q, pts0, hull0))) {
        return true;
    }
    const sets: [readonly P3[], readonly number[], readonly P3[], readonly number[]][] = [
        [pts0, hull0, pts1, hull1], [pts1, hull1, pts0, hull0]];
    for (const [pa, ha, pb, hb] of sets) {
        for (let t = 0; t + 2 < ha.length; t += 3) {
            for (let j = 0; j < 3; ++j) {
                const p = pa[ha[t + j]];
                const q = pa[ha[t + (j + 1) % 3]];
                for (let u = 0; u + 2 < hb.length; u += 3) {
                    if (segmentMeetsTriangle(p, q, pb[hb[u]], pb[hb[u + 1]], pb[hb[u + 2]])) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

function signedDistance(plane: Plane3, p: Vector): number {
    const n = plane.normal.values, x = p.values;
    return n[0] * x[0] + n[1] * x[1] + n[2] * x[2] - plane.constant;
}

function separationResidual(pts0: readonly Vector[], pts1: readonly Vector[],
    plane: Plane3): number {
    let maxNeg0 = 0, maxPos0 = 0, maxNeg1 = 0, maxPos1 = 0;
    for (const p of pts0) {
        const s = signedDistance(plane, p);
        maxPos0 = Math.max(maxPos0, s);
        maxNeg0 = Math.max(maxNeg0, -s);
    }
    for (const p of pts1) {
        const s = signedDistance(plane, p);
        maxPos1 = Math.max(maxPos1, s);
        maxNeg1 = Math.max(maxNeg1, -s);
    }
    return Math.min(Math.max(maxPos0, maxNeg1), Math.max(maxNeg0, maxPos1));
}

function maxAbsCoordinate(pts0: readonly Vector[], pts1: readonly Vector[]): number {
    let scale = 1;
    for (const p of [...pts0, ...pts1]) {
        for (const x of p.values) {
            scale = Math.max(scale, Math.abs(x));
        }
    }
    return scale;
}

function separationReference(pts0: readonly Vector[], pts1: readonly Vector[],
    separated: boolean, plane: Plane3): boolean {
    if (separated) {
        const tol = 1e-12 * maxAbsCoordinate(pts0, pts1);
        const n = plane.normal.values;
        const unit = n[0] * n[0] + n[1] * n[1] + n[2] * n[2] - 1;
        return Math.abs(unit) <= 1e-12 && separationResidual(pts0, pts1, plane) <= tol;
    }
    const ch0 = new ConvexHull3();
    ch0.compute(pts0);
    const ch1 = new ConvexHull3();
    ch1.compute(pts1);
    if (ch0.getDimension() !== 3 || ch1.getDimension() !== 3) {
        return true;
    }
    const flat = exactDyadic([...pts0, ...pts1].flatMap(p => [...p.values]));
    const exact: P3[] = [];
    for (let i = 0; i < flat.length; i += 3) {
        exact.push(flat.slice(i, i + 3));
    }
    return closedHullsMeet(exact.slice(0, pts0.length), ch0.getHull(),
        exact.slice(pts0.length), ch1.getHull());
}

// ---- TriangulateCDT input and output --------------------------------------

function readTree(io: OracleIO): PolygonTree {
    const node = new PolygonTree();
    const numIndices = io.integer();
    for (let i = 0; i < numIndices; ++i) {
        node.polygon.push(io.integer());
    }
    const numChildren = io.integer();
    for (let c = 0; c < numChildren; ++c) {
        node.child.push(readTree(io));
    }
    return node;
}

function readCdtInput(io: OracleIO): { pool: Vector[]; root: PolygonTree } {
    const n = io.integer();
    const pool: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        pool.push(io.vec(2));
    }
    return { pool, root: readTree(io) };
}

type Tri = [number, number, number];

function emitTriangles(io: OracleIO, tris: readonly Tri[]): void {
    io.outInt(tris.length);
    for (const t of tris) {
        io.outInt(t[0]);
        io.outInt(t[1]);
        io.outInt(t[2]);
    }
}

function emitNodeIndices(io: OracleIO, indices: readonly number[]): void {
    io.outInt(indices.length);
    for (const i of indices) {
        io.outInt(i);
    }
}

function compareTri(a: Tri, b: Tri): number {
    return a[0] !== b[0] ? a[0] - b[0] : (a[1] !== b[1] ? a[1] - b[1] : a[2] - b[2]);
}

function emitSortedTriangles(io: OracleIO, tris: readonly Tri[]): void {
    emitTriangles(io, [...tris].sort(compareTri));
}

function emitTreeEx(io: OracleIO, tree: PolygonTreeEx): void {
    io.outInt(tree.nodes.length);
    for (const node of tree.nodes) {
        io.outInt(node.self);
        io.outInt(node.chirality);
        io.outInt(node.parent === PolygonTreeEx.INVALID ? -1 : node.parent);
        io.outInt(node.minChild);
        io.outInt(node.supChild);
        io.outInt(node.polygon.length);
        for (const v of node.polygon) {
            io.outInt(v);
        }
        emitTriangles(io, node.triangulation);
    }
    emitTriangles(io, tree.interiorTriangles);
    emitNodeIndices(io, tree.interiorNodeIndices);
    emitTriangles(io, tree.exteriorTriangles);
    emitNodeIndices(io, tree.exteriorNodeIndices);
    emitTriangles(io, tree.insideTriangles);
    emitNodeIndices(io, tree.insideNodeIndices);
    emitSortedTriangles(io, tree.outsideTriangles);
    emitSortedTriangles(io, tree.allTriangles);
}

// ---- exact 2D reference (TriangulateCDT) ----------------------------------

type P2 = readonly bigint[];

function twiceArea(pool: readonly P2[], polygon: readonly number[]): bigint {
    let sum = 0n;
    const n = polygon.length;
    for (let i = 0, j = n - 1; i < n; j = i++) {
        const a = pool[polygon[j]], b = pool[polygon[i]];
        sum += a[0] * b[1] - b[0] * a[1];
    }
    return sum;
}

// Winding number about q = a + b + c (three times a triangle centroid) of
// the polygon scaled by 3.
function winding3(pool: readonly P2[], polygon: readonly number[], q: P2): number {
    let winding = 0;
    const n = polygon.length;
    for (let i = 0, j = n - 1; i < n; j = i++) {
        const pa = pool[polygon[j]], pb = pool[polygon[i]];
        const ax = 3n * pa[0], ay = 3n * pa[1], bx = 3n * pb[0], by = 3n * pb[1];
        if (ay <= q[1]) {
            if (q[1] < by && orient2(ax, ay, bx, by, q[0], q[1]) > 0) {
                ++winding;
            }
        }
        else if (by <= q[1] && orient2(ax, ay, bx, by, q[0], q[1]) < 0) {
            --winding;
        }
    }
    return winding;
}

function canonicalKey(t: readonly number[]): string {
    const k = t[0] <= t[1] && t[0] <= t[2] ? 0 : (t[1] <= t[2] ? 1 : 2);
    return `${t[k]},${t[(k + 1) % 3]},${t[(k + 2) % 3]}`;
}

// The same four checks as TreeReference in the C++ file.
function treeReference(poolD: readonly Vector[], tree: PolygonTreeEx): boolean {
    const flat = exactDyadic(poolD.flatMap(p => [...p.values]));
    const pool: P2[] = [];
    for (let i = 0; i < flat.length; i += 2) {
        pool.push(flat.slice(i, i + 2));
    }
    const sum3 = (t: readonly number[]): P2 => [
        pool[t[0]][0] + pool[t[1]][0] + pool[t[2]][0],
        pool[t[0]][1] + pool[t[1]][1] + pool[t[2]][1]];

    const allEdges = new Set<string>();
    const all = new Map<string, number>();
    for (const t of tree.allTriangles) {
        const key = canonicalKey(t);
        all.set(key, (all.get(key) ?? 0) + 1);
        for (let j = 0; j < 3; ++j) {
            const u = t[j], v = t[(j + 1) % 3];
            allEdges.add(`${Math.min(u, v)},${Math.max(u, v)}`);
        }
    }

    for (const node of tree.nodes) {
        let expected = twiceArea(pool, node.polygon);
        for (let c = node.minChild; c < node.supChild; ++c) {
            expected += twiceArea(pool, tree.nodes[c].polygon);
        }
        let total = 0n;
        for (const t of node.triangulation) {
            const a = pool[t[0]], b = pool[t[1]], c = pool[t[2]];
            if (orient2(a[0], a[1], b[0], b[1], c[0], c[1]) !== node.chirality) {
                return false;
            }
            total += (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            const q = sum3(t);
            if (winding3(pool, node.polygon, q) !== node.chirality) {
                return false;
            }
            for (let ci = node.minChild; ci < node.supChild; ++ci) {
                if (winding3(pool, tree.nodes[ci].polygon, q) !== 0) {
                    return false;
                }
            }
        }
        if (total !== expected) {
            return false;
        }
        const n = node.polygon.length;
        for (let i = 0, j = n - 1; i < n; j = i++) {
            const u = node.polygon[j], v = node.polygon[i];
            if (!allEdges.has(`${Math.min(u, v)},${Math.max(u, v)}`)) {
                return false;
            }
        }
    }

    const parts = new Map<string, number>();
    const addPart = (t: readonly number[]): void => {
        const key = canonicalKey(t);
        parts.set(key, (parts.get(key) ?? 0) + 1);
    };
    for (const t of tree.outsideTriangles) {
        if (winding3(pool, tree.nodes[0].polygon, sum3(t)) !== 0) {
            return false;
        }
        addPart(t);
    }

    let ie = 0, ee = 0;
    for (let k = 0; k < tree.insideTriangles.length; ++k) {
        const nIndex = tree.insideNodeIndices[k];
        const t = tree.insideTriangles[k];
        const same = (s: Tri | undefined): boolean =>
            s !== undefined && s[0] === t[0] && s[1] === t[1] && s[2] === t[2];
        if (tree.nodes[nIndex].chirality > 0) {
            if (!same(tree.interiorTriangles[ie]) || tree.interiorNodeIndices[ie] !== nIndex) {
                return false;
            }
            ++ie;
            addPart(t);
        }
        else {
            if (!same(tree.exteriorTriangles[ee]) || tree.exteriorNodeIndices[ee] !== nIndex) {
                return false;
            }
            ++ee;
            addPart([t[0], t[2], t[1]]);
        }
    }
    if (ie !== tree.interiorTriangles.length || ee !== tree.exteriorTriangles.length) {
        return false;
    }
    if (parts.size !== all.size) {
        return false;
    }
    for (const [key, count] of all) {
        if (parts.get(key) !== count) {
            return false;
        }
    }
    return true;
}

// Stages recorded by the C++ trace (TraceSeparate3).
const STAGE_FACE0 = 1;
const STAGE_FACE1 = 2;

describe('oracle: v13-compgeom', () => {
    const family = new OracleFamily('v13-compgeom');

    // ---- SeparatePoints3 ----------------------------------------------------

    family.case('SeparatePoints3.compute', (io) => {
        const n0 = io.integer();
        const n1 = io.integer();
        const pts0 = points3(io, n0);
        const pts1 = points3(io, n1);
        const stage = io.integer();
        const determined = io.boolean();
        const r = new SeparatePoints3().compute(pts0, pts1);
        io.outBool(r.separated);
        if (r.separated && determined) {
            io.outVec(r.separatingPlane.normal);
            io.outReal(r.separatingPlane.constant);
            if (stage === STAGE_FACE0 || stage === STAGE_FACE1) {
                io.outVec(r.separatingPlane.origin);
            }
        }
        const ok = separationReference(pts0, pts1, r.separated, r.separatingPlane);
        if (!ok) {
            referenceFailures.push(`SeparatePoints3.compute record ${io.index}`);
        }
        io.outBool(ok);
    }, { exact: true });

    // #348 finding 1 (fixed): upstream reports a separation for overlapping
    // clouds because a vertex of its own candidate face classifies as
    // strictly positive by round-off.
    family.case('SeparatePoints3.deviation.roundoff', (io) => {
        const n0 = io.integer();
        const n1 = io.integer();
        const pts0 = points3(io, n0);
        const pts1 = points3(io, n1);
        const r = new SeparatePoints3().compute(pts0, pts1);
        io.outBool(r.separated);
    }, { exact: true, deviation: '#348 SeparatePoints3 finding 1 (exact side predicate)' });

    // #348 finding 2 (fixed): upstream's edge-edge plane keeps a stale origin.
    family.case('SeparatePoints3.deviation.edgeAxisOrigin', (io) => {
        const n0 = io.integer();
        const n1 = io.integer();
        const pts0 = points3(io, n0);
        const pts1 = points3(io, n1);
        const r = new SeparatePoints3().compute(pts0, pts1);
        io.outBool(r.separated);
        io.outVec(r.separatingPlane.normal);
        io.outReal(r.separatingPlane.constant);
        io.outVec(r.separatingPlane.origin);
    }, { exact: true, deviation: '#348 SeparatePoints3 finding 2 (consistent plane origin)' });

    // ---- TriangulateCDT -----------------------------------------------------

    family.case('TriangulateCDT.compute', (io) => {
        const { pool, root } = readCdtInput(io);
        const tree = new TriangulateCDT().compute(pool, root);
        emitTreeEx(io, tree);
        const ok = treeReference(pool, tree);
        if (!ok) {
            referenceFailures.push(`TriangulateCDT.compute record ${io.index}`);
        }
        io.outBool(ok);
    }, { exact: true });

    family.case('TriangulateCDT.compute.invalidThrows', (io) => {
        const { pool, root } = readCdtInput(io);
        const tree = new TriangulateCDT().compute(pool, root);
        emitTreeEx(io, tree);
    }, { exact: true });

    // Upstream throws "Unexpected condition." when a hole shares an edge
    // with the outer polygon on the convex hull of the referenced points;
    // the port preserves it (new suspect, oracle/reports/v13-compgeom.md).
    family.case('TriangulateCDT.compute.hullSharedEdgeThrows', (io) => {
        const { pool, root } = readCdtInput(io);
        const tree = new TriangulateCDT().compute(pool, root);
        emitTreeEx(io, tree);
    }, { exact: true });

    family.finish();

    it('every independent reference check holds on the port output', () => {
        expect(referenceFailures).toEqual([]);
    });
});
