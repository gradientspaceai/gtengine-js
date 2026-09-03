import { describe, it, expect } from 'vitest';
import { Delaunay3, Delaunay3SearchInfo } from '../src/Delaunay3.js';
import { Vector } from '../src/Vector.js';
import { SWInterval } from '../src/SWInterval.js';

const v3 = (x: number, y: number, z: number): Vector => Vector.fromArray([x, y, z]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The signed volume determinant Dot(V1-V0, Cross(V2-V0, V3-V0)). It is
// positive for the GTE ordered-tetrahedron convention.
function orient3d(a: Vector, b: Vector, c: Vector, d: Vector): number {
    const p0 = a.values, p1 = b.values, p2 = c.values, p3 = d.values;
    const ux = p1[0] - p0[0], uy = p1[1] - p0[1], uz = p1[2] - p0[2];
    const vx = p2[0] - p0[0], vy = p2[1] - p0[1], vz = p2[2] - p0[2];
    const wx = p3[0] - p0[0], wy = p3[1] - p0[1], wz = p3[2] - p0[2];
    return ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx)
        + uz * (vx * wy - vy * wx);
}

// The circumsphere of a nondegenerate tetrahedron, solved as a 3x3 linear
// system for the center. Returns null if the tetrahedron is (numerically)
// degenerate.
function circumsphere(a: Vector, b: Vector, c: Vector, d: Vector):
    { center: number[]; r2: number } | null {
    const p0 = a.values;
    const rows: number[][] = [];
    for (const q of [b, c, d]) {
        const p = q.values;
        const ex = p[0] - p0[0], ey = p[1] - p0[1], ez = p[2] - p0[2];
        rows.push([ex, ey, ez, 0.5 * (ex * ex + ey * ey + ez * ez)]);
    }
    const det3 = (m: number[][], c0: number, c1: number, c2: number): number =>
        m[0][c0] * (m[1][c1] * m[2][c2] - m[1][c2] * m[2][c1])
        - m[0][c1] * (m[1][c0] * m[2][c2] - m[1][c2] * m[2][c0])
        + m[0][c2] * (m[1][c0] * m[2][c1] - m[1][c1] * m[2][c0]);
    const det = det3(rows, 0, 1, 2);
    if (det === 0) {
        return null;
    }
    // Cramer's rule with the right-hand side in column j.
    const withRhs = (j: number): number[][] => rows.map(r => {
        const s = r.slice(0, 3);
        s[j] = r[3];
        return s;
    });
    const cx = det3(withRhs(0), 0, 1, 2) / det;
    const cy = det3(withRhs(1), 0, 1, 2) / det;
    const cz = det3(withRhs(2), 0, 1, 2) / det;
    const center = [p0[0] + cx, p0[1] + cy, p0[2] + cz];
    const r2 = cx * cx + cy * cy + cz * cz;
    return { center, r2 };
}

function dist2(center: number[], p: Vector): number {
    const q = p.values;
    const dx = q[0] - center[0], dy = q[1] - center[1], dz = q[2] - center[2];
    return dx * dx + dy * dy + dz * dz;
}

// Structural checks applied to every successful tetrahedralization.
function checkTetrahedralization(del: Delaunay3, points: Vector[]): void {
    const indices = del.getIndices();
    const adjacencies = del.getAdjacencies();
    const numTetrahedra = del.getNumTetrahedra();
    expect(indices.length).toBe(4 * numTetrahedra);
    expect(adjacencies.length).toBe(4 * numTetrahedra);
    expect(numTetrahedra).toBeGreaterThan(0);

    // Every tetrahedron is positively oriented and nondegenerate.
    for (let t = 0; t < numTetrahedra; ++t) {
        const i0 = indices[4 * t], i1 = indices[4 * t + 1];
        const i2 = indices[4 * t + 2], i3 = indices[4 * t + 3];
        expect(new Set([i0, i1, i2, i3]).size).toBe(4);
        expect(orient3d(points[i0], points[i1], points[i2], points[i3]))
            .toBeGreaterThan(0);
    }

    // The face opposite V[j], as listed in the file comment.
    const opposite = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

    // Each face of the mesh is shared by at most two tetrahedra, and the
    // adjacency relation is symmetric and consistent with the shared faces.
    const faceOwners = new Map<string, number>();
    for (let t = 0; t < numTetrahedra; ++t) {
        for (let j = 0; j < 4; ++j) {
            const f = opposite[j].map(k => indices[4 * t + k]);
            const key = f.slice().sort((p, q) => p - q).join(',');
            faceOwners.set(key, (faceOwners.get(key) ?? 0) + 1);
            expect(faceOwners.get(key)).toBeLessThanOrEqual(2);

            const adj = adjacencies[4 * t + j];
            if (adj === -1) {
                continue;
            }
            expect(adj).toBeGreaterThanOrEqual(0);
            expect(adj).toBeLessThan(numTetrahedra);
            expect(adj).not.toBe(t);

            // The adjacent tetrahedron shares exactly this face and points
            // back at t through that face.
            let found = false;
            for (let k = 0; k < 4; ++k) {
                const g = opposite[k].map(m => indices[4 * adj + m]);
                const gkey = g.slice().sort((p, q) => p - q).join(',');
                if (gkey === key) {
                    found = true;
                    expect(adjacencies[4 * adj + k]).toBe(t);
                }
            }
            expect(found).toBe(true);
        }
    }

    // Every face with only one owner must be a hull face (adjacency -1), and
    // vice versa.
    let numBoundary = 0;
    for (const count of faceOwners.values()) {
        if (count === 1) {
            ++numBoundary;
        }
    }
    const numMinusOne = adjacencies.filter(a => a === -1).length;
    expect(numMinusOne).toBe(numBoundary);

    // The accessors agree with the compact arrays.
    for (let t = 0; t < numTetrahedra; ++t) {
        expect(del.getTetrahedronIndices(t)).toEqual(
            [indices[4 * t], indices[4 * t + 1], indices[4 * t + 2], indices[4 * t + 3]]);
        expect(del.getTetrahedronAdjacencies(t)).toEqual(
            [adjacencies[4 * t], adjacencies[4 * t + 1], adjacencies[4 * t + 2],
                adjacencies[4 * t + 3]]);
    }
    expect(del.getTetrahedronIndices(numTetrahedra)).toBeNull();
    expect(del.getTetrahedronAdjacencies(-1)).toBeNull();
}

// The Delaunay (empty circumsphere) property: no input vertex lies strictly
// inside the circumsphere of a tetrahedron. 'scale' is the coordinate scale
// used to build the relative tolerance.
function checkEmptyCircumsphere(del: Delaunay3, points: Vector[], scale: number): void {
    const indices = del.getIndices();
    const numTetrahedra = del.getNumTetrahedra();
    for (let t = 0; t < numTetrahedra; ++t) {
        const i0 = indices[4 * t], i1 = indices[4 * t + 1];
        const i2 = indices[4 * t + 2], i3 = indices[4 * t + 3];
        const sphere = circumsphere(points[i0], points[i1], points[i2], points[i3]);
        expect(sphere).not.toBeNull();
        const s = sphere as { center: number[]; r2: number };
        const tolerance = 1e-8 * Math.max(s.r2, scale * scale);
        for (let i = 0; i < points.length; ++i) {
            if (i === i0 || i === i1 || i === i2 || i === i3) {
                continue;
            }
            expect(dist2(s.center, points[i])).toBeGreaterThan(s.r2 - tolerance);
        }
    }
}

// A subclass that counts how many predicate evaluations fall through the
// interval-arithmetic fast path into the exact rational fallback. The
// interval expression trees are replicated here so the count is measured
// independently of the class internals.
class CountingDelaunay3 extends Delaunay3 {
    exactPlane = 0;
    exactCircumsphere = 0;

    private point(index: number): readonly number[] {
        return this.mVertices[index].values;
    }

    protected override toPlane(pIndex: number, v0Index: number, v1Index: number,
        v2Index: number): number {
        if (pIndex >= 0) {
            const inP = this.point(pIndex);
            const inV0 = this.point(v0Index);
            const inV1 = this.point(v1Index);
            const inV2 = this.point(v2Index);
            const x0 = SWInterval.sub(inP[0], inV0[0]);
            const y0 = SWInterval.sub(inP[1], inV0[1]);
            const z0 = SWInterval.sub(inP[2], inV0[2]);
            const x1 = SWInterval.sub(inV1[0], inV0[0]);
            const y1 = SWInterval.sub(inV1[1], inV0[1]);
            const z1 = SWInterval.sub(inV1[2], inV0[2]);
            const x2 = SWInterval.sub(inV2[0], inV0[0]);
            const y2 = SWInterval.sub(inV2[1], inV0[1]);
            const z2 = SWInterval.sub(inV2[2], inV0[2]);
            const c0 = y1.mul(z2).sub(y2.mul(z1));
            const c1 = y2.mul(z0).sub(y0.mul(z2));
            const c2 = y0.mul(z1).sub(y1.mul(z0));
            const det = x0.mul(c0).add(x1.mul(c1)).add(x2.mul(c2));
            if (!(det.get(0) > 0) && !(det.get(1) < 0)) {
                ++this.exactPlane;
            }
        }
        return super.toPlane(pIndex, v0Index, v1Index, v2Index);
    }

    protected override toCircumsphere(pIndex: number, v0Index: number,
        v1Index: number, v2Index: number, v3Index: number): number {
        if (pIndex >= 0) {
            const inP = this.point(pIndex);
            const vs = [this.point(v0Index), this.point(v1Index),
                this.point(v2Index), this.point(v3Index)];
            const x: SWInterval[] = [], y: SWInterval[] = [];
            const z: SWInterval[] = [], w: SWInterval[] = [];
            for (let k = 0; k < 4; ++k) {
                const xk = SWInterval.sub(vs[k][0], inP[0]);
                const yk = SWInterval.sub(vs[k][1], inP[1]);
                const zk = SWInterval.sub(vs[k][2], inP[2]);
                const s0 = SWInterval.add(vs[k][0], inP[0]);
                const s1 = SWInterval.add(vs[k][1], inP[1]);
                const s2 = SWInterval.add(vs[k][2], inP[2]);
                x.push(xk);
                y.push(yk);
                z.push(zk);
                w.push(s0.mul(xk).add(s1.mul(yk)).add(s2.mul(zk)));
            }
            const u = (i: number, j: number): SWInterval =>
                x[i].mul(y[j]).sub(x[j].mul(y[i]));
            const vv = (i: number, j: number): SWInterval =>
                z[i].mul(w[j]).sub(z[j].mul(w[i]));
            const det = u(0, 1).mul(vv(2, 3))
                .sub(u(0, 2).mul(vv(1, 3)))
                .add(u(0, 3).mul(vv(1, 2)))
                .add(u(1, 2).mul(vv(0, 3)))
                .sub(u(1, 3).mul(vv(0, 2)))
                .add(u(2, 3).mul(vv(0, 1)));
            if (!(det.get(0) > 0) && !(det.get(1) < 0)) {
                ++this.exactCircumsphere;
            }
        }
        return super.toCircumsphere(pIndex, v0Index, v1Index, v2Index, v3Index);
    }
}

function randomCloud(rnd: () => number, n: number, scale = 1): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        points.push(v3(scale * (2 * rnd() - 1), scale * (2 * rnd() - 1),
            scale * (2 * rnd() - 1)));
    }
    return points;
}

describe('Delaunay3', () => {
    it('tetrahedralizes a single tetrahedron', () => {
        const points = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        const del = new Delaunay3();
        expect(del.compute(points)).toBe(true);
        expect(del.getDimension()).toBe(3);
        expect(del.getNumTetrahedra()).toBe(1);
        expect(del.getNumVertices()).toBe(4);
        expect(del.getNumUniqueVertices()).toBe(4);
        expect(del.getVertices()).toBe(points);
        expect(Array.from(del.getDuplicates())).toEqual([0, 1, 2, 3]);
        checkTetrahedralization(del, points);

        // The one tetrahedron has four hull faces and no adjacencies.
        expect(Array.from(del.getAdjacencies())).toEqual([-1, -1, -1, -1]);
        const hull = del.getHull();
        expect(hull.length).toBe(12);
        // Each hull triangle is counterclockwise viewed from outside, i.e.
        // the fourth (omitted) vertex is on the negative side.
        const indices = del.getIndices();
        const all = new Set(indices);
        for (let f = 0; f < 4; ++f) {
            const a = points[hull[3 * f]], b = points[hull[3 * f + 1]];
            const c = points[hull[3 * f + 2]];
            const rest = Array.from(all).filter(
                i => i !== hull[3 * f] && i !== hull[3 * f + 1] && i !== hull[3 * f + 2]);
            expect(rest.length).toBe(1);
            expect(orient3d(a, b, c, points[rest[0]])).toBeLessThan(0);
        }
    });

    it('tetrahedralizes the corners of a cube into six tetrahedra', () => {
        const points = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0),
            v3(0, 0, 1), v3(1, 0, 1), v3(0, 1, 1), v3(1, 1, 1)];
        const del = new Delaunay3();
        expect(del.compute(points)).toBe(true);
        // A cube has volume 1 and every Delaunay tetrahedron of the eight
        // corners has volume 1/6, so there are exactly six of them.
        expect(del.getNumTetrahedra()).toBe(6);
        checkTetrahedralization(del, points);

        let volume = 0;
        const indices = del.getIndices();
        for (let t = 0; t < 6; ++t) {
            volume += orient3d(points[indices[4 * t]], points[indices[4 * t + 1]],
                points[indices[4 * t + 2]], points[indices[4 * t + 3]]) / 6;
        }
        expect(volume).toBeCloseTo(1, 12);

        // The hull of a cube is 12 triangles (two per square face).
        expect(del.getHull().length).toBe(36);
    });

    it('tetrahedralizes a 3x3x3 lattice with the correct total volume', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                for (let k = 0; k < 3; ++k) {
                    points.push(v3(i, j, k));
                }
            }
        }
        const del = new Delaunay3();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumUniqueVertices()).toBe(27);
        checkTetrahedralization(del, points);

        let volume = 0;
        const indices = del.getIndices();
        for (let t = 0; t < del.getNumTetrahedra(); ++t) {
            volume += orient3d(points[indices[4 * t]], points[indices[4 * t + 1]],
                points[indices[4 * t + 2]], points[indices[4 * t + 3]]) / 6;
        }
        // The convex hull of the lattice is the 2x2x2 cube.
        expect(volume).toBeCloseTo(8, 10);

        // Every input vertex is used by some tetrahedron.
        expect(new Set(indices).size).toBe(27);
    });

    it('satisfies the empty-circumsphere property on random clouds', () => {
        const rnd = makeRandom(20260901);
        for (const n of [8, 16, 24, 32, 40]) {
            const points = randomCloud(rnd, n);
            const del = new Delaunay3();
            expect(del.compute(points)).toBe(true);
            checkTetrahedralization(del, points);
            checkEmptyCircumsphere(del, points, 1);

            // Every vertex appears in the tetrahedralization.
            expect(new Set(del.getIndices()).size).toBe(n);
        }
    });

    it('satisfies the empty-circumsphere property on a 60-point cloud', () => {
        const rnd = makeRandom(777);
        const points = randomCloud(rnd, 60, 10);
        const del = new Delaunay3();
        expect(del.compute(points)).toBe(true);
        checkTetrahedralization(del, points);
        checkEmptyCircumsphere(del, points, 10);
    });

    it('is deterministic and reusable across calls', () => {
        const rnd = makeRandom(4242);
        const points = randomCloud(rnd, 30);
        const other = [v3(0, 0, 0), v3(2, 0, 0), v3(0, 2, 0), v3(0, 0, 2), v3(1, 1, 1)];

        const del0 = new Delaunay3();
        expect(del0.compute(points)).toBe(true);
        const indices0 = Array.from(del0.getIndices());
        const adjacencies0 = Array.from(del0.getAdjacencies());

        const del1 = new Delaunay3();
        expect(del1.compute(points)).toBe(true);
        expect(Array.from(del1.getIndices())).toEqual(indices0);
        expect(Array.from(del1.getAdjacencies())).toEqual(adjacencies0);

        // A second call on the same object resets all state.
        expect(del0.compute(other)).toBe(true);
        expect(del0.getNumVertices()).toBe(5);
        expect(del0.getNumTetrahedra()).toBeGreaterThan(0);
        expect(del0.compute(points)).toBe(true);
        expect(Array.from(del0.getIndices())).toEqual(indices0);
        expect(Array.from(del0.getAdjacencies())).toEqual(adjacencies0);
    });

    it('reports duplicate vertices', () => {
        const points = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1),
            v3(1, 0, 0), v3(0.25, 0.25, 0.25), v3(0, 0, 0)];
        const del = new Delaunay3();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumVertices()).toBe(7);
        expect(del.getNumUniqueVertices()).toBe(5);
        expect(Array.from(del.getDuplicates())).toEqual([0, 1, 2, 3, 1, 5, 0]);
        checkTetrahedralization(del, points);
        // The interior point splits the tetrahedron into four.
        expect(del.getNumTetrahedra()).toBe(4);
    });

    it('handles the degenerate dimensions 0, 1 and 2', () => {
        // Dimension 0: all points identical.
        const same = [v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3)];
        const d0 = new Delaunay3();
        expect(d0.compute(same)).toBe(false);
        expect(d0.getDimension()).toBe(0);
        expect(d0.getNumVertices()).toBe(4);
        expect(d0.getNumTetrahedra()).toBe(0);
        expect(d0.getIndices().length).toBe(0);
        expect(Array.from(d0.getLine().origin.values)).toEqual([1, 2, 3]);
        expect(() => d0.getHull()).toThrow(/dimension must be 3/);

        // Dimension 1: collinear points. The intrinsics epsilon is zero, so
        // the points must be exactly collinear after the normalization of
        // the line direction; an axis-aligned set is exact.
        const collinear = [v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0), v3(-1, 0, 0),
            v3(5, 0, 0)];
        const d1 = new Delaunay3();
        expect(d1.compute(collinear)).toBe(false);
        expect(d1.getDimension()).toBe(1);
        expect(d1.getNumVertices()).toBe(5);
        expect(d1.getNumTetrahedra()).toBe(0);
        const dir = d1.getLine().direction.values;
        expect(Math.hypot(dir[0], dir[1], dir[2])).toBeCloseTo(1, 12);
        expect(Math.abs(Math.abs(dir[0]) - 1)).toBeLessThan(1e-15);
        expect(dir[1]).toBe(0);
        expect(dir[2]).toBe(0);

        // Dimension 2: coplanar points (the plane z = 0).
        const coplanar = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0),
            v3(2, 3, 0)];
        const d2 = new Delaunay3();
        expect(d2.compute(coplanar)).toBe(false);
        expect(d2.getDimension()).toBe(2);
        const normal = d2.getPlane().normal.values;
        expect(Math.abs(normal[0])).toBeLessThan(1e-15);
        expect(Math.abs(normal[1])).toBeLessThan(1e-15);
        expect(Math.abs(Math.abs(normal[2]) - 1)).toBeLessThan(1e-15);
        expect(() => d2.getTetrahedronIndices(0)).not.toThrow();
        expect(d2.getTetrahedronIndices(0)).toBeNull();

        // A near-coplanar set is 3D, since the intrinsics epsilon is zero.
        const nearCoplanar = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0),
            v3(1, 1, 1e-14), v3(0.5, 0.5, -1e-14)];
        const d3 = new Delaunay3();
        expect(d3.compute(nearCoplanar)).toBe(true);
        expect(d3.getDimension()).toBe(3);
        checkTetrahedralization(d3, nearCoplanar);
    });

    it('rejects invalid inputs', () => {
        const del = new Delaunay3();
        expect(() => del.compute([])).toThrow(/Invalid argument/);
        expect(() => del.compute([Vector.fromArray([0, 0])])).toThrow(/3D vertices/);
        expect(() => del.getContainingTetrahedron(v3(0, 0, 0), new Delaunay3SearchInfo()))
            .toThrow(/Invalid dimension/);
    });

    it('searches for the tetrahedron containing a point', () => {
        const rnd = makeRandom(31337);
        const points = randomCloud(rnd, 30);
        const del = new Delaunay3();
        expect(del.compute(points)).toBe(true);
        const indices = del.getIndices();

        // The centroid of each tetrahedron is found in a tetrahedron, and the
        // reported barycentric coordinates of the search result are valid.
        const info = new Delaunay3SearchInfo();
        for (let t = 0; t < del.getNumTetrahedra(); ++t) {
            const i0 = indices[4 * t], i1 = indices[4 * t + 1];
            const i2 = indices[4 * t + 2], i3 = indices[4 * t + 3];
            const p = v3(
                (points[i0].values[0] + points[i1].values[0] + points[i2].values[0]
                    + points[i3].values[0]) / 4,
                (points[i0].values[1] + points[i1].values[1] + points[i2].values[1]
                    + points[i3].values[1]) / 4,
                (points[i0].values[2] + points[i1].values[2] + points[i2].values[2]
                    + points[i3].values[2]) / 4);
            const found = del.getContainingTetrahedron(p, info);
            expect(found).toBe(t);
            expect(info.numPath).toBeGreaterThan(0);
            expect(info.path[info.numPath - 1]).toBe(t);
            expect(info.finalTetrahedron).toBe(t);
            // Warm-starting from the previous result also finds the point.
            info.initialTetrahedron = found;
        }

        // A point far outside the hull is not contained.
        const outside = new Delaunay3SearchInfo();
        expect(del.getContainingTetrahedron(v3(100, 100, 100), outside)).toBe(-1);
        expect(outside.numPath).toBeGreaterThan(0);

        // An out-of-range initial tetrahedron is reset to zero.
        const reset = new Delaunay3SearchInfo();
        reset.initialTetrahedron = 10000;
        del.getContainingTetrahedron(v3(100, 100, 100), reset);
        expect(reset.initialTetrahedron).toBe(0);
    });

    it('exercises the exact rational fallback on degenerate configurations', () => {
        // The corners of a cube are cospherical and many of the toPlane
        // queries are exactly coplanar, so the interval fast path cannot
        // resolve the sign and the exact predicates must decide.
        const cube = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0),
            v3(0, 0, 1), v3(1, 0, 1), v3(0, 1, 1), v3(1, 1, 1)];
        const counting = new CountingDelaunay3();
        expect(counting.compute(cube)).toBe(true);
        expect(counting.exactPlane).toBeGreaterThan(0);
        expect(counting.exactCircumsphere).toBeGreaterThan(0);
        expect(counting.getNumTetrahedra()).toBe(6);
        checkTetrahedralization(counting, cube);

        // The same result is produced by the plain class, so counting does
        // not perturb the algorithm.
        const plain = new Delaunay3();
        expect(plain.compute(cube)).toBe(true);
        expect(Array.from(counting.getIndices())).toEqual(Array.from(plain.getIndices()));

        // Nearly cospherical points: 12 points on a sphere of radius 1
        // perturbed by about one ulp. The interval widths overlap zero, so
        // the exact fallback again decides.
        const near: Vector[] = [];
        const rnd = makeRandom(99);
        for (let i = 0; i < 12; ++i) {
            const u = 2 * rnd() - 1;
            const phi = 2 * Math.PI * rnd();
            const s = Math.sqrt(1 - u * u);
            const r = 1 + (i % 2 === 0 ? 1 : -1) * 4 * Number.EPSILON;
            near.push(v3(r * s * Math.cos(phi), r * s * Math.sin(phi), r * u));
        }
        const counting2 = new CountingDelaunay3();
        expect(counting2.compute(near)).toBe(true);
        expect(counting2.exactPlane + counting2.exactCircumsphere).toBeGreaterThan(0);
        checkTetrahedralization(counting2, near);
    });

    it('produces a closed, consistently oriented hull', () => {
        const rnd = makeRandom(5150);
        const points = randomCloud(rnd, 40);
        const del = new Delaunay3();
        expect(del.compute(points)).toBe(true);
        const hull = del.getHull();
        expect(hull.length % 3).toBe(0);
        const numFaces = hull.length / 3;
        expect(numFaces).toBeGreaterThan(3);

        // Every directed hull edge occurs exactly once, which means the hull
        // surface is closed and consistently oriented.
        const directed = new Set<string>();
        for (let f = 0; f < numFaces; ++f) {
            const a = hull[3 * f], b = hull[3 * f + 1], c = hull[3 * f + 2];
            expect(new Set([a, b, c]).size).toBe(3);
            for (const [p, q] of [[a, b], [b, c], [c, a]]) {
                const key = `${p},${q}`;
                expect(directed.has(key)).toBe(false);
                directed.add(key);
            }
        }
        for (const key of directed) {
            const [p, q] = key.split(',');
            expect(directed.has(`${q},${p}`)).toBe(true);
        }

        // Every hull face is outward-oriented: all points are on the
        // non-positive side of its plane.
        for (let f = 0; f < numFaces; ++f) {
            const a = points[hull[3 * f]].values;
            const b = points[hull[3 * f + 1]].values;
            const c = points[hull[3 * f + 2]].values;
            const ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
            const vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
            const nx = uy * vz - uz * vy;
            const ny = uz * vx - ux * vz;
            const nz = ux * vy - uy * vx;
            const scale = Math.hypot(nx, ny, nz);
            for (const p of points) {
                const w = p.values;
                const d = nx * (w[0] - a[0]) + ny * (w[1] - a[1]) + nz * (w[2] - a[2]);
                expect(d / scale).toBeLessThan(1e-10);
            }
        }
    });
});
