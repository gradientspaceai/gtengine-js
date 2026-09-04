import { describe, it, expect } from 'vitest';
import { ConvexHull3 } from '../src/ConvexHull3.js';
import { Vector } from '../src/Vector.js';
import { fc, check, latticeVector } from './helpers/arbitraries.js';
import { orient3 } from './helpers/exact.js';

const v3 = (x: number, y: number, z: number): Vector =>
    Vector.fromArray([x, y, z]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function sub(a: Vector, b: Vector): [number, number, number] {
    return [a.values[0] - b.values[0], a.values[1] - b.values[1],
        a.values[2] - b.values[2]];
}

function cross(a: readonly number[], b: readonly number[]): [number, number, number] {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]];
}

function dot(a: readonly number[], b: readonly number[]): number {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Every input point must be on the negative side of (or on) every face plane,
// where the face normal is the outward normal implied by the counterclockwise
// face ordering.
function checkConvexAndContains(points: readonly Vector[],
    hull: readonly number[], tolerance: number): void {
    expect(hull.length % 3).toBe(0);
    expect(hull.length).toBeGreaterThan(0);
    for (let f = 0; f < hull.length; f += 3) {
        const v0 = points[hull[f]];
        const v1 = points[hull[f + 1]];
        const v2 = points[hull[f + 2]];
        const normal = cross(sub(v1, v0), sub(v2, v0));
        const normalLength = Math.sqrt(dot(normal, normal));
        expect(normalLength).toBeGreaterThan(0);
        for (const p of points) {
            const signedDistance = dot(normal, sub(p, v0)) / normalLength;
            expect(signedDistance).toBeLessThanOrEqual(tolerance);
        }
    }
}

// The hull mesh is closed and oriented, so counting the distinct undirected
// edges must give E = 3*T/2 and Euler's formula V - E + T = 2 must hold.
function checkEuler(numVertices: number, hull: readonly number[]): void {
    const numTriangles = hull.length / 3;
    const edges = new Set<string>();
    for (let f = 0; f < hull.length; f += 3) {
        for (let i = 0; i < 3; ++i) {
            const a = hull[f + i];
            const b = hull[f + (i + 1) % 3];
            edges.add(a < b ? `${a},${b}` : `${b},${a}`);
        }
    }
    expect(edges.size).toBe(3 * numTriangles / 2);
    expect(numVertices - edges.size + numTriangles).toBe(2);
}

describe('ConvexHull3', () => {
    it('rejects invalid input', () => {
        const ch = new ConvexHull3();
        expect(() => ch.compute([])).toThrow();
        expect(() => ch.compute([Vector.fromArray([0, 0])])).toThrow();
    });

    it('computes a 0-dimensional hull for coincident points', () => {
        const points = [v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3)];
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(0);
        expect(ch.getHull().length).toBe(1);
        expect(points[ch.getHull()[0]].equals(v3(1, 2, 3))).toBe(true);
        expect(ch.getVertices().length).toBe(1);
    });

    it('computes a 1-dimensional hull for collinear points', () => {
        // Points on the line (t, 2t, 3t), given out of order and with a
        // duplicate.
        const ts = [3, -1, 0, 5, 2, 5];
        const points = ts.map(t => v3(t, 2 * t, 3 * t));
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(1);
        const hull = ch.getHull();
        expect(hull.length).toBe(2);
        const endpointTs = [ts[hull[0]], ts[hull[1]]].sort((a, b) => a - b);
        expect(endpointTs).toEqual([-1, 5]);
        expect(ch.getVertices().length).toBe(2);
    });

    it('computes a 2-dimensional hull for coplanar points', () => {
        // A unit square in the plane z = 4, plus interior and edge-interior
        // points that must not be hull vertices.
        const points = [
            v3(0, 0, 4), v3(1, 0, 4), v3(1, 1, 4), v3(0, 1, 4),
            v3(0.5, 0.5, 4), v3(0.5, 0, 4), v3(0.25, 0.75, 4)
        ];
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(2);
        const hull = ch.getHull();
        expect(hull.length).toBe(4);
        expect([...hull].sort((a, b) => a - b)).toEqual([0, 1, 2, 3]);

        // The polygon must be ordered (consecutive corners of the square).
        for (let i = 0; i < 4; ++i) {
            const a = points[hull[i]];
            const b = points[hull[(i + 1) % 4]];
            const d = sub(b, a);
            expect(Math.abs(dot(d, d) - 1)).toBeLessThan(1e-12);
        }
        expect(ch.getVertices().length).toBe(4);
    });

    it('computes the hull of the corners of a cube', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            points.push(v3(i & 1, (i >> 1) & 1, (i >> 2) & 1));
        }
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        expect(ch.getVertices().length).toBe(8);

        const hull = ch.getHull();
        // Each square face of the cube is split into 2 triangles.
        expect(hull.length).toBe(3 * 12);
        checkConvexAndContains(points, hull, 1e-12);
        checkEuler(8, hull);

        // The hull mesh is a closed, oriented manifold.
        const mesh = ch.getHullMesh();
        expect(mesh.getNumTriangles()).toBe(12);
        expect(mesh.getNumEdges()).toBe(18);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(true);
        expect(mesh.getNumVertices()).toBe(8);
    });

    it('ignores interior points and duplicates', () => {
        // Cube corners, the cube center (interior), a face center (on the
        // boundary but not a vertex) and a repeated corner.
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            points.push(v3(i & 1, (i >> 1) & 1, (i >> 2) & 1));
        }
        points.push(v3(0.5, 0.5, 0.5));
        points.push(v3(0.5, 0.5, 0));
        points.push(v3(1, 1, 1));
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        const vertices = [...ch.getVertices()];
        // All 8 distinct corners are hull vertices.
        for (let i = 0; i < 8; ++i) {
            expect(vertices).toContain(i);
        }
        // The cube center (index 8) is strictly interior, so it is not a hull
        // vertex. The repeated corner (index 10) is dropped by the uniqueness
        // filter in favor of the first occurrence (index 7).
        expect(vertices).not.toContain(8);
        expect(vertices).not.toContain(10);
        // The face center (index 9) lies on the boundary; upstream keeps such
        // a point as a vertex of coplanar hull triangles, which the header
        // documents as expected behavior.
        expect(vertices.length).toBe(9);
        checkConvexAndContains(points, ch.getHull(), 1e-12);
        checkEuler(vertices.length, ch.getHull());
    });

    it('computes the hull of a regular octahedron', () => {
        const points = [
            v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
            v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)
        ];
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        expect(ch.getVertices().length).toBe(6);
        const hull = ch.getHull();
        expect(hull.length).toBe(3 * 8);
        checkConvexAndContains(points, hull, 1e-12);
        checkEuler(6, hull);
    });

    it('computes hulls of random clouds', () => {
        const random = makeRandom(20260902);
        for (let trial = 0; trial < 6; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 30; ++i) {
                points.push(v3(2 * random() - 1, 2 * random() - 1, 2 * random() - 1));
            }
            const ch = new ConvexHull3();
            ch.compute(points);
            expect(ch.getDimension()).toBe(3);

            const hull = ch.getHull();
            checkConvexAndContains(points, hull, 1e-9);
            checkEuler(ch.getVertices().length, hull);

            // The hull vertices are exactly the indices appearing in the
            // triangle list.
            const used = new Set<number>(hull);
            expect([...used].sort((a, b) => a - b))
                .toEqual([...ch.getVertices()].sort((a, b) => a - b));

            const mesh = ch.getHullMesh();
            expect(mesh.isClosed()).toBe(true);
            expect(mesh.isOriented()).toBe(true);
        }
    });

    it('computes hulls of random points on a sphere (all points are vertices)', () => {
        const random = makeRandom(777);
        const points: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            const z = 2 * random() - 1;
            const theta = 2 * Math.PI * random();
            const r = Math.sqrt(1 - z * z);
            points.push(v3(r * Math.cos(theta), r * Math.sin(theta), z));
        }
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        // Points on a sphere are in convex position, so all are hull vertices.
        expect(ch.getVertices().length).toBe(points.length);
        checkConvexAndContains(points, ch.getHull(), 1e-9);
        checkEuler(points.length, ch.getHull());
    });

    it('is exact on a near-degenerate set (interval fallback to rational)', () => {
        // A grid of points on the plane z = 0 with one point lifted by a
        // quantity so tiny that the interval predicate cannot resolve the
        // sign of the determinant. The exact fallback must classify it as
        // above the plane, making the hull 3-dimensional.
        const points: Vector[] = [];
        for (let x = 0; x <= 2; ++x) {
            for (let y = 0; y <= 2; ++y) {
                points.push(v3(x, y, 0));
            }
        }
        const lift = Math.pow(2, -60);
        points.push(v3(1, 1, lift));

        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        // The lifted point is a hull vertex even though it is only 2^-60
        // above the plane of the others.
        expect([...ch.getVertices()]).toContain(points.length - 1);
        checkConvexAndContains(points, ch.getHull(), 1e-12);
        checkEuler(ch.getVertices().length, ch.getHull());

        // Without the lift the same set is coplanar and the hull is 2D.
        const flat = points.slice(0, points.length - 1);
        const chFlat = new ConvexHull3();
        chFlat.compute(flat);
        expect(chFlat.getDimension()).toBe(2);
        expect(chFlat.getHull().length).toBe(4);
    });

    it('reuses the functor across data sets', () => {
        const ch = new ConvexHull3();
        ch.compute([v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)]);
        expect(ch.getDimension()).toBe(3);
        expect(ch.getHull().length).toBe(3 * 4);

        ch.compute([v3(5, 5, 5), v3(5, 5, 5)]);
        expect(ch.getDimension()).toBe(0);
        expect(ch.getHull().length).toBe(1);
        expect(ch.getVertices().length).toBe(1);
        expect(ch.getHullMesh().getNumTriangles()).toBe(0);
    });
});

// ---------------------------------------------------------------------------
// Verification (V11): property-based cross-checks against exact bigint
// oracles. The generators are integer lattices, so every predicate below is
// evaluated exactly and the assertions are deterministic.
// ---------------------------------------------------------------------------

type BigPoint = [bigint, bigint, bigint];

const toBig = (p: Vector): BigPoint =>
    [BigInt(p.values[0]), BigInt(p.values[1]), BigInt(p.values[2])];

const bigEquals = (a: BigPoint, b: BigPoint): boolean =>
    a[0] === b[0] && a[1] === b[1] && a[2] === b[2];

const bigSub = (a: BigPoint, b: BigPoint): BigPoint =>
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

const bigCross = (a: BigPoint, b: BigPoint): BigPoint =>
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]];

const bigDot = (a: BigPoint, b: BigPoint): bigint =>
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

const isZero = (a: BigPoint): boolean =>
    a[0] === 0n && a[1] === 0n && a[2] === 0n;

function uniquePoints(points: readonly Vector[]): BigPoint[] {
    const uniq: BigPoint[] = [];
    for (const p of points) {
        const b = toBig(p);
        if (!uniq.some(q => bigEquals(q, b))) { uniq.push(b); }
    }
    return uniq;
}

/**
 * The intrinsic dimension of an integer point set, computed exactly: 0 for a
 * single point, 1 when all points are collinear, 2 when all are coplanar and
 * 3 otherwise. This is the independent oracle for getDimension().
 */
function exactDimension(points: readonly Vector[]): number {
    const uniq = uniquePoints(points);
    if (uniq.length === 1) { return 0; }
    const a = uniq[0], b = uniq[1];
    let c: BigPoint | null = null;
    for (let i = 2; i < uniq.length; ++i) {
        if (!isZero(bigCross(bigSub(b, a), bigSub(uniq[i], a)))) {
            c = uniq[i];
            break;
        }
    }
    if (c === null) { return 1; }
    for (const d of uniq) {
        if (orient3(a, b, c, d) !== 0) { return 3; }
    }
    return 2;
}

/** Lexicographic order on integer points, matching the port's point sort. */
function lexLess(a: BigPoint, b: BigPoint): boolean {
    for (let i = 0; i < 3; ++i) {
        if (a[i] !== b[i]) { return a[i] < b[i]; }
    }
    return false;
}

/** Canonical string for a face, invariant under cyclic rotation. */
function faceKey(a: BigPoint, b: BigPoint, c: BigPoint): string {
    const s = [a, b, c].map(p => p.join(',')) as [string, string, string];
    const rotations = [
        s[0] + '|' + s[1] + '|' + s[2],
        s[1] + '|' + s[2] + '|' + s[0],
        s[2] + '|' + s[0] + '|' + s[1]
    ];
    rotations.sort();
    return rotations[0];
}

function faceKeys(points: readonly Vector[], hull: readonly number[]): string[] {
    const keys: string[] = [];
    for (let f = 0; f < hull.length; f += 3) {
        keys.push(faceKey(toBig(points[hull[f]]), toBig(points[hull[f + 1]]),
            toBig(points[hull[f + 2]])));
    }
    keys.sort();
    return keys;
}

const latticeCloud3 = (count: number, range: number): fc.Arbitrary<Vector[]> =>
    fc.array(latticeVector(3, -range, range),
        { minLength: count, maxLength: count });

describe('ConvexHull3 verification', () => {
    it('classifies the intrinsic dimension exactly', () => {
        // A tiny lattice range makes duplicate, collinear and coplanar draws
        // common, so all four branches of the classification are exercised.
        // ConvexHull3 uses exact colocated/colinear/toPlane predicates rather
        // than the hardcoded epsilon of Delaunay2/Delaunay3 (issue #391), so
        // there is no misclassification to fix here.
        // Degenerate families are generated explicitly: random lattice
        // clouds are almost always 3-dimensional.
        const repeated = fc.tuple(latticeVector(3, -6, 6),
            fc.integer({ min: 1, max: 5 }))
            .map(([p, n]) => new Array<Vector>(n).fill(p));
        const collinear = fc.tuple(latticeVector(3, -6, 6),
            latticeVector(3, -3, 3), fc.array(fc.integer({ min: -4, max: 4 }),
                { minLength: 2, maxLength: 6 }))
            .map(([o, d, ts]) => ts.map(t => Vector.fromArray([
                o.values[0] + t * d.values[0], o.values[1] + t * d.values[1],
                o.values[2] + t * d.values[2]])));
        const coplanar = fc.tuple(fc.integer({ min: -3, max: 3 }),
            fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 }),
            fc.array(fc.tuple(fc.integer({ min: -4, max: 4 }),
                fc.integer({ min: -4, max: 4 })),
                { minLength: 2, maxLength: 7 }))
            .map(([a, b, c, xy]) => xy.map(([x, y]) =>
                Vector.fromArray([x, y, a * x + b * y + c])));
        check(fc.oneof(repeated, collinear, coplanar, latticeCloud3(5, 1),
            latticeCloud3(7, 2), latticeCloud3(9, 4)), points => {
            const ch = new ConvexHull3();
            ch.compute(points);
            expect(ch.getDimension()).toBe(exactDimension(points));
            return true;
        }, 200);
    }, 30000);

    it('reports the lexicographic extremes for a 1-dimensional hull', () => {
        check(fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -6, 6),
            fc.array(fc.integer({ min: -6, max: 6 }),
                { minLength: 1, maxLength: 6 })),
            ([origin, direction, ts]) => {
                if (isZero(toBig(direction))) { return true; }
                const points = ts.map(t => Vector.fromArray([
                    origin.values[0] + t * direction.values[0],
                    origin.values[1] + t * direction.values[1],
                    origin.values[2] + t * direction.values[2]]));
                const ch = new ConvexHull3();
                ch.compute(points);
                if (ch.getDimension() !== 1) {
                    // All the sample parameters coincided: a single point.
                    expect(ch.getDimension()).toBe(0);
                    return true;
                }
                const hull = ch.getHull();
                expect(hull.length).toBe(2);
                const uniq = uniquePoints(points);
                let lo = uniq[0], hi = uniq[0];
                for (const p of uniq) {
                    if (lexLess(p, lo)) { lo = p; }
                    if (lexLess(hi, p)) { hi = p; }
                }
                expect(toBig(points[hull[0]])).toEqual(lo);
                expect(toBig(points[hull[1]])).toEqual(hi);
                expect(ch.getVertices().slice()).toEqual(hull.slice());
                return true;
            }, 200);
    }, 30000);

    it('reports a convex ordered polygon for a 2-dimensional hull', () => {
        // Points on the plane z = a*x + b*y + c are exactly coplanar for
        // integer coefficients.
        check(fc.tuple(fc.integer({ min: -3, max: 3 }),
            fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 }),
            fc.array(fc.tuple(fc.integer({ min: -5, max: 5 }),
                fc.integer({ min: -5, max: 5 })),
                { minLength: 4, maxLength: 9 })),
            ([a, b, c, xy]) => {
                const points = xy.map(([x, y]) =>
                    Vector.fromArray([x, y, a * x + b * y + c]));
                const ch = new ConvexHull3();
                ch.compute(points);
                if (ch.getDimension() !== 2) { return true; }

                const hull = ch.getHull();
                expect(hull.length).toBeGreaterThanOrEqual(3);
                // The polygon vertices are distinct hull points and the
                // vertex list agrees with the hull list.
                expect(ch.getVertices().slice()).toEqual(hull.slice());
                const poly = hull.map(h => toBig(points[h]));
                for (let i = 0; i < poly.length; ++i) {
                    for (let j = i + 1; j < poly.length; ++j) {
                        expect(bigEquals(poly[i], poly[j])).toBe(false);
                    }
                }

                // The polygon is counterclockwise when viewed from the side
                // its own normal points to, and every input point is inside
                // or on the boundary.
                const normal = bigCross(bigSub(poly[1], poly[0]),
                    bigSub(poly[2], poly[0]));
                expect(isZero(normal)).toBe(false);
                const uniq = uniquePoints(points);
                for (let i = 0; i < poly.length; ++i) {
                    const p0 = poly[i];
                    const p1 = poly[(i + 1) % poly.length];
                    const edge = bigSub(p1, p0);
                    for (const q of uniq) {
                        const side = bigDot(normal, bigCross(edge, bigSub(q, p0)));
                        expect(side >= 0n).toBe(true);
                    }
                    // Strict convexity: the next vertex is strictly left of
                    // the current edge, so no three consecutive polygon
                    // vertices are collinear.
                    const p2 = poly[(i + 2) % poly.length];
                    expect(bigDot(normal, bigCross(edge, bigSub(p2, p0))) > 0n)
                        .toBe(true);
                }
                return true;
            }, 150);
    }, 30000);

    it('is a closed oriented manifold that contains every input point', () => {
        check(fc.oneof(latticeCloud3(8, 3), latticeCloud3(12, 5),
            latticeCloud3(16, 8)), points => {
            const ch = new ConvexHull3();
            ch.compute(points);
            if (ch.getDimension() !== 3) { return true; }

            const hull = ch.getHull();
            expect(hull.length % 3).toBe(0);
            const numTriangles = hull.length / 3;
            const uniq = uniquePoints(points);

            // Every face is nondegenerate and supports the point set: with
            // counterclockwise-outward ordering, orient3 <= 0 for every
            // input point (exact).
            for (let f = 0; f < hull.length; f += 3) {
                const a = toBig(points[hull[f]]);
                const b = toBig(points[hull[f + 1]]);
                const c = toBig(points[hull[f + 2]]);
                expect(isZero(bigCross(bigSub(b, a), bigSub(c, a)))).toBe(false);
                for (const q of uniq) {
                    expect(orient3(a, b, c, q) <= 0).toBe(true);
                }
            }

            // Every directed edge occurs exactly once and its reverse occurs
            // exactly once: the surface is closed, oriented and manifold.
            const directed = new Set<string>();
            for (let f = 0; f < hull.length; f += 3) {
                for (let i = 0; i < 3; ++i) {
                    const key = hull[f + i] + '->' + hull[f + (i + 1) % 3];
                    expect(directed.has(key)).toBe(false);
                    directed.add(key);
                }
            }
            for (const key of directed) {
                const [u, v] = key.split('->');
                expect(directed.has(v + '->' + u)).toBe(true);
            }

            // Euler's formula with E = 3T/2 (the header claims E = T/2, an
            // upstream documentation bug corrected in the port).
            const numEdges = directed.size / 2;
            expect(numEdges).toBe(3 * numTriangles / 2);
            const vertexSet = new Set<number>(hull);
            expect(ch.getVertices().slice().sort((p, q) => p - q))
                .toEqual(Array.from(vertexSet).sort((p, q) => p - q));
            expect(vertexSet.size - numEdges + numTriangles).toBe(2);
            return true;
        }, 120);
    }, 30000);

    it('does not depend on the order of the input points', () => {
        // The port sorts the points lexicographically before building the
        // hull, so a permutation of the input must produce the same faces
        // (compared by coordinates, since duplicates may pick a different
        // representative index).
        check(fc.tuple(latticeCloud3(11, 4), fc.array(fc.nat(1000),
            { minLength: 11, maxLength: 11 })), ([points, keys]) => {
            const ch = new ConvexHull3();
            ch.compute(points);
            const expected = faceKeys(points, ch.getHull());
            const dimension = ch.getDimension();

            const order = points.map((p, i) => i);
            order.sort((i, j) => (keys[i] - keys[j]) || (i - j));
            const permuted = order.map(i => points[i]);
            const ch2 = new ConvexHull3();
            ch2.compute(permuted);
            expect(ch2.getDimension()).toBe(dimension);
            if (dimension === 3) {
                expect(faceKeys(permuted, ch2.getHull())).toEqual(expected);
            }
            return true;
        }, 120);
    }, 30000);

    it('is idempotent on its own hull vertices', () => {
        // Recomputing the hull from just the hull vertices must produce the
        // same faces, which pins down that no extreme point was dropped.
        check(fc.oneof(latticeCloud3(10, 4), latticeCloud3(14, 7)), points => {
            const ch = new ConvexHull3();
            ch.compute(points);
            if (ch.getDimension() !== 3) { return true; }
            const expected = faceKeys(points, ch.getHull());

            const reduced = ch.getVertices().map(i => points[i]);
            const ch2 = new ConvexHull3();
            ch2.compute(reduced);
            expect(ch2.getDimension()).toBe(3);
            expect(faceKeys(reduced, ch2.getHull())).toEqual(expected);
            return true;
        }, 120);
    }, 30000);

    it('reuses the functor without leaking state between data sets', () => {
        check(fc.tuple(latticeCloud3(9, 4), latticeCloud3(9, 4)),
            ([first, second]) => {
                const fresh = new ConvexHull3();
                fresh.compute(second);

                const reused = new ConvexHull3();
                reused.compute(first);
                reused.compute(second);

                expect(reused.getDimension()).toBe(fresh.getDimension());
                expect(reused.getHull().slice()).toEqual(fresh.getHull().slice());
                expect(reused.getVertices().slice())
                    .toEqual(fresh.getVertices().slice());
                return true;
            }, 120);
    }, 30000);
});
