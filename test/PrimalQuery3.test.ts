import { describe, it, expect } from 'vitest';
import { PrimalQuery3 } from '../src/PrimalQuery3.js';
import { Vector } from '../src/Vector.js';

const v3 = (x: number, y: number, z: number): Vector => Vector.fromArray([x, y, z]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function sign(value: number): number {
    return (value > 0 ? +1 : (value < 0 ? -1 : 0));
}

// Dot(P-V0, Cross(V1-V0, V2-V0)), computed independently of the query.
function tripleProduct(P: Vector, v0: Vector, v1: Vector, v2: Vector): number {
    const a = [P.values[0] - v0.values[0], P.values[1] - v0.values[1],
        P.values[2] - v0.values[2]];
    const u = [v1.values[0] - v0.values[0], v1.values[1] - v0.values[1],
        v1.values[2] - v0.values[2]];
    const w = [v2.values[0] - v0.values[0], v2.values[1] - v0.values[1],
        v2.values[2] - v0.values[2]];
    const cross = [
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0]
    ];
    return a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2];
}

describe('PrimalQuery3', () => {
    it('has a default constructor with no vertices', () => {
        const query = new PrimalQuery3();
        expect(query.getNumVertices()).toBe(0);
        expect(query.getVertices()).toEqual([]);
    });

    it('supports set() member access', () => {
        const vertices = [v3(0, 0, 0), v3(1, 0, 0)];
        const query = new PrimalQuery3();
        query.set(vertices.length, vertices);
        expect(query.getNumVertices()).toBe(2);
        expect(query.getVertices()).toBe(vertices);
    });

    describe('toPlane', () => {
        // The plane z = 0 with V0 = (0,0,0), V1 = (1,0,0), V2 = (0,1,0), so
        // N = Cross(V1-V0, V2-V0) = (0,0,1).
        const vertices = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0),
            v3(3, -4, 5),   // 3: positive side
            v3(3, -4, -5),  // 4: negative side
            v3(7, 9, 0)     // 5: on the plane
        ];
        const query = new PrimalQuery3(vertices.length, vertices);

        it('classifies points by the side of the plane', () => {
            expect(query.toPlane(3, 0, 1, 2)).toBe(+1);
            expect(query.toPlane(4, 0, 1, 2)).toBe(-1);
            expect(query.toPlane(5, 0, 1, 2)).toBe(0);
            for (const i of [0, 1, 2]) {
                expect(query.toPlane(i, 0, 1, 2)).toBe(0);
            }
        });

        it('accepts a test point as well as a vertex index', () => {
            expect(query.toPlane(v3(0, 0, 2), 0, 1, 2)).toBe(+1);
            expect(query.toPlane(v3(0, 0, -2), 0, 1, 2)).toBe(-1);
        });

        it('flips sign when the plane orientation is reversed', () => {
            expect(query.toPlane(3, 0, 2, 1)).toBe(-1);
            expect(query.toPlane(4, 0, 2, 1)).toBe(+1);
        });

        it('agrees with the sign of Dot(P-V0, Cross(V1-V0, V2-V0))', () => {
            const random = makeRandom(80808);
            for (let trial = 0; trial < 400; ++trial) {
                const coordinate = (): number => Math.floor(9 * random()) - 4;
                const points = [
                    v3(coordinate(), coordinate(), coordinate()),
                    v3(coordinate(), coordinate(), coordinate()),
                    v3(coordinate(), coordinate(), coordinate()),
                    v3(coordinate(), coordinate(), coordinate())
                ];
                const query2 = new PrimalQuery3(4, points);
                expect(query2.toPlane(3, 0, 1, 2))
                    .toBe(sign(tripleProduct(points[3], points[0], points[1], points[2])));
            }
        });
    });

    describe('toTetrahedron', () => {
        // The tetrahedron V0 = (0,0,0), V1 = (1,0,0), V2 = (0,1,0),
        // V3 = (0,0,1) is ordered as PrimalQuery3 expects: an interior point
        // classifies as -1.
        const vertices = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1),
            v3(0.1, 0.1, 0.1),      // 4: interior
            v3(5, 5, 5),            // 5: exterior
            v3(0.25, 0.25, 0),      // 6: on the face z = 0
            v3(-0.1, 0.1, 0.1),     // 7: exterior, x < 0
            v3(0.5, 0.25, 0.25)     // 8: on the slanted face x+y+z = 1
        ];
        const query = new PrimalQuery3(vertices.length, vertices);

        it('returns -1 inside, +1 outside and 0 on the boundary', () => {
            expect(query.toTetrahedron(4, 0, 1, 2, 3)).toBe(-1);
            expect(query.toTetrahedron(5, 0, 1, 2, 3)).toBe(+1);
            expect(query.toTetrahedron(7, 0, 1, 2, 3)).toBe(+1);
            expect(query.toTetrahedron(6, 0, 1, 2, 3)).toBe(0);
            expect(query.toTetrahedron(8, 0, 1, 2, 3)).toBe(0);
        });

        it('returns 0 for the tetrahedron vertices themselves', () => {
            for (const i of [0, 1, 2, 3]) {
                expect(query.toTetrahedron(i, 0, 1, 2, 3)).toBe(0);
            }
        });

        it('agrees with a barycentric containment test (randomized)', () => {
            const random = makeRandom(20260901);
            let numInside = 0;
            let numOutside = 0;
            let numOn = 0;
            for (let trial = 0; trial < 600; ++trial) {
                // Test points on a coarse grid inside and around the
                // reference tetrahedron; the coordinates are exact multiples
                // of 1/8 so the boundary cases are hit exactly.
                const eighth = (): number => (Math.floor(9 * random()) - 2) / 8;
                const p = v3(eighth(), eighth(), eighth());
                const x = p.values[0];
                const y = p.values[1];
                const z = p.values[2];
                // Barycentric coordinates of the reference tetrahedron.
                const bary = [1 - x - y - z, x, y, z];
                let expected: number;
                if (bary.some(b => b < 0)) {
                    expected = +1;
                    ++numOutside;
                }
                else if (bary.some(b => b === 0)) {
                    expected = 0;
                    ++numOn;
                }
                else {
                    expected = -1;
                    ++numInside;
                }
                expect(query.toTetrahedron(p, 0, 1, 2, 3)).toBe(expected);
            }
            expect(numInside).toBeGreaterThan(0);
            expect(numOutside).toBeGreaterThan(0);
            expect(numOn).toBeGreaterThan(0);
        });
    });

    describe('toCircumsphere', () => {
        // The reference tetrahedron has circumcenter (0.5,0.5,0.5) and
        // circumradius sqrt(3)/2.
        const vertices = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1),
            v3(0.5, 0.5, 0.5),  // 4: the circumcenter, strictly inside
            v3(5, 5, 5),        // 5: far outside
            v3(1, 1, 1),        // 6: exactly on the circumsphere
            v3(1, 1, 0)         // 7: exactly on the circumsphere
        ];
        const query = new PrimalQuery3(vertices.length, vertices);

        it('returns -1 inside, +1 outside and 0 on the circumsphere', () => {
            expect(query.toCircumsphere(4, 0, 1, 2, 3)).toBe(-1);
            expect(query.toCircumsphere(5, 0, 1, 2, 3)).toBe(+1);
            expect(query.toCircumsphere(6, 0, 1, 2, 3)).toBe(0);
            expect(query.toCircumsphere(7, 0, 1, 2, 3)).toBe(0);
        });

        it('returns 0 for the tetrahedron vertices themselves', () => {
            for (const i of [0, 1, 2, 3]) {
                expect(query.toCircumsphere(i, 0, 1, 2, 3)).toBe(0);
            }
        });

        it('agrees with an explicit circumcenter/radius test (randomized)', () => {
            const random = makeRandom(1234567);
            let numInside = 0;
            let numOutside = 0;
            for (let trial = 0; trial < 400; ++trial) {
                const coordinate = (): number => Math.floor(9 * random()) - 4;
                const t = [
                    v3(coordinate(), coordinate(), coordinate()),
                    v3(coordinate(), coordinate(), coordinate()),
                    v3(coordinate(), coordinate(), coordinate()),
                    v3(coordinate(), coordinate(), coordinate())
                ];
                // The query documents the TetrahedronKey ordering, for which
                // Dot(V3-V0, Cross(V1-V0, V2-V0)) > 0.
                const volume = tripleProduct(t[3], t[0], t[1], t[2]);
                if (volume <= 0) {
                    continue;
                }

                // Solve for the circumcenter: |C-Vi|^2 equal for all i gives
                // three linear equations.
                const rows: number[][] = [];
                for (let i = 1; i < 4; ++i) {
                    const dx = t[i].values[0] - t[0].values[0];
                    const dy = t[i].values[1] - t[0].values[1];
                    const dz = t[i].values[2] - t[0].values[2];
                    const rhs = 0.5 * (dx * (t[i].values[0] + t[0].values[0])
                        + dy * (t[i].values[1] + t[0].values[1])
                        + dz * (t[i].values[2] + t[0].values[2]));
                    rows.push([dx, dy, dz, rhs]);
                }
                const center = solve3x3(rows);
                if (center === null) {
                    continue;
                }
                const radius = Math.hypot(center[0] - t[0].values[0],
                    center[1] - t[0].values[1], center[2] - t[0].values[2]);

                const p = v3(coordinate(), coordinate(), coordinate());
                const distance = Math.hypot(p.values[0] - center[0],
                    p.values[1] - center[1], p.values[2] - center[2]);
                if (Math.abs(distance - radius) < 1e-8 * Math.max(1, radius)) {
                    // Skip the on-sphere cases, where the floating-point
                    // reference computation is itself ambiguous.
                    continue;
                }

                const q = new PrimalQuery3(4, t);
                const expected = (distance < radius ? -1 : +1);
                if (expected < 0) {
                    ++numInside;
                }
                else {
                    ++numOutside;
                }
                expect(q.toCircumsphere(p, 0, 1, 2, 3)).toBe(expected);
            }
            expect(numInside).toBeGreaterThan(0);
            expect(numOutside).toBeGreaterThan(0);
        });
    });
});

// Gaussian elimination with partial pivoting for a 3x3 system given as three
// rows [a, b, c, rhs]. Returns null when the matrix is (numerically) singular.
function solve3x3(rows: number[][]): number[] | null {
    const m = rows.map(row => row.slice());
    for (let col = 0; col < 3; ++col) {
        let pivot = col;
        for (let r = col + 1; r < 3; ++r) {
            if (Math.abs(m[r][col]) > Math.abs(m[pivot][col])) {
                pivot = r;
            }
        }
        if (Math.abs(m[pivot][col]) < 1e-12) {
            return null;
        }
        [m[col], m[pivot]] = [m[pivot], m[col]];
        for (let r = 0; r < 3; ++r) {
            if (r === col) {
                continue;
            }
            const factor = m[r][col] / m[col][col];
            for (let c = col; c < 4; ++c) {
                m[r][c] -= factor * m[col][c];
            }
        }
    }
    return [m[0][3] / m[0][0], m[1][3] / m[1][1], m[2][3] / m[2][2]];
}
