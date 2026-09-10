import { describe, it, expect } from 'vitest';
import { computeMassProperties } from '../src/PolyhedralMassProperties.js';
import { Vector } from '../src/Vector.js';
import { Matrix, mulMatrix, multiplyAB, multiplyABT } from '../src/Matrix.js';
import {
    check, expectClose, fc, scaled as scaledArb
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A box [x0,x1] x [y0,y1] x [z0,z1] with outward-pointing triangle normals.
function makeBox(x0: number, y0: number, z0: number,
    x1: number, y1: number, z1: number): { vertices: Vector[], indices: number[] } {
    const vertices = [
        v3(x0, y0, z0), v3(x1, y0, z0), v3(x1, y1, z0), v3(x0, y1, z0),
        v3(x0, y0, z1), v3(x1, y0, z1), v3(x1, y1, z1), v3(x0, y1, z1)
    ];
    const indices = [
        0, 3, 2, 0, 2, 1,   // z = z0, normal -z
        4, 5, 6, 4, 6, 7,   // z = z1, normal +z
        0, 1, 5, 0, 5, 4,   // y = y0, normal -y
        3, 7, 6, 3, 6, 2,   // y = y1, normal +y
        0, 4, 7, 0, 7, 3,   // x = x0, normal -x
        1, 2, 6, 1, 6, 5    // x = x1, normal +x
    ];
    return { vertices, indices };
}

// The tetrahedron with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1), scaled by
// 's' and with outward-pointing triangle normals.
function makeTetrahedron(s: number): { vertices: Vector[], indices: number[] } {
    const vertices = [v3(0, 0, 0), v3(s, 0, 0), v3(0, s, 0), v3(0, 0, s)];
    const indices = [
        0, 2, 1,   // z = 0, normal -z
        0, 1, 3,   // y = 0, normal -y
        0, 3, 2,   // x = 0, normal -x
        1, 2, 3    // slanted face, normal (1,1,1)/sqrt(3)
    ];
    return { vertices, indices };
}

function translate(vertices: Vector[], t: Vector): Vector[] {
    return vertices.map(v => v3(v.values[0] + t.values[0], v.values[1] + t.values[1],
        v.values[2] + t.values[2]));
}

function expectMatrixClose(M: Matrix, expected: number[][], tol: number): void {
    for (let r = 0; r < 3; ++r) {
        for (let c = 0; c < 3; ++c) {
            expect(M.get(r, c)).toBeCloseTo(expected[r][c], tol);
        }
    }
}

describe('PolyhedralMassProperties', () => {
    it('computes the unit cube [0,1]^3 mass, center and world inertia', () => {
        const { vertices, indices } = makeBox(0, 0, 0, 1, 1, 1);
        const result = computeMassProperties(vertices, 12, indices, false);

        expect(result.mass).toBeCloseTo(1, 12);
        expect(result.center.values[0]).toBeCloseTo(0.5, 12);
        expect(result.center.values[1]).toBeCloseTo(0.5, 12);
        expect(result.center.values[2]).toBeCloseTo(0.5, 12);

        // Relative to the world origin: I_xx = integral of y^2 + z^2 = 2/3,
        // I_xy = -integral of x*y = -1/4.
        expectMatrixClose(result.inertia, [
            [2 / 3, -1 / 4, -1 / 4],
            [-1 / 4, 2 / 3, -1 / 4],
            [-1 / 4, -1 / 4, 2 / 3]
        ], 12);
    });

    it('computes the unit cube body-coordinate inertia (diagonal, m/6)', () => {
        const { vertices, indices } = makeBox(0, 0, 0, 1, 1, 1);
        const result = computeMassProperties(vertices, 12, indices, true);

        // For a cube of side L and mass m, I = m*(L^2 + L^2)/12 = m/6 and the
        // products of inertia vanish.
        expectMatrixClose(result.inertia, [
            [1 / 6, 0, 0],
            [0, 1 / 6, 0],
            [0, 0, 1 / 6]
        ], 12);
    });

    it('computes a 2x2x2 box centered at the origin', () => {
        const { vertices, indices } = makeBox(-1, -1, -1, 1, 1, 1);
        const result = computeMassProperties(vertices, 12, indices, true);

        expect(result.mass).toBeCloseTo(8, 12);
        expect(result.center.values[0]).toBeCloseTo(0, 12);
        expect(result.center.values[1]).toBeCloseTo(0, 12);
        expect(result.center.values[2]).toBeCloseTo(0, 12);
        // m*(w^2 + h^2)/12 = 8*(4 + 4)/12 = 16/3.
        expectMatrixClose(result.inertia, [
            [16 / 3, 0, 0],
            [0, 16 / 3, 0],
            [0, 0, 16 / 3]
        ], 12);
    });

    it('computes a non-cubic box against the analytic rigid-body formula', () => {
        const a = 2, b = 3, c = 5;
        const { vertices, indices } = makeBox(0, 0, 0, a, b, c);
        const result = computeMassProperties(vertices, 12, indices, true);
        const m = a * b * c;

        expect(result.mass).toBeCloseTo(m, 10);
        expect(result.center.values[0]).toBeCloseTo(a / 2, 12);
        expect(result.center.values[1]).toBeCloseTo(b / 2, 12);
        expect(result.center.values[2]).toBeCloseTo(c / 2, 12);
        expectMatrixClose(result.inertia, [
            [m * (b * b + c * c) / 12, 0, 0],
            [0, m * (c * c + a * a) / 12, 0],
            [0, 0, m * (a * a + b * b) / 12]
        ], 9);
    });

    it('computes the corner tetrahedron against analytic integrals', () => {
        const { vertices, indices } = makeTetrahedron(1);
        const world = computeMassProperties(vertices, 4, indices, false);

        expect(world.mass).toBeCloseTo(1 / 6, 12);
        for (let i = 0; i < 3; ++i) {
            expect(world.center.values[i]).toBeCloseTo(0.25, 12);
        }

        // integral of x^2 = integral of y^2 = integral of z^2 = 1/60, so the
        // diagonal is 1/30. integral of x*y = 1/120, so the off-diagonal
        // entries are -1/120.
        expectMatrixClose(world.inertia, [
            [1 / 30, -1 / 120, -1 / 120],
            [-1 / 120, 1 / 30, -1 / 120],
            [-1 / 120, -1 / 120, 1 / 30]
        ], 12);

        const body = computeMassProperties(vertices, 4, indices, true);
        // 1/30 - (1/6)*(1/16 + 1/16) = 1/80 and -1/120 + (1/6)*(1/16) = 1/480.
        expectMatrixClose(body.inertia, [
            [1 / 80, 1 / 480, 1 / 480],
            [1 / 480, 1 / 80, 1 / 480],
            [1 / 480, 1 / 480, 1 / 80]
        ], 12);
    });

    it('scales mass as s^3 and inertia as s^5 under a uniform scale', () => {
        const s = 3;
        const unit = computeMassProperties(makeTetrahedron(1).vertices, 4,
            makeTetrahedron(1).indices, true);
        const scaled = computeMassProperties(makeTetrahedron(s).vertices, 4,
            makeTetrahedron(s).indices, true);

        expect(scaled.mass).toBeCloseTo(unit.mass * s * s * s, 10);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(scaled.inertia.get(r, c)).toBeCloseTo(
                    unit.inertia.get(r, c) * Math.pow(s, 5), 9);
            }
        }
    });

    it('is translation invariant in body coordinates and shifts the center', () => {
        const t = v3(-7.25, 3.5, 11);
        const base = makeBox(0, 0, 0, 1, 2, 3);
        const shifted = { vertices: translate(base.vertices, t), indices: base.indices };

        const b0 = computeMassProperties(base.vertices, 12, base.indices, true);
        const b1 = computeMassProperties(shifted.vertices, 12, shifted.indices, true);

        expect(b1.mass).toBeCloseTo(b0.mass, 10);
        for (let i = 0; i < 3; ++i) {
            expect(b1.center.values[i]).toBeCloseTo(b0.center.values[i] + t.values[i], 10);
        }
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(b1.inertia.get(r, c)).toBeCloseTo(b0.inertia.get(r, c), 8);
            }
        }
    });

    it('produces a symmetric inertia tensor and negative mass for inverted winding', () => {
        const { vertices, indices } = makeBox(0, 0, 0, 1, 2, 3);
        const outward = computeMassProperties(vertices, 12, indices, false);
        expect(outward.inertia.get(0, 1)).toBe(outward.inertia.get(1, 0));
        expect(outward.inertia.get(0, 2)).toBe(outward.inertia.get(2, 0));
        expect(outward.inertia.get(1, 2)).toBe(outward.inertia.get(2, 1));

        // Reversing the winding of every triangle negates every normal, hence
        // negates every integral.
        const reversed: number[] = [];
        for (let t = 0; t < 12; ++t) {
            reversed.push(indices[3 * t], indices[3 * t + 2], indices[3 * t + 1]);
        }
        const inward = computeMassProperties(vertices, 12, reversed, false);
        expect(inward.mass).toBeCloseTo(-outward.mass, 10);
        // The center is unchanged because both numerator and denominator flip.
        for (let i = 0; i < 3; ++i) {
            expect(inward.center.values[i]).toBeCloseTo(outward.center.values[i], 10);
        }
    });

    it('matches a Monte Carlo estimate of the integrals over a random tetrahedron',
        () => {
            // A random-but-fixed tetrahedron with outward normals; the winding
            // is chosen so that the mass is positive.
            const p0 = v3(0.3, -0.7, 0.2);
            const p1 = v3(1.9, 0.1, -0.4);
            const p2 = v3(-0.5, 1.6, 0.3);
            const p3 = v3(0.1, 0.2, 2.3);
            const vertices = [p0, p1, p2, p3];
            const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
            const result = computeMassProperties(vertices, 4, indices, false);

            // Exact volume: |det(p1-p0, p2-p0, p3-p0)| / 6.
            const e1 = [p1.values[0] - p0.values[0], p1.values[1] - p0.values[1],
                p1.values[2] - p0.values[2]];
            const e2 = [p2.values[0] - p0.values[0], p2.values[1] - p0.values[1],
                p2.values[2] - p0.values[2]];
            const e3 = [p3.values[0] - p0.values[0], p3.values[1] - p0.values[1],
                p3.values[2] - p0.values[2]];
            const det =
                e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) -
                e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) +
                e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
            expect(result.mass).toBeGreaterThan(0);
            expect(result.mass).toBeCloseTo(Math.abs(det) / 6, 12);

            // The centroid of a tetrahedron is the average of its vertices.
            for (let i = 0; i < 3; ++i) {
                const avg = (p0.values[i] + p1.values[i] + p2.values[i] + p3.values[i]) / 4;
                expect(result.center.values[i]).toBeCloseTo(avg, 12);
            }

            // Monte Carlo over the tetrahedron using barycentric sampling.
            let seed = 0x9e3779b9;
            const rand = (): number => {
                seed = (seed + 0x6d2b79f5) >>> 0;
                let t = seed;
                t = Math.imul(t ^ (t >>> 15), t | 1);
                t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
                return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
            };
            const n = 200000;
            let sx = 0, sy = 0, sz = 0, syy = 0, szz = 0, sxy = 0;
            for (let k = 0; k < n; ++k) {
                let a = rand(), b = rand(), c = rand();
                // Sample uniformly in the tetrahedron (fold the cube).
                if (a + b > 1) { a = 1 - a; b = 1 - b; }
                if (b + c > 1) { const t = c; c = 1 - a - b; b = 1 - t; }
                else if (a + b + c > 1) { const t = c; c = a + b + c - 1; a = 1 - b - t; }
                const d = 1 - a - b - c;
                const x = d * p0.values[0] + a * p1.values[0] + b * p2.values[0] + c * p3.values[0];
                const y = d * p0.values[1] + a * p1.values[1] + b * p2.values[1] + c * p3.values[1];
                const z = d * p0.values[2] + a * p1.values[2] + b * p2.values[2] + c * p3.values[2];
                sx += x; sy += y; sz += z;
                syy += y * y; szz += z * z; sxy += x * y;
            }
            const vol = result.mass;
            expect(result.center.values[0]).toBeCloseTo(sx / n, 2);
            expect(result.center.values[1]).toBeCloseTo(sy / n, 2);
            expect(result.center.values[2]).toBeCloseTo(sz / n, 2);
            expect(result.inertia.get(0, 0)).toBeCloseTo(vol * (syy + szz) / n, 1);
            expect(result.inertia.get(0, 1)).toBeCloseTo(-vol * sxy / n, 1);
        });
});

// ---------------------------------------------------------------------------
// Verification (V41): the divergence-theorem integrals are cross-checked
// against the closed-form moments of a tetrahedron, which are independent of
// the algorithm under test.
// ---------------------------------------------------------------------------

const pmpCoord = () => scaledArb(-3, 3, 512);

function det3(a: number[], b: number[], c: number[]): number {
    return a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);
}

// Four points spanning a tetrahedron of comfortable volume, ordered so that
// det(p1-p0, p2-p0, p3-p0) > 0.
const tetrahedronPoints = (): fc.Arbitrary<number[][]> =>
    fc.array(fc.array(pmpCoord(), { minLength: 3, maxLength: 3 }),
        { minLength: 4, maxLength: 4 })
        .filter(p => {
            const e = (k: number) => [p[k][0] - p[0][0], p[k][1] - p[0][1],
                p[k][2] - p[0][2]];
            return Math.abs(det3(e(1), e(2), e(3))) > 0.5;
        })
        .map(p => {
            const e = (k: number) => [p[k][0] - p[0][0], p[k][1] - p[0][1],
                p[k][2] - p[0][2]];
            return det3(e(1), e(2), e(3)) > 0
                ? p : [p[0], p[2], p[1], p[3]];
        });

// The four outward-facing triangles of a positively oriented tetrahedron
// p0, p1, p2, p3 (see the makeTetrahedron helper above for the canonical
// instance).
const TETRAHEDRON_FACES = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];

function tetraMesh(p: number[][]): { vertices: Vector[], indices: number[] } {
    return {
        vertices: p.map(a => v3(a[0], a[1], a[2])),
        indices: TETRAHEDRON_FACES.slice()
    };
}

// The exact moments of a solid tetrahedron with vertices p0..p3. Using
// barycentric coordinates, integral(lambda_a * lambda_b) = V*(1 + delta_ab)/20
// and x = sum_a lambda_a p_a, hence
//   integral(x_i) = V * (sum_a p_a,i) / 4
//   integral(x_i x_j) = (V/20) * [ (sum_a p_a,i)(sum_b p_b,j)
//                                  + sum_a p_a,i p_a,j ].
function tetraMoments(p: number[][]):
    { volume: number, first: number[], second: number[][] } {
    const e = (k: number) => [p[k][0] - p[0][0], p[k][1] - p[0][1],
        p[k][2] - p[0][2]];
    const volume = det3(e(1), e(2), e(3)) / 6;
    const s = [0, 0, 0];
    for (const q of p) { for (let i = 0; i < 3; ++i) { s[i] += q[i]; } }
    const first = s.map(x => (volume * x) / 4);
    const second = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; ++i) {
        for (let j = 0; j < 3; ++j) {
            let sameVertex = 0;
            for (const q of p) { sameVertex += q[i] * q[j]; }
            second[i][j] = (volume / 20) * (s[i] * s[j] + sameVertex);
        }
    }
    return { volume, first, second };
}

// Recover the six second moments from the world-frame inertia tensor:
// I(0,0) = int(y^2 + z^2), I(0,1) = -int(x*y), and so on.
function secondMomentsFromInertia(J: Matrix): number[][] {
    const d0 = J.get(0, 0), d1 = J.get(1, 1), d2 = J.get(2, 2);
    const xx = (d1 + d2 - d0) / 2;
    const yy = (d0 + d2 - d1) / 2;
    const zz = (d0 + d1 - d2) / 2;
    return [
        [xx, -J.get(0, 1), -J.get(0, 2)],
        [-J.get(0, 1), yy, -J.get(1, 2)],
        [-J.get(0, 2), -J.get(1, 2), zz]
    ];
}

function maxMatrixDiff(A: Matrix, B: Matrix): number {
    let m = 0;
    for (let i = 0; i < A.numElements; ++i) {
        m = Math.max(m, Math.abs(A.values[i] - B.values[i]));
    }
    return m;
}

describe('PolyhedralMassProperties verification', () => {
    it('reproduces the closed-form volume, centroid and second moments of a ' +
        'tetrahedron', () => {
            check(tetrahedronPoints(), (p) => {
                const { vertices, indices } = tetraMesh(p);
                const result = computeMassProperties(vertices, 4, indices,
                    false);
                const exact = tetraMoments(p);

                expectClose(result.mass, exact.volume, 1e-12, 1e-11);
                for (let i = 0; i < 3; ++i) {
                    expectClose(result.center.values[i],
                        exact.first[i] / exact.volume, 1e-11, 1e-10);
                }
                const moments = secondMomentsFromInertia(result.inertia);
                for (let i = 0; i < 3; ++i) {
                    for (let j = 0; j < 3; ++j) {
                        // The tensor entries are differences of O(|p|^5)
                        // integrals, so the tolerance is scale relative.
                        expectClose(moments[i][j], exact.second[i][j],
                            1e-10, 1e-9);
                    }
                }
            }, 150);
        });

    it('the bodyCoords tensor is the world tensor of the centered mesh', () => {
        check(tetrahedronPoints(), (p) => {
            const { vertices, indices } = tetraMesh(p);
            const body = computeMassProperties(vertices, 4, indices, true);
            const centered = vertices.map(v => v3(
                v.values[0] - body.center.values[0],
                v.values[1] - body.center.values[1],
                v.values[2] - body.center.values[2]));
            const world = computeMassProperties(centered, 4, indices, false);
            // The parallel-axis correction and an explicit recentering must
            // agree; the entries are O(|p|^5) so the tolerance is relative.
            let scale = 1;
            for (let i = 0; i < 9; ++i) {
                scale = Math.max(scale, Math.abs(world.inertia.values[i]));
            }
            expect(maxMatrixDiff(body.inertia, world.inertia))
                .toBeLessThan(1e-9 * scale);
        }, 150);
    });

    it('is additive over a disjoint union of tetrahedra', () => {
        check(fc.tuple(tetrahedronPoints(), tetrahedronPoints(),
            fc.array(scaledArb(20, 40, 32), { minLength: 3, maxLength: 3 })),
            ([pa, pbRaw, offset]) => {
                // Move the second tetrahedron far away so the two solids are
                // disjoint; the integrals are additive regardless, but this
                // keeps the union a genuine polyhedron.
                const pb = pbRaw.map(q => [q[0] + offset[0], q[1] + offset[1],
                    q[2] + offset[2]]);
                const a = tetraMesh(pa);
                const b = tetraMesh(pb);
                const union = {
                    vertices: a.vertices.concat(b.vertices),
                    indices: a.indices.concat(b.indices.map(i => i + 4))
                };
                const ra = computeMassProperties(a.vertices, 4, a.indices,
                    false);
                const rb = computeMassProperties(b.vertices, 4, b.indices,
                    false);
                const ru = computeMassProperties(union.vertices, 8,
                    union.indices, false);

                expectClose(ru.mass, ra.mass + rb.mass, 1e-11, 1e-10);
                let scale = 1;
                for (let i = 0; i < 9; ++i) {
                    scale = Math.max(scale, Math.abs(ru.inertia.values[i]));
                }
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        expect(Math.abs(ru.inertia.get(r, c)
                            - (ra.inertia.get(r, c) + rb.inertia.get(r, c))))
                            .toBeLessThan(1e-10 * scale);
                    }
                }
                for (let i = 0; i < 3; ++i) {
                    expectClose(ru.center.values[i],
                        (ra.mass * ra.center.values[i]
                            + rb.mass * rb.center.values[i])
                        / (ra.mass + rb.mass), 1e-10, 1e-9);
                }
            }, 100);
    });

    it('negates the mass and the world tensor when the winding is reversed, ' +
        'leaving the center of mass unchanged', () => {
            check(tetrahedronPoints(), (p) => {
                const { vertices, indices } = tetraMesh(p);
                const forward = computeMassProperties(vertices, 4, indices,
                    false);
                const reversed: number[] = [];
                for (let t = 0; t < 4; ++t) {
                    reversed.push(indices[3 * t], indices[3 * t + 2],
                        indices[3 * t + 1]);
                }
                const backward = computeMassProperties(vertices, 4, reversed,
                    false);
                expectClose(backward.mass, -forward.mass, 1e-13, 1e-12);
                let scale = 1;
                for (let i = 0; i < 9; ++i) {
                    scale = Math.max(scale,
                        Math.abs(forward.inertia.values[i]));
                }
                for (let i = 0; i < 9; ++i) {
                    expect(Math.abs(backward.inertia.values[i]
                        + forward.inertia.values[i]))
                        .toBeLessThan(1e-11 * scale);
                }
                for (let i = 0; i < 3; ++i) {
                    // center = (integral x)/mass, and both numerator and
                    // denominator flip sign.
                    expectClose(backward.center.values[i],
                        forward.center.values[i], 1e-11, 1e-10);
                }
            }, 100);
        });

    it('transforms covariantly under a rotation of the mesh', () => {
        const angle = () => scaledArb(-Math.PI, Math.PI, 512);
        check(fc.tuple(tetrahedronPoints(), angle(), angle(), angle()),
            ([p, a0, a1, a2]) => {
                const c0 = Math.cos(a0), s0 = Math.sin(a0);
                const c1 = Math.cos(a1), s1 = Math.sin(a1);
                const c2 = Math.cos(a2), s2 = Math.sin(a2);
                const R = Matrix.fromArray(3, 3, [
                    c0, -s0, 0, s0, c0, 0, 0, 0, 1]);
                const Ry = Matrix.fromArray(3, 3, [
                    c1, 0, s1, 0, 1, 0, -s1, 0, c1]);
                const Rx = Matrix.fromArray(3, 3, [
                    1, 0, 0, 0, c2, -s2, 0, s2, c2]);
                const rot = multiplyAB(multiplyAB(R, Ry), Rx);

                const base = tetraMesh(p);
                const rotatedVertices = base.vertices.map(v =>
                    mulMatrix(rot, v));
                const a = computeMassProperties(base.vertices, 4,
                    base.indices, true);
                const b = computeMassProperties(rotatedVertices, 4,
                    base.indices, true);

                expectClose(b.mass, a.mass, 1e-12, 1e-11);
                const rotatedCenter = mulMatrix(rot, a.center);
                for (let i = 0; i < 3; ++i) {
                    expectClose(b.center.values[i], rotatedCenter.values[i],
                        1e-10, 1e-9);
                }
                const expected = multiplyABT(multiplyAB(rot, a.inertia), rot);
                let scale = 1;
                for (let i = 0; i < 9; ++i) {
                    scale = Math.max(scale, Math.abs(expected.values[i]));
                }
                expect(maxMatrixDiff(b.inertia, expected))
                    .toBeLessThan(1e-9 * scale);
            }, 100);
    });

    it('reproduces the closed form for an axis-aligned box', () => {
        const arb = fc.tuple(
            fc.array(scaledArb(-3, 3, 128), { minLength: 3, maxLength: 3 }),
            fc.array(scaledArb(0.2, 4, 128), { minLength: 3, maxLength: 3 }));
        check(arb, ([lo, extent]) => {
            const { vertices, indices } = makeBox(lo[0], lo[1], lo[2],
                lo[0] + extent[0], lo[1] + extent[1], lo[2] + extent[2]);
            const result = computeMassProperties(vertices, 12, indices, true);
            const m = extent[0] * extent[1] * extent[2];
            expectClose(result.mass, m, 1e-12, 1e-11);
            for (let i = 0; i < 3; ++i) {
                expectClose(result.center.values[i], lo[i] + extent[i] / 2,
                    1e-11, 1e-10);
            }
            // The body-frame tensor of a box is diagonal with entries
            // m*(b^2 + c^2)/12.
            for (let r = 0; r < 3; ++r) {
                for (let c = 0; c < 3; ++c) {
                    const expected = (r === c)
                        ? (m * (extent[(r + 1) % 3] ** 2
                            + extent[(r + 2) % 3] ** 2)) / 12
                        : 0;
                    expectClose(result.inertia.get(r, c), expected,
                        1e-9 * Math.max(1, m * 25), 1e-9);
                }
            }
        }, 150);
    });
});
