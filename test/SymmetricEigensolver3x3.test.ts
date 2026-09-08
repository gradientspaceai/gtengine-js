import { describe, it, expect } from 'vitest';
import {
    SymmetricEigensolver3x3,
    NISymmetricEigensolver3x3,
    SortEigenstuff,
    type EigenTriple,
    type EigenBasis3
} from '../src/SymmetricEigensolver3x3.js';
import { check, expectClose, fc, wellScaled } from './helpers/arbitraries.js';

// A simple deterministic pseudorandom generator so test runs are repeatable.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

interface Sym3 {
    a00: number; a01: number; a02: number;
    a11: number; a12: number; a22: number;
}

function randomSym3(random: () => number): Sym3 {
    return {
        a00: 4 * random() - 2,
        a01: 4 * random() - 2,
        a02: 4 * random() - 2,
        a11: 4 * random() - 2,
        a12: 4 * random() - 2,
        a22: 4 * random() - 2
    };
}

function matVec(m: Sym3, v: EigenTriple): EigenTriple {
    return [
        m.a00 * v[0] + m.a01 * v[1] + m.a02 * v[2],
        m.a01 * v[0] + m.a11 * v[1] + m.a12 * v[2],
        m.a02 * v[0] + m.a12 * v[1] + m.a22 * v[2]
    ];
}

function dot(u: EigenTriple, v: EigenTriple): number {
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

function det3(evecs: EigenBasis3): number {
    return evecs[0][0] * (evecs[1][1] * evecs[2][2] - evecs[1][2] * evecs[2][1])
        - evecs[0][1] * (evecs[1][0] * evecs[2][2] - evecs[1][2] * evecs[2][0])
        + evecs[0][2] * (evecs[1][0] * evecs[2][1] - evecs[1][1] * evecs[2][0]);
}

// Verify A*v = lambda*v for each pair, orthonormality of the basis and,
// when sorted, the eigenvalue ordering and right-handedness.
function checkEigensystem(m: Sym3, evals: EigenTriple, evecs: EigenBasis3,
    sortType: number, tol: number): void {
    for (let i = 0; i < 3; ++i) {
        const av = matVec(m, evecs[i]);
        for (let j = 0; j < 3; ++j) {
            expect(Math.abs(av[j] - evals[i] * evecs[i][j])).toBeLessThanOrEqual(tol);
        }
    }

    for (let i = 0; i < 3; ++i) {
        for (let j = i; j < 3; ++j) {
            const expected = (i === j ? 1 : 0);
            expect(Math.abs(dot(evecs[i], evecs[j]) - expected)).toBeLessThanOrEqual(1e-12);
        }
    }

    if (sortType === 1) {
        expect(evals[0]).toBeLessThanOrEqual(evals[1]);
        expect(evals[1]).toBeLessThanOrEqual(evals[2]);
    }
    else if (sortType === -1) {
        expect(evals[0]).toBeGreaterThanOrEqual(evals[1]);
        expect(evals[1]).toBeGreaterThanOrEqual(evals[2]);
    }

    if (sortType !== 0) {
        // {evecs[0], evecs[1], evecs[2]} must be right-handed when sorted.
        expect(Math.abs(det3(evecs) - 1)).toBeLessThanOrEqual(1e-12);
    }
}

// Build A = R * diag(d) * R^T from a rotation R specified by axis and angle.
function fromEigensystem(d: EigenTriple, axis: EigenTriple, angle: number): Sym3 {
    const len = Math.sqrt(dot(axis, axis));
    const u: EigenTriple = [axis[0] / len, axis[1] / len, axis[2] / len];
    const c = Math.cos(angle), s = Math.sin(angle), omc = 1 - c;
    // Rotation matrix R (rows).
    const r: number[][] = [
        [c + u[0] * u[0] * omc, u[0] * u[1] * omc - u[2] * s, u[0] * u[2] * omc + u[1] * s],
        [u[1] * u[0] * omc + u[2] * s, c + u[1] * u[1] * omc, u[1] * u[2] * omc - u[0] * s],
        [u[2] * u[0] * omc - u[1] * s, u[2] * u[1] * omc + u[0] * s, c + u[2] * u[2] * omc]
    ];
    // A[i][j] = sum_k R[i][k] * d[k] * R[j][k].
    const a = (i: number, j: number): number =>
        r[i][0] * d[0] * r[j][0] + r[i][1] * d[1] * r[j][1] + r[i][2] * d[2] * r[j][2];
    return { a00: a(0, 0), a01: a(0, 1), a02: a(0, 2), a11: a(1, 1), a12: a(1, 2), a22: a(2, 2) };
}

describe('SymmetricEigensolver3x3 (iterative)', () => {
    const solver = new SymmetricEigensolver3x3();

    it('solves randomized symmetric matrices for all sort types', () => {
        const random = makeRandom(3030);
        for (let trial = 0; trial < 100; ++trial) {
            const m = randomSym3(random);
            for (const sortType of [-1, 0, 1]) {
                for (const aggressive of [false, true]) {
                    const { iterations, evals, evecs } =
                        solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, aggressive, sortType);
                    expect(iterations).toBeGreaterThanOrEqual(0);
                    checkEigensystem(m, evals, evecs, sortType, 1e-12);
                }
            }
        }
    });

    it('reproduces known eigenvalues of a constructed matrix', () => {
        const d: EigenTriple = [-2, 1, 4];
        const m = fromEigensystem(d, [1, 2, 3], 0.7);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, false, +1);
        expect(evals[0]).toBeCloseTo(-2, 12);
        expect(evals[1]).toBeCloseTo(1, 12);
        expect(evals[2]).toBeCloseTo(4, 12);
        checkEigensystem(m, evals, evecs, +1, 1e-12);
    });

    it('handles a diagonal matrix', () => {
        const m: Sym3 = { a00: 3, a01: 0, a02: 0, a11: -5, a12: 0, a22: 1 };
        const { evals, evecs } = solver.solve(3, 0, 0, -5, 0, 1, false, +1);
        expect(evals[0]).toBeCloseTo(-5, 14);
        expect(evals[1]).toBeCloseTo(1, 14);
        expect(evals[2]).toBeCloseTo(3, 14);
        checkEigensystem(m, evals, evecs, +1, 1e-13);
    });

    it('handles the zero matrix', () => {
        const m: Sym3 = { a00: 0, a01: 0, a02: 0, a11: 0, a12: 0, a22: 0 };
        const { evals, evecs } = solver.solve(0, 0, 0, 0, 0, 0, false, +1);
        expect(evals).toEqual([0, 0, 0]);
        checkEigensystem(m, evals, evecs, +1, 1e-14);
    });

    it('handles repeated eigenvalues', () => {
        // A with eigenvalues {2, 2, 5}.
        const m = fromEigensystem([2, 2, 5], [1, -1, 2], 1.1);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, false, +1);
        expect(evals[0]).toBeCloseTo(2, 12);
        expect(evals[1]).toBeCloseTo(2, 12);
        expect(evals[2]).toBeCloseTo(5, 12);
        checkEigensystem(m, evals, evecs, +1, 1e-12);
    });

    it('handles a triple eigenvalue (multiple of identity)', () => {
        const m: Sym3 = { a00: 7, a01: 0, a02: 0, a11: 7, a12: 0, a22: 7 };
        const { evals, evecs } = solver.solve(7, 0, 0, 7, 0, 7, false, +1);
        expect(evals[0]).toBeCloseTo(7, 14);
        expect(evals[1]).toBeCloseTo(7, 14);
        expect(evals[2]).toBeCloseTo(7, 14);
        checkEigensystem(m, evals, evecs, +1, 1e-13);
    });
});

describe('NISymmetricEigensolver3x3 (noniterative)', () => {
    const solver = new NISymmetricEigensolver3x3();

    it('solves randomized symmetric matrices for all sort types', () => {
        const random = makeRandom(4040);
        for (let trial = 0; trial < 100; ++trial) {
            const m = randomSym3(random);
            for (const sortType of [-1, 0, 1]) {
                const { evals, evecs } =
                    solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, sortType);
                checkEigensystem(m, evals, evecs, sortType, 1e-10);
            }
        }
    });

    it('agrees with the iterative solver on random matrices', () => {
        const random = makeRandom(5050);
        const iterative = new SymmetricEigensolver3x3();
        for (let trial = 0; trial < 50; ++trial) {
            const m = randomSym3(random);
            const ni = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, +1);
            const it2 = iterative.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, false, +1);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(ni.evals[i] - it2.evals[i])).toBeLessThanOrEqual(1e-10);
            }
        }
    });

    it('reproduces known eigenvalues of a constructed matrix', () => {
        const m = fromEigensystem([-3, 0.5, 6], [2, 1, -1], -0.4);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, +1);
        expect(evals[0]).toBeCloseTo(-3, 11);
        expect(evals[1]).toBeCloseTo(0.5, 11);
        expect(evals[2]).toBeCloseTo(6, 11);
        checkEigensystem(m, evals, evecs, +1, 1e-11);
    });

    it('handles a diagonal matrix', () => {
        const m: Sym3 = { a00: 3, a01: 0, a02: 0, a11: -5, a12: 0, a22: 1 };
        const { evals, evecs } = solver.solve(3, 0, 0, -5, 0, 1, +1);
        expect(evals[0]).toBeCloseTo(-5, 14);
        expect(evals[1]).toBeCloseTo(1, 14);
        expect(evals[2]).toBeCloseTo(3, 14);
        checkEigensystem(m, evals, evecs, +1, 1e-13);
    });

    it('handles the zero matrix', () => {
        const m: Sym3 = { a00: 0, a01: 0, a02: 0, a11: 0, a12: 0, a22: 0 };
        const { evals, evecs } = solver.solve(0, 0, 0, 0, 0, 0, +1);
        expect(evals).toEqual([0, 0, 0]);
        checkEigensystem(m, evals, evecs, +1, 1e-14);
    });

    it('handles repeated eigenvalues', () => {
        // Repeated eigenvalues are computed to about sqrt(epsilon) accuracy
        // by the noniterative trigonometric formulas; the eigen equation
        // residual is of the same order.
        const m = fromEigensystem([2, 2, 5], [1, -1, 2], 1.1);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, +1);
        expect(evals[0]).toBeCloseTo(2, 7);
        expect(evals[1]).toBeCloseTo(2, 7);
        expect(evals[2]).toBeCloseTo(5, 11);
        checkEigensystem(m, evals, evecs, +1, 1e-7);
    });

    it('handles decreasing sort order', () => {
        const m = fromEigensystem([-1, 2, 9], [0, 1, 1], 0.3);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, -1);
        expect(evals[0]).toBeCloseTo(9, 11);
        expect(evals[1]).toBeCloseTo(2, 11);
        expect(evals[2]).toBeCloseTo(-1, 11);
        checkEigensystem(m, evals, evecs, -1, 1e-11);
    });
});

describe('SortEigenstuff', () => {
    it('sorts eigenvalues increasingly and keeps eigenvector pairing', () => {
        const evals: EigenTriple = [3, 1, 2];
        const evecs: EigenBasis3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        new SortEigenstuff().sort(+1, true, evals, evecs);
        expect(evals).toEqual([1, 2, 3]);
        // Eigenvector rows follow their eigenvalues (up to sign of row 2,
        // which may be negated to keep the basis right-handed).
        expect(evecs[0]).toEqual([0, 1, 0]);
        expect(evecs[1]).toEqual([0, 0, 1]);
        expect(Math.abs(evecs[2][0])).toBe(1);
        expect(Math.abs(det3(evecs) - 1)).toBeLessThanOrEqual(1e-15);
    });

    it('sorts eigenvalues decreasingly and preserves right-handedness', () => {
        const evals: EigenTriple = [3, 1, 2];
        const evecs: EigenBasis3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        new SortEigenstuff().sort(-1, true, evals, evecs);
        expect(evals).toEqual([3, 2, 1]);
        expect(Math.abs(det3(evecs) - 1)).toBeLessThanOrEqual(1e-15);
    });

    it('with sortType 0 only enforces right-handedness', () => {
        const evals: EigenTriple = [3, 1, 2];
        const evecs: EigenBasis3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        new SortEigenstuff().sort(0, false, evals, evecs);
        expect(evals).toEqual([3, 1, 2]);
        expect(evecs[2]).toEqual([-0, -0, -1]);
    });
});

// ---------------------------------------------------------------------------
// Verification (group V38).
// ---------------------------------------------------------------------------
describe('SymmetricEigensolver3x3 verification', () => {
    const solver = new SymmetricEigensolver3x3();
    const niSolver = new NISymmetricEigensolver3x3();

    // wellScaled snaps |x| < 1e-3 to exactly zero, so no entry is subnormal
    // and the squares formed inside the solver stay in the normal range.
    const sym3 = fc.tuple(wellScaled(-10, 10), wellScaled(-10, 10),
        wellScaled(-10, 10), wellScaled(-10, 10), wellScaled(-10, 10),
        wellScaled(-10, 10)).map(([a00, a01, a02, a11, a12, a22]): Sym3 =>
            ({ a00, a01, a02, a11, a12, a22 }));

    function normOf(m: Sym3): number {
        return Math.max(Math.abs(m.a00), Math.abs(m.a01), Math.abs(m.a02),
            Math.abs(m.a11), Math.abs(m.a12), Math.abs(m.a22), 1);
    }

    function traceOf(m: Sym3): number {
        return m.a00 + m.a11 + m.a22;
    }

    function detOf(m: Sym3): number {
        return m.a00 * (m.a11 * m.a22 - m.a12 * m.a12)
            - m.a01 * (m.a01 * m.a22 - m.a12 * m.a02)
            + m.a02 * (m.a01 * m.a12 - m.a11 * m.a02);
    }

    it('the iterative solver produces an eigensystem for every flag combination',
        () => {
            check(fc.tuple(sym3, fc.constantFrom(-1, 0, 1), fc.boolean()),
                ([m, sortType, aggressive]) => {
                    const { iterations, evals, evecs } = solver.solve(m.a00,
                        m.a01, m.a02, m.a11, m.a12, m.a22, aggressive, sortType);
                    expect(iterations).toBeGreaterThanOrEqual(0);
                    checkEigensystem(m, evals, evecs, sortType,
                        1e-12 * normOf(m));
                });
        });

    it('the iterative solver preserves the trace and the determinant', () => {
        check(sym3, m => {
            const { evals } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12,
                m.a22, false, +1);
            const norm = normOf(m);
            expectClose(evals[0] + evals[1] + evals[2], traceOf(m),
                1e-12 * norm, 1e-12);
            expectClose(evals[0] * evals[1] * evals[2], detOf(m),
                1e-10 * norm * norm * norm, 1e-10);
        });
    });

    it('the noniterative solver agrees with the iterative one', () => {
        check(fc.tuple(sym3, fc.constantFrom(-1, 0, 1)), ([m, sortType]) => {
            const it3 = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22,
                false, sortType);
            const ni = niSolver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22,
                sortType);
            checkEigensystem(m, ni.evals, ni.evecs, sortType, 1e-7 * normOf(m));
            // With sortType 0 neither solver promises an order, so compare
            // the eigenvalues as sets. The two algorithms are independent (a
            // QR-style iteration versus the closed-form cubic of the PDF), so
            // they agree only to the conditioning of the eigenvalues: the
            // closed form evaluates acos of a value clamped to [-1,1], which
            // costs about half the mantissa for nearly equal eigenvalues.
            const aIt = it3.evals.slice().sort((x, y) => x - y);
            const aNi = ni.evals.slice().sort((x, y) => x - y);
            for (let i = 0; i < 3; ++i) {
                expectClose(aNi[i], aIt[i], 1e-7 * normOf(m), 1e-7);
            }
        });
    });

    it('both solvers report the diagonal entries of a diagonal matrix', () => {
        check(fc.tuple(wellScaled(-10, 10), wellScaled(-10, 10),
            wellScaled(-10, 10)), ([d0, d1, d2]) => {
                const sorted = [d0, d1, d2].slice().sort((x, y) => x - y);
                // The noniterative solver takes its 'norm == 0' branch here
                // and returns the diagonal entries with no arithmetic at all.
                const ni = niSolver.solve(d0, 0, 0, d1, 0, d2, +1);
                for (let i = 0; i < 3; ++i) {
                    // Not bit-exact: upstream divides the matrix by its
                    // largest absolute entry and multiplies the eigenvalues
                    // back afterwards, and d/max*max need not round-trip.
                    expectClose(ni.evals[i], sorted[i], 0, 4 * Number.EPSILON);
                }
                expectClose(Math.abs(det3(ni.evecs)), 1, 1e-15, 0);

                const itr = solver.solve(d0, 0, 0, d1, 0, d2, false, +1);
                for (let i = 0; i < 3; ++i) {
                    expectClose(itr.evals[i], sorted[i], 1e-13, 1e-13);
                }
            });
    });

    it('the zero matrix gives zero eigenvalues and an orthonormal basis', () => {
        for (const sortType of [-1, 0, 1]) {
            const ni = niSolver.solve(0, 0, 0, 0, 0, 0, sortType);
            expect(ni.evals).toEqual([0, 0, 0]);
            expect(ni.evecs).toEqual([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);

            const it3 = solver.solve(0, 0, 0, 0, 0, 0, false, sortType);
            for (let i = 0; i < 3; ++i) {
                expect(it3.evals[i] + 0).toBe(0);
                expectClose(Math.sqrt(dot(it3.evecs[i], it3.evecs[i])), 1,
                    1e-15, 0);
            }
            if (sortType !== 0) {
                expectClose(det3(it3.evecs), 1, 1e-14, 0);
            }
        }
    });

    it('SortEigenstuff permutes eigenpairs together and keeps a right-handed basis',
        () => {
            // Cross-check against an independent sort of (value, vector) pairs.
            check(fc.tuple(fc.tuple(wellScaled(-10, 10), wellScaled(-10, 10),
                wellScaled(-10, 10)), fc.constantFrom(-1, 0, 1), fc.boolean()),
                ([[d0, d1, d2], sortType, flip]) => {
                    const evals: EigenTriple = [d0, d1, d2];
                    // An orthonormal basis of known handedness. 'isRotation'
                    // tells SortEigenstuff whether the incoming basis is a
                    // rotation; it must be the truth about that basis, since
                    // the routine uses it to decide whether to negate a row.
                    const evecs: EigenBasis3 = flip
                        ? [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]
                        : [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
                    new SortEigenstuff().sort(sortType, !flip, evals, evecs);

                    const expected = [d0, d1, d2].slice();
                    if (sortType > 0) { expected.sort((x, y) => x - y); }
                    else if (sortType < 0) { expected.sort((x, y) => y - x); }
                    for (let i = 0; i < 3; ++i) {
                        expect(evals[i] + 0).toBe(expected[i] + 0);
                    }
                    // Each row is still a signed axis and the basis stays
                    // right-handed, which is what the isRotation bookkeeping
                    // is for.
                    for (let i = 0; i < 3; ++i) {
                        expectClose(dot(evecs[i], evecs[i]), 1, 0, 0);
                    }
                    expectClose(det3(evecs), 1, 0, 0);
                    // Row i is +-e_j where d[j] is the eigenvalue now in slot i.
                    for (let i = 0; i < 3; ++i) {
                        const j = evecs[i].findIndex(x => x !== 0);
                        expect([d0, d1, d2][j] + 0).toBe(evals[i] + 0);
                    }
                });
        });

    it('regression (#379): the iterative solver is scale equivariant over the whole exponent range',
        () => {
            // Upstream's GetCosSin computes sqrt(u*u + v*v) with no rescaling,
            // unlike the 2x2 sibling. Before the port fix the squares overflow
            // above 2^511 (the returned (c,s) is (0,0), so the eigenvectors
            // have length 0 and the eigenvalues carry no correct digits) and
            // underflow below 2^-511 (the eigenvalues are wrong by tens of
            // percent). A*v = lambda*v is scale equivariant, so t*A must have
            // the same eigenvectors and t times the eigenvalues for every t.
            const B: Sym3 = {
                a00: 3, a01: 1, a02: 0.5, a11: 2, a12: -0.25, a22: -1
            };
            const base = solver.solve(B.a00, B.a01, B.a02, B.a11, B.a12, B.a22,
                false, +1);
            for (const t of [1e-300, 1e-200, 1e-170, 1e-160, 1e-155, 1e-100, 1,
                1e100, 1e155, 1e160, 1e200, 1e300]) {
                const r = solver.solve(t * B.a00, t * B.a01, t * B.a02,
                    t * B.a11, t * B.a12, t * B.a22, false, +1);
                for (let i = 0; i < 3; ++i) {
                    // Unit-length eigenvectors (length 0 before the fix).
                    expectClose(Math.sqrt(dot(r.evecs[i], r.evecs[i])), 1,
                        1e-14, 0);
                    // Eigenvalues scale by exactly t, to full precision.
                    expectClose(r.evals[i] / t, base.evals[i], 0, 1e-12);
                }
            }
        });

    it('regression (#379): a covariance-scale matrix keeps unit eigenvectors',
        () => {
            // The reported case: entries near 1e-160, where the squares formed
            // by GetCosSin are subnormal. Every eigenpair must still satisfy
            // A*v = lambda*v to full relative precision.
            check(fc.tuple(sym3, fc.constantFrom(1e-170, 1e-160, 1e170, 1e160)),
                ([m, t]) => {
                    const scaled: Sym3 = {
                        a00: t * m.a00, a01: t * m.a01, a02: t * m.a02,
                        a11: t * m.a11, a12: t * m.a12, a22: t * m.a22
                    };
                    const r = solver.solve(scaled.a00, scaled.a01, scaled.a02,
                        scaled.a11, scaled.a12, scaled.a22, false, +1);
                    const norm = normOf(m);
                    for (let i = 0; i < 3; ++i) {
                        expectClose(Math.sqrt(dot(r.evecs[i], r.evecs[i])), 1,
                            1e-14, 0);
                        const av = matVec(scaled, r.evecs[i]);
                        for (let j = 0; j < 3; ++j) {
                            expect(Math.abs(av[j] - r.evals[i] * r.evecs[i][j]))
                                .toBeLessThanOrEqual(1e-12 * t * norm);
                        }
                    }
                }, 100);
        });
});
