import { describe, it, expect } from 'vitest';
import { SymmetricEigensolver2x2 } from '../src/SymmetricEigensolver2x2.js';
import { check, expectClose, fc, wellScaled } from './helpers/arbitraries.js';

// A simple deterministic pseudorandom generator so test runs are repeatable.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        // Numerical Recipes LCG constants.
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function checkEigenpair(a00: number, a01: number, a11: number,
    lambda: number, v: [number, number], tol: number): void {
    // A*v = lambda*v
    const r0 = a00 * v[0] + a01 * v[1] - lambda * v[0];
    const r1 = a01 * v[0] + a11 * v[1] - lambda * v[1];
    expect(Math.abs(r0)).toBeLessThanOrEqual(tol);
    expect(Math.abs(r1)).toBeLessThanOrEqual(tol);
}

describe('SymmetricEigensolver2x2', () => {
    const solver = new SymmetricEigensolver2x2();

    it('solves randomized symmetric matrices (eigen equation, orthonormality)', () => {
        const random = makeRandom(1010);
        for (let trial = 0; trial < 200; ++trial) {
            const a00 = 4 * random() - 2;
            const a01 = 4 * random() - 2;
            const a11 = 4 * random() - 2;
            const scale = Math.max(Math.abs(a00), Math.abs(a01), Math.abs(a11), 1);
            const tol = 1e-13 * scale;

            for (const sortType of [-1, 0, 1]) {
                const { evals, evecs } = solver.solve(a00, a01, a11, sortType);

                checkEigenpair(a00, a01, a11, evals[0], evecs[0], tol);
                checkEigenpair(a00, a01, a11, evals[1], evecs[1], tol);

                // Orthonormality of the eigenvectors.
                const len0 = Math.hypot(evecs[0][0], evecs[0][1]);
                const len1 = Math.hypot(evecs[1][0], evecs[1][1]);
                const dot01 = evecs[0][0] * evecs[1][0] + evecs[0][1] * evecs[1][1];
                expect(Math.abs(len0 - 1)).toBeLessThanOrEqual(1e-14);
                expect(Math.abs(len1 - 1)).toBeLessThanOrEqual(1e-14);
                expect(Math.abs(dot01)).toBeLessThanOrEqual(1e-14);

                // Right-handedness: det[evecs] = +1.
                const det = evecs[0][0] * evecs[1][1] - evecs[0][1] * evecs[1][0];
                expect(Math.abs(det - 1)).toBeLessThanOrEqual(1e-14);
            }
        }
    });

    it('orders the eigenvalues according to sortType', () => {
        const random = makeRandom(2020);
        for (let trial = 0; trial < 100; ++trial) {
            const a00 = 4 * random() - 2;
            const a01 = 4 * random() - 2;
            const a11 = 4 * random() - 2;

            const inc = solver.solve(a00, a01, a11, +1);
            expect(inc.evals[0]).toBeLessThanOrEqual(inc.evals[1]);

            const dec = solver.solve(a00, a01, a11, -1);
            expect(dec.evals[0]).toBeGreaterThanOrEqual(dec.evals[1]);
        }
    });

    it('matches the analytic eigenvalues of a known matrix', () => {
        // A = [[2, 1], [1, 2]] has eigenvalues 1 and 3.
        const { evals, evecs } = solver.solve(2, 1, 2, +1);
        expect(evals[0]).toBeCloseTo(1, 14);
        expect(evals[1]).toBeCloseTo(3, 14);
        // Eigenvector for eigenvalue 1 is parallel to (1,-1).
        expect(Math.abs(evecs[0][0] + evecs[0][1])).toBeLessThanOrEqual(1e-14);
        // Eigenvector for eigenvalue 3 is parallel to (1,1).
        expect(Math.abs(evecs[1][0] - evecs[1][1])).toBeLessThanOrEqual(1e-14);
    });

    it('handles a diagonal matrix', () => {
        const { evals, evecs } = solver.solve(5, 0, -7, +1);
        expect(evals[0]).toBeCloseTo(-7, 14);
        expect(evals[1]).toBeCloseTo(5, 14);
        checkEigenpair(5, 0, -7, evals[0], evecs[0], 1e-13);
        checkEigenpair(5, 0, -7, evals[1], evecs[1], 1e-13);
    });

    it('handles the zero matrix', () => {
        const { evals, evecs } = solver.solve(0, 0, 0, +1);
        expect(evals[0]).toBe(0);
        expect(evals[1]).toBe(0);
        const det = evecs[0][0] * evecs[1][1] - evecs[0][1] * evecs[1][0];
        expect(Math.abs(det - 1)).toBeLessThanOrEqual(1e-14);
    });

    it('handles repeated eigenvalues (multiple of identity)', () => {
        const { evals, evecs } = solver.solve(3, 0, 3, +1);
        expect(evals[0]).toBeCloseTo(3, 14);
        expect(evals[1]).toBeCloseTo(3, 14);
        const dot01 = evecs[0][0] * evecs[1][0] + evecs[0][1] * evecs[1][1];
        expect(Math.abs(dot01)).toBeLessThanOrEqual(1e-14);
    });
});

// ---------------------------------------------------------------------------
// Verification (group V38).
// ---------------------------------------------------------------------------
describe('SymmetricEigensolver2x2 verification', () => {
    const solver = new SymmetricEigensolver2x2();

    // The analytic eigenvalues of [[a00,a01],[a01,a11]].
    function analytic(a00: number, a01: number, a11: number): [number, number] {
        const half = 0.5 * (a00 + a11);
        const disc = Math.sqrt(0.25 * (a00 - a11) * (a00 - a11) + a01 * a01);
        return [half - disc, half + disc];
    }

    const sym2 = fc.tuple(wellScaled(-10, 10), wellScaled(-10, 10),
        wellScaled(-10, 10));

    it('produces an orthonormal right-handed eigenbasis for every sortType',
        () => {
            check(fc.tuple(sym2, fc.constantFrom(-1, 0, 1)),
                ([[a00, a01, a11], sortType]) => {
                    const { evals, evecs } = solver.solve(a00, a01, a11, sortType);
                    const norm = Math.max(Math.abs(a00), Math.abs(a01),
                        Math.abs(a11), 1);
                    for (let i = 0; i < 2; ++i) {
                        // A*v = lambda*v.
                        const r0 = a00 * evecs[i][0] + a01 * evecs[i][1]
                            - evals[i] * evecs[i][0];
                        const r1 = a01 * evecs[i][0] + a11 * evecs[i][1]
                            - evals[i] * evecs[i][1];
                        expect(Math.abs(r0)).toBeLessThanOrEqual(1e-14 * norm);
                        expect(Math.abs(r1)).toBeLessThanOrEqual(1e-14 * norm);
                        // Unit length.
                        expectClose(Math.hypot(evecs[i][0], evecs[i][1]), 1,
                            1e-14, 0);
                    }
                    // Orthogonal and right handed (det = +1).
                    expectClose(evecs[0][0] * evecs[1][0]
                        + evecs[0][1] * evecs[1][1], 0, 1e-14, 0);
                    expectClose(evecs[0][0] * evecs[1][1]
                        - evecs[0][1] * evecs[1][0], 1, 1e-14, 0);
                });
        });

    it('matches the closed-form eigenvalues and preserves trace and determinant',
        () => {
            check(sym2, ([a00, a01, a11]) => {
                const { evals } = solver.solve(a00, a01, a11, +1);
                const [lo, hi] = analytic(a00, a01, a11);
                const norm = Math.max(Math.abs(a00), Math.abs(a01),
                    Math.abs(a11), 1);
                expectClose(evals[0], lo, 1e-13 * norm, 1e-13);
                expectClose(evals[1], hi, 1e-13 * norm, 1e-13);
                // The similarity transform preserves trace and determinant.
                expectClose(evals[0] + evals[1], a00 + a11, 1e-13 * norm, 1e-13);
                expectClose(evals[0] * evals[1], a00 * a11 - a01 * a01,
                    1e-12 * norm * norm, 1e-12);
            });
        });

    it('orders the eigenvalues and the eigenvectors according to sortType',
        () => {
            check(sym2, ([a00, a01, a11]) => {
                const inc = solver.solve(a00, a01, a11, +1);
                const dec = solver.solve(a00, a01, a11, -1);
                const none = solver.solve(a00, a01, a11, 0);
                expect(inc.evals[0]).toBeLessThanOrEqual(inc.evals[1]);
                expect(dec.evals[0]).toBeGreaterThanOrEqual(dec.evals[1]);
                // The two sorted orders are reverses of one another and both
                // are permutations of the unsorted pair.
                expectClose(inc.evals[0], dec.evals[1], 0, 0);
                expectClose(inc.evals[1], dec.evals[0], 0, 0);
                const unsorted = [none.evals[0], none.evals[1]].sort(
                    (x, y) => x - y);
                expectClose(inc.evals[0], unsorted[0], 0, 0);
                expectClose(inc.evals[1], unsorted[1], 0, 0);
            });
        });

    it('is equivariant under scaling of the matrix over the full exponent range',
        () => {
            // Upstream normalizes (c2,s2) by max(|c2|,|s2|) before the sqrt
            // precisely so that the solver survives extreme scales. This is
            // the property that SymmetricEigensolver3x3's GetCosSin loses by
            // omitting the same rescaling (see issue #379); pin it here for
            // the sibling that has it.
            const scales = [1e-300, 1e-200, 1e-100, 1, 1e100, 1e200, 1e300];
            check(sym2, ([a00, a01, a11]) => {
                const norm = Math.max(Math.abs(a00), Math.abs(a01),
                    Math.abs(a11));
                if (norm === 0) {
                    return;
                }
                const base = solver.solve(a00, a01, a11, +1);
                for (const t of scales) {
                    const r = solver.solve(t * a00, t * a01, t * a11, +1);
                    for (let i = 0; i < 2; ++i) {
                        expectClose(Math.hypot(r.evecs[i][0], r.evecs[i][1]), 1,
                            1e-14, 0);
                        // An eigenvalue much smaller than the matrix norm
                        // comes from a cancelling sum of entries, so the
                        // tolerance is relative to the norm, not to the
                        // eigenvalue.
                        expectClose(r.evals[i] / t, base.evals[i],
                            1e-13 * norm, 1e-13);
                    }
                }
            }, 100);
        });

    it('handles diagonal matrices exactly', () => {
        check(fc.tuple(wellScaled(-10, 10), wellScaled(-10, 10)),
            ([a00, a11]) => {
                const { evals, evecs } = solver.solve(a00, 0, a11, +1);
                const lo = Math.min(a00, a11);
                const hi = Math.max(a00, a11);
                // A diagonal matrix has a01 = 0, so c2 = (a00-a11)/2 and
                // s2 = 0; the eigenvalues are the diagonal entries exactly.
                expect(evals[0] + 0).toBe(lo + 0);
                expect(evals[1] + 0).toBe(hi + 0);
                // The eigenvectors are the axes (up to sign and order).
                for (let i = 0; i < 2; ++i) {
                    expectClose(Math.abs(evecs[i][0]) + Math.abs(evecs[i][1]), 1,
                        1e-15, 0);
                }
            });
    });

    it('returns the identity basis for a multiple of the identity', () => {
        check(wellScaled(-10, 10), a => {
            const { evals, evecs } = solver.solve(a, 0, a, +1);
            expect(evals[0] + 0).toBe(a + 0);
            expect(evals[1] + 0).toBe(a + 0);
            // (c2,s2) is (0,0) here, so upstream falls back to (-1,0), giving
            // s = 1 and c = 0: the basis is {(0,1), (-1,0)}.
            expect(evecs[0][0] + 0).toBe(0);
            expect(evecs[0][1]).toBe(1);
            expect(evecs[1][0]).toBe(-1);
            expect(evecs[1][1] + 0).toBe(0);
        });
    });
});
