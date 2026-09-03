import { describe, expect, it } from 'vitest';
import { ApprPolynomial3 } from '../src/ApprPolynomial3.js';
import { ApprQuery } from '../src/ApprQuery.js';

// Solve the (small, well-conditioned) system A*X = B by Gauss-Jordan
// elimination with partial pivoting; an independent check of the fit.
function solve(A: number[][], B: number[]): number[] {
    const n = B.length;
    const M = A.map((row, i) => row.concat([B[i]]));
    for (let c = 0; c < n; ++c) {
        let piv = c;
        for (let r = c + 1; r < n; ++r) {
            if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) {
                piv = r;
            }
        }
        const t = M[c]; M[c] = M[piv]; M[piv] = t;
        for (let r = 0; r < n; ++r) {
            if (r !== c) {
                const f = M[r][c] / M[c][c];
                for (let k = c; k <= n; ++k) {
                    M[r][k] -= f * M[c][k];
                }
            }
        }
    }
    return M.map((row, i) => row[n] / row[i]);
}

// w = sum_{i,j} c[i + (xDegree+1)*j] * x^i * y^j
function poly(c: readonly number[], xDegree: number, yDegree: number,
    x: number, y: number): number {
    let w = 0;
    for (let j = 0; j <= yDegree; ++j) {
        for (let i = 0; i <= xDegree; ++i) {
            w += c[i + (xDegree + 1) * j] * Math.pow(x, i) * Math.pow(y, j);
        }
    }
    return w;
}

describe('ApprPolynomial3', () => {
    it('initializes the parameters to zero and reports the sizes', () => {
        const fitter = new ApprPolynomial3(2, 1);
        expect(fitter.getParameters()).toEqual([0, 0, 0, 0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(6);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getYDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('recovers a degree (2,1) polynomial exactly from data sampled on it', () => {
        // c[i + 3*j] for i in [0,2], j in [0,1].
        const c = [1, -2, 0.5, 3, 0.25, -1.5];
        const observations: number[][] = [];
        for (let a = 0; a <= 4; ++a) {
            for (let b = 0; b <= 4; ++b) {
                const x = -1 + 0.5 * a;
                const y = -2 + 1 * b;
                observations.push([x, y, poly(c, 2, 1, x, y)]);
            }
        }

        const fitter = new ApprPolynomial3(2, 1);
        expect(fitter.fit(observations)).toBe(true);

        const p = fitter.getParameters();
        expect(p.length).toBe(6);
        for (let i = 0; i < 6; ++i) {
            expect(p[i]).toBeCloseTo(c[i], 8);
        }

        for (const obs of observations) {
            expect(fitter.error(obs)).toBeLessThan(1e-8);
        }

        // The evaluate() Horner scheme agrees with the direct sum, both
        // inside and outside the sampled domain.
        for (const [x, y] of [[0.3, -1.2], [-0.9, 1.7], [5, 5], [-4, 3]]) {
            expect(fitter.evaluate(x, y)).toBeCloseTo(poly(c, 2, 1, x, y), 6);
        }
    });

    it('records the x- and y-domains of the fitted samples', () => {
        const fitter = new ApprPolynomial3(1, 1);
        fitter.fit([[-1, 4, 0], [2, -3, 1], [0.5, 7, 2], [3, 1, 3]]);
        expect(fitter.getXDomain()).toEqual([-1, 3]);
        expect(fitter.getYDomain()).toEqual([-3, 7]);
    });

    it('is a plane fit when both degrees are (1,0)', () => {
        // w = c0 + c1*x reduces to a 1D line fit in x, independent of y.
        const fitter = new ApprPolynomial3(1, 0);
        expect(fitter.fit([[0, 5, 1], [1, -3, 3], [2, 8, 5], [3, 0, 7]])).toBe(true);
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(1, 10);
        expect(p[1]).toBeCloseTo(2, 10);
    });

    it('matches the normal equations for an overdetermined bilinear fit', () => {
        const observations = [
            [0, 0, 1], [1, 0, 2.1], [0, 1, 3.2], [1, 1, 4.4],
            [2, 0, 3.0], [2, 1, 5.4], [0, 2, 5.1], [1, 2, 6.6]
        ];

        // Terms (1, x, y, x*y) in the parameter order i + 2*j.
        const term = (x: number, y: number, k: number): number => {
            const i = k % 2;
            const j = (k - i) / 2;
            return Math.pow(x, i) * Math.pow(y, j);
        };
        const A: number[][] = [];
        const B: number[] = [];
        for (let r = 0; r < 4; ++r) {
            A.push(new Array<number>(4).fill(0));
            B.push(0);
            for (const [x, y, w] of observations) {
                B[r] += w * term(x, y, r);
                for (let c = 0; c < 4; ++c) {
                    A[r][c] += term(x, y, r) * term(x, y, c);
                }
            }
        }
        const expected = solve(A, B);

        const fitter = new ApprPolynomial3(1, 1);
        expect(fitter.fit(observations)).toBe(true);
        const p = fitter.getParameters();
        for (let i = 0; i < 4; ++i) {
            expect(p[i]).toBeCloseTo(expected[i], 8);
        }

        // The residuals are orthogonal to each basis term.
        for (let r = 0; r < 4; ++r) {
            let sum = 0;
            for (const [x, y, w] of observations) {
                sum += (fitter.evaluate(x, y) - w) * term(x, y, r);
            }
            expect(sum).toBeCloseTo(0, 8);
        }
    });

    it('fails when the system is singular', () => {
        // Three samples cannot determine four bilinear coefficients.
        const fitter = new ApprPolynomial3(1, 1);
        expect(fitter.fit([[0, 0, 1], [1, 0, 2], [0, 1, 3]])).toBe(false);
        expect(fitter.getParameters()).toEqual([0, 0, 0, 0]);
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomial3(0, 0);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0, 1]], [7])).toBe(false);
            expect(fitter.getParameters()).toEqual([0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const source = new ApprPolynomial3(1, 1);
        source.fit([[0, 0, 1], [1, 0, 2], [0, 1, 3], [1, 1, 5]]);
        const target = new ApprPolynomial3(1, 1);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.getYDomain()).toEqual(source.getYDomain());
        expect(target.evaluate(0.5, 0.5)).toBe(source.evaluate(0.5, 0.5));

        const copied = target.getParameters().slice();
        source.fit([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 1]]);
        expect(target.getParameters()).toEqual(copied);
    });
});
