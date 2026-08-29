// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SymmetricEigensolver2x2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream 'eval' cannot be used as an identifier in strict-mode
// JavaScript, so the eigenvalue/eigenvector names are 'evals' and 'evecs'
// throughout the eigensolver ports. The output reference parameters become a
// returned object literal per PORTING.md, and 'operator()' becomes 'solve'.

// evals[i] is the eigenvalue associated with the (row) eigenvector evecs[i].
export interface SymmetricEigensolver2x2Result {
    evals: [number, number];
    evecs: [[number, number], [number, number]];
}

export class SymmetricEigensolver2x2 {
    // The input matrix must be symmetric, so only the unique elements must
    // be specified: a00, a01, and a11.
    //
    // The order of the eigenvalues is specified by sortType: -1 (decreasing),
    // 0 (no sorting) or +1 (increasing). When sorted, the eigenvectors are
    // ordered accordingly, and {evecs[0], evecs[1]} is guaranteed to be a
    // right-handed orthonormal set.
    solve(a00: number, a01: number, a11: number, sortType: number): SymmetricEigensolver2x2Result {
        // Normalize (c2,s2) robustly, avoiding floating-point overflow in
        // the sqrt call.
        const zero = 0, one = 1, half = 0.5;
        let c2 = half * (a00 - a11);
        let s2 = a01;
        const maxAbsComp = Math.max(Math.abs(c2), Math.abs(s2));
        if (maxAbsComp > zero) {
            c2 /= maxAbsComp;  // in [-1,1]
            s2 /= maxAbsComp;  // in [-1,1]
            const length = Math.sqrt(c2 * c2 + s2 * s2);
            c2 /= length;
            s2 /= length;
            if (c2 > zero) {
                c2 = -c2;
                s2 = -s2;
            }
        }
        else {
            c2 = -one;
            s2 = zero;
        }

        const s = Math.sqrt(half * (one - c2));  // >= 1/sqrt(2)
        const c = half * s2 / s;

        const csqr = c * c, ssqr = s * s, mid = s2 * a01;
        const diagonal: [number, number] = [
            csqr * a00 + mid + ssqr * a11,
            csqr * a11 - mid + ssqr * a00
        ];

        if (sortType === 0 || sortType * diagonal[0] <= sortType * diagonal[1]) {
            return {
                evals: [diagonal[0], diagonal[1]],
                evecs: [[c, s], [-s, c]]
            };
        }
        else {
            return {
                evals: [diagonal[1], diagonal[0]],
                evecs: [[s, -c], [c, s]]
            };
        }
    }
}
