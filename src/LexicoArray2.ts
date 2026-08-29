// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) LexicoArray2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A class to provide 2D array access to a flat array that conforms to
// row-major order (rowMajor = true) or column-major order (rowMajor = false).
//
// Port notes: upstream is a template with a compile-time bool RowMajor
// parameter and optional compile-time dimensions, wrapping a Real* it does
// not own. The port collapses the four specializations into one class with a
// runtime rowMajor flag and runtime dimensions, wrapping (aliasing, not
// copying) a caller-owned number[]. 'operator()(r, c)' is ported as
// 'get(r, c)' and 'set(r, c, value)'.
export class LexicoArray2 {
    private mRowMajor: boolean;
    private mNumRows: number;
    private mNumCols: number;
    private mMatrix: number[];

    constructor(rowMajor: boolean, numRows: number, numCols: number, matrix: number[]) {
        this.mRowMajor = rowMajor;
        this.mNumRows = numRows;
        this.mNumCols = numCols;
        this.mMatrix = matrix;
    }

    isRowMajor(): boolean {
        return this.mRowMajor;
    }

    getNumRows(): number {
        return this.mNumRows;
    }

    getNumCols(): number {
        return this.mNumCols;
    }

    get(r: number, c: number): number {
        return this.mMatrix[this.index(r, c)];
    }

    set(r: number, c: number, value: number): void {
        this.mMatrix[this.index(r, c)] = value;
    }

    private index(r: number, c: number): number {
        return this.mRowMajor ? c + this.mNumCols * r : r + this.mNumRows * c;
    }
}
