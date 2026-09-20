// Replays oracle/cpp/cases/v01-algebra.cpp (verify group 1, algebra).
// Keep the two files in the same order.
//
// Every case is arithmetic only (+ - * / sqrt fabs and comparisons), so every
// case is declared exact: the port must be bit-identical to the MSVC build of
// upstream GTE, signed zeros included.
import { describe } from 'vitest';
import { GMatrix } from '../../src/GMatrix.js';
import { GVector } from '../../src/GVector.js';
import {
    Matrix, addMatrix, determinant, divMatrix, hliftMatrix, hprojectMatrix,
    inverse, l1Norm, l2Norm, lInfinityNorm, makeDiagonal, mulMatrix,
    multiplyAB, multiplyABT, multiplyATB, multiplyATBT, multiplyDM,
    multiplyMD, negateMatrix, outerProduct, subMatrix, transpose
} from '../../src/Matrix.js';
import {
    Vector, add, compDiv, compMul, computeExtremes, div, dot, getOrthogonal,
    hlift, hproject, length, lift, mul, negate, normalize, orthonormalize,
    project, sub
} from '../../src/Vector.js';
import {
    IntrinsicsVector2, computeBarycentrics2, computeOrthogonalComplement2,
    dotPerp, perp, unitPerp
} from '../../src/Vector2.js';
import {
    IntrinsicsVector3, computeBarycentrics3, computeOrthogonalComplement3,
    cross, dotCross, fastComputeOrthogonalComplement, unitCross
} from '../../src/Vector3.js';
import {
    computeOrthogonalComplement4, dotHyperCross, hyperCross, unitHyperCross
} from '../../src/Vector4.js';
import { OracleFamily, type OracleIO } from './harness.js';

function mat(io: OracleIO, numRows: number, numCols: number): Matrix {
    return Matrix.fromArray(numRows, numCols, io.reals(numRows * numCols));
}

function gmat(io: OracleIO, numRows: number, numCols: number): GMatrix {
    return GMatrix.fromArray(numRows, numCols, io.reals(numRows * numCols));
}

function gvec(io: OracleIO, size: number): GVector {
    return GVector.fromArray(io.reals(size));
}

// The (R, C) shapes of the three-way shape switch in the C++ cases.
const SHAPES: ReadonlyArray<readonly [number, number]> =
    [[2, 3], [3, 3], [4, 2]];
// The (R, K, C) shapes of Matrix.multiply.
const MUL_SHAPES: ReadonlyArray<readonly [number, number, number]> =
    [[2, 3, 4], [3, 3, 3], [4, 2, 3]];

describe('oracle: v01-algebra', () => {
    const family = new OracleFamily('v01-algebra');

    // ---- Vector.h ----------------------------------------------------------

    family.case('Vector.arithmetic', (io) => {
        const v0 = io.vec(4);
        const v1 = io.vec(4);
        const s = io.real();
        io.outVec(negate(v0));
        io.outVec(v0);                  // upstream unary operator+ is identity
        io.outVec(add(v0, v1));
        io.outVec(sub(v0, v1));
        io.outVec(mul(v0, s));
        io.outVec(mul(s, v0));
        io.outVec(div(v0, s));
        io.outVec(compMul(v0, v1));
        io.outVec(compDiv(v0, v1));
    }, { exact: true });

    family.case('Vector.dot', (io) => {
        const v0 = io.vec(4);
        const v1 = io.vec(4);
        io.outReal(dot(v0, v1));
    }, { exact: true });

    family.case('Vector.dot.signedZero', (io) => {
        const v0 = io.vec(4);
        const v1 = io.vec(4);
        io.outReal(dot(v0, v1));
        io.outReal(dot(v1, v0));
    }, { exact: true });

    family.case('Vector.length', (io) => {
        const v = io.vec(4);
        io.outReal(length(v, false));
        io.outReal(length(v, true));
    }, { exact: true });

    family.case('Vector.normalize', (io) => {
        const v = io.vec(4);
        const plain = v.clone();
        const robust = v.clone();
        const lengthPlain = normalize(plain, false);
        const lengthRobust = normalize(robust, true);
        io.outReal(lengthPlain);
        io.outVec(plain);
        io.outReal(lengthRobust);
        io.outVec(robust);
    }, { exact: true });

    family.case('Vector.orthonormalize', (io) => {
        const numInputs = io.integer();
        const robust = io.boolean();
        const v = [io.vec(4), io.vec(4), io.vec(4), io.vec(4)];
        io.outReal(orthonormalize(numInputs, v, robust));
        for (const vi of v) {
            io.outVec(vi);
        }
    }, { exact: true });

    family.case('Vector.getOrthogonal', (io) => {
        const n = io.integer();
        const unitLength = io.boolean();
        const v = io.vec(n);
        io.outVec(getOrthogonal(v, unitLength));
    }, { exact: true });

    family.case('Vector.computeExtremes', (io) => {
        const numVectors = io.integer();
        const v: Vector[] = [];
        for (let i = 0; i < numVectors; ++i) {
            v.push(io.vec(3));
        }
        const result = computeExtremes(v);
        io.outBool(result !== null);
        io.outVec(result!.vmin);
        io.outVec(result!.vmax);
    }, { exact: true });

    family.case('Vector.liftProject', (io) => {
        const inject = io.integer();
        const reject = io.integer();
        const last = io.real();
        const v = io.vec(4);
        io.outVec(hlift(v, last));
        io.outVec(hproject(v));
        io.outVec(lift(v, inject, last));
        io.outVec(project(v, reject));
    }, { exact: true });

    family.case('Vector.special', (io) => {
        const d = io.integer();
        const v0 = io.vec(4);
        const v1 = io.vec(4);
        const zero = new Vector(4);
        const ones = new Vector(4);
        const unit = new Vector(4);
        zero.makeZero();
        ones.makeOnes();
        unit.makeUnit(d);
        io.outVec(zero);
        io.outVec(ones);
        io.outVec(unit);
        io.outVec(Vector.zero(4));
        io.outVec(Vector.ones(4));
        io.outVec(Vector.unit(4, d));
        io.outInt(v0.size);
        io.outBool(v0.equals(v1));
        io.outBool(v0.notEquals(v1));
        io.outBool(v0.lessThan(v1));
        io.outBool(v0.lessThanOrEqual(v1));
        io.outBool(v0.greaterThan(v1));
        io.outBool(v0.greaterThanOrEqual(v1));
        io.outBool(v0.equals(v0));
        io.outBool(v0.lessThanOrEqual(v0));
    }, { exact: true });

    // ---- Vector2.h ---------------------------------------------------------

    family.case('Vector2.perp', (io) => {
        const robust = io.boolean();
        const v0 = io.vec(2);
        const v1 = io.vec(2);
        io.outVec(perp(v0));
        io.outVec(unitPerp(v0, robust));
        io.outReal(dotPerp(v0, v1));
        io.outReal(dotPerp(v1, v0));
    }, { exact: true });

    family.case('Vector2.computeOrthogonalComplement', (io) => {
        const numInputs = io.integer();
        const robust = io.boolean();
        const v = [io.vec(2), io.vec(2)];
        io.outReal(computeOrthogonalComplement2(numInputs, v, robust));
        io.outVec(v[0]);
        io.outVec(v[1]);
    }, { exact: true });

    family.case('Vector2.computeBarycentrics', (io) => {
        const epsilon = io.real();
        const v0 = io.vec(2);
        const v1 = io.vec(2);
        const v2 = io.vec(2);
        const p = io.vec(2);
        const r = computeBarycentrics2(p, v0, v1, v2, epsilon);
        io.outBool(r.valid);
        io.outReals(r.bary);
    }, { exact: true });

    family.case('Vector2.intrinsicsVector2', (io) => {
        const numVectors = io.integer();
        const epsilon = io.real();
        const v: Vector[] = [];
        for (let i = 0; i < numVectors; ++i) {
            v.push(io.vec(2));
        }
        const intr = new IntrinsicsVector2(v, epsilon);
        io.outReal(intr.epsilon);
        io.outInt(intr.dimension);
        io.outReal(intr.min[0]);
        io.outReal(intr.min[1]);
        io.outReal(intr.max[0]);
        io.outReal(intr.max[1]);
        io.outReal(intr.maxRange);
        io.outVec(intr.origin);
        io.outVec(intr.direction[0]);
        io.outVec(intr.direction[1]);
        io.outInt(intr.extreme[0]);
        io.outInt(intr.extreme[1]);
        io.outInt(intr.extreme[2]);
        io.outBool(intr.extremeCCW);
    }, { exact: true });

    // ---- Vector3.h ---------------------------------------------------------

    family.case('Vector3.cross', (io) => {
        const robust = io.boolean();
        const v0 = io.vec(3);
        const v1 = io.vec(3);
        const v2 = io.vec(3);
        io.outVec(cross(v0, v1));
        io.outVec(unitCross(v0, v1, robust));
        io.outReal(dotCross(v0, v1, v2));
        io.outReal(dotCross(v2, v0, v1));
    }, { exact: true });

    family.case('Vector3.cross.4d', (io) => {
        const robust = io.boolean();
        io.boolean();                   // 'affine'; only shapes the inputs
        const v0 = io.vec(4);
        const v1 = io.vec(4);
        const v2 = io.vec(4);
        io.outVec(cross(v0, v1));
        io.outVec(unitCross(v0, v1, robust));
        io.outReal(dotCross(v0, v1, v2));
    }, { exact: true });

    family.case('Vector3.computeOrthogonalComplement', (io) => {
        const numInputs = io.integer();
        const robust = io.boolean();
        const v = [io.vec(3), io.vec(3), io.vec(3)];
        io.outReal(computeOrthogonalComplement3(numInputs, v, robust));
        for (const vi of v) {
            io.outVec(vi);
        }
    }, { exact: true });

    family.case('Vector3.fastComputeOrthogonalComplement', (io) => {
        const v2 = io.vec(3);
        const r = fastComputeOrthogonalComplement(v2);
        io.outVec(r.v0);
        io.outVec(r.v1);
    }, { exact: true });

    family.case('Vector3.computeBarycentrics', (io) => {
        const epsilon = io.real();
        const v0 = io.vec(3);
        const v1 = io.vec(3);
        const v2 = io.vec(3);
        const v3 = io.vec(3);
        const p = io.vec(3);
        const r = computeBarycentrics3(p, v0, v1, v2, v3, epsilon);
        io.outBool(r.valid);
        io.outReals(r.bary);
    }, { exact: true });

    family.case('Vector3.intrinsicsVector3', (io) => {
        const numVectors = io.integer();
        const epsilon = io.real();
        const v: Vector[] = [];
        for (let i = 0; i < numVectors; ++i) {
            v.push(io.vec(3));
        }
        const intr = new IntrinsicsVector3(v, epsilon);
        io.outReal(intr.epsilon);
        io.outInt(intr.dimension);
        io.outReals(intr.min);
        io.outReals(intr.max);
        io.outReal(intr.maxRange);
        io.outVec(intr.origin);
        io.outVec(intr.direction[0]);
        io.outVec(intr.direction[1]);
        io.outVec(intr.direction[2]);
        for (let k = 0; k < 4; ++k) {
            io.outInt(intr.extreme[k]);
        }
        io.outBool(intr.extremeCCW);
    }, { exact: true });

    // ---- Vector4.h ---------------------------------------------------------

    family.case('Vector4.hyperCross', (io) => {
        const robust = io.boolean();
        const v0 = io.vec(4);
        const v1 = io.vec(4);
        const v2 = io.vec(4);
        const v3 = io.vec(4);
        io.outVec(hyperCross(v0, v1, v2));
        io.outVec(unitHyperCross(v0, v1, v2, robust));
        io.outReal(dotHyperCross(v0, v1, v2, v3));
        io.outReal(dotHyperCross(v3, v0, v1, v2));
    }, { exact: true });

    family.case('Vector4.computeOrthogonalComplement', (io) => {
        const numInputs = io.integer();
        const robust = io.boolean();
        const v = [io.vec(4), io.vec(4), io.vec(4), io.vec(4)];
        io.outReal(computeOrthogonalComplement4(numInputs, v, robust));
        for (const vi of v) {
            io.outVec(vi);
        }
    }, { exact: true });

    family.case('Vector4.computeOrthogonalComplement.zeroMiddle', (io) => {
        const robust = io.boolean();
        const v = [io.vec(4), io.vec(4), io.vec(4), io.vec(4)];
        io.outReal(computeOrthogonalComplement4(1, v, robust));
        for (const vi of v) {
            io.outVec(vi);
        }
    }, { exact: true });

    // ---- GVector.h ---------------------------------------------------------

    family.case('GVector.arithmetic', (io) => {
        const size = io.integer();
        const v0 = gvec(io, size);
        const v1 = gvec(io, size);
        const s = io.real();
        io.outGVec(negate(v0));
        io.outGVec(v0);                 // upstream unary operator+ is identity
        io.outGVec(add(v0, v1));
        io.outGVec(sub(v0, v1));
        io.outGVec(mul(v0, s));
        io.outGVec(mul(s, v0));
        io.outGVec(div(v0, s));
    }, { exact: true });

    family.case('GVector.dot', (io) => {
        const size = io.integer();
        const v0 = gvec(io, size);
        const v1 = gvec(io, size);
        io.outReal(dot(v0, v1));
    }, { exact: true });

    family.case('GVector.dot.signedZero', (io) => {
        const size = io.integer();
        const v0 = gvec(io, size);
        const v1 = gvec(io, size);
        io.outReal(dot(v0, v1));
        io.outReal(dot(v1, v0));
    }, { exact: true });

    family.case('GVector.lengthNormalize', (io) => {
        const size = io.integer();
        const v = gvec(io, size);
        const plain = v.clone();
        const robust = v.clone();
        io.outReal(length(v, false));
        io.outReal(length(v, true));
        const lengthPlain = normalize(plain, false);
        const lengthRobust = normalize(robust, true);
        io.outReal(lengthPlain);
        io.outGVec(plain);
        io.outReal(lengthRobust);
        io.outGVec(robust);
    }, { exact: true });

    family.case('GVector.orthonormalize', (io) => {
        const size = io.integer();
        const numInputs = io.integer();
        const robust = io.boolean();
        const v: Vector[] = [];
        for (let i = 0; i < size; ++i) {
            v.push(gvec(io, size));
        }
        io.outReal(orthonormalize(numInputs, v, robust));
        for (const vi of v) {
            io.outGVec(vi);
        }
    }, { exact: true });

    family.case('GVector.computeExtremes', (io) => {
        const size = io.integer();
        const numVectors = io.integer();
        const v: Vector[] = [];
        for (let i = 0; i < numVectors; ++i) {
            v.push(gvec(io, size));
        }
        const result = computeExtremes(v);
        io.outBool(result !== null);
        io.outGVec(result!.vmin);
        io.outGVec(result!.vmax);
    }, { exact: true });

    family.case('GVector.liftProject', (io) => {
        const size = io.integer();
        const inject = io.integer();
        const reject = io.integer();
        const last = io.real();
        const v = gvec(io, size);
        io.outGVec(hlift(v, last));
        io.outGVec(hproject(v));
        io.outGVec(lift(v, inject, last));
        io.outGVec(project(v, reject));
    }, { exact: true });

    family.case('GVector.special', (io) => {
        const size0 = io.integer();
        const size1 = io.integer();
        const d = io.integer();
        const newSize = io.integer();
        const v0 = gvec(io, size0);
        const v1 = gvec(io, size1);
        const zero = new GVector(size0);
        const unit = new GVector(size0);
        zero.makeZero();
        unit.makeUnit(d);
        io.outGVec(zero);
        io.outGVec(unit);
        io.outGVec(GVector.zero(size0));
        io.outGVec(GVector.unit(size0, d));
        io.outInt(v0.size);
        io.outBool(v0.equals(v1));
        io.outBool(v0.notEquals(v1));
        io.outBool(v0.lessThan(v1));
        io.outBool(v0.lessThanOrEqual(v1));
        io.outBool(v0.greaterThan(v1));
        io.outBool(v0.greaterThanOrEqual(v1));
        const resized = v0.clone();
        resized.setSize(newSize);
        io.outGVec(resized);
    }, { exact: true });

    // ---- Matrix.h ----------------------------------------------------------

    family.case('Matrix.arithmetic', (io) => {
        const [r, c] = SHAPES[io.integer()];
        const m0 = mat(io, r, c);
        const m1 = mat(io, r, c);
        const s = io.real();
        io.outMat(negateMatrix(m0));
        io.outMat(m0);                  // upstream unary operator+ is identity
        io.outMat(addMatrix(m0, m1));
        io.outMat(subMatrix(m0, m1));
        io.outMat(mulMatrix(m0, s));
        io.outMat(mulMatrix(s, m0));
        io.outMat(divMatrix(m0, s));
        io.outReal(l1Norm(m0));
        io.outReal(l2Norm(m0));
        io.outReal(lInfinityNorm(m0));
    }, { exact: true });

    family.case('Matrix.multiply', (io) => {
        const [r, k, c] = MUL_SHAPES[io.integer()];
        const A = mat(io, r, k);
        const B = mat(io, k, c);
        const Bt = mat(io, c, k);
        const At = mat(io, k, r);
        const v = io.vec(k);
        const w = io.vec(r);
        io.outMat(mulMatrix(A, B));
        io.outMat(multiplyAB(A, B));
        io.outMat(multiplyABT(A, Bt));
        io.outMat(multiplyATB(At, B));
        io.outMat(multiplyATBT(At, Bt));
        io.outVec(mulMatrix(A, v));
        io.outVec(mulMatrix(w, A));
        io.outMat(transpose(A));
    }, { exact: true });

    family.case('Matrix.multiplyDiagonal', (io) => {
        const [r, c] = SHAPES[io.integer()];
        const M = mat(io, r, c);
        const dc = io.vec(c);
        const dr = io.vec(r);
        io.outMat(multiplyMD(M, dc));
        io.outMat(multiplyDM(dr, M));
        io.outMat(outerProduct(dr, dc));
        const diag = new Matrix(c, c);
        makeDiagonal(dc, diag);
        io.outMat(diag);
    }, { exact: true });

    family.case('Matrix.inverse', (io) => {
        const n = io.integer();
        const M = mat(io, n, n);
        const r = inverse(M);
        io.outBool(r.invertible);
        io.outMat(r.inverse);
        io.outReal(determinant(M));
        io.outMat(inverse(M).inverse);
    }, { exact: true });

    family.case('Matrix.access', (io) => {
        const [numRows, numCols] = SHAPES[io.integer()];
        const r = io.integer();
        const c = io.integer();
        const ur = io.integer();
        const uc = io.integer();
        const flat = io.integer();
        const M = mat(io, numRows, numCols);
        const M2 = mat(io, numRows, numCols);
        const row = io.vec(numCols);
        const col = io.vec(numRows);
        const A = M.clone();
        A.setRow(r, row);
        A.setCol(c, col);
        io.outMat(A);
        io.outVec(M.getRow(r));
        io.outVec(M.getCol(c));
        io.outReal(M.get(r, c));
        io.outReal(M.getFlat(flat));
        const zero = new Matrix(numRows, numCols);
        const unit = new Matrix(numRows, numCols);
        const identity = new Matrix(numRows, numCols);
        zero.makeZero();
        unit.makeUnit(ur, uc);
        identity.makeIdentity();
        io.outMat(zero);
        io.outMat(unit);
        io.outMat(identity);
        io.outMat(Matrix.zero(numRows, numCols));
        io.outMat(Matrix.unit(numRows, numCols, ur, uc));
        io.outMat(Matrix.identity(numRows, numCols));
        io.outBool(M.equals(M2));
        io.outBool(M.notEquals(M2));
        io.outBool(M.lessThan(M2));
        io.outBool(M.lessThanOrEqual(M2));
        io.outBool(M.greaterThan(M2));
        io.outBool(M.greaterThanOrEqual(M2));
        io.outBool(M.equals(M));
    }, { exact: true });

    family.case('Matrix.hliftProject', (io) => {
        const n = io.integer();
        const M = mat(io, n, n);
        io.outMat(hliftMatrix(M));
        io.outMat(hprojectMatrix(M));
    }, { exact: true });

    // ---- GMatrix.h ---------------------------------------------------------

    family.case('GMatrix.arithmetic', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const m0 = gmat(io, numRows, numCols);
        const m1 = gmat(io, numRows, numCols);
        const s = io.real();
        io.outGMat(negateMatrix(m0));
        io.outGMat(m0);                 // upstream unary operator+ is identity
        io.outGMat(addMatrix(m0, m1));
        io.outGMat(subMatrix(m0, m1));
        io.outGMat(mulMatrix(m0, s));
        io.outGMat(mulMatrix(s, m0));
        io.outGMat(divMatrix(m0, s));
        io.outReal(l1Norm(m0));
        io.outReal(l2Norm(m0));
        io.outReal(lInfinityNorm(m0));
    }, { exact: true });

    family.case('GMatrix.multiply', (io) => {
        const r = io.integer();
        const k = io.integer();
        const c = io.integer();
        const A = gmat(io, r, k);
        const B = gmat(io, k, c);
        const Bt = gmat(io, c, k);
        const At = gmat(io, k, r);
        const v = gvec(io, k);
        const w = gvec(io, r);
        io.outGMat(mulMatrix(A, B));
        io.outGMat(multiplyAB(A, B));
        io.outGMat(multiplyABT(A, Bt));
        io.outGMat(multiplyATB(At, B));
        io.outGMat(multiplyATBT(At, Bt));
        io.outGVec(mulMatrix(A, v));
        io.outGVec(mulMatrix(w, A));
        io.outGMat(transpose(A));
    }, { exact: true });

    family.case('GMatrix.multiplyDiagonal', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const M = gmat(io, numRows, numCols);
        const dc = gvec(io, numCols);
        const dr = gvec(io, numRows);
        io.outGMat(multiplyMD(M, dc));
        io.outGMat(multiplyDM(dr, M));
        io.outGMat(outerProduct(dr, dc));
        const diag = new GMatrix(numRows, numCols);
        makeDiagonal(dc, diag);
        io.outGMat(diag);
    }, { exact: true });

    family.case('GMatrix.inverse', (io) => {
        const n = io.integer();
        const M = gmat(io, n, n);
        const r = inverse(M);
        io.outBool(r.invertible);
        io.outGMat(GMatrix.fromMatrix(r.inverse));
        io.outReal(determinant(M));
        io.outGMat(GMatrix.fromMatrix(inverse(M).inverse));
    }, { exact: true });

    family.case('GMatrix.access', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const r = io.integer();
        const c = io.integer();
        const flat = io.integer();
        const newRows = io.integer();
        const newCols = io.integer();
        const M = gmat(io, numRows, numCols);
        const M2 = gmat(io, numRows, numCols);
        const row = gvec(io, numCols);
        const col = gvec(io, numRows);
        const A = M.clone();
        A.setRow(r, row);
        A.setCol(c, col);
        io.outGMat(A);
        io.outGVec(M.getRow(r));
        io.outGVec(M.getCol(c));
        io.outReal(M.get(r, c));
        io.outReal(M.getFlat(flat));
        io.outInt(M.getNumRows());
        io.outInt(M.getNumCols());
        io.outInt(M.getNumElements());
        const zero = new GMatrix(numRows, numCols);
        const unit = new GMatrix(numRows, numCols);
        const identity = new GMatrix(numRows, numCols);
        zero.makeZero();
        unit.makeUnit(r, c);
        identity.makeIdentity();
        io.outGMat(zero);
        io.outGMat(unit);
        io.outGMat(identity);
        io.outGMat(GMatrix.zero(numRows, numCols));
        io.outGMat(GMatrix.unit(numRows, numCols, r, c));
        io.outGMat(GMatrix.identity(numRows, numCols));
        io.outBool(M.equals(M2));
        io.outBool(M.notEquals(M2));
        io.outBool(M.lessThan(M2));
        io.outBool(M.lessThanOrEqual(M2));
        io.outBool(M.greaterThan(M2));
        io.outBool(M.greaterThanOrEqual(M2));
        io.outBool(M.equals(M));
        const resized = M.clone();
        resized.setSize(newRows, newCols);
        io.outGMat(resized);
    }, { exact: true });

    family.case('GMatrix.access.mismatchedSizes', (io) => {
        const rows0 = io.integer();
        const cols0 = io.integer();
        const rows1 = io.integer();
        const cols1 = io.integer();
        const M0 = gmat(io, rows0, cols0);
        const M1 = gmat(io, rows1, cols1);
        io.outBool(M0.equals(M1));
        io.outBool(M0.notEquals(M1));
        io.outBool(M0.lessThan(M1));
        io.outBool(M0.lessThanOrEqual(M1));
        io.outBool(M0.greaterThan(M1));
        io.outBool(M0.greaterThanOrEqual(M1));
    }, { exact: true });

    family.case('GMatrix.access.invalidIndex', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const r = io.integer();
        const c = io.integer();
        const M = gmat(io, numRows, numCols);
        io.outReal(M.get(r, c));
    }, { exact: true });

    family.finish();
});
