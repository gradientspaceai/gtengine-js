// Replays oracle/cpp/cases/v43-primitives.cpp. Keep the two files in the
// same order.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { AlignedBoxBV } from '../../src/AlignedBoxBV.js';
import { AlignedBoxTreeOfPoints } from '../../src/AlignedBoxTreeOfPoints.js';
import { AlignedBoxTreeOfSegments } from '../../src/AlignedBoxTreeOfSegments.js';
import { AlignedBoxTreeOfTriangles } from '../../src/AlignedBoxTreeOfTriangles.js';
import { BVTree, BVTreeNode } from '../../src/BVTree.js';
import type { BVTreeOfTrianglesIntersection } from '../../src/BVTreeOfTriangles.js';
import { Cone } from '../../src/Cone.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { IndexAttribute } from '../../src/IndexAttribute.js';
import { Matrix } from '../../src/Matrix.js';
import { Mesh, MeshDescription, MeshTopology } from '../../src/Mesh.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { OrientedBoxBV } from '../../src/OrientedBoxBV.js';
import { OrientedBoxTreeOfPoints } from '../../src/OrientedBoxTreeOfPoints.js';
import { OrientedBoxTreeOfSegments } from '../../src/OrientedBoxTreeOfSegments.js';
import { OrientedBoxTreeOfTriangles } from '../../src/OrientedBoxTreeOfTriangles.js';
import { ParametricSurface } from '../../src/ParametricSurface.js';
import { Polygon2 } from '../../src/Polygon2.js';
import { Ray } from '../../src/Ray.js';
import { Rectangle } from '../../src/Rectangle.js';
import { RectangleManager } from '../../src/RectangleManager.js';
import { RectangleMesh } from '../../src/RectangleMesh.js';
import { RectanglePatchMesh } from '../../src/RectanglePatchMesh.js';
import { Tetrahedron3 } from '../../src/Tetrahedron3.js';
import { Vector, add, sub, mul, normalize } from '../../src/Vector.js';
import { VertexAttribute } from '../../src/VertexAttribute.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The ports of the C++ helpers in the case file. Each reads exactly the
// doubles its C++ counterpart recorded.
function frame3(io: OracleIO): Vector[] {
    return [io.vec(3), io.vec(3), io.vec(3)];
}

function frame2(io: OracleIO): Vector[] {
    return [io.vec(2), io.vec(2)];
}

function ellipsoid(io: OracleIO, n: number): Hyperellipsoid {
    const center = io.vec(n);
    const axis = (n === 2 ? frame2(io) : frame3(io));
    const extent = io.vec(n);
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

function symmetric(n: number, upper: readonly number[]): Matrix {
    const A = Matrix.zero(n, n);
    let k = 0;
    for (let r = 0; r < n; ++r) {
        for (let c = r; c < n; ++c, ++k) {
            A.set(r, c, upper[k]);
            A.set(c, r, upper[k]);
        }
    }
    return A;
}

describe('oracle: v43-primitives', () => {
    const family = new OracleFamily('v43-primitives');

    // ------------------------------------------------ Hyperellipsoid

    for (const n of [2, 3]) {
        family.case(`Hyperellipsoid.getM.${n}d`, (io) => {
            const E = ellipsoid(io, n);
            io.outMat(E.getM());
            io.outMat(E.getMInverse());
        }, { exact: true });

        family.case(`Hyperellipsoid.toCoefficients.${n}d`, (io) => {
            const E = ellipsoid(io, n);
            const abc = E.toCoefficientsABC();
            const coeff = E.toCoefficients();
            io.outMat(abc.A);
            io.outVec(abc.B);
            io.outReal(abc.C);
            io.outReals(coeff);
        }, { exact: true });

        family.case(`Hyperellipsoid.fromCoefficients.${n}d`, (io) => {
            const E = ellipsoid(io, n);
            const coeff = E.toCoefficients();
            const F = new Hyperellipsoid(n);
            const valid = F.fromCoefficients(coeff);
            io.outBool(valid);
            if (valid) {
                io.outVec(F.center);
                for (let d = 0; d < n; ++d) { io.outVec(F.axis[d]); }
                io.outVec(F.extent);
            }
        }, { exact: true });
    }

    family.case('Hyperellipsoid.fromCoefficientsABC.2d', (io) => {
        io.boolean();  // 'definite', used only by the generator
        const a00 = io.real();
        const a01 = io.real();
        const a11 = io.real();
        const B = io.vec(2);
        const C = io.real();
        const A = symmetric(2, [a00, a01, a11]);
        const F = new Hyperellipsoid(2);
        const valid = F.fromCoefficientsABC(A, B, C);
        io.outBool(valid);
        if (valid) {
            io.outVec(F.center);
            io.outVec(F.axis[0]);
            io.outVec(F.axis[1]);
            io.outVec(F.extent);
        }
    }, { exact: true });

    family.case('Hyperellipsoid.fromCoefficientsABC.3d', (io) => {
        io.boolean();
        const a00 = io.real();
        const a02 = io.real();
        const a11 = io.real();
        const a12 = io.real();
        const a22 = io.real();
        const B = io.vec(3);
        const C = io.real();
        const a01 = io.real();
        const A = symmetric(3, [a00, a01, a02, a11, a12, a22]);
        const F = new Hyperellipsoid(3);
        const valid = F.fromCoefficientsABC(A, B, C);
        io.outBool(valid);
        if (valid) {
            io.outVec(F.center);
            io.outVec(F.axis[0]);
            io.outVec(F.axis[1]);
            io.outVec(F.axis[2]);
            io.outVec(F.extent);
        }
    }, { exact: true });

    // Issue #80: upstream's SymmetricEigensolver stores the reflection
    // parameter of a degenerate Householder step, which flips the sign of one
    // row of the eigenvector matrix. An axis-aligned ellipsoid produces a
    // diagonal M, which is exactly that case.
    family.case('Hyperellipsoid.fromCoefficients.decoupledDeviation.3d', (io) => {
        const E = ellipsoid(io, 3);
        const coeff = E.toCoefficients();
        const F = new Hyperellipsoid(3);
        const valid = F.fromCoefficients(coeff);
        io.outBool(valid);
        if (valid) {
            io.outVec(F.center);
            for (let d = 0; d < 3; ++d) { io.outVec(F.axis[d]); }
            io.outVec(F.extent);
        }
    }, { exact: true, deviation: '#80' });

    family.case('Hyperellipsoid.compare.2d', (io) => {
        const E0 = Hyperellipsoid.fromCenterAxisExtent(io.vec(2),
            [io.vec(2), io.vec(2)], io.vec(2));
        const E1 = Hyperellipsoid.fromCenterAxisExtent(io.vec(2),
            [io.vec(2), io.vec(2)], io.vec(2));
        io.outBool(E0.equals(E1));
        io.outBool(E0.notEquals(E1));
        io.outBool(E0.lessThan(E1));
        io.outBool(E0.lessThanOrEqual(E1));
        io.outBool(E0.greaterThan(E1));
        io.outBool(E0.greaterThanOrEqual(E1));
    }, { exact: true });

    // ---------------------------------------------------- Hyperplane

    family.case('Hyperplane.construct.3d', (io) => {
        const normal = io.vec(3);
        const constant = io.real();
        const origin = io.vec(3);
        const defaultPlane = new Hyperplane(3);
        const fromConstant = Hyperplane.fromNormalConstant(normal, constant);
        const fromOrigin = Hyperplane.fromNormalOrigin(normal, origin);
        for (const plane of [defaultPlane, fromConstant, fromOrigin]) {
            io.outVec(plane.normal);
            io.outVec(plane.origin);
            io.outReal(plane.constant);
        }
    }, { exact: true });

    family.case('Hyperplane.fromPoints.3d', (io) => {
        const plane = Hyperplane.fromPoints([io.vec(3), io.vec(3), io.vec(3)]);
        io.outVec(plane.normal);
        io.outVec(plane.origin);
        io.outReal(plane.constant);
    }, { exact: true });

    family.case('Hyperplane.compare.3d', (io) => {
        const h0 = new Hyperplane(3);
        h0.normal = io.vec(3);
        h0.origin = io.vec(3);
        h0.constant = io.real();
        const h1 = new Hyperplane(3);
        h1.normal = io.vec(3);
        h1.origin = io.vec(3);
        h1.constant = io.real();
        io.outBool(h0.equals(h1));
        io.outBool(h0.notEquals(h1));
        io.outBool(h0.lessThan(h1));
        io.outBool(h0.lessThanOrEqual(h1));
        io.outBool(h0.greaterThan(h1));
        io.outBool(h0.greaterThanOrEqual(h1));
    }, { exact: true });

    // Upstream cannot construct a hyperplane from points for any N != 3: the
    // SVD constructor rejects N = 2 and Solve rejects the multiplier -1 for
    // N >= 4, so every C++ record throws. See docs/UPSTREAM-FINDINGS.md,
    // issue #217.
    for (const n of [2, 4]) {
        family.case(`Hyperplane.fromPoints.deviation.${n}d`, (io) => {
            const p: Vector[] = [];
            for (let i = 0; i < n; ++i) { p.push(io.vec(n)); }
            const plane = Hyperplane.fromPoints(p);
            io.outVec(plane.normal);
            io.outVec(plane.origin);
            io.outReal(plane.constant);
        }, { exact: true, deviation: '#217' });
    }

    family.finish();
});
