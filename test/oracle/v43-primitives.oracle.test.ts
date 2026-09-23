// Replays oracle/cpp/cases/v43-primitives.cpp. Keep the two files in the
// same order.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { AlignedBoxBV } from '../../src/AlignedBoxBV.js';
import { AlignedBoxTreeOfPoints } from '../../src/AlignedBoxTreeOfPoints.js';
import { AlignedBoxTreeOfSegments } from '../../src/AlignedBoxTreeOfSegments.js';
import { AlignedBoxTreeOfTriangles } from '../../src/AlignedBoxTreeOfTriangles.js';
import { BVTree, BVTreeNode } from '../../src/BVTree.js';
import {
    intersectLineTriangle, intersectRayTriangle, intersectSegmentTriangle,
    type BVTreeOfTrianglesIntersection
} from '../../src/BVTreeOfTriangles.js';
import { Cone } from '../../src/Cone.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { IndexAttribute } from '../../src/IndexAttribute.js';
import { Matrix } from '../../src/Matrix.js';
import { MeshDescription, MeshTopology } from '../../src/Mesh.js';
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
import { Triangle } from '../../src/Triangle.js';
import { Vector, add } from '../../src/Vector.js';
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

    // ---------------------------------------------------------- Cone

    // cos/sin/tan come from the MSVC runtime in C++ and from V8 here, so the
    // derived members carry the default tolerance; the angle itself is
    // arithmetic and is held to bit-identity.
    family.case('Cone.setAngle', (io) => {
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        const cone = Cone.fromRayAngle(ray, io.real());
        io.outRealExact(cone.angle);
        io.outReal(cone.cosAngle);
        io.outReal(cone.sinAngle);
        io.outReal(cone.tanAngle);
        io.outReal(cone.cosAngleSqr);
        io.outReal(cone.sinAngleSqr);
        io.outReal(cone.invSinAngle);
    });

    family.case('Cone.heights', (io) => {
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        const angle = io.real();
        const mode = io.integer();
        const h0 = io.real();
        const h1 = io.real();
        const probe = io.reals(4);
        const cone = Cone.fromRayAngle(ray, angle);
        if (mode === 1) {
            cone.makeInfiniteTruncatedCone(h0);
        } else if (mode === 2) {
            cone.makeFiniteCone(h1 > 0 ? h1 : 1);
        } else if (mode === 3) {
            cone.makeConeFrustum(h0, h0 + (h1 > 0 ? h1 : 1));
        } else if (mode === 4) {
            cone.makeInfiniteCone();
        }
        io.outReal(cone.getMinHeight());
        io.outReal(cone.getMaxHeight());
        io.outBool(cone.isFinite());
        io.outBool(cone.isInfinite());
        for (let i = 0; i < 4; ++i) {
            io.outBool(cone.heightInRange(probe[i]));
            io.outBool(cone.heightLessThanMin(probe[i]));
            io.outBool(cone.heightGreaterThanMax(probe[i]));
        }
    }, { exact: true });

    family.case('Cone.heights.throwParity', (io) => {
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        const angle = io.real();
        const mode = io.integer();
        const h0 = io.real();
        const h1 = io.real();
        const cone = Cone.fromRayAngle(ray, angle);
        if (mode === 0) {
            cone.makeInfiniteTruncatedCone(h0);
        } else if (mode === 1) {
            cone.makeFiniteCone(h1);
        } else {
            cone.makeConeFrustum(h0, h1);
        }
        io.outReal(cone.getMinHeight());
        io.outReal(cone.getMaxHeight());
    }, { exact: true });

    family.case('Cone.compare', (io) => {
        const cone0 = Cone.fromRayAngleMinMaxHeight(
            Ray.fromOriginDirection(io.vec(3), io.vec(3)), io.real(),
            io.real(), io.real());
        const cone1 = Cone.fromRayAngleMinMaxHeight(
            Ray.fromOriginDirection(io.vec(3), io.vec(3)), io.real(),
            io.real(), io.real());
        io.outBool(cone0.equals(cone1));
        io.outBool(cone0.notEquals(cone1));
        io.outBool(cone0.lessThan(cone1));
        io.outBool(cone0.lessThanOrEqual(cone1));
        io.outBool(cone0.greaterThan(cone1));
        io.outBool(cone0.greaterThanOrEqual(cone1));
    }, { exact: true });

    // std::cos/std::sin/std::tan decide the mesh vertices, so they carry the
    // default tolerance. The mesh structure is integer output and must match.
    family.case('Cone.createMesh', (io) => {
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        const angle = io.real();
        const hMin = io.real();
        const hMax = io.real();
        const inscribed = io.boolean();
        const cone = Cone.fromRayAngleMinMaxHeight(ray, angle, hMin, hMax);
        const { vertices, indices } = cone.createMesh(3, inscribed);
        io.outInt(vertices.length);
        io.outInt(indices.length);
        for (const i of indices) { io.outInt(i); }
        for (const v of vertices) { io.outVec(v); }
    });

    // ------------------------------------------------------ Polygon2

    family.case('Polygon2.queries', (io) => {
        const n = io.integer();
        io.integer();  // 'mode', used only by the generator
        const pool: Vector[] = [];
        for (let i = 0; i <= n; ++i) { pool.push(io.vec(2)); }
        const indices: number[] = [];
        for (let i = 0; i < n; ++i) { indices.push(io.integer()); }
        const polygon = new Polygon2(pool, indices, true);
        io.outBool(polygon.isValid());
        io.outInt(polygon.getVertices().length);
        for (const v of polygon.getVertices()) { io.outInt(v); }
        io.outInt(polygon.getIndices().length);
        for (const v of polygon.getIndices()) { io.outInt(v); }
        io.outBool(polygon.counterClockwise());
        io.outVec(polygon.computeVertexAverage());
        io.outReal(polygon.computePerimeterLength());
        io.outReal(polygon.computeArea());
        io.outBool(polygon.isSimple());
        io.outBool(polygon.isConvex());
    }, { exact: true });

    family.case('Polygon2.queries.clockwise', (io) => {
        const n = io.integer();
        const pool: Vector[] = [];
        for (let i = 0; i < n; ++i) { pool.push(io.vec(2)); }
        const indices: number[] = [];
        for (let i = 0; i < n; ++i) { indices.push(io.integer()); }
        const polygon = new Polygon2(pool, indices, false);
        io.outBool(polygon.isValid());
        io.outVec(polygon.computeVertexAverage());
        io.outReal(polygon.computePerimeterLength());
        io.outReal(polygon.computeArea());
        io.outBool(polygon.isSimple());
        io.outBool(polygon.isConvex());
    }, { exact: true });

    // -------------------------------------------------- Tetrahedron3

    const tetrahedron = (io: OracleIO): Tetrahedron3 =>
        Tetrahedron3.fromVertices(io.vec(3), io.vec(3), io.vec(3), io.vec(3));

    family.case('Tetrahedron3.normals', (io) => {
        const t = tetrahedron(io);
        for (let i = 0; i < 4; ++i) { io.outVec(t.computeFaceNormal(i)); }
        for (let i = 0; i < 6; ++i) { io.outVec(t.computeEdgeNormal(i)); }
        for (let i = 0; i < 4; ++i) { io.outVec(t.computeVertexNormal(i)); }
    }, { exact: true });

    family.case('Tetrahedron3.computeCentroid', (io) => {
        io.outVec(tetrahedron(io).computeCentroid());
    }, { exact: true });

    family.case('Tetrahedron3.getPlanes', (io) => {
        const plane = tetrahedron(io).getPlanes();
        for (let i = 0; i < 4; ++i) {
            io.outVec(plane[i].normal);
            io.outReal(plane[i].constant);
        }
    }, { exact: true });

    // Issue #268: GetPlanes never writes Plane3::origin, so upstream leaves
    // the caller's buffer untouched and the returned planes violate the
    // Hyperplane invariant Dot(normal, origin) == constant. The port builds
    // each plane from its normal and constant.
    family.case('Tetrahedron3.getPlanes.deviation', (io) => {
        const t = tetrahedron(io);
        io.vec(3);   // the stale normal, used only by the C++ generator
        io.real();   // the stale constant
        const plane = t.getPlanes();
        for (let i = 0; i < 4; ++i) { io.outVec(plane[i].origin); }
    }, { exact: true, deviation: '#268' });

    family.case('Tetrahedron3.tables', (io) => {
        const face = io.integer();
        const edge = io.integer();
        const vertex = io.integer();
        for (const i of Tetrahedron3.getFaceIndices(face)) { io.outInt(i); }
        for (const i of Tetrahedron3.getAllFaceIndices()) { io.outInt(i); }
        for (const i of Tetrahedron3.getEdgeIndices(edge)) { io.outInt(i); }
        for (const i of Tetrahedron3.getAllEdgeIndices()) { io.outInt(i); }
        for (const i of Tetrahedron3.getEdgeAugmented(edge)) { io.outInt(i); }
        for (const i of Tetrahedron3.getVertexAugmented(vertex)) { io.outInt(i); }
    }, { exact: true });

    family.case('Tetrahedron3.compare', (io) => {
        const t0 = tetrahedron(io);
        const t1 = tetrahedron(io);
        io.outBool(t0.equals(t1));
        io.outBool(t0.notEquals(t1));
        io.outBool(t0.lessThan(t1));
        io.outBool(t0.lessThanOrEqual(t1));
        io.outBool(t0.greaterThan(t1));
        io.outBool(t0.greaterThanOrEqual(t1));
    }, { exact: true });

    family.case('Tetrahedron3.defaultConstruct', (io) => {
        const v = [io.vec(3), io.vec(3), io.vec(3), io.vec(3)];
        const canonical = new Tetrahedron3();
        const fromFour = Tetrahedron3.fromVertices(v[0], v[1], v[2], v[3]);
        for (let i = 0; i < 4; ++i) { io.outVec(canonical.v[i]); }
        for (let i = 0; i < 4; ++i) { io.outVec(fromFour.v[i]); }
        io.outVec(canonical.computeCentroid());
    }, { exact: true });

    // ---------------------------------------------- RectangleManager

    function drawRectangles(io: OracleIO, n: number): AlignedBox[] {
        const rectangles: AlignedBox[] = [];
        for (let i = 0; i < n; ++i) {
            const x0 = io.real();
            const y0 = io.real();
            const w = io.real();
            const h = io.real();
            rectangles.push(AlignedBox.fromMinMax(
                Vector.fromArray([x0, y0]),
                Vector.fromArray([x0 + w, y0 + h])));
        }
        return rectangles;
    }

    // getOverlap() returns the keys in the order the upstream std::set
    // iterates them, lexicographically by (V[0], V[1]).
    function outOverlap(io: OracleIO, manager: RectangleManager): void {
        const overlap = manager.getOverlap();
        io.outInt(overlap.length);
        for (const e of overlap) {
            io.outInt(e.V[0]);
            io.outInt(e.V[1]);
        }
    }

    family.case('RectangleManager.initialize', (io) => {
        const n = io.integer();
        const rectangles = drawRectangles(io, n);
        const manager = new RectangleManager(rectangles);
        outOverlap(io, manager);
        for (let i = 0; i < n; ++i) {
            const box = manager.getRectangle(i);
            io.outVec(box.min);
            io.outVec(box.max);
        }
    }, { exact: true });

    family.case('RectangleManager.update', (io) => {
        const n = io.integer();
        const rectangles = drawRectangles(io, n);
        const manager = new RectangleManager(rectangles);
        outOverlap(io, manager);
        for (let pass = 0; pass < 2; ++pass) {
            for (let i = 0; i < n; ++i) {
                const dx = io.real();
                const dy = io.real();
                const box = manager.getRectangle(i);
                box.min.values[0] += dx;
                box.max.values[0] += dx;
                box.min.values[1] += dy;
                box.max.values[1] += dy;
                manager.setRectangle(i, box);
            }
            manager.update();
            outOverlap(io, manager);
        }
    }, { exact: true });

    // ------------------------------------- RectangleMesh and patches

    function rectangle3(io: OracleIO): Rectangle {
        const center = io.vec(3);
        const axis = frame3(io);
        const extent = io.vec(2);
        return Rectangle.fromCenterAxisExtent(center, [axis[0], axis[1]], extent);
    }

    function channel(numVertices: number, numComponents: number):
        { data: Float64Array, attribute: (s: string) => VertexAttribute } {
        const data = new Float64Array(numVertices * numComponents);
        return {
            data,
            attribute: (semantic: string) =>
                new VertexAttribute(semantic, data, 8 * numComponents)
        };
    }

    function outBuffer(io: OracleIO, data: Float64Array): void {
        for (const x of data) { io.outReal(x); }
    }

    family.case('RectangleMesh.construct', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const rectangle = rectangle3(io);
        const description = new MeshDescription(MeshTopology.RECTANGLE,
            numRows, numCols);
        const positions = channel(description.numVertices, 3);
        const normals = channel(description.numVertices, 3);
        const indices = new Uint32Array(3 * description.numTriangles);
        description.vertexAttributes.push(positions.attribute('position'));
        description.vertexAttributes.push(normals.attribute('normal'));
        description.indexAttribute = new IndexAttribute(indices, 4);
        const mesh = new RectangleMesh(description, rectangle);
        io.outInt(mesh.getDescription().numVertices);
        io.outInt(mesh.getDescription().numTriangles);
        io.outBool(mesh.getDescription().constructed);
        outBuffer(io, positions.data);
        outBuffer(io, normals.data);
        for (const i of indices) { io.outInt(i); }
        io.outVec(mesh.getRectangle().center);
    }, { exact: true });

    family.case('RectangleMesh.tcoords', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const rectangle = rectangle3(io);
        const wantCCW = io.boolean();
        const description = new MeshDescription(MeshTopology.RECTANGLE,
            numRows, numCols);
        description.wantCCW = wantCCW;
        const positions = channel(description.numVertices, 3);
        const tcoords = channel(description.numVertices, 2);
        const indices = new Uint32Array(3 * description.numTriangles);
        description.vertexAttributes.push(positions.attribute('position'));
        description.vertexAttributes.push(tcoords.attribute('tcoord'));
        description.indexAttribute = new IndexAttribute(indices, 4);
        void new RectangleMesh(description, rectangle);
        outBuffer(io, tcoords.data);
        outBuffer(io, positions.data);
        for (const i of indices) { io.outInt(i); }
    }, { exact: true });

    // Issue #268: RectangleMesh::InitializeFrame hardcodes tangent = (1,0,0)
    // and bitangent = (0,1,0) regardless of the rectangle's axes, so the
    // tangent, bitangent, dpdu and dpdv channels are not in the rectangle's
    // plane. The port uses the rectangle's own orthonormal axes. The normal
    // channel is unaffected and is compared in the main case above.
    family.case('RectangleMesh.frame.deviation', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const rectangle = rectangle3(io);
        const description = new MeshDescription(MeshTopology.RECTANGLE,
            numRows, numCols);
        description.wantDynamicTangentSpaceUpdate = true;
        const positions = channel(description.numVertices, 3);
        const names = ['normal', 'tangent', 'bitangent', 'dpdu', 'dpdv'];
        const channels = names.map(() => channel(description.numVertices, 3));
        const indices = new Uint32Array(3 * description.numTriangles);
        description.vertexAttributes.push(positions.attribute('position'));
        names.forEach((name, i) => {
            description.vertexAttributes.push(channels[i].attribute(name));
        });
        description.indexAttribute = new IndexAttribute(indices, 4);
        const mesh = new RectangleMesh(description, rectangle);
        io.outBool(mesh.getDescription().allowUpdateFrame);
        io.outBool(mesh.getDescription().hasTangentSpaceVectors);
        for (const c of channels) { outBuffer(io, c.data); }
    }, { exact: true, deviation: '#268' });

    // The port of the QuadraticGraphSurface defined in the case file. The
    // expression grouping matches the C++ verbatim.
    class QuadraticGraphSurface extends ParametricSurface {
        private readonly mA: readonly number[];

        constructor(umin: number, umax: number, vmin: number, vmax: number,
            a: readonly number[]) {
            super(3, umin, umax, vmin, vmax, true);
            this.mA = a;
            this.mConstructed = true;
        }

        evaluate(u: number, v: number, order: number, jet: Vector[]): void {
            const a = this.mA;
            const z = ((((a[0] + a[1] * u) + a[2] * v) + a[3] * (u * u))
                + a[4] * (u * v)) + a[5] * (v * v);
            jet[0].values[0] = u;
            jet[0].values[1] = v;
            jet[0].values[2] = z;
            if (order >= 1) {
                const zu = (a[1] + (a[3] + a[3]) * u) + a[4] * v;
                const zv = (a[2] + a[4] * u) + (a[5] + a[5]) * v;
                jet[1].values[0] = 1; jet[1].values[1] = 0; jet[1].values[2] = zu;
                jet[2].values[0] = 0; jet[2].values[1] = 1; jet[2].values[2] = zv;
            }
            if (order >= 2) {
                jet[3].values[0] = 0; jet[3].values[1] = 0;
                jet[3].values[2] = a[3] + a[3];
                jet[4].values[0] = 0; jet[4].values[1] = 0;
                jet[4].values[2] = a[4];
                jet[5].values[0] = 0; jet[5].values[1] = 0;
                jet[5].values[2] = a[5] + a[5];
            }
        }
    }

    function graphSurface(io: OracleIO): QuadraticGraphSurface {
        const umin = io.real();
        const umax = io.real();
        const vmin = io.real();
        const vmax = io.real();
        const a = io.reals(6);
        return new QuadraticGraphSurface(umin, umax, vmin, vmax, a);
    }

    family.case('RectanglePatchMesh.construct', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const surface = graphSurface(io);
        const description = new MeshDescription(MeshTopology.RECTANGLE,
            numRows, numCols);
        const positions = channel(description.numVertices, 3);
        const normals = channel(description.numVertices, 3);
        const tcoords = channel(description.numVertices, 2);
        const indices = new Uint32Array(3 * description.numTriangles);
        description.vertexAttributes.push(positions.attribute('position'));
        description.vertexAttributes.push(normals.attribute('normal'));
        description.vertexAttributes.push(tcoords.attribute('tcoord'));
        description.indexAttribute = new IndexAttribute(indices, 4);
        void new RectanglePatchMesh(description, surface);
        outBuffer(io, tcoords.data);
        outBuffer(io, positions.data);
        outBuffer(io, normals.data);
        for (const i of indices) { io.outInt(i); }
    }, { exact: true });

    family.case('RectanglePatchMesh.frame', (io) => {
        const numRows = io.integer();
        const numCols = io.integer();
        const surface = graphSurface(io);
        const description = new MeshDescription(MeshTopology.RECTANGLE,
            numRows, numCols);
        description.wantDynamicTangentSpaceUpdate = true;
        const positions = channel(description.numVertices, 3);
        const names = ['normal', 'tangent', 'bitangent', 'dpdu', 'dpdv'];
        const channels = names.map(() => channel(description.numVertices, 3));
        const indices = new Uint32Array(3 * description.numTriangles);
        description.vertexAttributes.push(positions.attribute('position'));
        names.forEach((name, i) => {
            description.vertexAttributes.push(channels[i].attribute(name));
        });
        description.indexAttribute = new IndexAttribute(indices, 4);
        const mesh = new RectanglePatchMesh(description, surface);
        io.outBool(mesh.getDescription().allowUpdateFrame);
        io.outBool(mesh.getDescription().hasTangentSpaceVectors);
        for (const c of channels) { outBuffer(io, c.data); }
    }, { exact: true });

    // ------------------------------- AlignedBoxBV and OrientedBoxBV

    family.case('AlignedBoxBV.queries', (io) => {
        const useDefault = io.boolean();
        const boxMin = io.vec(3);
        const boxExtent = io.vec(3);
        const P = io.vec(3);
        const D = io.vec(3);
        const R = io.vec(3);
        const bv = new AlignedBoxBV();
        if (!useDefault) {
            bv.box.min = boxMin;
            bv.box.max = add(boxMin, boxExtent);
        }
        const axis = bv.getSplittingAxis();
        io.outVec(bv.box.min);
        io.outVec(bv.box.max);
        io.outVec(axis.origin);
        io.outVec(axis.direction);
        io.outBool(AlignedBoxBV.intersectLine(P, D, bv));
        io.outBool(AlignedBoxBV.intersectRay(P, D, bv));
        io.outBool(AlignedBoxBV.intersectSegment(P, R, bv));
    }, { exact: true });

    family.case('OrientedBoxBV.queries', (io) => {
        const useDefault = io.boolean();
        const center = io.vec(3);
        const frame = frame3(io);
        const extent = io.vec(3);
        const P = io.vec(3);
        const D = io.vec(3);
        const R = io.vec(3);
        const bv = new OrientedBoxBV();
        if (!useDefault) {
            bv.box = OrientedBox.fromCenterAxisExtent(center, frame, extent);
        }
        const axis = bv.getSplittingAxis();
        io.outVec(bv.box.center);
        io.outVec(bv.box.axis[0]);
        io.outVec(bv.box.axis[1]);
        io.outVec(bv.box.axis[2]);
        io.outVec(bv.box.extent);
        io.outVec(axis.origin);
        io.outVec(axis.direction);
        io.outBool(OrientedBoxBV.intersectLine(P, D, bv));
        io.outBool(OrientedBoxBV.intersectRay(P, D, bv));
        io.outBool(OrientedBoxBV.intersectSegment(P, R, bv));
    }, { exact: true });

    // --------------------------------------- Bounding-volume trees

    function nodeIndexOut(i: number): number {
        return i === BVTreeNode.invalid ? -1 : i;
    }

    function vertices(io: OracleIO, count: number): Vector[] {
        const v: Vector[] = [];
        for (let i = 0; i < count; ++i) { v.push(io.vec(3)); }
        return v;
    }

    function linear(io: OracleIO): { P: Vector, D: Vector, R: Vector } {
        return { P: io.vec(3), D: io.vec(3), R: io.vec(3) };
    }

    function treeHeight(io: OracleIO): number {
        const height = io.integer();
        return height < 0 ? BVTree.fullHeight : height;
    }

    type AnyTree = AlignedBoxTreeOfPoints | AlignedBoxTreeOfSegments
        | AlignedBoxTreeOfTriangles | OrientedBoxTreeOfPoints
        | OrientedBoxTreeOfSegments | OrientedBoxTreeOfTriangles;

    function outTreeStructure(io: OracleIO, tree: AnyTree): void {
        io.outInt(tree.getHeight());
        const partition = tree.getPartition();
        io.outInt(partition.length);
        for (const p of partition) { io.outInt(p); }
        const nodes = tree.getNodes();
        io.outInt(nodes.length);
        for (const n of nodes) {
            io.outInt(nodeIndexOut(n.minIndex));
            io.outInt(nodeIndexOut(n.maxIndex));
            io.outInt(nodeIndexOut(n.leftChild));
            io.outInt(nodeIndexOut(n.rightChild));
        }
    }

    function outAlignedBoxes(io: OracleIO,
        tree: AlignedBoxTreeOfPoints | AlignedBoxTreeOfSegments
            | AlignedBoxTreeOfTriangles): void {
        for (const n of tree.getNodes()) {
            io.outVec(n.boundingVolume.box.min);
            io.outVec(n.boundingVolume.box.max);
        }
    }

    // 'skipLeafExtent' omits the extents of leaf nodes, which upstream's
    // OrientedBoxTreeOfTriangles collapses along the wrong axis (issue #343).
    function outOrientedBoxes(io: OracleIO,
        tree: OrientedBoxTreeOfPoints | OrientedBoxTreeOfSegments
            | OrientedBoxTreeOfTriangles, skipLeafExtent: boolean): void {
        for (const n of tree.getNodes()) {
            const box = n.boundingVolume.box;
            io.outVec(box.center);
            io.outVec(box.axis[0]);
            io.outVec(box.axis[1]);
            io.outVec(box.axis[2]);
            const isLeaf = (n.leftChild === BVTreeNode.invalid);
            if (!(skipLeafExtent && isLeaf)) { io.outVec(box.extent); }
        }
    }

    function outTreeQueries(io: OracleIO,
        tree: AlignedBoxTreeOfPoints | AlignedBoxTreeOfSegments
            | OrientedBoxTreeOfPoints | OrientedBoxTreeOfSegments,
        L: { P: Vector, D: Vector, R: Vector }): void {
        for (const [queryType, Q] of [[BVTree.LINE_QUERY, L.D],
            [BVTree.RAY_QUERY, L.D], [BVTree.SEGMENT_QUERY, L.R]] as
            [number, Vector][]) {
            const nodeIndices = tree.execute(queryType, L.P, Q);
            io.outInt(nodeIndices.length);
            for (const i of nodeIndices) { io.outInt(i); }
        }
    }

    family.case('AlignedBoxTreeOfPoints.create', (io) => {
        const count = io.integer();
        io.integer();  // 'mode', used only by the generator
        const height = treeHeight(io);
        const v = vertices(io, count);
        const L = linear(io);
        const tree = new AlignedBoxTreeOfPoints();
        tree.create(v, height);
        outTreeStructure(io, tree);
        outAlignedBoxes(io, tree);
        outTreeQueries(io, tree, L);
    }, { exact: true });

    family.case('OrientedBoxTreeOfPoints.create', (io) => {
        const count = io.integer();
        io.integer();
        const height = treeHeight(io);
        const v = vertices(io, count);
        const L = linear(io);
        const tree = new OrientedBoxTreeOfPoints();
        tree.create(v, height);
        outTreeStructure(io, tree);
        outOrientedBoxes(io, tree, false);
        outTreeQueries(io, tree, L);
    }, { exact: true });

    function segments(io: OracleIO, count: number): [number, number][] {
        const s: [number, number][] = [];
        for (let i = 0; i < count; ++i) { s.push([io.integer(), io.integer()]); }
        return s;
    }

    family.case('AlignedBoxTreeOfSegments.create', (io) => {
        const numVertices = io.integer();
        const numSegments = io.integer();
        io.integer();
        const height = treeHeight(io);
        const v = vertices(io, numVertices);
        const s = segments(io, numSegments);
        const L = linear(io);
        const tree = new AlignedBoxTreeOfSegments();
        tree.createFromSegments(v, s, height);
        outTreeStructure(io, tree);
        outAlignedBoxes(io, tree);
        outTreeQueries(io, tree, L);
        for (const c of tree.getCentroids()) { io.outVec(c); }
    }, { exact: true });

    family.case('OrientedBoxTreeOfSegments.create', (io) => {
        const numVertices = io.integer();
        const numSegments = io.integer();
        io.integer();
        const height = treeHeight(io);
        const v = vertices(io, numVertices);
        const s = segments(io, numSegments);
        const L = linear(io);
        const tree = new OrientedBoxTreeOfSegments();
        tree.createFromSegments(v, s, height);
        outTreeStructure(io, tree);
        outOrientedBoxes(io, tree, false);
        outTreeQueries(io, tree, L);
        for (const c of tree.getCentroids()) { io.outVec(c); }
    }, { exact: true });

    function triangles(io: OracleIO, count: number): [number, number, number][] {
        const t: [number, number, number][] = [];
        for (let i = 0; i < count; ++i) {
            t.push([io.integer(), io.integer(), io.integer()]);
        }
        return t;
    }

    // The port of the case file's HasCoincidentParameters probe: true when
    // two triangles are hit at bit-identical parameters, which upstream's
    // std::set<Intersection> collapses (issue #167).
    function hasCoincidentParameters(queryType: number, A: Vector, B: Vector,
        v: readonly Vector[], t: readonly [number, number, number][]): boolean {
        const query = [intersectLineTriangle, intersectRayTriangle,
            intersectSegmentTriangle][queryType];
        const params: number[] = [];
        for (const tri of t) {
            const triangle = Triangle.fromVertices(v[tri[0]], v[tri[1]], v[tri[2]]);
            const r = query(A, B, triangle);
            if (r.intersect) { params.push(r.parameter); }
        }
        for (let i = 0; i < params.length; ++i) {
            for (let j = i + 1; j < params.length; ++j) {
                if (params[i] === params[j]) { return true; }
            }
        }
        return false;
    }

    function outIntersections(io: OracleIO,
        intersections: readonly BVTreeOfTrianglesIntersection[]): void {
        io.outInt(intersections.length);
        for (const it of intersections) {
            io.outInt(it.triangleIndex);
            io.outReal(it.parameter);
            io.outVec(it.point);
        }
    }

    function outTriangleTreeQueries(io: OracleIO,
        tree: AlignedBoxTreeOfTriangles | OrientedBoxTreeOfTriangles,
        L: { P: Vector, D: Vector, R: Vector }, v: readonly Vector[],
        t: readonly [number, number, number][]): void {
        for (let queryType = 0; queryType < 3; ++queryType) {
            const B = (queryType === 2 ? L.R : L.D);
            const result = tree.execute(queryType, L.P, B);
            io.outInt(result.nodeIndices.length);
            for (const i of result.nodeIndices) { io.outInt(i); }
            const coincident = hasCoincidentParameters(queryType, L.P, B, v, t);
            io.outBool(coincident);
            if (!coincident) { outIntersections(io, result.intersections); }
        }
    }

    family.case('AlignedBoxTreeOfTriangles.create', (io) => {
        const numVertices = io.integer();
        const numTriangles = io.integer();
        io.integer();
        const height = treeHeight(io);
        const v = vertices(io, numVertices);
        const t = triangles(io, numTriangles);
        const L = linear(io);
        const tree = new AlignedBoxTreeOfTriangles();
        tree.createFromTriangles(v, t, height);
        outTreeStructure(io, tree);
        outAlignedBoxes(io, tree);
        outTriangleTreeQueries(io, tree, L, v, t);
        for (const c of tree.getCentroids()) { io.outVec(c); }
    }, { exact: true });

    family.case('OrientedBoxTreeOfTriangles.create', (io) => {
        const numVertices = io.integer();
        const numTriangles = io.integer();
        io.integer();
        const height = treeHeight(io);
        const v = vertices(io, numVertices);
        const t = triangles(io, numTriangles);
        const L = linear(io);
        const tree = new OrientedBoxTreeOfTriangles();
        tree.createFromTriangles(v, t, height);
        outTreeStructure(io, tree);
        outOrientedBoxes(io, tree, true);
        outTriangleTreeQueries(io, tree, L, v, t);
        for (const c of tree.getCentroids()) { io.outVec(c); }
    }, { exact: true });

    // Issue #343: the smallest-extent scan of
    // OrientedBoxTreeOfTriangles::ComputeLeafBoundingVolume ends with
    // 'absExtent > minAbsExtent', so upstream zeroes the largest extent and
    // the leaf box collapses along its longest axis. The port uses '<'.
    family.case('OrientedBoxTreeOfTriangles.leafExtent.deviation', (io) => {
        const numVertices = io.integer();
        const numTriangles = io.integer();
        io.integer();
        const v = vertices(io, numVertices);
        const t = triangles(io, numTriangles);
        const tree = new OrientedBoxTreeOfTriangles();
        tree.createFromTriangles(v, t);
        for (const n of tree.getNodes()) {
            if (n.leftChild === BVTreeNode.invalid) {
                io.outVec(n.boundingVolume.box.extent);
            }
        }
    }, { exact: true, deviation: '#343' });

    // Issue #167: BVTreeOfTriangles::Execute collects hits in a
    // std::set<Intersection> ordered by parameter alone, so two triangles hit
    // at the same parameter are set-equivalent and all but one are dropped.
    // The port orders by (parameter, triangleIndex) and keeps both.
    family.case('BVTreeOfTriangles.coincident.deviation', (io) => {
        const numVertices = io.integer();
        const v = vertices(io, numVertices);
        const i0 = io.integer();
        const i1 = io.integer();
        const i2 = io.integer();
        const t: [number, number, number][] = [[i0, i1, i2], [i0, i1, i2]];
        const P = io.vec(3);
        const d = io.vec(3);
        const tree = new AlignedBoxTreeOfTriangles();
        tree.createFromTriangles(v, t);
        const result = tree.execute(BVTree.LINE_QUERY, P, d);
        outIntersections(io, result.intersections);
    }, { exact: true, deviation: '#167' });

    family.finish();
});
