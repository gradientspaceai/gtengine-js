// Replays oracle/cpp/cases/v22-distance.cpp (verify group 22, distance).
// Keep the two files in the same order.
//
// Every generator in the C++ file records the same number of doubles in every
// mode, so these replays read the inputs without knowing which mode produced
// them.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { CanonicalBox } from '../../src/CanonicalBox.js';
import { Circle3 } from '../../src/Circle3.js';
import { Cone } from '../../src/Cone.js';
import { DistAlignedBox3OrientedBox3 } from '../../src/DistAlignedBox3OrientedBox3.js';
import { DistCircle3Circle3 } from '../../src/DistCircle3Circle3.js';
import { DistOrientedBox3Cone3 } from '../../src/DistOrientedBox3Cone3.js';
import { DistOrientedBox3OrientedBox3 } from '../../src/DistOrientedBox3OrientedBox3.js';
import { DistPlane3AlignedBox3 } from '../../src/DistPlane3AlignedBox3.js';
import { DistPlane3OrientedBox3 } from '../../src/DistPlane3OrientedBox3.js';
import { DistPoint3Parallelepiped3 } from '../../src/DistPoint3Parallelepiped3.js';
import { DistPoint3Tetrahedron3 } from '../../src/DistPoint3Tetrahedron3.js';
import { DistRectangle3AlignedBox3 } from '../../src/DistRectangle3AlignedBox3.js';
import { DistRectangle3CanonicalBox3 } from '../../src/DistRectangle3CanonicalBox3.js';
import { DistRectangle3OrientedBox3 } from '../../src/DistRectangle3OrientedBox3.js';
import { DistRectangle3Rectangle3 } from '../../src/DistRectangle3Rectangle3.js';
import { DistTetrahedron3Tetrahedron3 } from '../../src/DistTetrahedron3Tetrahedron3.js';
import { DistTriangle3AlignedBox3 } from '../../src/DistTriangle3AlignedBox3.js';
import { DistTriangle3CanonicalBox3 } from '../../src/DistTriangle3CanonicalBox3.js';
import { DistTriangle3OrientedBox3 } from '../../src/DistTriangle3OrientedBox3.js';
import { DistTriangle3Rectangle3 } from '../../src/DistTriangle3Rectangle3.js';
import { DistTriangle3Triangle3 } from '../../src/DistTriangle3Triangle3.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { Matrix } from '../../src/Matrix.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Parallelepiped3 } from '../../src/Parallelepiped3.js';
import { Ray } from '../../src/Ray.js';
import { Rectangle } from '../../src/Rectangle.js';
import { Tetrahedron3 } from '../../src/Tetrahedron3.js';
import { Triangle } from '../../src/Triangle.js';
import { Vector, add } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The C++ AlignedBoxInputs records the minimum corner and the size; the
// maximum corner is derived.
function alignedBox(io: OracleIO): AlignedBox {
    const lo = io.vec(3);
    const size = io.vec(3);
    return AlignedBox.fromMinMax(lo, add(lo, size));
}

// The C++ OrientedBoxInputs records the center, the three axes and the extent.
function orientedBox(io: OracleIO): OrientedBox {
    const center = io.vec(3);
    const axis: Vector[] = [io.vec(3), io.vec(3), io.vec(3)];
    return OrientedBox.fromCenterAxisExtent(center, axis, io.vec(3));
}

// The C++ RectangleInputs records the center, the two axes and the extent.
function rectangle3(io: OracleIO): Rectangle {
    const center = io.vec(3);
    const axis: Vector[] = [io.vec(3), io.vec(3)];
    return Rectangle.fromCenterAxisExtent(center, axis, io.vec(2));
}

// The C++ TriangleInputs records the three vertices.
function triangle3(io: OracleIO): Triangle {
    const v0 = io.vec(3);
    const v1 = io.vec(3);
    return Triangle.fromVertices(v0, v1, io.vec(3));
}

// The C++ TetrahedronInputs records the four vertices.
function tetrahedron3(io: OracleIO): Tetrahedron3 {
    const v0 = io.vec(3);
    const v1 = io.vec(3);
    const v2 = io.vec(3);
    return Tetrahedron3.fromVertices(v0, v1, v2, io.vec(3));
}

function emitTwoClosest(io: OracleIO, r: {
    distance: number, sqrDistance: number, closest: readonly Vector[]
}): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

function emitCartesian(io: OracleIO, r: {
    distance: number, sqrDistance: number, cartesian: readonly number[],
    closest: readonly Vector[]
}): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.cartesian[0]);
    io.outReal(r.cartesian[1]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

function emitBarycentric3(io: OracleIO, r: {
    distance: number, sqrDistance: number, barycentric: readonly number[],
    closest: readonly Vector[]
}): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.barycentric[0]);
    io.outReal(r.barycentric[1]);
    io.outReal(r.barycentric[2]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

function emitBarycentric4Pair(io: OracleIO, r: {
    distance: number, sqrDistance: number, barycentric0: readonly number[],
    barycentric1: readonly number[], closest: readonly Vector[]
}): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    for (let i = 0; i < 4; ++i) { io.outReal(r.barycentric0[i]); }
    for (let i = 0; i < 4; ++i) { io.outReal(r.barycentric1[i]); }
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

// The C++ MakeConeInputs records the ray origin and direction, the angle,
// Cone::SetAngle's tanAngle, the minimum height and the increment from the
// minimum to the maximum height. tanAngle is replayed into the cone rather
// than recomputed: it is the only libm value the query reads from the cone,
// std::tan and Math.tan belong to Cone.h (a header of another group), and
// the point of this case is the box-cone query itself.
function cone3(io: OracleIO, finite: boolean): Cone {
    const origin = io.vec(3);
    const direction = io.vec(3);
    const angle = io.real();
    const tanAngle = io.real();
    const minHeight = io.real();
    const maxHeight = minHeight + io.real();
    const ray = Ray.fromOriginDirection(origin, direction);
    const cone = finite
        ? Cone.fromRayAngleMinMaxHeight(ray, angle, minHeight, maxHeight)
        : Cone.fromRayAngleMinHeight(ray, angle, minHeight);
    cone.tanAngle = tanAngle;
    return cone;
}

describe('oracle: v22-distance', () => {
    const family = new OracleFamily('v22-distance');

    // ---- DistPlane3AlignedBox3.h, DistPlane3OrientedBox3.h ----
    family.case('DistPlane3AlignedBox3.compute', (io) => {
        const box = alignedBox(io);
        const normal = io.vec(3);
        const plane = Hyperplane.fromNormalOrigin(normal, io.vec(3));
        const r = new DistPlane3AlignedBox3().compute(plane, box);
        emitTwoClosest(io, r);
    }, { exact: true });

    family.case('DistPlane3OrientedBox3.compute', (io) => {
        const box = orientedBox(io);
        const normal = io.vec(3);
        const plane = Hyperplane.fromNormalOrigin(normal, io.vec(3));
        const r = new DistPlane3OrientedBox3().compute(plane, box);
        emitTwoClosest(io, r);
    }, { exact: true });

    // ---- DistRectangle3*Box3.h ----
    family.case('DistRectangle3CanonicalBox3.compute', (io) => {
        const box = CanonicalBox.fromExtent(io.vec(3));
        const rectangle = rectangle3(io);
        const r = new DistRectangle3CanonicalBox3().compute(rectangle, box);
        emitCartesian(io, r);
    }, { exact: true });

    family.case('DistRectangle3AlignedBox3.compute', (io) => {
        const box = alignedBox(io);
        const rectangle = rectangle3(io);
        const r = new DistRectangle3AlignedBox3().compute(rectangle, box);
        emitCartesian(io, r);
    }, { exact: true });

    family.case('DistRectangle3OrientedBox3.compute', (io) => {
        const box = orientedBox(io);
        const rectangle = rectangle3(io);
        const r = new DistRectangle3OrientedBox3().compute(rectangle, box);
        emitCartesian(io, r);
    }, { exact: true });

    // ---- DistTriangle3*Box3.h ----
    family.case('DistTriangle3CanonicalBox3.compute', (io) => {
        const box = CanonicalBox.fromExtent(io.vec(3));
        const triangle = triangle3(io);
        const r = new DistTriangle3CanonicalBox3().compute(triangle, box);
        emitBarycentric3(io, r);
    }, { exact: true });

    family.case('DistTriangle3AlignedBox3.compute', (io) => {
        const box = alignedBox(io);
        const triangle = triangle3(io);
        const r = new DistTriangle3AlignedBox3().compute(triangle, box);
        emitBarycentric3(io, r);
    }, { exact: true });

    family.case('DistTriangle3OrientedBox3.compute', (io) => {
        const box = orientedBox(io);
        const triangle = triangle3(io);
        const r = new DistTriangle3OrientedBox3().compute(triangle, box);
        emitBarycentric3(io, r);
    }, { exact: true });

    // ---- DistRectangle3Rectangle3.h, DistTriangle3Rectangle3.h,
    //      DistTriangle3Triangle3.h ----
    family.case('DistRectangle3Rectangle3.compute', (io) => {
        const rectangle0 = rectangle3(io);
        const rectangle1 = rectangle3(io);
        const r = new DistRectangle3Rectangle3().compute(rectangle0,
            rectangle1);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.cartesian0[0]);
        io.outReal(r.cartesian0[1]);
        io.outReal(r.cartesian1[0]);
        io.outReal(r.cartesian1[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('DistTriangle3Rectangle3.compute', (io) => {
        const triangle = triangle3(io);
        const rectangle = rectangle3(io);
        const r = new DistTriangle3Rectangle3().compute(triangle, rectangle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.barycentric[0]);
        io.outReal(r.barycentric[1]);
        io.outReal(r.barycentric[2]);
        io.outReal(r.cartesian[0]);
        io.outReal(r.cartesian[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('DistTriangle3Triangle3.compute', (io) => {
        const triangle0 = triangle3(io);
        const triangle1 = triangle3(io);
        const r = new DistTriangle3Triangle3().compute(triangle0, triangle1);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.barycentric0[0]);
        io.outReal(r.barycentric0[1]);
        io.outReal(r.barycentric0[2]);
        io.outReal(r.barycentric1[0]);
        io.outReal(r.barycentric1[1]);
        io.outReal(r.barycentric1[2]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    // ---- DistOrientedBox3OrientedBox3.h, DistAlignedBox3OrientedBox3.h ----
    family.case('DistOrientedBox3OrientedBox3.compute', (io) => {
        const box0 = orientedBox(io);
        const box1 = orientedBox(io);
        const r = new DistOrientedBox3OrientedBox3().compute(box0, box1);
        emitTwoClosest(io, r);
    }, { exact: true });

    family.case('DistAlignedBox3OrientedBox3.compute', (io) => {
        const box0 = alignedBox(io);
        const box1 = orientedBox(io);
        const r = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        emitTwoClosest(io, r);
    }, { exact: true });

    // ---- DistPoint3Tetrahedron3.h, DistTetrahedron3Tetrahedron3.h ----
    family.case('DistPoint3Tetrahedron3.compute', (io) => {
        const tetra = tetrahedron3(io);
        const point = io.vec(3);
        const r = new DistPoint3Tetrahedron3().compute(point, tetra);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        for (let i = 0; i < 4; ++i) { io.outReal(r.barycentric[i]); }
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('DistTetrahedron3Tetrahedron3.compute', (io) => {
        const tetra0 = tetrahedron3(io);
        const tetra1 = tetrahedron3(io);
        const r = new DistTetrahedron3Tetrahedron3().compute(tetra0, tetra1);
        emitBarycentric4Pair(io, r);
    }, { exact: true });

    family.case('DistTetrahedron3Tetrahedron3.compute.nested', (io) => {
        const tetra0 = tetrahedron3(io);
        const tetra1 = tetrahedron3(io);
        const r = new DistTetrahedron3Tetrahedron3().compute(tetra0, tetra1);
        emitBarycentric4Pair(io, r);
    }, { exact: true });

    // ---- DistPoint3Parallelepiped3.h ----
    for (const name of ['DistPoint3Parallelepiped3.compute',
        'DistPoint3Parallelepiped3.ctor.leftHanded']) {
        family.case(name, (io) => {
            const center = io.vec(3);
            const axis: Vector[] = [io.vec(3), io.vec(3), io.vec(3)];
            const point = io.vec(3);
            const ppd = Parallelepiped3.fromCenterAxis(center, axis);
            const r = new DistPoint3Parallelepiped3().compute(point, ppd);
            emitTwoClosest(io, r);
        }, { exact: true });
    }

    family.case('DistPoint3Parallelepiped3.getMinimizer', (io) => {
        const A = Matrix.fromArray(3, 3, io.reals(9));
        const Z = io.vec(3);
        const K = new DistPoint3Parallelepiped3().getMinimizer(A, Z);
        io.outVec(K);
    }, { exact: true });

    // ---- DistOrientedBox3Cone3.h ----
    family.case('DistOrientedBox3Cone3.compute.control', (io) => {
        const cone = cone3(io, true);
        const box = orientedBox(io);
        const r = new DistOrientedBox3Cone3().compute(box, cone, {
            maxSubdivisions: 1, maxBisections: 1, epsilon: 10, tolerance: 0
        });
        io.outReal(r.distance);
        io.outVec(r.boxClosestPoint);
        io.outVec(r.coneClosestPoint);
    }, { exact: true });

    family.case('DistOrientedBox3Cone3.compute.deviation', (io) => {
        const cone = cone3(io, true);
        const box = orientedBox(io);
        const r = new DistOrientedBox3Cone3().compute(box, cone);
        io.outReal(r.distance);
        io.outVec(r.boxClosestPoint);
        io.outVec(r.coneClosestPoint);
    }, { deviation: 'UPSTREAM-FINDINGS Minimize1.h item 1, issue #298' });

    family.case('DistOrientedBox3Cone3.compute.infinite', (io) => {
        const cone = cone3(io, false);
        const box = orientedBox(io);
        const r = new DistOrientedBox3Cone3().compute(box, cone);
        io.outReal(r.distance);
        io.outVec(r.boxClosestPoint);
        io.outVec(r.coneClosestPoint);
    }, { deviation: 'UPSTREAM-FINDINGS DistOrientedBox3Cone3.h item 2, issue #298' });

    // ---- DistCircle3Circle3.h ----
    // numClosestPairs, distance, sqrDistance, equidistant and then the pairs.
    // 'firstPairOnly' emits pair 0 alone and omits the count.
    function circle3Circle3(io: OracleIO, firstPairOnly: boolean): void {
        const center0 = io.vec(3);
        const normal0 = io.vec(3);
        const circle0 = Circle3.fromCenterNormalRadius(center0, normal0,
            io.real());
        const center1 = io.vec(3);
        const normal1 = io.vec(3);
        const circle1 = Circle3.fromCenterNormalRadius(center1, normal1,
            io.real());
        const r = new DistCircle3Circle3().compute(circle0, circle1);
        if (!firstPairOnly) { io.outInt(r.numClosestPairs); }
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outBool(r.equidistant);
        const n = firstPairOnly ? 1 : r.numClosestPairs;
        for (let j = 0; j < n; ++j) {
            io.outVec(r.circle0Closest[j]);
            io.outVec(r.circle1Closest[j]);
        }
    }

    // The generator keeps only the configurations on which upstream's own
    // distance agrees with an independent reference to 1e-12 relative, so a
    // difference larger than that band is a finding rather than the
    // conditioning of the port's documented fix of issue #331. The closest
    // point coordinates are more sensitive than the distance the generator
    // screens on; the measured maxima over 2000 deep-run records are 7.6e-12
    // here and 2.3e-12 for the identity-rotation case, so the tolerance is
    // set an order of magnitude above them. See
    // oracle/reports/v22-distance.md.
    family.case('DistCircle3Circle3.compute', (io) => {
        circle3Circle3(io, false);
    }, { tol: 1e-10 });

    family.case('DistCircle3Circle3.compute.parallelPlanes', (io) => {
        circle3Circle3(io, false);
    }, { exact: true });

    // No libm call other than sqrt is on this path, so the residual is the
    // port's sn = +/-sqrt(1 - cs^2) alone; see oracle/reports/v22-distance.md.
    family.case('DistCircle3Circle3.compute.identityRotation', (io) => {
        circle3Circle3(io, false);
    }, { tol: 1e-10 });

    for (const name of ['DistCircle3Circle3.compute.doubleRoot.deviation',
        'DistCircle3Circle3.compute.antiParallel.deviation',
        'DistCircle3Circle3.compute.coaxial.deviation']) {
        family.case(name, (io) => {
            circle3Circle3(io, true);
        }, { deviation: 'UPSTREAM-FINDINGS DistCircle3Circle3.h, issues #331 and #431' });
    }

    family.finish();
});
