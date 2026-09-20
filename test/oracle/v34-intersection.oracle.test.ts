// Replays oracle/cpp/cases/v34-intersection.cpp (verify group 34,
// intersection). Keep the two files in the same order.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { Arc2 } from '../../src/Arc2.js';
import { Capsule } from '../../src/Capsule.js';
import { Circle3 } from '../../src/Circle3.js';
import { Cone } from '../../src/Cone.js';
import { ConvexMesh3 } from '../../src/ConvexMesh3.js';
import type { ConvexMesh3Triangle } from '../../src/ConvexMesh3.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Ray } from '../../src/Ray.js';
import { Segment } from '../../src/Segment.js';
import { SegmentMesh } from '../../src/SegmentMesh.js';
import { Vector } from '../../src/Vector.js';
import { QFNumber } from '../../src/QFNumber.js';
import {
    AreaEllipse2Ellipse2
} from '../../src/IntrAreaEllipse2Ellipse2.js';
import {
    IntrAlignedBox3Cone3TI,
    intrAlignedBox3Cone3ComputeBoxHeightInterval,
    intrAlignedBox3Cone3ConeAxisIntersectsBox,
    intrAlignedBox3Cone3HasPointInsideCone
} from '../../src/IntrAlignedBox3Cone3.js';
import { IntrConvexMesh3Plane3FI } from '../../src/IntrConvexMesh3Plane3.js';
import type { IntrConvexMesh3Plane3FIResult } from '../../src/IntrConvexMesh3Plane3.js';
import {
    IntrLine3Cone3FI, defaultIntrLine3Cone3FIResult,
    intrLine3Cone3FIDoQuery
} from '../../src/IntrLine3Cone3.js';
import type { IntrLine3Cone3FIResult } from '../../src/IntrLine3Cone3.js';
import {
    IntrLine3Plane3FI, IntrLine3Plane3TI,
    defaultIntrLine3Plane3FIResult, intrLine3Plane3FIDoQuery
} from '../../src/IntrLine3Plane3.js';
import {
    IntrOrientedBox2Circle2FI, IntrOrientedBox2Circle2TI
} from '../../src/IntrOrientedBox2Circle2.js';
import { IntrOrientedBox2Cone2TI } from '../../src/IntrOrientedBox2Cone2.js';
import {
    IntrOrientedBox3Sphere3FI, IntrOrientedBox3Sphere3TI
} from '../../src/IntrOrientedBox3Sphere3.js';
import { IntrPlane3Capsule3TI } from '../../src/IntrPlane3Capsule3.js';
import {
    IntrPlane3Circle3FI, IntrPlane3Circle3TI,
    intrPlane3Circle3InfinitePoints
} from '../../src/IntrPlane3Circle3.js';
import type { IntrPlane3Circle3FIResult } from '../../src/IntrPlane3Circle3.js';
import {
    IntrPlane3Cylinder3FI, IntrPlane3Cylinder3TI
} from '../../src/IntrPlane3Cylinder3.js';
import type { IntrPlane3Cylinder3FIResult } from '../../src/IntrPlane3Cylinder3.js';
import { IntrPlane3Ellipsoid3TI } from '../../src/IntrPlane3Ellipsoid3.js';
import { IntrPlane3OrientedBox3TI } from '../../src/IntrPlane3OrientedBox3.js';
import {
    IntrPlane3Sphere3FI, IntrPlane3Sphere3TI
} from '../../src/IntrPlane3Sphere3.js';
import type { IntrPlane3Sphere3FIResult } from '../../src/IntrPlane3Sphere3.js';
import { IntrRay2Arc2FI, IntrRay2Arc2TI } from '../../src/IntrRay2Arc2.js';
import { IntrRay2SegmentMesh2FI } from '../../src/IntrRay2SegmentMesh2.js';
import {
    IntrRay3Ellipsoid3FI, IntrRay3Ellipsoid3TI,
    defaultIntrRay3Ellipsoid3FIResult, intrRay3Ellipsoid3FIDoQuery
} from '../../src/IntrRay3Ellipsoid3.js';
import {
    IntrSegment2Arc2FI, IntrSegment2Arc2TI
} from '../../src/IntrSegment2Arc2.js';
import { IntrSegment2SegmentMesh2FI } from '../../src/IntrSegment2SegmentMesh2.js';
import {
    IntrSegment3Ellipsoid3FI, IntrSegment3Ellipsoid3TI,
    defaultIntrSegment3Ellipsoid3FIResult, intrSegment3Ellipsoid3FIDoQuery
} from '../../src/IntrSegment3Ellipsoid3.js';
import type { IntrLine3Ellipsoid3FIResult } from '../../src/IntrLine3Ellipsoid3.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- input readers, mirroring the C++ generators ----

function segment(io: OracleIO, n: number): Segment {
    const p0 = io.vec(n);
    const p1 = io.vec(n);
    return Segment.fromEndpoints(p0, p1);
}

function sphere(io: OracleIO, n: number): Hypersphere {
    const center = io.vec(n);
    const radius = io.real();
    return Hypersphere.fromCenterRadius(center, radius);
}

// The C++ Pln3 records the normal and then a point on the plane.
function plane3(io: OracleIO): Hyperplane {
    const normal = io.vec(3);
    const origin = io.vec(3);
    return Hyperplane.fromNormalOrigin(normal, origin);
}

function orientedBox(io: OracleIO, n: number): OrientedBox {
    const center = io.vec(n);
    const axis: Vector[] = [];
    for (let i = 0; i < n; ++i) { axis.push(io.vec(n)); }
    const extent = io.vec(n);
    return OrientedBox.fromCenterAxisExtent(center, axis, extent);
}

function ellipsoid(io: OracleIO): Hyperellipsoid {
    const center = io.vec(3);
    const axis = [io.vec(3), io.vec(3), io.vec(3)];
    const extent = io.vec(3);
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

function ellipse2(io: OracleIO): Hyperellipsoid {
    const center = io.vec(2);
    const axis = [io.vec(2), io.vec(2)];
    const extent = io.vec(2);
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

// The C++ Cn<N> records the ray, then the angle and the six values Cone's
// SetAngle derives from it with cos/sin/tan, then the height pair with -1
// standing for an infinite maximum height. The replay assigns the recorded
// trigonometry instead of recomputing it, so that the C math library never
// enters the compared computation.
function cone(io: OracleIO, n: number): Cone {
    const origin = io.vec(n);
    const direction = io.vec(n);
    const c = new Cone(n);
    c.ray = Ray.fromOriginDirection(origin, direction);
    c.angle = io.real();
    c.cosAngle = io.real();
    c.sinAngle = io.real();
    c.tanAngle = io.real();
    c.cosAngleSqr = io.real();
    c.sinAngleSqr = io.real();
    c.invSinAngle = io.real();
    const minHeight = io.real();
    const maxHeight = io.real();
    if (maxHeight < 0) { c.makeInfiniteTruncatedCone(minHeight); }
    else { c.makeConeFrustum(minHeight, maxHeight); }
    return c;
}

function circle3(io: OracleIO): Circle3 {
    const center = io.vec(3);
    const normal = io.vec(3);
    const radius = io.real();
    return Circle3.fromCenterNormalRadius(center, normal, radius);
}

function cylinder3(io: OracleIO): Cylinder3 {
    const origin = io.vec(3);
    const direction = io.vec(3);
    const radius = io.real();
    const height = io.real();
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(origin, direction), radius, height);
}

function arc2(io: OracleIO): Arc2 {
    const center = io.vec(2);
    const radius = io.real();
    const end0 = io.vec(2);
    const end1 = io.vec(2);
    return Arc2.fromCenterRadiusEnds(center, radius, end0, end1);
}

describe('oracle: v34-intersection', () => {
    const family = new OracleFamily('v34-intersection');

    // ============================ IntrLine3Plane3 ========================

    family.case('IntrLine3Plane3.test', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const p = plane3(io);
        const r = new IntrLine3Plane3TI().test(
            Line.fromOriginDirection(origin, direction), p);
        io.outBool(r.intersect);
    }, { exact: true });

    const line3Plane3Find = (io: OracleIO): void => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const p = plane3(io);
        const r = new IntrLine3Plane3FI().find(
            Line.fromOriginDirection(origin, direction), p);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter);
        io.outVec(r.point);
    };

    family.case('IntrLine3Plane3.find', line3Plane3Find, { exact: true });

    family.case('IntrLine3Plane3.doQuery.fi', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const p = plane3(io);
        const r = defaultIntrLine3Plane3FIResult();
        intrLine3Plane3FIDoQuery(origin, direction, p, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter);
        io.outVec(r.point);
    }, { exact: true });

    // The C++ case records the plane normal, then the plane origin, then the
    // line direction and the line origin.
    family.case('IntrLine3Plane3.find.parallel', (io) => {
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const direction = io.vec(3);
        const lineOrigin = io.vec(3);
        const p = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        const r = new IntrLine3Plane3FI().find(
            Line.fromOriginDirection(lineOrigin, direction), p);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter);
        io.outVec(r.point);
    }, { exact: true });

    // ========================== IntrPlane3Capsule3 =======================

    family.case('IntrPlane3Capsule3.test', (io) => {
        const p = plane3(io);
        const sg = segment(io, 3);
        const radius = io.real();
        const capsule = Capsule.fromSegmentRadius(sg, radius);
        const r = new IntrPlane3Capsule3TI().test(p, capsule);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrPlane3Capsule3.test.touching', (io) => {
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const radius = io.real();
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const p = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        const capsule = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), radius);
        const r = new IntrPlane3Capsule3TI().test(p, capsule);
        io.outBool(r.intersect);
    }, { exact: true });

    // ========================= IntrPlane3Ellipsoid3 ======================

    family.case('IntrPlane3Ellipsoid3.test', (io) => {
        const p = plane3(io);
        const e = ellipsoid(io);
        const r = new IntrPlane3Ellipsoid3TI().test(p, e);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrPlane3Ellipsoid3.test.tangent', (io) => {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3), io.vec(3)];
        const extent = io.vec(3);
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const e = Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
        const p = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        const r = new IntrPlane3Ellipsoid3TI().test(p, e);
        io.outBool(r.intersect);
    }, { exact: true });

    // ======================== IntrPlane3OrientedBox3 =====================

    family.case('IntrPlane3OrientedBox3.test', (io) => {
        const p = plane3(io);
        const box = orientedBox(io, 3);
        const r = new IntrPlane3OrientedBox3TI().test(p, box);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrPlane3OrientedBox3.test.touching', (io) => {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3), io.vec(3)];
        const extent = io.vec(3);
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const box = OrientedBox.fromCenterAxisExtent(center, axis, extent);
        const p = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        const r = new IntrPlane3OrientedBox3TI().test(p, box);
        io.outBool(r.intersect);
    }, { exact: true });

    // ========================== IntrPlane3Sphere3 ========================

    family.case('IntrPlane3Sphere3.test', (io) => {
        const p = plane3(io);
        const s = sphere(io, 3);
        const r = new IntrPlane3Sphere3TI().test(p, s);
        io.outBool(r.intersect);
    }, { exact: true });

    const outPlane3Sphere3 = (io: OracleIO,
        r: IntrPlane3Sphere3FIResult): void => {
        io.outBool(r.intersect);
        io.outBool(r.isCircle);
        io.outVec(r.circle.center);
        io.outVec(r.circle.normal);
        io.outReal(r.circle.radius);
        io.outVec(r.point);
    };

    family.case('IntrPlane3Sphere3.find', (io) => {
        const p = plane3(io);
        const s = sphere(io, 3);
        outPlane3Sphere3(io, new IntrPlane3Sphere3FI().find(p, s));
    }, { exact: true });

    family.case('IntrPlane3Sphere3.find.tangent', (io) => {
        const center = io.vec(3);
        const radius = io.real();
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const s = Hypersphere.fromCenterRadius(center, radius);
        const p = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        outPlane3Sphere3(io, new IntrPlane3Sphere3FI().find(p, s));
    }, { exact: true });

    // ======================== IntrOrientedBox2Circle2 ====================

    family.case('IntrOrientedBox2Circle2.test', (io) => {
        const box = orientedBox(io, 2);
        const circle = sphere(io, 2);
        const r = new IntrOrientedBox2Circle2TI().test(box, circle);
        io.outBool(r.intersect);
    }, { exact: true });

    const oBox2Circle2Find = (io: OracleIO): void => {
        const box = orientedBox(io, 2);
        const boxVelocity = io.vec(2);
        const circle = sphere(io, 2);
        const circleVelocity = io.vec(2);
        const r = new IntrOrientedBox2Circle2FI().find(box, boxVelocity,
            circle, circleVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    };

    family.case('IntrOrientedBox2Circle2.find', oBox2Circle2Find, { exact: true });
    family.case('IntrOrientedBox2Circle2.find.aimed', oBox2Circle2Find, { exact: true });

    // ========================= IntrOrientedBox2Cone2 =====================

    family.case('IntrOrientedBox2Cone2.test', (io) => {
        const box = orientedBox(io, 2);
        const c = cone(io, 2);
        const r = new IntrOrientedBox2Cone2TI().test(box, c);
        io.outBool(r.intersect);
    }, { exact: true });

    // The C++ case records the box, then the cone direction, the cone origin
    // and the seven angle-derived values (the cone has no height pair here,
    // because the 2D query does not use one).
    family.case('IntrOrientedBox2Cone2.test.aimed', (io) => {
        const box = orientedBox(io, 2);
        const direction = io.vec(2);
        const origin = io.vec(2);
        const c = new Cone(2);
        c.ray = Ray.fromOriginDirection(origin, direction);
        c.angle = io.real();
        c.cosAngle = io.real();
        c.sinAngle = io.real();
        c.tanAngle = io.real();
        c.cosAngleSqr = io.real();
        c.sinAngleSqr = io.real();
        c.invSinAngle = io.real();
        const r = new IntrOrientedBox2Cone2TI().test(box, c);
        io.outBool(r.intersect);
    }, { exact: true });

    // ======================== IntrOrientedBox3Sphere3 ====================

    family.case('IntrOrientedBox3Sphere3.test', (io) => {
        const box = orientedBox(io, 3);
        const s = sphere(io, 3);
        const r = new IntrOrientedBox3Sphere3TI().test(box, s);
        io.outBool(r.intersect);
    }, { exact: true });

    const oBox3Sphere3Find = (io: OracleIO): void => {
        const box = orientedBox(io, 3);
        const boxVelocity = io.vec(3);
        const s = sphere(io, 3);
        const sphereVelocity = io.vec(3);
        const r = new IntrOrientedBox3Sphere3FI().find(box, boxVelocity, s,
            sphereVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    };

    family.case('IntrOrientedBox3Sphere3.find', oBox3Sphere3Find, { exact: true });
    family.case('IntrOrientedBox3Sphere3.find.aimed', oBox3Sphere3Find, { exact: true });

    // The port probes every candidate piece of the Minkowski sum and keeps
    // the earliest contact; upstream accepts the first probe.
    // docs/UPSTREAM-FINDINGS.md IntrAlignedBox3Sphere3.h, issues #458, #465.
    family.case('IntrOrientedBox3Sphere3.find.probeDeviation', oBox3Sphere3Find,
        { deviation: 'UPSTREAM-FINDINGS IntrAlignedBox3Sphere3.h DoQuery; #458, #465' });

    // ========================== IntrPlane3Circle3 ========================

    family.case('IntrPlane3Circle3.test', (io) => {
        const p = plane3(io);
        const c = circle3(io);
        const r = new IntrPlane3Circle3TI().test(p, c);
        io.outBool(r.intersect);
    }, { exact: true });

    // The C++ side emits the code 3 for the SIZE_MAX "the whole circle"
    // sentinel, which is not representable as a double.
    const outPlane3Circle3 = (io: OracleIO,
        r: IntrPlane3Circle3FIResult): void => {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections === intrPlane3Circle3InfinitePoints
            ? 3 : r.numIntersections);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
        io.outVec(r.circle.center);
        io.outVec(r.circle.normal);
        io.outReal(r.circle.radius);
    };

    family.case('IntrPlane3Circle3.find', (io) => {
        const p = plane3(io);
        const c = circle3(io);
        outPlane3Circle3(io, new IntrPlane3Circle3FI().find(p, c));
    }, { exact: true });

    family.case('IntrPlane3Circle3.find.parallel', (io) => {
        const planeNormal = io.vec(3);
        const planeOrigin = io.vec(3);
        const circleNormal = io.vec(3);
        const circleCenter = io.vec(3);
        const radius = io.real();
        const p = Hyperplane.fromNormalOrigin(planeNormal, planeOrigin);
        const c = Circle3.fromCenterNormalRadius(circleCenter, circleNormal,
            radius);
        outPlane3Circle3(io, new IntrPlane3Circle3FI().find(p, c));
    }, { exact: true });

    family.case('IntrPlane3Circle3.find.tangent', (io) => {
        const planeNormal = io.vec(3);
        const circleCenter = io.vec(3);
        const circleNormal = io.vec(3);
        const radius = io.real();
        const planeOrigin = io.vec(3);
        const p = Hyperplane.fromNormalOrigin(planeNormal, planeOrigin);
        const c = Circle3.fromCenterNormalRadius(circleCenter, circleNormal,
            radius);
        outPlane3Circle3(io, new IntrPlane3Circle3FI().find(p, c));
    }, { exact: true });

    // ========================= IntrPlane3Cylinder3 =======================

    family.case('IntrPlane3Cylinder3.test', (io) => {
        const p = plane3(io);
        const cyl = cylinder3(io);
        const r = new IntrPlane3Cylinder3TI().test(p, cyl);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrPlane3Cylinder3.test.infinite', (io) => {
        const p = plane3(io);
        const cyl = cylinder3(io);
        const r = new IntrPlane3Cylinder3TI().test(p, cyl);
        io.outBool(r.intersect);
    }, { exact: true });

    const outPlane3Cylinder3 = (io: OracleIO,
        r: IntrPlane3Cylinder3FIResult): void => {
        io.outBool(r.intersect);
        io.outInt(r.type);
        for (let i = 0; i < 2; ++i) {
            io.outVec(r.line[i].origin);
            io.outVec(r.line[i].direction);
        }
        io.outVec(r.ellipse.center);
        io.outVec(r.ellipse.normal);
        io.outVec(r.ellipse.axis[0]);
        io.outVec(r.ellipse.axis[1]);
        io.outReal(r.ellipse.extent.values[0]);
        io.outReal(r.ellipse.extent.values[1]);
        for (let i = 0; i < 2; ++i) {
            io.outVec(r.trimLine[i].origin);
            io.outVec(r.trimLine[i].direction);
        }
    };

    family.case('IntrPlane3Cylinder3.find', (io) => {
        const p = plane3(io);
        const cyl = cylinder3(io);
        outPlane3Cylinder3(io, new IntrPlane3Cylinder3FI().find(p, cyl));
    }, { exact: true });

    family.case('IntrPlane3Cylinder3.find.parallel', (io) => {
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const center = io.vec(3);
        const direction = io.vec(3);
        const radius = io.real();
        const height = io.real();
        const p = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        const cyl = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(center, direction), radius, height);
        outPlane3Cylinder3(io, new IntrPlane3Cylinder3FI().find(p, cyl));
    }, { exact: true });

    // ==================== IntrRay2Arc2, IntrSegment2Arc2 =================

    interface ArcResult {
        intersect: boolean;
        numIntersections: number;
        parameter: [number, number];
        point: [Vector, Vector];
    }

    const outArcResult = (io: OracleIO, r: ArcResult): void => {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    };

    family.case('IntrRay2Arc2.test', (io) => {
        const arc = arc2(io);
        const origin = io.vec(2);
        const direction = io.vec(2);
        const r = new IntrRay2Arc2TI().test(
            Ray.fromOriginDirection(origin, direction), arc);
        io.outBool(r.intersect);
    }, { exact: true });

    const ray2Arc2Find = (io: OracleIO): void => {
        const arc = arc2(io);
        const origin = io.vec(2);
        const direction = io.vec(2);
        outArcResult(io, new IntrRay2Arc2FI().find(
            Ray.fromOriginDirection(origin, direction), arc));
    };

    family.case('IntrRay2Arc2.find', ray2Arc2Find, { exact: true });

    // Upstream clips against the solid disk, so a ray origin inside the disk
    // is reported as a point "on the arc" although it is not even on the
    // circle. docs/UPSTREAM-FINDINGS.md IntrRay2Arc2.h, issue #304.
    family.case('IntrRay2Arc2.find.insideDeviation', ray2Arc2Find,
        { deviation: 'UPSTREAM-FINDINGS IntrRay2Arc2.h operator(); #304' });

    family.case('IntrSegment2Arc2.test', (io) => {
        const arc = arc2(io);
        const sg = segment(io, 2);
        const r = new IntrSegment2Arc2TI().test(sg, arc);
        io.outBool(r.intersect);
    }, { exact: true });

    const segment2Arc2Find = (io: OracleIO): void => {
        const arc = arc2(io);
        const sg = segment(io, 2);
        outArcResult(io, new IntrSegment2Arc2FI().find(sg, arc));
    };

    family.case('IntrSegment2Arc2.find', segment2Arc2Find, { exact: true });

    family.case('IntrSegment2Arc2.find.insideDeviation', segment2Arc2Find,
        { deviation: 'UPSTREAM-FINDINGS IntrSegment2Arc2.h operator(); #304' });

    // ============ IntrRay2SegmentMesh2, IntrSegment2SegmentMesh2 =========
    //
    // The delegated line/mesh query orders its output with std::sort on the
    // line parameter alone. std::sort is not stable and Array.prototype.sort
    // is, so both sides re-sort by the total key
    // (parameter, indexPair, meshSegmentParameter, point) before emitting.

    interface MeshHit {
        parameter: number;
        meshSegmentParameter: number;
        indexPair: [number, number];
        point: Vector;
    }

    const meshHitLess = (a: MeshHit, b: MeshHit): number => {
        if (a.parameter !== b.parameter) { return a.parameter < b.parameter ? -1 : 1; }
        if (a.indexPair[0] !== b.indexPair[0]) { return a.indexPair[0] - b.indexPair[0]; }
        if (a.indexPair[1] !== b.indexPair[1]) { return a.indexPair[1] - b.indexPair[1]; }
        if (a.meshSegmentParameter !== b.meshSegmentParameter) {
            return a.meshSegmentParameter < b.meshSegmentParameter ? -1 : 1;
        }
        if (a.point.values[0] !== b.point.values[0]) {
            return a.point.values[0] < b.point.values[0] ? -1 : 1;
        }
        if (a.point.values[1] !== b.point.values[1]) {
            return a.point.values[1] < b.point.values[1] ? -1 : 1;
        }
        return 0;
    };

    const outMeshHits = (io: OracleIO, hits: MeshHit[]): void => {
        hits.sort(meshHitLess);
        io.outInt(hits.length);
        for (const hit of hits) {
            io.outInt(hit.indexPair[0]);
            io.outInt(hit.indexPair[1]);
            io.outReal(hit.parameter);
            io.outReal(hit.meshSegmentParameter);
            io.outVec(hit.point);
        }
    };

    const readMesh2 = (io: OracleIO): SegmentMesh => {
        const numVertices = io.integer();
        const vertices: Vector[] = [];
        for (let i = 0; i < numVertices; ++i) { vertices.push(io.vec(2)); }
        const isOpen = io.boolean();
        return SegmentMesh.fromContiguous(vertices, isOpen);
    };

    const ray2MeshFind = (io: OracleIO): void => {
        const mesh = readMesh2(io);
        const origin = io.vec(2);
        const direction = io.vec(2);
        const r = new IntrRay2SegmentMesh2FI().find(
            Ray.fromOriginDirection(origin, direction), mesh);
        outMeshHits(io, r.intersections.map((object) => ({
            parameter: object.rayParameter,
            meshSegmentParameter: object.meshSegmentParameter,
            indexPair: object.indexPair,
            point: object.point
        })));
    };

    family.case('IntrRay2SegmentMesh2.find', ray2MeshFind, { exact: true });
    family.case('IntrRay2SegmentMesh2.find.throughVertex', ray2MeshFind,
        { exact: true });

    const segment2MeshFind = (io: OracleIO): void => {
        const mesh = readMesh2(io);
        const sg = segment(io, 2);
        const r = new IntrSegment2SegmentMesh2FI().find(sg, mesh);
        outMeshHits(io, r.intersections.map((object) => ({
            parameter: object.segmentParameter,
            meshSegmentParameter: object.meshSegmentParameter,
            indexPair: object.indexPair,
            point: object.point
        })));
    };

    family.case('IntrSegment2SegmentMesh2.find', segment2MeshFind,
        { exact: true });
    family.case('IntrSegment2SegmentMesh2.find.throughVertex', segment2MeshFind,
        { exact: true });

    // ============ IntrRay3Ellipsoid3, IntrSegment3Ellipsoid3 =============

    const outEllipsoidResult = (io: OracleIO,
        r: IntrLine3Ellipsoid3FIResult, withPoints: boolean): void => {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (withPoints) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    };

    family.case('IntrRay3Ellipsoid3.test', (io) => {
        const e = ellipsoid(io);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = new IntrRay3Ellipsoid3TI().test(
            Ray.fromOriginDirection(origin, direction), e);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay3Ellipsoid3.find', (io) => {
        const e = ellipsoid(io);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = new IntrRay3Ellipsoid3FI().find(
            Ray.fromOriginDirection(origin, direction), e);
        outEllipsoidResult(io, r, true);
    }, { exact: true });

    family.case('IntrRay3Ellipsoid3.doQuery.fi', (io) => {
        const e = ellipsoid(io);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = defaultIntrRay3Ellipsoid3FIResult();
        intrRay3Ellipsoid3FIDoQuery(origin, direction, e, r);
        outEllipsoidResult(io, r, false);
    }, { exact: true });

    family.case('IntrRay3Ellipsoid3.find.tangent', (io) => {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3), io.vec(3)];
        const extent = io.vec(3);
        const e = Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = new IntrRay3Ellipsoid3FI().find(
            Ray.fromOriginDirection(origin, direction), e);
        outEllipsoidResult(io, r, true);
    }, { exact: true });

    const segment3Ellipsoid3Test = (io: OracleIO): void => {
        const e = ellipsoid(io);
        const sg = segment(io, 3);
        const r = new IntrSegment3Ellipsoid3TI().test(sg, e);
        io.outBool(r.intersect);
    };

    family.case('IntrSegment3Ellipsoid3.test', segment3Ellipsoid3Test,
        { exact: true });

    // A segment strictly inside the solid ellipsoid: upstream reports no
    // intersection, contradicting its own FI query.
    // docs/UPSTREAM-FINDINGS.md IntrSegment3Ellipsoid3.h TIQuery, issue #304.
    family.case('IntrSegment3Ellipsoid3.test.containedDeviation',
        segment3Ellipsoid3Test,
        { deviation: 'UPSTREAM-FINDINGS IntrSegment3Ellipsoid3.h TIQuery; #304' });

    family.case('IntrSegment3Ellipsoid3.find', (io) => {
        const e = ellipsoid(io);
        const sg = segment(io, 3);
        const r = new IntrSegment3Ellipsoid3FI().find(sg, e);
        outEllipsoidResult(io, r, true);
    }, { exact: true });

    family.case('IntrSegment3Ellipsoid3.doQuery.fi', (io) => {
        const e = ellipsoid(io);
        const segOrigin = io.vec(3);
        const segDirection = io.vec(3);
        const segExtent = io.real();
        const r = defaultIntrSegment3Ellipsoid3FIResult();
        intrSegment3Ellipsoid3FIDoQuery(segOrigin, segDirection, segExtent, e, r);
        outEllipsoidResult(io, r, false);
    }, { exact: true });

    family.case('IntrSegment3Ellipsoid3.find.tangent', (io) => {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3), io.vec(3)];
        const extent = io.vec(3);
        const e = Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
        const sg = segment(io, 3);
        const r = new IntrSegment3Ellipsoid3FI().find(sg, e);
        outEllipsoidResult(io, r, true);
    }, { exact: true });

    // ============================= IntrLine3Cone3 ========================

    const outQFN = (io: OracleIO, q: QFNumber): void => {
        io.outReal(q.x[0] as number);
        io.outReal(q.x[1] as number);
        io.outReal(q.d);
    };

    const outLine3Cone3 = (io: OracleIO, r: IntrLine3Cone3FIResult,
        withPoints: boolean): void => {
        io.outBool(r.intersect);
        io.outInt(r.type);
        outQFN(io, r.t[0]);
        outQFN(io, r.t[1]);
        if (withPoints) {
            for (let i = 0; i < 2; ++i) {
                for (let k = 0; k < 3; ++k) { outQFN(io, r.P[i][k]); }
            }
        }
    };

    const line3Cone3Find = (io: OracleIO): void => {
        const c = cone(io, 3);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = new IntrLine3Cone3FI().find(
            Line.fromOriginDirection(origin, direction), c);
        outLine3Cone3(io, r, true);
    };

    family.case('IntrLine3Cone3.find', line3Cone3Find, { exact: true });

    family.case('IntrLine3Cone3.doQuery.fi', (io) => {
        const c = cone(io, 3);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = defaultIntrLine3Cone3FIResult();
        intrLine3Cone3FIDoQuery(origin, direction, c, r);
        outLine3Cone3(io, r, false);
    }, { exact: true });

    // A line through the cone vertex makes the quadratic have a double root,
    // so upstream's case analysis is decided by a discriminant with no
    // significant digits. docs/UPSTREAM-FINDINGS.md IntrLine3Cone3.h,
    // issues #304 and #465.
    family.case('IntrLine3Cone3.find.vertexDeviation', line3Cone3Find,
        { deviation: 'UPSTREAM-FINDINGS IntrLine3Cone3.h through-vertex; #304, #465' });

    // ========================== IntrAlignedBox3Cone3 =====================

    const alignedBox3 = (io: OracleIO): AlignedBox => {
        const min = io.vec(3);
        const max = io.vec(3);
        return AlignedBox.fromMinMax(min, max);
    };

    family.case('IntrAlignedBox3Cone3.test', (io) => {
        const c = cone(io, 3);
        const box = alignedBox3(io);
        const r = new IntrAlignedBox3Cone3TI().test(box, c);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrAlignedBox3Cone3.test.straddle', (io) => {
        const c = cone(io, 3);
        const box = alignedBox3(io);
        const r = new IntrAlignedBox3Cone3TI().test(box, c);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrAlignedBox3Cone3.computeBoxHeightInterval', (io) => {
        const box = alignedBox3(io);
        const c = cone(io, 3);
        const r = intrAlignedBox3Cone3ComputeBoxHeightInterval(box, c);
        io.outReal(r.boxMinHeight);
        io.outReal(r.boxMaxHeight);
    }, { exact: true });

    family.case('IntrAlignedBox3Cone3.coneAxisIntersectsBox', (io) => {
        const c = cone(io, 3);
        const box = alignedBox3(io);
        io.outBool(intrAlignedBox3Cone3ConeAxisIntersectsBox(box, c));
    }, { exact: true });

    family.case('IntrAlignedBox3Cone3.hasPointInsideCone', (io) => {
        const c = cone(io, 3);
        const P0 = io.vec(3);
        const P1 = io.vec(3);
        io.outBool(intrAlignedBox3Cone3HasPointInsideCone(P0, P1, c));
    }, { exact: true });

    // Upstream's BoxFullyInConeSlab never clears the adjacency matrix, so a
    // later clipping query on the same object silently drops candidate edges.
    // The three queries run on ONE query object.
    // docs/UPSTREAM-FINDINGS.md IntrAlignedBox3Cone3.h, issue #301.
    family.case('IntrAlignedBox3Cone3.test.staleAdjacencyDeviation', (io) => {
        const c = cone(io, 3);
        const boxA = alignedBox3(io);
        const boxB = alignedBox3(io);
        const query = new IntrAlignedBox3Cone3TI();
        io.outBool(query.test(boxA, c).intersect);
        io.outBool(query.test(boxB, c).intersect);
        io.outBool(query.test(boxA, c).intersect);
    }, { deviation: 'UPSTREAM-FINDINGS IntrAlignedBox3Cone3.h BoxFullyInConeSlab; #301' });

    // ======================== IntrAreaEllipse2Ellipse2 ===================
    //
    // The compared computation calls atan2, atan, sin and cos, and the
    // underlying ellipse/ellipse find-intersection query calls
    // RootsPolynomial::SolveQuartic (pow, cos, atan2), so the area is
    // compared with a tolerance. The configuration and the intersection point
    // count are discrete and compared exactly.

    const areaEllipsesCompute = (io: OracleIO): void => {
        const e0 = ellipse2(io);
        const e1 = ellipse2(io);
        const r = new AreaEllipse2Ellipse2().compute(e0, e1);
        io.outInt(r.configuration);
        io.outInt(r.findResult.numPoints === Number.MAX_SAFE_INTEGER
            ? -1 : r.findResult.numPoints);
        io.outReal(r.area);
    };

    family.case('IntrAreaEllipse2Ellipse2.compute', areaEllipsesCompute,
        { tol: 1e-8 });

    // Upstream never assigns mZero, mOne, mTwo, mPi and mTwoPi, so every area
    // it computes is garbage. docs/UPSTREAM-FINDINGS.md
    // IntrAreaEllipse2Ellipse2.h, issue #301.
    family.case('IntrAreaEllipse2Ellipse2.compute.uninitializedDeviation',
        areaEllipsesCompute,
        { deviation: 'UPSTREAM-FINDINGS IntrAreaEllipse2Ellipse2.h uninitialized members; #301' });

    // ======================== IntrConvexMesh3Plane3 ======================

    const Q = IntrConvexMesh3Plane3FI;

    const outMesh3 = (io: OracleIO, mesh: ConvexMesh3): void => {
        io.outInt(mesh.configuration);
        io.outInt(mesh.vertices.length);
        for (const vertex of mesh.vertices) { io.outVec(vertex); }
        io.outInt(mesh.triangles.length);
        for (const triangle of mesh.triangles) {
            io.outInt(triangle[0]);
            io.outInt(triangle[1]);
            io.outInt(triangle[2]);
        }
    };

    const outMesh3Result = (io: OracleIO,
        r: IntrConvexMesh3Plane3FIResult, withPolyhedra: boolean): void => {
        io.outInt(r.configuration);
        io.outInt(r.requested);
        outMesh3(io, r.intersectionMesh);
        io.outInt(r.intersectionPolygon.length);
        for (const vertex of r.intersectionPolygon) { io.outVec(vertex); }
        if (withPolyhedra) {
            outMesh3(io, r.positivePolyhedron);
            outMesh3(io, r.negativePolyhedron);
        }
    };

    const tetraMesh = (io: OracleIO): ConvexMesh3 => {
        const mesh = new ConvexMesh3();
        mesh.configuration = ConvexMesh3.CFG_POLYHEDRON;
        for (let i = 0; i < 4; ++i) { mesh.vertices.push(io.vec(3)); }
        mesh.triangles = ([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]] as
            ConvexMesh3Triangle[]);
        return mesh;
    };

    const boxMesh = (io: OracleIO): ConvexMesh3 => {
        const mesh = new ConvexMesh3();
        mesh.configuration = ConvexMesh3.CFG_POLYHEDRON;
        for (let i = 0; i < 8; ++i) { mesh.vertices.push(io.vec(3)); }
        mesh.triangles = ([
            [0, 2, 3], [0, 3, 1],
            [4, 5, 7], [4, 7, 6],
            [0, 1, 5], [0, 5, 4],
            [2, 6, 7], [2, 7, 3],
            [0, 4, 6], [0, 6, 2],
            [1, 3, 7], [1, 7, 5]
        ] as ConvexMesh3Triangle[]);
        return mesh;
    };

    family.case('IntrConvexMesh3Plane3.find.tetrahedron', (io) => {
        const mesh = tetraMesh(io);
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const plane = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        const r = new Q().find(mesh, plane, Q.REQ_ALL);
        outMesh3Result(io, r, true);
    }, { exact: true });

    family.case('IntrConvexMesh3Plane3.find.configurationOnly', (io) => {
        const mesh = tetraMesh(io);
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const plane = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        const r = new Q().find(mesh, plane, Q.REQ_CONFIGURATION_ONLY);
        outMesh3Result(io, r, true);
    }, { exact: true });

    family.case('IntrConvexMesh3Plane3.find.tangent', (io) => {
        const mesh = tetraMesh(io);
        const normal = io.vec(3);
        const through = io.vec(3);
        const plane = Hyperplane.fromNormalOrigin(normal, through);
        const r = new Q().find(mesh, plane, Q.REQ_ALL);
        outMesh3Result(io, r, true);
    }, { exact: true });

    family.case('IntrConvexMesh3Plane3.find.box', (io) => {
        const mesh = boxMesh(io);
        const normal = io.vec(3);
        const planeOrigin = io.vec(3);
        const plane = Hyperplane.fromNormalOrigin(normal, planeOrigin);
        const r = new Q().find(mesh, plane, Q.REQ_INTR_BOTH);
        outMesh3Result(io, r, false);
    }, { exact: true });

    // Upstream reads the predecessor map of the coplanar face boundary in
    // index order instead of traversing it as a cycle.
    // docs/UPSTREAM-FINDINGS.md IntrConvexMesh3Plane3.h, issue #301.
    family.case('IntrConvexMesh3Plane3.find.coplanarFaceDeviation', (io) => {
        const mesh = boxMesh(io);
        const normal = io.vec(3);
        const through = io.vec(3);
        const plane = Hyperplane.fromNormalOrigin(normal, through);
        const r = new Q().find(mesh, plane, Q.REQ_INTR_BOTH);
        outMesh3Result(io, r, false);
    }, { deviation: 'UPSTREAM-FINDINGS IntrConvexMesh3Plane3.h GetIntersectionPolygon; #301' });

    family.finish();
});
