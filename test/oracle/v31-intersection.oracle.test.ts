// Replays oracle/cpp/cases/v31-intersection.cpp (verify group 31,
// intersection). Keep the two files in the same order.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { Arc2 } from '../../src/Arc2.js';
import { CanonicalBox } from '../../src/CanonicalBox.js';
import { Capsule } from '../../src/Capsule.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { Frustum3 } from '../../src/Frustum3.js';
import { Halfspace } from '../../src/Halfspace.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Ray } from '../../src/Ray.js';
import { Rectangle } from '../../src/Rectangle.js';
import { Sector2 } from '../../src/Sector2.js';
import { Segment } from '../../src/Segment.js';
import { Torus3 } from '../../src/Torus3.js';
import { Triangle } from '../../src/Triangle.js';
import { Vector, add, mul, negate } from '../../src/Vector.js';
import { IntrArc2Arc2FI } from '../../src/IntrArc2Arc2.js';
import { IntrCanonicalBox3Cylinder3TI } from '../../src/IntrCanonicalBox3Cylinder3.js';
import { IntrCapsule3Capsule3TI } from '../../src/IntrCapsule3Capsule3.js';
import { IntrCircle2Arc2FI } from '../../src/IntrCircle2Arc2.js';
import { IntrCylinder3Cylinder3TI } from '../../src/IntrCylinder3Cylinder3.js';
import { IntrDisk2Sector2TI } from '../../src/IntrDisk2Sector2.js';
import { IntrHalfspace3Capsule3TI } from '../../src/IntrHalfspace3Capsule3.js';
import { IntrHalfspace3Cylinder3TI } from '../../src/IntrHalfspace3Cylinder3.js';
import {
    IntrLine2Circle2FI, IntrLine2Circle2TI,
    defaultIntrLine2Circle2FIResult, intrLine2Circle2FIDoQuery
} from '../../src/IntrLine2Circle2.js';
import {
    IntrLine2OrientedBox2FI, IntrLine2OrientedBox2TI
} from '../../src/IntrLine2OrientedBox2.js';
import { IntrLine2Ray2FI, IntrLine2Ray2TI } from '../../src/IntrLine2Ray2.js';
import { IntrLine2Segment2FI, IntrLine2Segment2TI } from '../../src/IntrLine2Segment2.js';
import {
    IntrLine3Capsule3FI, IntrLine3Capsule3TI,
    defaultIntrLine3Capsule3FIResult, intrLine3Capsule3FIDoQuery
} from '../../src/IntrLine3Capsule3.js';
import {
    IntrLine3Cylinder3FI,
    defaultIntrLine3Cylinder3FIResult, intrLine3Cylinder3FIDoQuery
} from '../../src/IntrLine3Cylinder3.js';
import {
    IntrLine3OrientedBox3FI, IntrLine3OrientedBox3TI
} from '../../src/IntrLine3OrientedBox3.js';
import { IntrLine3Torus3FI } from '../../src/IntrLine3Torus3.js';
import { IntrOrientedBox2Sector2TI } from '../../src/IntrOrientedBox2Sector2.js';
import { IntrOrientedBox3Frustum3TI } from '../../src/IntrOrientedBox3Frustum3.js';
import {
    IntrRay2AlignedBox2FI, IntrRay2AlignedBox2TI,
    intrRay2AlignedBox2FIDoQuery, intrRay2AlignedBox2TIDoQuery
} from '../../src/IntrRay2AlignedBox2.js';
import {
    defaultIntrLine2AlignedBox2FIResult, defaultIntrLine2AlignedBox2TIResult
} from '../../src/IntrLine2AlignedBox2.js';
import {
    defaultIntrLine3AlignedBox3FIResult, defaultIntrLine3AlignedBox3TIResult
} from '../../src/IntrLine3AlignedBox3.js';
import { defaultIntrLine2Triangle2FIResult } from '../../src/IntrLine2Triangle2.js';
import { IntrRay2Ray2FI, IntrRay2Ray2TI } from '../../src/IntrRay2Ray2.js';
import { IntrRay2Segment2FI, IntrRay2Segment2TI } from '../../src/IntrRay2Segment2.js';
import {
    IntrRay2Triangle2FI, IntrRay2Triangle2TI, intrRay2Triangle2FIDoQuery
} from '../../src/IntrRay2Triangle2.js';
import {
    IntrRay3AlignedBox3FI, IntrRay3AlignedBox3TI,
    intrRay3AlignedBox3FIDoQuery, intrRay3AlignedBox3TIDoQuery
} from '../../src/IntrRay3AlignedBox3.js';
import { IntrRay3Rectangle3FI, IntrRay3Rectangle3TI } from '../../src/IntrRay3Rectangle3.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- input readers, mirroring the C++ generators ----

function box(io: OracleIO, n: number): AlignedBox {
    const min = io.vec(n);
    const max = io.vec(n);
    return AlignedBox.fromMinMax(min, max);
}

function obox(io: OracleIO, n: number): OrientedBox {
    const center = io.vec(n);
    const axis: Vector[] = [];
    for (let i = 0; i < n; ++i) { axis.push(io.vec(n)); }
    const extent = io.vec(n);
    return OrientedBox.fromCenterAxisExtent(center, axis, extent);
}

// The C++ generator records a full orthonormal frame; the rectangle uses only
// the first two axes.
function rectangle3(io: OracleIO): Rectangle {
    const center = io.vec(3);
    const axis0 = io.vec(3);
    const axis1 = io.vec(3);
    io.vec(3);
    const extent = io.vec(2);
    return Rectangle.fromCenterAxisExtent(center, [axis0, axis1], extent);
}

function triangle(io: OracleIO, n: number): Triangle {
    const v0 = io.vec(n);
    const v1 = io.vec(n);
    const v2 = io.vec(n);
    return Triangle.fromVertices(v0, v1, v2);
}

function line(io: OracleIO, n: number): Line {
    const origin = io.vec(n);
    const direction = io.vec(n);
    return Line.fromOriginDirection(origin, direction);
}

function ray(io: OracleIO, n: number): Ray {
    const origin = io.vec(n);
    const direction = io.vec(n);
    return Ray.fromOriginDirection(origin, direction);
}

function segment(io: OracleIO, n: number): Segment {
    const p0 = io.vec(n);
    const p1 = io.vec(n);
    return Segment.fromEndpoints(p0, p1);
}

function circle(io: OracleIO): Hypersphere {
    const center = io.vec(2);
    const radius = io.real();
    return Hypersphere.fromCenterRadius(center, radius);
}

function capsule(io: OracleIO): Capsule {
    const p0 = io.vec(3);
    const p1 = io.vec(3);
    const radius = io.real();
    return Capsule.fromSegmentRadius(Segment.fromEndpoints(p0, p1), radius);
}

function cylinder(io: OracleIO): Cylinder3 {
    const origin = io.vec(3);
    const direction = io.vec(3);
    const radius = io.real();
    const height = io.real();
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(origin, direction), radius, height);
}

// The C++ generator records cosAngle and sinAngle directly, so no libm call
// enters the compared computation. 'angle' is unused by both queries and keeps
// its default value on either side.
function sector(io: OracleIO): Sector2 {
    const s = new Sector2();
    s.vertex = io.vec(2);
    s.radius = io.real();
    s.direction = io.vec(2);
    s.cosAngle = io.real();
    s.sinAngle = io.real();
    return s;
}

function arc(io: OracleIO): Arc2 {
    const center = io.vec(2);
    const radius = io.real();
    const end0 = io.vec(2);
    const end1 = io.vec(2);
    return Arc2.fromCenterRadiusEnds(center, radius, end0, end1);
}

function outArc(io: OracleIO, a: Arc2): void {
    io.outVec(a.center);
    io.outReal(a.radius);
    io.outVec(a.end[0]);
    io.outVec(a.end[1]);
}

describe('oracle: v31-intersection', () => {
    const family = new OracleFamily('v31-intersection');

    // ============================ IntrLine2Circle2 =======================

    family.case('IntrLine2Circle2.test', (io) => {
        const ln = line(io, 2);
        const c = circle(io);
        const r = new IntrLine2Circle2TI().test(ln, c);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine2Circle2.find', (io) => {
        const ln = line(io, 2);
        const c = circle(io);
        const r = new IntrLine2Circle2FI().find(ln, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    family.case('IntrLine2Circle2.doQuery.fi', (io) => {
        const lineOrigin = io.vec(2);
        const lineDirection = io.vec(2);
        const c = circle(io);
        const r = defaultIntrLine2Circle2FIResult();
        intrLine2Circle2FIDoQuery(lineOrigin, lineDirection, c, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrLine2Circle2.find.tangent', (io) => {
        const c = circle(io);
        const origin = io.vec(2);
        const direction = io.vec(2);
        const r = new IntrLine2Circle2FI().find(
            Line.fromOriginDirection(origin, direction), c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    // ========================= IntrLine2OrientedBox2 =====================

    family.case('IntrLine2OrientedBox2.test', (io) => {
        const ln = line(io, 2);
        const b = obox(io, 2);
        const r = new IntrLine2OrientedBox2TI().test(ln, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine2OrientedBox2.find', (io) => {
        const ln = line(io, 2);
        const b = obox(io, 2);
        const r = new IntrLine2OrientedBox2FI().find(ln, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    // ========================= IntrLine3OrientedBox3 =====================

    family.case('IntrLine3OrientedBox3.test', (io) => {
        const ln = line(io, 3);
        const b = obox(io, 3);
        const r = new IntrLine3OrientedBox3TI().test(ln, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3OrientedBox3.find', (io) => {
        const ln = line(io, 3);
        const b = obox(io, 3);
        const r = new IntrLine3OrientedBox3FI().find(ln, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    // ============================= IntrLine2Ray2 =========================

    family.case('IntrLine2Ray2.test', (io) => {
        const ln = line(io, 2);
        const ry = ray(io, 2);
        const r = new IntrLine2Ray2TI().test(ln, ry);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    family.case('IntrLine2Ray2.find', (io) => {
        const ln = line(io, 2);
        const ry = ray(io, 2);
        const r = new IntrLine2Ray2FI().find(ln, ry);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.lineParameter[0]);
        io.outReal(r.lineParameter[1]);
        io.outReal(r.rayParameter[0]);
        io.outReal(r.rayParameter[1]);
        io.outVec(r.point);
    }, { exact: true });

    // =========================== IntrLine2Segment2 =======================

    family.case('IntrLine2Segment2.test', (io) => {
        const ln = line(io, 2);
        const sg = segment(io, 2);
        const r = new IntrLine2Segment2TI().test(ln, sg);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    family.case('IntrLine2Segment2.find', (io) => {
        const ln = line(io, 2);
        const sg = segment(io, 2);
        const r = new IntrLine2Segment2FI().find(ln, sg);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.lineParameter[0]);
        io.outReal(r.lineParameter[1]);
        io.outReal(r.segmentParameter[0]);
        io.outReal(r.segmentParameter[1]);
        io.outVec(r.point);
    }, { exact: true });

    // ============================== IntrRay2Ray2 =========================

    family.case('IntrRay2Ray2.test', (io) => {
        const ray0 = ray(io, 2);
        const ray1 = ray(io, 2);
        const r = new IntrRay2Ray2TI().test(ray0, ray1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    const ray2Ray2Find = (io: OracleIO, ray0: Ray, ray1: Ray): void => {
        const r = new IntrRay2Ray2FI().find(ray0, ray1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.ray0Parameter[0]);
        io.outReal(r.ray0Parameter[1]);
        io.outReal(r.ray1Parameter[0]);
        io.outReal(r.ray1Parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    };

    family.case('IntrRay2Ray2.find', (io) => {
        const ray0 = ray(io, 2);
        const ray1 = ray(io, 2);
        ray2Ray2Find(io, ray0, ray1);
    }, { exact: true });

    family.case('IntrRay2Ray2.find.collinear', (io) => {
        const base = io.vec(2);
        const direction = io.vec(2);
        const t0 = io.real();
        const t1 = io.real();
        const flip = io.boolean();
        const origin0 = add(base, mul(t0, direction));
        const origin1 = add(base, mul(t1, direction));
        const direction1 = flip ? negate(direction) : direction;
        ray2Ray2Find(io, Ray.fromOriginDirection(origin0, direction),
            Ray.fromOriginDirection(origin1, direction1));
    }, { exact: true });

    family.case('IntrRay2Ray2.test.collinear', (io) => {
        const base = io.vec(2);
        const direction = io.vec(2);
        const t0 = io.real();
        const t1 = io.real();
        const flip = io.boolean();
        const origin0 = add(base, mul(t0, direction));
        const origin1 = add(base, mul(t1, direction));
        const direction1 = flip ? negate(direction) : direction;
        const r = new IntrRay2Ray2TI().test(
            Ray.fromOriginDirection(origin0, direction),
            Ray.fromOriginDirection(origin1, direction1));
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    // ============================ IntrRay2Segment2 =======================

    family.case('IntrRay2Segment2.test', (io) => {
        const ry = ray(io, 2);
        const sg = segment(io, 2);
        const r = new IntrRay2Segment2TI().test(ry, sg);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    const ray2Segment2Find = (io: OracleIO, ry: Ray, sg: Segment): void => {
        const r = new IntrRay2Segment2FI().find(ry, sg);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.rayParameter[0]);
        io.outReal(r.rayParameter[1]);
        io.outReal(r.segmentParameter[0]);
        io.outReal(r.segmentParameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    };

    family.case('IntrRay2Segment2.find', (io) => {
        const ry = ray(io, 2);
        const sg = segment(io, 2);
        ray2Segment2Find(io, ry, sg);
    }, { exact: true });

    family.case('IntrRay2Segment2.find.collinear', (io) => {
        const base = io.vec(2);
        const direction = io.vec(2);
        const t0 = io.real();
        const t1 = io.real();
        const t2 = io.real();
        const origin = add(base, mul(t0, direction));
        const p0 = add(base, mul(t1, direction));
        const p1 = add(base, mul(t2, direction));
        ray2Segment2Find(io, Ray.fromOriginDirection(origin, direction),
            Segment.fromEndpoints(p0, p1));
    }, { exact: true });

    family.case('IntrRay2Segment2.test.collinear', (io) => {
        const base = io.vec(2);
        const direction = io.vec(2);
        const t0 = io.real();
        const t1 = io.real();
        const t2 = io.real();
        const origin = add(base, mul(t0, direction));
        const p0 = add(base, mul(t1, direction));
        const p1 = add(base, mul(t2, direction));
        const r = new IntrRay2Segment2TI().test(
            Ray.fromOriginDirection(origin, direction),
            Segment.fromEndpoints(p0, p1));
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    // =========================== IntrRay2AlignedBox2 =====================

    family.case('IntrRay2AlignedBox2.test', (io) => {
        const ry = ray(io, 2);
        const b = box(io, 2);
        const r = new IntrRay2AlignedBox2TI().test(ry, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay2AlignedBox2.find', (io) => {
        const ry = ray(io, 2);
        const b = box(io, 2);
        const r = new IntrRay2AlignedBox2FI().find(ry, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    family.case('IntrRay2AlignedBox2.doQuery.ti', (io) => {
        const rayOrigin = io.vec(2);
        const rayDirection = io.vec(2);
        const boxExtent = io.vec(2);
        const r = defaultIntrLine2AlignedBox2TIResult();
        intrRay2AlignedBox2TIDoQuery(rayOrigin, rayDirection, boxExtent, r);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay2AlignedBox2.doQuery.fi', (io) => {
        const rayOrigin = io.vec(2);
        const rayDirection = io.vec(2);
        const boxExtent = io.vec(2);
        const r = defaultIntrLine2AlignedBox2FIResult();
        intrRay2AlignedBox2FIDoQuery(rayOrigin, rayDirection, boxExtent, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    // =========================== IntrRay3AlignedBox3 =====================

    family.case('IntrRay3AlignedBox3.test', (io) => {
        const ry = ray(io, 3);
        const b = box(io, 3);
        const r = new IntrRay3AlignedBox3TI().test(ry, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay3AlignedBox3.find', (io) => {
        const ry = ray(io, 3);
        const b = box(io, 3);
        const r = new IntrRay3AlignedBox3FI().find(ry, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrRay3AlignedBox3.doQuery.ti', (io) => {
        const rayOrigin = io.vec(3);
        const rayDirection = io.vec(3);
        const boxExtent = io.vec(3);
        const r = defaultIntrLine3AlignedBox3TIResult();
        intrRay3AlignedBox3TIDoQuery(rayOrigin, rayDirection, boxExtent, r);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay3AlignedBox3.doQuery.fi', (io) => {
        const rayOrigin = io.vec(3);
        const rayDirection = io.vec(3);
        const boxExtent = io.vec(3);
        const r = defaultIntrLine3AlignedBox3FIResult();
        intrRay3AlignedBox3FIDoQuery(rayOrigin, rayDirection, boxExtent, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    // ============================ IntrRay2Triangle2 ======================

    family.case('IntrRay2Triangle2.test', (io) => {
        const ry = ray(io, 2);
        const t = triangle(io, 2);
        const r = new IntrRay2Triangle2TI().test(ry, t);
        io.outBool(r.intersect);
    }, { exact: true });

    const ray2Triangle2Find = (io: OracleIO, ry: Ray, t: Triangle): void => {
        const r = new IntrRay2Triangle2FI().find(ry, t);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    };

    family.case('IntrRay2Triangle2.find', (io) => {
        const ry = ray(io, 2);
        const t = triangle(io, 2);
        ray2Triangle2Find(io, ry, t);
    }, { exact: true });

    family.case('IntrRay2Triangle2.doQuery.fi', (io) => {
        const origin = io.vec(2);
        const direction = io.vec(2);
        const t = triangle(io, 2);
        const r = defaultIntrLine2Triangle2FIResult();
        intrRay2Triangle2FIDoQuery(origin, direction, t, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrRay2Triangle2.find.throughPoint', (io) => {
        const t = triangle(io, 2);
        const direction = io.vec(2);
        const origin = io.vec(2);
        ray2Triangle2Find(io, Ray.fromOriginDirection(origin, direction), t);
    }, { exact: true });

    // =========================== IntrRay3Rectangle3 ======================

    const ray3Rectangle3Find = (io: OracleIO, ry: Ray, rect: Rectangle): void => {
        const r = new IntrRay3Rectangle3FI().find(ry, rect);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.rectCoord[0]);
            io.outReal(r.rectCoord[1]);
            io.outReal(r.rectCoord[2]);
            io.outVec(r.point);
        }
    };

    family.case('IntrRay3Rectangle3.test', (io) => {
        const ry = ray(io, 3);
        const rect = rectangle3(io);
        const r = new IntrRay3Rectangle3TI().test(ry, rect);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay3Rectangle3.find', (io) => {
        const ry = ray(io, 3);
        const rect = rectangle3(io);
        ray3Rectangle3Find(io, ry, rect);
    }, { exact: true });

    family.case('IntrRay3Rectangle3.find.throughPoint', (io) => {
        const rect = rectangle3(io);
        const direction = io.vec(3);
        const origin = io.vec(3);
        ray3Rectangle3Find(io, Ray.fromOriginDirection(origin, direction), rect);
    }, { exact: true });

    // ================================ arcs ===============================

    family.case('IntrArc2Arc2.find', (io) => {
        // The lattice centre shared by the two cocircular arcs of mode 0.
        io.vec(2);
        const arc0 = arc(io);
        const arc1 = arc(io);
        const r = new IntrArc2Arc2FI().find(arc0, arc1);
        io.outBool(r.intersect);
        io.outInt(r.configuration);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
        outArc(io, r.arc[0]);
        outArc(io, r.arc[1]);
    }, { exact: true });

    family.case('IntrCircle2Arc2.find', (io) => {
        io.vec(2);
        const a = arc(io);
        const c = circle(io);
        const r = new IntrCircle2Arc2FI().find(c, a);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
        outArc(io, r.arc);
    }, { exact: true });

    // ============================ IntrDisk2Sector2 =======================

    family.case('IntrDisk2Sector2.test', (io) => {
        const disk = circle(io);
        const s = sector(io);
        const r = new IntrDisk2Sector2TI().test(disk, s);
        io.outBool(r.intersect);
    }, { exact: true });

    // ========================= IntrOrientedBox2Sector2 ===================

    const orientedBox2Sector2 = (io: OracleIO): void => {
        const b = obox(io, 2);
        const s = sector(io);
        const r = new IntrOrientedBox2Sector2TI().test(b, s);
        io.outBool(r.intersect);
    };

    family.case('IntrOrientedBox2Sector2.test', orientedBox2Sector2,
        { exact: true });

    // The port keeps the working polygon when IntrHalfspace2Polygon2 reports
    // "no clipping necessary" with an empty polygon; upstream discards it.
    family.case('IntrOrientedBox2Sector2.test.clipDeviation',
        orientedBox2Sector2,
        { deviation: 'docs/UPSTREAM-FINDINGS.md IntrOrientedBox2Sector2.h '
            + 'boundary clipping; issue #200' });

    // =========================== IntrCapsule3Capsule3 ====================

    family.case('IntrCapsule3Capsule3.test', (io) => {
        const capsule0 = capsule(io);
        const capsule1 = capsule(io);
        const r = new IntrCapsule3Capsule3TI().test(capsule0, capsule1);
        io.outBool(r.intersect);
    }, { exact: true });

    // ========================== IntrHalfspace3Capsule3 ===================

    family.case('IntrHalfspace3Capsule3.test', (io) => {
        const c = capsule(io);
        const normal = io.vec(3);
        const constant = io.real();
        const r = new IntrHalfspace3Capsule3TI().test(
            Halfspace.fromNormalConstant(normal, constant), c);
        io.outBool(r.intersect);
    }, { exact: true });

    // ========================= IntrHalfspace3Cylinder3 ===================

    const halfspace3Cylinder3 = (io: OracleIO): void => {
        const origin = io.vec(3);
        const axisDirection = io.vec(3);
        const radius = io.real();
        const height = io.real();
        const normal = io.vec(3);
        const constant = io.real();
        const cyl = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin, axisDirection), radius, height);
        const r = new IntrHalfspace3Cylinder3TI().test(
            Halfspace.fromNormalConstant(normal, constant), cyl);
        io.outBool(r.intersect);
    };

    family.case('IntrHalfspace3Cylinder3.test.perpendicular',
        halfspace3Cylinder3, { exact: true });

    // The port computes sqrt(max(0, 1 - Dot(N,W)^2)); upstream's
    // sqrt(max(1, ...)) is always 1.
    family.case('IntrHalfspace3Cylinder3.test.rootDeviation',
        halfspace3Cylinder3,
        { deviation: 'docs/UPSTREAM-FINDINGS.md IntrHalfspace3Cylinder3.h '
            + 'root computation; issue #197' });

    // ============================ IntrLine3Capsule3 ======================

    family.case('IntrLine3Capsule3.test', (io) => {
        const ln = line(io, 3);
        const c = capsule(io);
        const r = new IntrLine3Capsule3TI().test(ln, c);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3Capsule3.find', (io) => {
        const ln = line(io, 3);
        const c = capsule(io);
        const r = new IntrLine3Capsule3FI().find(ln, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    const line3Capsule3DoQuery = (io: OracleIO): void => {
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const c = capsule(io);
        const r = defaultIntrLine3Capsule3FIResult();
        intrLine3Capsule3FIDoQuery(lineOrigin, lineDirection, c, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    };

    family.case('IntrLine3Capsule3.doQuery.fi', line3Capsule3DoQuery,
        { exact: true });

    // The two deliberate fixes of gtengine-js issue #461: the capsule-cap
    // 'intersect' flag and the cap-junction double count. The C++ generator
    // records the capsule before the line, so the reader order differs from
    // the doQuery case above.
    const line3Capsule3Deviation = (io: OracleIO): void => {
        io.vec(3);
        const c = capsule(io);
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const r = defaultIntrLine3Capsule3FIResult();
        intrLine3Capsule3FIDoQuery(lineOrigin, lineDirection, c, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    };

    family.case('IntrLine3Capsule3.doQuery.capTangentDeviation',
        line3Capsule3Deviation,
        { deviation: 'docs/UPSTREAM-FINDINGS.md IntrLine3Capsule3.h hemisphere '
            + 'roots; issue #461 item 2' });

    family.case('IntrLine3Capsule3.doQuery.junctionDeviation',
        line3Capsule3Deviation,
        { deviation: 'docs/UPSTREAM-FINDINGS.md IntrLine3Capsule3.h '
            + 'cap-junction plane; issue #461 item 3' });

    // =========================== IntrLine3Cylinder3 ======================

    const line3Cylinder3Find = (io: OracleIO, ln: Line, cyl: Cylinder3): void => {
        const r = new IntrLine3Cylinder3FI().find(ln, cyl);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    };

    family.case('IntrLine3Cylinder3.find', (io) => {
        const ln = line(io, 3);
        const cyl = cylinder(io);
        line3Cylinder3Find(io, ln, cyl);
    }, { exact: true });

    family.case('IntrLine3Cylinder3.doQuery.fi', (io) => {
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const cyl = cylinder(io);
        const r = defaultIntrLine3Cylinder3FIResult();
        intrLine3Cylinder3FIDoQuery(lineOrigin, lineDirection, cyl, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrLine3Cylinder3.find.axisAligned', (io) => {
        const axisOrigin = io.vec(3);
        const axisDirection = io.vec(3);
        const radius = io.real();
        const height = io.real();
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const cyl = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(axisOrigin, axisDirection), radius, height);
        line3Cylinder3Find(io,
            Line.fromOriginDirection(lineOrigin, lineDirection), cyl);
    }, { exact: true });

    // ============================= IntrLine3Torus3 =======================

    // The only case in the family whose compared path calls the C math
    // library, so it is compared with a tolerance instead of bit equality.
    // RootsPolynomial::SolveQuartic reduces the quartic to a depressed quartic
    // and solves a resolvent cubic with cbrt and inverse-trigonometric
    // functions, so a 1 ulp libm difference is amplified through the
    // back-substitutions; Torus3::GetParameters then adds atan2. In a
    // 2000-record deep run the worst scaled error over all outputs was 1.7e-12
    // (one record; 1999 records were within 1e-12 and most were bit
    // identical), so the tolerance is 1e-9, which is still four orders of
    // magnitude below the smallest root separation seen.
    family.case('IntrLine3Torus3.find', (io) => {
        const center = io.vec(3);
        const normal = io.vec(3);
        const direction0 = io.vec(3);
        const direction1 = io.vec(3);
        const radius1 = io.real();
        // The C++ records the increment, not radius0 itself.
        const radius0 = radius1 + io.real();
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const torus = Torus3.fromCenterFrameRadii(center, direction0,
            direction1, normal, radius0, radius1);
        const r = new IntrLine3Torus3FI().find(
            Line.fromOriginDirection(lineOrigin, lineDirection), torus);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.lineParameter[i]);
            io.outVec(r.point[i]);
            io.outReal(r.torusParameter[i][0]);
            io.outReal(r.torusParameter[i][1]);
        }
    }, { tol: 1e-9 });

    // ========================= IntrCylinder3Cylinder3 ====================

    // The hemisphere sampling loop calls cos and sin, so the default scaled
    // tolerance applies.
    family.case('IntrCylinder3Cylinder3.test', (io) => {
        const numTheta = io.integer();
        const numPhi = io.integer();
        const cylinder0 = cylinder(io);
        const cylinder1 = cylinder(io);
        const r = new IntrCylinder3Cylinder3TI(1, numTheta, numPhi)
            .test(cylinder0, cylinder1);
        io.outBool(r.separated);
        io.outVec(r.separatingDirection);
    });

    // Parallel axes take the branch that never calls cos or sin.
    family.case('IntrCylinder3Cylinder3.test.parallel', (io) => {
        const origin0 = io.vec(3);
        const direction = io.vec(3);
        const radius0 = io.real();
        const height0 = io.real();
        const origin1 = io.vec(3);
        const flip = io.boolean();
        const radius1 = io.real();
        const height1 = io.real();
        const direction1 = flip ? negate(direction) : direction;
        const cylinder0 = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin0, direction), radius0, height0);
        const cylinder1 = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin1, direction1), radius1, height1);
        const r = new IntrCylinder3Cylinder3TI(1, 4, 3)
            .test(cylinder0, cylinder1);
        io.outBool(r.separated);
        io.outVec(r.separatingDirection);
    }, { exact: true });

    // The port asserts that both cylinders are finite; upstream feeds the
    // height = -1 sentinel into the finite formulas.
    family.case('IntrCylinder3Cylinder3.test.infiniteDeviation', (io) => {
        const cylinder0 = cylinder(io);
        const origin1 = io.vec(3);
        const direction1 = io.vec(3);
        const radius1 = io.real();
        const height1 = io.real();
        const cylinder1 = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin1, direction1), radius1, height1);
        const r = new IntrCylinder3Cylinder3TI(1, 4, 3)
            .test(cylinder0, cylinder1);
        io.outBool(r.separated);
        io.outVec(r.separatingDirection);
    }, { deviation: 'docs/UPSTREAM-FINDINGS.md IntrCylinder3Cylinder3.h '
        + 'missing IsFinite guard; issue #197' });

    // ======================= IntrCanonicalBox3Cylinder3 ==================

    const canonicalBox3Cylinder3 = (io: OracleIO): void => {
        const extent = io.vec(3);
        const cyl = cylinder(io);
        const r = new IntrCanonicalBox3Cylinder3TI().test(
            CanonicalBox.fromExtent(extent), cyl);
        io.outBool(r.intersect);
    };

    family.case('IntrCanonicalBox3Cylinder3.test', canonicalBox3Cylinder3,
        { exact: true });
    family.case('IntrCanonicalBox3Cylinder3.test.oneZero',
        canonicalBox3Cylinder3, { exact: true });
    family.case('IntrCanonicalBox3Cylinder3.test.twoZeros',
        canonicalBox3Cylinder3, { exact: true });

    // The port corrects the sign typo in the (U1,-D) block of DoQueryNoZeros.
    family.case('IntrCanonicalBox3Cylinder3.test.edgeTypoDeviation',
        canonicalBox3Cylinder3,
        { deviation: 'docs/UPSTREAM-FINDINGS.md IntrCanonicalBox3Cylinder3.h '
            + 'DoQueryNoZeros (U1,-D) sign typo; issue #197' });

    // ======================== IntrOrientedBox3Frustum3 ===================

    family.case('IntrOrientedBox3Frustum3.test', (io) => {
        const frustumOrigin = io.vec(3);
        const dVector = io.vec(3);
        const uVector = io.vec(3);
        const rVector = io.vec(3);
        const dMin = io.real();
        // The C++ records the increment, not dMax itself.
        const dMax = dMin + io.real();
        const uBound = io.real();
        const rBound = io.real();
        const center = io.vec(3);
        const axis: Vector[] = [io.vec(3), io.vec(3), io.vec(3)];
        const extent = io.vec(3);
        const frustum = Frustum3.fromParameters(frustumOrigin, dVector,
            uVector, rVector, dMin, dMax, uBound, rBound);
        const b = OrientedBox.fromCenterAxisExtent(center, axis, extent);
        const r = new IntrOrientedBox3Frustum3TI().test(b, frustum);
        io.outBool(r.intersect);
    }, { exact: true });

    family.finish();
});
