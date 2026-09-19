// Replays oracle/cpp/cases/v33-intersection.cpp (verify group 33,
// intersection). Keep the two files in the same order.
import { describe } from 'vitest';
import { Arc2 } from '../../src/Arc2.js';
import { Capsule } from '../../src/Capsule.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { Frustum3 } from '../../src/Frustum3.js';
import { Halfspace } from '../../src/Halfspace.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Ray } from '../../src/Ray.js';
import { Segment } from '../../src/Segment.js';
import { SegmentMesh } from '../../src/SegmentMesh.js';
import { Triangle } from '../../src/Triangle.js';
import { Vector } from '../../src/Vector.js';

import { IntrEllipsoid3Ellipsoid3TI } from '../../src/IntrEllipsoid3Ellipsoid3.js';
import { IntrHalfspace3Ellipsoid3TI } from '../../src/IntrHalfspace3Ellipsoid3.js';
import { IntrLine2Arc2FI, IntrLine2Arc2TI } from '../../src/IntrLine2Arc2.js';
import { IntrLine2SegmentMesh2FI } from '../../src/IntrLine2SegmentMesh2.js';
import {
    IntrLine3Ellipsoid3FI, IntrLine3Ellipsoid3TI,
    defaultIntrLine3Ellipsoid3FIResult, intrLine3Ellipsoid3FIDoQuery
} from '../../src/IntrLine3Ellipsoid3.js';
import { IntrOrientedBox3Cylinder3TI } from '../../src/IntrOrientedBox3Cylinder3.js';
import { IntrPlane3Plane3FI, IntrPlane3Plane3TI } from '../../src/IntrPlane3Plane3.js';
import {
    IntrPlane3Triangle3FI, IntrPlane3Triangle3TI
} from '../../src/IntrPlane3Triangle3.js';
import {
    IntrRay2Circle2FI, IntrRay2Circle2TI,
    defaultIntrRay2Circle2FIResult, intrRay2Circle2FIDoQuery
} from '../../src/IntrRay2Circle2.js';
import {
    IntrRay2OrientedBox2FI, IntrRay2OrientedBox2TI
} from '../../src/IntrRay2OrientedBox2.js';
import {
    IntrRay3Capsule3FI, IntrRay3Capsule3TI,
    defaultIntrRay3Capsule3FIResult, intrRay3Capsule3FIDoQuery
} from '../../src/IntrRay3Capsule3.js';
import {
    IntrRay3Cylinder3FI,
    defaultIntrRay3Cylinder3FIResult, intrRay3Cylinder3FIDoQuery
} from '../../src/IntrRay3Cylinder3.js';
import {
    IntrRay3OrientedBox3FI, IntrRay3OrientedBox3TI
} from '../../src/IntrRay3OrientedBox3.js';
import {
    IntrSegment2Circle2FI, IntrSegment2Circle2TI,
    defaultIntrSegment2Circle2FIResult, intrSegment2Circle2FIDoQuery
} from '../../src/IntrSegment2Circle2.js';
import {
    IntrSegment2OrientedBox2FI, IntrSegment2OrientedBox2TI
} from '../../src/IntrSegment2OrientedBox2.js';
import {
    IntrSegment3Capsule3FI, IntrSegment3Capsule3TI,
    defaultIntrSegment3Capsule3FIResult, intrSegment3Capsule3FIDoQuery
} from '../../src/IntrSegment3Capsule3.js';
import {
    IntrSegment3Cylinder3FI,
    defaultIntrSegment3Cylinder3FIResult, intrSegment3Cylinder3FIDoQuery
} from '../../src/IntrSegment3Cylinder3.js';
import {
    IntrSegment3OrientedBox3FI, IntrSegment3OrientedBox3TI
} from '../../src/IntrSegment3OrientedBox3.js';
import { IntrSphere3Frustum3TI } from '../../src/IntrSphere3Frustum3.js';
import type {
    IntrLine2SegmentMesh2FIResult
} from '../../src/IntrLine2SegmentMesh2.js';
import type {
    IntrSegment2OrientedBox2FIResult
} from '../../src/IntrSegment2OrientedBox2.js';
import type { IntrPlane3Plane3FIResult } from '../../src/IntrPlane3Plane3.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- input readers, mirroring the C++ generators ----

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

function orientedBox(io: OracleIO, n: number): OrientedBox {
    const center = io.vec(n);
    const axis: Vector[] = [];
    for (let i = 0; i < n; ++i) { axis.push(io.vec(n)); }
    const extent = io.vec(n);
    return OrientedBox.fromCenterAxisExtent(center, axis, extent);
}

function triangle(io: OracleIO, n: number): Triangle {
    const v0 = io.vec(n);
    const v1 = io.vec(n);
    const v2 = io.vec(n);
    return Triangle.fromVertices(v0, v1, v2);
}

function plane(io: OracleIO): Hyperplane {
    const normal = io.vec(3);
    const constant = io.real();
    return Hyperplane.fromNormalConstant(normal, constant);
}

function ellipsoid(io: OracleIO): Hyperellipsoid {
    const center = io.vec(3);
    const axis0 = io.vec(3);
    const axis1 = io.vec(3);
    const axis2 = io.vec(3);
    const extent = io.vec(3);
    return Hyperellipsoid.fromCenterAxisExtent(center, [axis0, axis1, axis2],
        extent);
}

function cylinder(io: OracleIO): Cylinder3 {
    const origin = io.vec(3);
    const direction = io.vec(3);
    const radius = io.real();
    const height = io.real();
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(origin, direction), radius, height);
}

function arc(io: OracleIO): Arc2 {
    const center = io.vec(2);
    const radius = io.real();
    const end0 = io.vec(2);
    const end1 = io.vec(2);
    return Arc2.fromCenterRadiusEnds(center, radius, end0, end1);
}

describe('oracle: v33-intersection', () => {
    const family = new OracleFamily('v33-intersection');

    // ============================= IntrLine2Arc2 =========================

    family.case('IntrLine2Arc2.test', (io) => {
        const a = arc(io);
        const l = line(io, 2);
        const r = new IntrLine2Arc2TI().test(l, a);
        io.outBool(r.intersect);
    }, { exact: true });

    function line2Arc2Find(io: OracleIO): void {
        const a = arc(io);
        const l = line(io, 2);
        const r = new IntrLine2Arc2FI().find(l, a);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }

    family.case('IntrLine2Arc2.find', line2Arc2Find, { exact: true });
    family.case('IntrLine2Arc2.find.throughLatticePoint', line2Arc2Find,
        { exact: true });

    // ============================ IntrRay2Circle2 ========================

    family.case('IntrRay2Circle2.test', (io) => {
        const r0 = ray(io, 2);
        const c = circle(io);
        const r = new IntrRay2Circle2TI().test(r0, c);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay2Circle2.find', (io) => {
        const r0 = ray(io, 2);
        const c = circle(io);
        const r = new IntrRay2Circle2FI().find(r0, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    family.case('IntrRay2Circle2.doQuery.fi', (io) => {
        const rayOrigin = io.vec(2);
        const rayDirection = io.vec(2);
        const c = circle(io);
        const r = defaultIntrRay2Circle2FIResult();
        intrRay2Circle2FIDoQuery(rayOrigin, rayDirection, c, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrRay2Circle2.find.latticeTangent', (io) => {
        const c = circle(io);
        const r0 = ray(io, 2);
        const r = new IntrRay2Circle2FI().find(r0, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    // ========================== IntrSegment2Circle2 ======================

    family.case('IntrSegment2Circle2.test', (io) => {
        const s = segment(io, 2);
        const c = circle(io);
        const r = new IntrSegment2Circle2TI().test(s, c);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment2Circle2.find', (io) => {
        const s = segment(io, 2);
        const c = circle(io);
        const r = new IntrSegment2Circle2FI().find(s, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    family.case('IntrSegment2Circle2.doQuery.fi', (io) => {
        const segOrigin = io.vec(2);
        const segDirection = io.vec(2);
        const segExtent = io.real();
        const c = circle(io);
        const r = defaultIntrSegment2Circle2FIResult();
        intrSegment2Circle2FIDoQuery(segOrigin, segDirection, segExtent, c, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrSegment2Circle2.find.latticeTangent', (io) => {
        const c = circle(io);
        const s = segment(io, 2);
        const r = new IntrSegment2Circle2FI().find(s, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    // ========================= IntrRay2OrientedBox2 ======================

    family.case('IntrRay2OrientedBox2.test', (io) => {
        const r0 = ray(io, 2);
        const b = orientedBox(io, 2);
        const r = new IntrRay2OrientedBox2TI().test(r0, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay2OrientedBox2.find', (io) => {
        const r0 = ray(io, 2);
        const b = orientedBox(io, 2);
        const r = new IntrRay2OrientedBox2FI().find(r0, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    family.case('IntrRay2OrientedBox2.find.throughPoint', (io) => {
        const b = orientedBox(io, 2);
        const r0 = ray(io, 2);
        const r = new IntrRay2OrientedBox2FI().find(r0, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    // ========================= IntrRay3OrientedBox3 ======================

    function ray3Box3Find(io: OracleIO, boxFirst: boolean): void {
        const b = boxFirst ? orientedBox(io, 3) : undefined;
        const r0 = ray(io, 3);
        const box = b ?? orientedBox(io, 3);
        const r = new IntrRay3OrientedBox3FI().find(r0, box);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }

    family.case('IntrRay3OrientedBox3.test', (io) => {
        const r0 = ray(io, 3);
        const b = orientedBox(io, 3);
        const r = new IntrRay3OrientedBox3TI().test(r0, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay3OrientedBox3.find', (io) => { ray3Box3Find(io, false); },
        { exact: true });
    family.case('IntrRay3OrientedBox3.find.throughPoint',
        (io) => { ray3Box3Find(io, true); }, { exact: true });

    // ======================= IntrSegment3OrientedBox3 ====================

    function segment3Box3Find(io: OracleIO, boxFirst: boolean): void {
        const b = boxFirst ? orientedBox(io, 3) : undefined;
        const s = segment(io, 3);
        const box = b ?? orientedBox(io, 3);
        const r = new IntrSegment3OrientedBox3FI().find(s, box);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }

    family.case('IntrSegment3OrientedBox3.test', (io) => {
        const s = segment(io, 3);
        const b = orientedBox(io, 3);
        const r = new IntrSegment3OrientedBox3TI().test(s, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment3OrientedBox3.find',
        (io) => { segment3Box3Find(io, false); }, { exact: true });
    family.case('IntrSegment3OrientedBox3.find.throughPoint',
        (io) => { segment3Box3Find(io, true); }, { exact: true });

    // ======================= IntrSegment2OrientedBox2 ====================

    function outSeg2Box2FI(io: OracleIO, r: IntrSegment2OrientedBox2FIResult,
        withPoints: boolean): void {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outReal(r.cdeParameter[0]);
        io.outReal(r.cdeParameter[1]);
        if (withPoints) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }

    family.case('IntrSegment2OrientedBox2.test', (io) => {
        const s = segment(io, 2);
        const b = orientedBox(io, 2);
        const r = new IntrSegment2OrientedBox2TI().test(s, b);
        io.outBool(r.intersect);
    }, { exact: true });

    // 'point' is the field on which the port deliberately deviates from
    // upstream's frame-mixing expression; it is compared by find.originBox
    // (where upstream is sound) and by find.frameDeviation.
    family.case('IntrSegment2OrientedBox2.find', (io) => {
        const s = segment(io, 2);
        const b = orientedBox(io, 2);
        const r = new IntrSegment2OrientedBox2FI().find(s, b);
        outSeg2Box2FI(io, r, false);
    }, { exact: true });

    family.case('IntrSegment2OrientedBox2.find.throughPoint', (io) => {
        const b = orientedBox(io, 2);
        const s = segment(io, 2);
        const r = new IntrSegment2OrientedBox2FI().find(s, b);
        outSeg2Box2FI(io, r, false);
    }, { exact: true });

    family.case('IntrSegment2OrientedBox2.find.originBox', (io) => {
        const b = orientedBox(io, 2);
        const s = segment(io, 2);
        const r = new IntrSegment2OrientedBox2FI().find(s, b);
        outSeg2Box2FI(io, r, true);
    }, { exact: true });

    family.case('IntrSegment2OrientedBox2.find.degenerate', (io) => {
        const b = orientedBox(io, 2);
        const p = io.vec(2);
        const s = Segment.fromEndpoints(p, p);
        const r = new IntrSegment2OrientedBox2FI().find(s, b);
        outSeg2Box2FI(io, r, true);
    }, { exact: true });

    family.case('IntrSegment2OrientedBox2.find.frameDeviation', (io) => {
        const b = orientedBox(io, 2);
        const s = segment(io, 2);
        const r = new IntrSegment2OrientedBox2FI().find(s, b);
        outSeg2Box2FI(io, r, true);
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrSegment2OrientedBox2.h FIQuery '
            + '(box-frame vector added to the world-space box centre); issue #255'
    });

    // ============================ IntrPlane3Plane3 =======================

    function outPlane3Plane3FI(io: OracleIO, r: IntrPlane3Plane3FIResult): void {
        io.outBool(r.intersect);
        io.outBool(r.isLine);
        if (r.intersect) {
            if (r.isLine) {
                io.outVec(r.line.origin);
                io.outVec(r.line.direction);
            }
            else {
                io.outVec(r.plane.normal);
                io.outReal(r.plane.constant);
            }
        }
    }

    family.case('IntrPlane3Plane3.test', (io) => {
        const p0 = plane(io);
        const p1 = plane(io);
        const r = new IntrPlane3Plane3TI().test(p0, p1);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrPlane3Plane3.find', (io) => {
        const p0 = plane(io);
        const p1 = plane(io);
        outPlane3Plane3FI(io, new IntrPlane3Plane3FI().find(p0, p1));
    }, { exact: true });

    family.case('IntrPlane3Plane3.find.parallel', (io) => {
        const p0 = plane(io);
        const p1 = plane(io);
        outPlane3Plane3FI(io, new IntrPlane3Plane3FI().find(p0, p1));
        io.outBool(new IntrPlane3Plane3TI().test(p0, p1).intersect);
    }, { exact: true });

    // ========================== IntrPlane3Triangle3 ======================

    family.case('IntrPlane3Triangle3.test', (io) => {
        const p = plane(io);
        const t = triangle(io, 3);
        const r = new IntrPlane3Triangle3TI().test(p, t);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outBool(r.isInterior);
    }, { exact: true });

    family.case('IntrPlane3Triangle3.find', (io) => {
        const p = plane(io);
        const t = triangle(io, 3);
        const r = new IntrPlane3Triangle3FI().find(p, t);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outBool(r.isInterior);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
        io.outVec(r.point[2]);
    }, { exact: true });

    family.case('IntrPlane3Triangle3.find.onPlane', (io) => {
        const p = plane(io);
        const t = triangle(io, 3);
        const r = new IntrPlane3Triangle3FI().find(p, t);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outBool(r.isInterior);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
        io.outVec(r.point[2]);
        const ti = new IntrPlane3Triangle3TI().test(p, t);
        io.outBool(ti.intersect);
        io.outInt(ti.numIntersections);
        io.outBool(ti.isInterior);
    }, { exact: true });

    // ========================= IntrLine3Ellipsoid3 =======================

    family.case('IntrLine3Ellipsoid3.test', (io) => {
        const l = line(io, 3);
        const e = ellipsoid(io);
        const r = new IntrLine3Ellipsoid3TI().test(l, e);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3Ellipsoid3.find', (io) => {
        const l = line(io, 3);
        const e = ellipsoid(io);
        const r = new IntrLine3Ellipsoid3FI().find(l, e);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrLine3Ellipsoid3.doQuery.fi', (io) => {
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const e = ellipsoid(io);
        const r = defaultIntrLine3Ellipsoid3FIResult();
        intrLine3Ellipsoid3FIDoQuery(lineOrigin, lineDirection, e, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrLine3Ellipsoid3.find.latticeTangent', (io) => {
        const e = ellipsoid(io);
        const l = line(io, 3);
        const r = new IntrLine3Ellipsoid3FI().find(l, e);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    // ======================= IntrHalfspace3Ellipsoid3 ====================

    family.case('IntrHalfspace3Ellipsoid3.test', (io) => {
        const normal = io.vec(3);
        const constant = io.real();
        const h = Halfspace.fromNormalConstant(normal, constant);
        const e = ellipsoid(io);
        io.outBool(new IntrHalfspace3Ellipsoid3TI().test(h, e).intersect);
    }, { exact: true });

    family.case('IntrHalfspace3Ellipsoid3.test.tangent', (io) => {
        const e = ellipsoid(io);
        const normal = io.vec(3);
        const constant = io.real();
        const h = Halfspace.fromNormalConstant(normal, constant);
        io.outBool(new IntrHalfspace3Ellipsoid3TI().test(h, e).intersect);
    }, { exact: true });

    // ====================== IntrEllipsoid3Ellipsoid3 =====================
    //
    // The eigensolver and the bisection call the C math library, but both
    // outputs are discrete and compare exactly.

    function ellipsoidPairTest(io: OracleIO): void {
        const a = ellipsoid(io);
        const b = ellipsoid(io);
        const r = new IntrEllipsoid3Ellipsoid3TI().test(a, b);
        io.outBool(r.intersect);
        io.outInt(r.classification);
    }

    family.case('IntrEllipsoid3Ellipsoid3.test', ellipsoidPairTest,
        { exact: true });

    family.case('IntrEllipsoid3Ellipsoid3.test.concentric', (io) => {
        const center = io.vec(3);
        const aAxis = [io.vec(3), io.vec(3), io.vec(3)];
        const aExtent = io.vec(3);
        const bAxis = [io.vec(3), io.vec(3), io.vec(3)];
        const bExtent = io.vec(3);
        const a = Hyperellipsoid.fromCenterAxisExtent(center, aAxis, aExtent);
        const b = Hyperellipsoid.fromCenterAxisExtent(center, bAxis, bExtent);
        const r = new IntrEllipsoid3Ellipsoid3TI().test(a, b);
        io.outBool(r.intersect);
        io.outInt(r.classification);
    }, { exact: true });

    family.case('IntrEllipsoid3Ellipsoid3.test.equalEigenvalues',
        ellipsoidPairTest, { exact: true });

    family.case('IntrEllipsoid3Ellipsoid3.test.equalEigenvaluesDeviation',
        ellipsoidPairTest, {
        exact: true,
        deviation: 'src/IntrEllipsoid3Ellipsoid3.ts KNOWN UPSTREAM DEFECT: the '
            + "'d0 > d1 = d2' branch folds the coefficient of the distinct "
            + 'eigenvalue instead of the one that shares it'
    });

    // Both sides throw "Unexpected condition." on this pinned configuration.
    family.case('IntrEllipsoid3Ellipsoid3.test.bracketAssert',
        ellipsoidPairTest, { exact: true });

    // ====================== IntrOrientedBox3Cylinder3 ====================

    family.case('IntrOrientedBox3Cylinder3.test', (io) => {
        const b = orientedBox(io, 3);
        const c = cylinder(io);
        io.outBool(new IntrOrientedBox3Cylinder3TI().test(b, c).intersect);
    }, { exact: true });

    family.case('IntrOrientedBox3Cylinder3.test.edgeTypoDeviation', (io) => {
        const b = orientedBox(io, 3);
        const c = cylinder(io);
        io.outBool(new IntrOrientedBox3Cylinder3TI().test(b, c).intersect);
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrCanonicalBox3Cylinder3.h '
            + 'DoQueryNoZeros (U1,-D) sign typo, inherited by the wrapper; issue #197'
    });

    // Upstream asserts that the cylinder is finite, so both sides throw.
    family.case('IntrOrientedBox3Cylinder3.test.infinite', (io) => {
        const b = orientedBox(io, 3);
        const c = cylinder(io);
        io.outBool(new IntrOrientedBox3Cylinder3TI().test(b, c).intersect);
    }, { exact: true });

    // ========================== IntrRay3Cylinder3 ========================

    function ray3Cylinder3Find(io: OracleIO, cylinderFirst: boolean): void {
        const c0 = cylinderFirst ? cylinder(io) : undefined;
        const r0 = ray(io, 3);
        const c = c0 ?? cylinder(io);
        const r = new IntrRay3Cylinder3FI().find(r0, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }

    family.case('IntrRay3Cylinder3.find',
        (io) => { ray3Cylinder3Find(io, false); }, { exact: true });

    family.case('IntrRay3Cylinder3.doQuery.fi', (io) => {
        const rayOrigin = io.vec(3);
        const rayDirection = io.vec(3);
        const c = cylinder(io);
        const r = defaultIntrRay3Cylinder3FIResult();
        intrRay3Cylinder3FIDoQuery(rayOrigin, rayDirection, c, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrRay3Cylinder3.find.latticeAimed',
        (io) => { ray3Cylinder3Find(io, true); }, { exact: true });

    family.case('IntrRay3Cylinder3.find.infiniteDeviation',
        (io) => { ray3Cylinder3Find(io, false); }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md cylinder headers, missing IsFinite '
            + 'guard (the port asserts); issues #197, #206, #255'
    });

    // ======================== IntrSegment3Cylinder3 ======================

    function segment3Cylinder3Find(io: OracleIO, cylinderFirst: boolean): void {
        const c0 = cylinderFirst ? cylinder(io) : undefined;
        const s = segment(io, 3);
        const c = c0 ?? cylinder(io);
        const r = new IntrSegment3Cylinder3FI().find(s, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }

    family.case('IntrSegment3Cylinder3.find',
        (io) => { segment3Cylinder3Find(io, false); }, { exact: true });

    family.case('IntrSegment3Cylinder3.doQuery.fi', (io) => {
        const segOrigin = io.vec(3);
        const segDirection = io.vec(3);
        const segExtent = io.real();
        const c = cylinder(io);
        const r = defaultIntrSegment3Cylinder3FIResult();
        intrSegment3Cylinder3FIDoQuery(segOrigin, segDirection, segExtent, c, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrSegment3Cylinder3.find.latticeAimed',
        (io) => { segment3Cylinder3Find(io, true); }, { exact: true });

    // The documented upstream limitation, preserved: a zero-length segment
    // produces NaN parameters on both sides.
    family.case('IntrSegment3Cylinder3.find.degenerateSegment', (io) => {
        const c = cylinder(io);
        const p = io.vec(3);
        const s = Segment.fromEndpoints(p, p);
        const r = new IntrSegment3Cylinder3FI().find(s, c);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrSegment3Cylinder3.find.infiniteDeviation',
        (io) => { segment3Cylinder3Find(io, false); }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md cylinder headers, missing IsFinite '
            + 'guard (the port asserts); issues #197, #206, #255'
    });

    // =========================== IntrRay3Capsule3 ========================

    function ray3Capsule3Test(io: OracleIO): void {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const radius = io.real();
        const capsule = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), radius);
        const r = new IntrRay3Capsule3TI().test(
            Ray.fromOriginDirection(origin, direction), capsule);
        io.outBool(r.intersect);
    }

    family.case('IntrRay3Capsule3.test', ray3Capsule3Test, { exact: true });

    family.case('IntrRay3Capsule3.test.clampDeviation', ray3Capsule3Test, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md DistRaySegment.h, the missing '
            + 's0 >= 0 clamp in the parallel branch and regions 1 and 5; issue #126'
    });

    family.case('IntrRay3Capsule3.find', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const radius = io.real();
        const capsule = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), radius);
        const r = new IntrRay3Capsule3FI().find(
            Ray.fromOriginDirection(origin, direction), capsule);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrRay3Capsule3.doQuery.fi', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const radius = io.real();
        const capsule = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), radius);
        const r = defaultIntrRay3Capsule3FIResult();
        intrRay3Capsule3FIDoQuery(origin, direction, capsule, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrRay3Capsule3.doQuery.capTangentDeviation', (io) => {
        io.vec(3); // the capsule centre, implicit in the endpoints below
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const radius = io.real();
        const origin = io.vec(3);
        const direction = io.vec(3);
        const capsule = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), radius);
        const r = defaultIntrRay3Capsule3FIResult();
        intrRay3Capsule3FIDoQuery(origin, direction, capsule, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrLine3Capsule3.h, "intersect" is '
            + 'never set when only one hemisphere root is accepted; issue #461'
    });

    // ========================= IntrSegment3Capsule3 ======================

    family.case('IntrSegment3Capsule3.test', (io) => {
        const e0 = io.vec(3);
        const e1 = io.vec(3);
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const radius = io.real();
        const capsule = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), radius);
        const r = new IntrSegment3Capsule3TI().test(
            Segment.fromEndpoints(e0, e1), capsule);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment3Capsule3.find', (io) => {
        const e0 = io.vec(3);
        const e1 = io.vec(3);
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const radius = io.real();
        const capsule = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), radius);
        const r = new IntrSegment3Capsule3FI().find(
            Segment.fromEndpoints(e0, e1), capsule);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (r.intersect) {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    function segment3Capsule3DoQuery(io: OracleIO, centreFirst: boolean): void {
        if (centreFirst) { io.vec(3); }
        const segOrigin = centreFirst ? null : io.vec(3);
        const segDirection = centreFirst ? null : io.vec(3);
        const segExtent = centreFirst ? null : io.real();
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const radius = io.real();
        const origin = segOrigin ?? io.vec(3);
        const direction = segDirection ?? io.vec(3);
        const extent = segExtent ?? io.real();
        const capsule = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), radius);
        const r = defaultIntrSegment3Capsule3FIResult();
        intrSegment3Capsule3FIDoQuery(origin, direction, extent, capsule, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }

    family.case('IntrSegment3Capsule3.doQuery.fi',
        (io) => { segment3Capsule3DoQuery(io, false); }, { exact: true });

    family.case('IntrSegment3Capsule3.doQuery.junctionDeviation',
        (io) => { segment3Capsule3DoQuery(io, true); }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrLine3Capsule3.h, a cap-junction '
            + 'root accepted twice collapses the interval; issue #461'
    });

    // ======================== IntrLine2SegmentMesh2 ======================
    //
    // Upstream's std::sort is not stable, so records with equal line
    // parameters have an unspecified relative order. Both sides re-sort by a
    // total order before emitting; see the C++ case comment.

    function mesh2(vertices: Vector[], indices: [number, number][],
        topology: number): SegmentMesh {
        if (topology === 0) { return SegmentMesh.fromDisjoint(vertices); }
        if (topology === 1) { return SegmentMesh.fromContiguous(vertices, true); }
        if (topology === 2) { return SegmentMesh.fromContiguous(vertices, false); }
        return SegmentMesh.fromIndexed(vertices, indices, true);
    }

    function outMeshResult(io: OracleIO,
        r: IntrLine2SegmentMesh2FIResult): void {
        const hits = r.intersections.map((x) => ({
            lineParameter: x.lineParameter,
            meshParameter: x.meshSegmentParameter,
            px: x.point.get(0),
            py: x.point.get(1),
            i0: x.indexPair[0],
            i1: x.indexPair[1]
        }));
        hits.sort((a, b) => {
            if (a.lineParameter !== b.lineParameter) {
                return a.lineParameter < b.lineParameter ? -1 : 1;
            }
            if (a.i0 !== b.i0) { return a.i0 < b.i0 ? -1 : 1; }
            if (a.i1 !== b.i1) { return a.i1 < b.i1 ? -1 : 1; }
            if (a.meshParameter !== b.meshParameter) {
                return a.meshParameter < b.meshParameter ? -1 : 1;
            }
            if (a.px !== b.px) { return a.px < b.px ? -1 : 1; }
            if (a.py !== b.py) { return a.py < b.py ? -1 : 1; }
            return 0;
        });
        io.outInt(hits.length);
        for (const h of hits) {
            io.outInt(h.i0);
            io.outInt(h.i1);
            io.outReal(h.lineParameter);
            io.outReal(h.meshParameter);
            io.outReal(h.px);
            io.outReal(h.py);
        }
    }

    function meshCase(io: OracleIO): void {
        const topology = io.integer();
        const vertices: Vector[] = [];
        for (let i = 0; i < 6; ++i) { vertices.push(io.vec(2)); }
        const indices: [number, number][] = [];
        for (let i = 0; i < 4; ++i) {
            const a = io.integer();
            const b = io.integer();
            indices.push([a, b]);
        }
        const l = line(io, 2);
        const mesh = mesh2(vertices, indices, topology);
        const r = new IntrLine2SegmentMesh2FI().find(l, mesh);
        outMeshResult(io, r);
    }

    family.case('IntrLine2SegmentMesh2.find', meshCase, { exact: true });
    family.case('IntrLine2SegmentMesh2.find.collinear', meshCase,
        { exact: true });

    // ========================= IntrSphere3Frustum3 =======================

    function sphereFrustumCase(io: OracleIO): void {
        const origin = io.vec(3);
        const dVector = io.vec(3);
        const uVector = io.vec(3);
        const rVector = io.vec(3);
        const dMin = io.real();
        const dMax = io.real();
        const uBound = io.real();
        const rBound = io.real();
        const point = io.vec(3);
        const radius = io.real();
        const frustum = Frustum3.fromParameters(origin, dVector, uVector,
            rVector, dMin, dMax, uBound, rBound);
        const sphere = Hypersphere.fromCenterRadius(point, radius);
        io.outBool(new IntrSphere3Frustum3TI().test(sphere, frustum).intersect);
    }

    family.case('IntrSphere3Frustum3.test', sphereFrustumCase, { exact: true });

    family.case('IntrSphere3Frustum3.test.clampDeviation', sphereFrustumCase, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md DistPoint3Frustum3.h, the two '
            + 'unclamped far-edge assignments; issue #421'
    });

    family.finish();
});
