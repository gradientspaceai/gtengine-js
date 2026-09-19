// Replays oracle/cpp/cases/v32-intersection.cpp (verify group 32,
// intersection). Keep the two files in the same order.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { Matrix } from '../../src/Matrix.js';
import { Ray } from '../../src/Ray.js';
import { Rectangle } from '../../src/Rectangle.js';
import { Segment } from '../../src/Segment.js';
import { Triangle } from '../../src/Triangle.js';
import { Vector } from '../../src/Vector.js';
import { IntrAlignedBox2Circle2FI, IntrAlignedBox2Circle2TI } from '../../src/IntrAlignedBox2Circle2.js';
import { IntrAlignedBox3Cylinder3TI } from '../../src/IntrAlignedBox3Cylinder3.js';
import { IntrAlignedBox3Sphere3FI, IntrAlignedBox3Sphere3TI } from '../../src/IntrAlignedBox3Sphere3.js';
import {
    IntrConvexPolygonHyperplaneFI, IntrConvexPolygonHyperplaneTI
} from '../../src/IntrConvexPolygonHyperplane.js';
import {
    IntrEllipse2Ellipse2FI, IntrEllipse2Ellipse2TI,
    intrEllipse2Ellipse2InfinitePoints
} from '../../src/IntrEllipse2Ellipse2.js';
import {
    IntrRay3Sphere3FI, IntrRay3Sphere3TI,
    defaultIntrRay3Sphere3FIResult, intrRay3Sphere3FIDoQuery
} from '../../src/IntrRay3Sphere3.js';
import {
    IntrSegment2AlignedBox2FI, IntrSegment2AlignedBox2TI,
    defaultIntrSegment2AlignedBox2FIResult,
    defaultIntrSegment2AlignedBox2TIResult,
    intrSegment2AlignedBox2FIDoQuery, intrSegment2AlignedBox2TIDoQuery
} from '../../src/IntrSegment2AlignedBox2.js';
import {
    IntrSegment2Segment2FI, IntrSegment2Segment2TI
} from '../../src/IntrSegment2Segment2.js';
import {
    IntrSegment2Triangle2FI, IntrSegment2Triangle2TI,
    defaultIntrSegment2Triangle2FIResult, intrSegment2Triangle2FIDoQuery
} from '../../src/IntrSegment2Triangle2.js';
import {
    IntrSegment3AlignedBox3FI, IntrSegment3AlignedBox3TI,
    defaultIntrSegment3AlignedBox3FIResult,
    defaultIntrSegment3AlignedBox3TIResult,
    intrSegment3AlignedBox3FIDoQuery, intrSegment3AlignedBox3TIDoQuery
} from '../../src/IntrSegment3AlignedBox3.js';
import {
    IntrSegment3Rectangle3FI, IntrSegment3Rectangle3TI
} from '../../src/IntrSegment3Rectangle3.js';
import {
    IntrSegment3Sphere3FI, IntrSegment3Sphere3TI,
    defaultIntrSegment3Sphere3FIResult, intrSegment3Sphere3FIDoQuery
} from '../../src/IntrSegment3Sphere3.js';
import { IntrSphere3Sphere3FI, IntrSphere3Sphere3TI } from '../../src/IntrSphere3Sphere3.js';
import { IntrSphere3Triangle3FI } from '../../src/IntrSphere3Triangle3.js';
import { IntrTriangle3Cylinder3TI } from '../../src/IntrTriangle3Cylinder3.js';
import type { IntrSphere3Sphere3FIResult } from '../../src/IntrSphere3Sphere3.js';
import type { IntrSegment2AlignedBox2FIResult } from '../../src/IntrSegment2AlignedBox2.js';
import type { IntrSegment2Segment2FIResult } from '../../src/IntrSegment2Segment2.js';
import type { IntrEllipse2Ellipse2FIResult } from '../../src/IntrEllipse2Ellipse2.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- input readers, mirroring the C++ generators ----

function box(io: OracleIO, n: number): AlignedBox {
    const min = io.vec(n);
    const max = io.vec(n);
    return AlignedBox.fromMinMax(min, max);
}

function segment(io: OracleIO, n: number): Segment {
    const p0 = io.vec(n);
    const p1 = io.vec(n);
    return Segment.fromEndpoints(p0, p1);
}

function ray(io: OracleIO, n: number): Ray {
    const origin = io.vec(n);
    const direction = io.vec(n);
    return Ray.fromOriginDirection(origin, direction);
}

function sphere(io: OracleIO, n: number): Hypersphere {
    const center = io.vec(n);
    const radius = io.real();
    return Hypersphere.fromCenterRadius(center, radius);
}

function triangle(io: OracleIO, n: number): Triangle {
    const v0 = io.vec(n);
    const v1 = io.vec(n);
    const v2 = io.vec(n);
    return Triangle.fromVertices(v0, v1, v2);
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

function cylinder(io: OracleIO): Cylinder3 {
    const origin = io.vec(3);
    const direction = io.vec(3);
    const radius = io.real();
    const height = io.real();
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(origin, direction), radius, height);
}

function ellipse(io: OracleIO): Hyperellipsoid {
    const center = io.vec(2);
    const axis0 = io.vec(2);
    const axis1 = io.vec(2);
    const extent = io.vec(2);
    return Hyperellipsoid.fromCenterAxisExtent(center, [axis0, axis1], extent);
}

describe('oracle: v32-intersection', () => {
    const family = new OracleFamily('v32-intersection');

    // ============================ IntrRay3Sphere3 ========================

    family.case('IntrRay3Sphere3.test', (io) => {
        const r0 = ray(io, 3);
        const s = sphere(io, 3);
        const r = new IntrRay3Sphere3TI().test(r0, s);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay3Sphere3.find', (io) => {
        const r0 = ray(io, 3);
        const s = sphere(io, 3);
        const r = new IntrRay3Sphere3FI().find(r0, s);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrRay3Sphere3.doQuery.fi', (io) => {
        const rayOrigin = io.vec(3);
        const rayDirection = io.vec(3);
        const s = sphere(io, 3);
        const r = defaultIntrRay3Sphere3FIResult();
        intrRay3Sphere3FIDoQuery(rayOrigin, rayDirection, s, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrRay3Sphere3.find.tangent', (io) => {
        const s = sphere(io, 3);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = new IntrRay3Sphere3FI().find(
            Ray.fromOriginDirection(origin, direction), s);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    // ========================== IntrSegment3Sphere3 ======================

    family.case('IntrSegment3Sphere3.test', (io) => {
        const sg = segment(io, 3);
        const s = sphere(io, 3);
        const r = new IntrSegment3Sphere3TI().test(sg, s);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment3Sphere3.test.containedDeviation', (io) => {
        const s = sphere(io, 3);
        const sg = segment(io, 3);
        const r = new IntrSegment3Sphere3TI().test(sg, s);
        io.outBool(r.intersect);
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrSegment3Sphere3.h TIQuery '
            + '(a segment strictly inside the sphere reports no intersection); issue #203'
    });

    family.case('IntrSegment3Sphere3.find', (io) => {
        const sg = segment(io, 3);
        const s = sphere(io, 3);
        const r = new IntrSegment3Sphere3FI().find(sg, s);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrSegment3Sphere3.doQuery.fi', (io) => {
        const segOrigin = io.vec(3);
        const segDirection = io.vec(3);
        const segExtent = io.real();
        const s = sphere(io, 3);
        const r = defaultIntrSegment3Sphere3FIResult();
        intrSegment3Sphere3FIDoQuery(segOrigin, segDirection, segExtent, s, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrSegment3Sphere3.doQuery.fi.tangent', (io) => {
        const s = sphere(io, 3);
        const segOrigin = io.vec(3);
        const segDirection = io.vec(3);
        const segExtent = io.real();
        const r = defaultIntrSegment3Sphere3FIResult();
        intrSegment3Sphere3FIDoQuery(segOrigin, segDirection, segExtent, s, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    // ========================== IntrSphere3Sphere3 =======================

    function outSphereFI(io: OracleIO, r: IntrSphere3Sphere3FIResult): void {
        io.outBool(r.intersect);
        io.outInt(r.type);
        io.outVec(r.point);
        io.outVec(r.circle.center);
        io.outVec(r.circle.normal);
        io.outReal(r.circle.radius);
    }

    family.case('IntrSphere3Sphere3.test', (io) => {
        const s0 = sphere(io, 3);
        const s1 = sphere(io, 3);
        const r = new IntrSphere3Sphere3TI().test(s0, s1);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSphere3Sphere3.find', (io) => {
        const s0 = sphere(io, 3);
        const s1 = sphere(io, 3);
        outSphereFI(io, new IntrSphere3Sphere3FI().find(s0, s1));
    }, { exact: true });

    family.case('IntrSphere3Sphere3.find.tangent', (io) => {
        const s0 = sphere(io, 3);
        const s1 = sphere(io, 3);
        outSphereFI(io, new IntrSphere3Sphere3FI().find(s0, s1));
    }, { exact: true });

    family.case('IntrSphere3Sphere3.find.internalTangentDeviation', (io) => {
        const s0 = sphere(io, 3);
        const s1 = sphere(io, 3);
        outSphereFI(io, new IntrSphere3Sphere3FI().find(s0, s1));
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrSphere3Sphere3.h FIQuery type-4 '
            + 'branch (the internal-tangency contact point is the antipode); issue #203'
    });

    // ======================== IntrSegment2AlignedBox2 ====================

    function outSeg2Box2FI(io: OracleIO, r: IntrSegment2AlignedBox2FIResult): void {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outReal(r.cdeParameter[0]);
        io.outReal(r.cdeParameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }

    family.case('IntrSegment2AlignedBox2.test', (io) => {
        const sg = segment(io, 2);
        const b = box(io, 2);
        const r = new IntrSegment2AlignedBox2TI().test(sg, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment2AlignedBox2.find', (io) => {
        const sg = segment(io, 2);
        const b = box(io, 2);
        outSeg2Box2FI(io, new IntrSegment2AlignedBox2FI().find(sg, b));
    }, { exact: true });

    family.case('IntrSegment2AlignedBox2.find.degenerate', (io) => {
        const p = io.vec(2);
        const b = box(io, 2);
        const sg = Segment.fromEndpoints(p, p.clone());
        outSeg2Box2FI(io, new IntrSegment2AlignedBox2FI().find(sg, b));
    }, { exact: true });

    family.case('IntrSegment2AlignedBox2.doQuery.ti', (io) => {
        const segOrigin = io.vec(2);
        const segDirection = io.vec(2);
        const segExtent = io.real();
        const boxExtent = io.vec(2);
        const r = defaultIntrSegment2AlignedBox2TIResult();
        intrSegment2AlignedBox2TIDoQuery(segOrigin, segDirection, segExtent,
            boxExtent, r);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment2AlignedBox2.doQuery.fi', (io) => {
        const segOrigin = io.vec(2);
        const segDirection = io.vec(2);
        const segExtent = io.real();
        const boxExtent = io.vec(2);
        const r = defaultIntrSegment2AlignedBox2FIResult();
        intrSegment2AlignedBox2FIDoQuery(segOrigin, segDirection, segExtent,
            boxExtent, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    // ======================== IntrSegment3AlignedBox3 ====================

    family.case('IntrSegment3AlignedBox3.test', (io) => {
        const sg = segment(io, 3);
        const b = box(io, 3);
        const r = new IntrSegment3AlignedBox3TI().test(sg, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment3AlignedBox3.find', (io) => {
        const sg = segment(io, 3);
        const b = box(io, 3);
        const r = new IntrSegment3AlignedBox3FI().find(sg, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrSegment3AlignedBox3.doQuery.ti', (io) => {
        const segOrigin = io.vec(3);
        const segDirection = io.vec(3);
        const segExtent = io.real();
        const boxExtent = io.vec(3);
        const r = defaultIntrSegment3AlignedBox3TIResult();
        intrSegment3AlignedBox3TIDoQuery(segOrigin, segDirection, segExtent,
            boxExtent, r);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment3AlignedBox3.doQuery.fi', (io) => {
        const segOrigin = io.vec(3);
        const segDirection = io.vec(3);
        const segExtent = io.real();
        const boxExtent = io.vec(3);
        const r = defaultIntrSegment3AlignedBox3FIResult();
        intrSegment3AlignedBox3FIDoQuery(segOrigin, segDirection, segExtent,
            boxExtent, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrSegment3AlignedBox3.find.throughPoint', (io) => {
        const b = box(io, 3);
        const sg = segment(io, 3);
        const r = new IntrSegment3AlignedBox3FI().find(sg, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    // ========================= IntrSegment2Triangle2 =====================

    family.case('IntrSegment2Triangle2.test', (io) => {
        const sg = segment(io, 2);
        const t = triangle(io, 2);
        const r = new IntrSegment2Triangle2TI().test(sg, t);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment2Triangle2.find', (io) => {
        const sg = segment(io, 2);
        const t = triangle(io, 2);
        const r = new IntrSegment2Triangle2FI().find(sg, t);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrSegment2Triangle2.doQuery.fi', (io) => {
        const origin = io.vec(2);
        const direction = io.vec(2);
        const t = triangle(io, 2);
        const r = defaultIntrSegment2Triangle2FIResult();
        intrSegment2Triangle2FIDoQuery(origin, direction, t, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrSegment2Triangle2.find.throughPoint', (io) => {
        const t = triangle(io, 2);
        io.vec(2);  // the recorded lattice direction of the construction
        const sg = segment(io, 2);
        const r = new IntrSegment2Triangle2FI().find(sg, t);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    // ======================== IntrSegment3Rectangle3 =====================

    family.case('IntrSegment3Rectangle3.test', (io) => {
        const sg = segment(io, 3);
        const rect = rectangle3(io);
        const r = new IntrSegment3Rectangle3TI().test(sg, rect);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment3Rectangle3.find', (io) => {
        const sg = segment(io, 3);
        const rect = rectangle3(io);
        const r = new IntrSegment3Rectangle3FI().find(sg, rect);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.rectCoord[0]);
            io.outReal(r.rectCoord[1]);
            io.outReal(r.rectCoord[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    family.case('IntrSegment3Rectangle3.find.throughPoint', (io) => {
        const rect = rectangle3(io);
        io.vec(3);  // the recorded lattice direction of the construction
        const sg = segment(io, 3);
        const r = new IntrSegment3Rectangle3FI().find(sg, rect);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.rectCoord[0]);
            io.outReal(r.rectCoord[1]);
            io.outReal(r.rectCoord[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    // ========================= IntrSegment2Segment2 ======================

    function outSeg2Seg2FI(io: OracleIO, r: IntrSegment2Segment2FIResult): void {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.segment0Parameter[0]);
        io.outReal(r.segment0Parameter[1]);
        io.outReal(r.segment1Parameter[0]);
        io.outReal(r.segment1Parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }

    family.case('IntrSegment2Segment2.test', (io) => {
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        const r = new IntrSegment2Segment2TI().test(s0, s1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    family.case('IntrSegment2Segment2.testExact', (io) => {
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        const r = new IntrSegment2Segment2TI().testExact(s0, s1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    family.case('IntrSegment2Segment2.find', (io) => {
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        outSeg2Seg2FI(io, new IntrSegment2Segment2FI().find(s0, s1));
    }, { exact: true });

    family.case('IntrSegment2Segment2.findExact', (io) => {
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        outSeg2Seg2FI(io, new IntrSegment2Segment2FI().findExact(s0, s1));
    }, { exact: true });

    family.case('IntrSegment2Segment2.find.collinear', (io) => {
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        outSeg2Seg2FI(io, new IntrSegment2Segment2FI().find(s0, s1));
    }, { exact: true });

    family.case('IntrSegment2Segment2.findExact.collinear', (io) => {
        io.boolean();
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        outSeg2Seg2FI(io, new IntrSegment2Segment2FI().findExact(s0, s1));
    }, { exact: true });

    family.case('IntrSegment2Segment2.test.collinear', (io) => {
        io.boolean();
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        const r = new IntrSegment2Segment2TI().test(s0, s1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    family.case('IntrSegment2Segment2.testExact.collinear', (io) => {
        io.boolean();
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        const r = new IntrSegment2Segment2TI().testExact(s0, s1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    family.case('IntrSegment2Segment2.find.antiparallelDeviation', (io) => {
        const s0 = segment(io, 2);
        const s1 = segment(io, 2);
        outSeg2Seg2FI(io, new IntrSegment2Segment2FI().find(s0, s1));
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrSegment2Segment2.h FIQuery '
            + '(segment1Parameter = overlap - t has the wrong sign for antiparallel '
            + 'collinear segments); issue #458'
    });

    // ========================= IntrSphere3Triangle3 ======================

    family.case('IntrSphere3Triangle3.find', (io) => {
        const s = sphere(io, 3);
        const sphereVelocity = io.vec(3);
        const t = triangle(io, 3);
        const triangleVelocity = io.vec(3);
        const r = new IntrSphere3Triangle3FI().find(s, sphereVelocity, t,
            triangleVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }, { exact: true });

    family.case('IntrSphere3Triangle3.find.aimed', (io) => {
        const t = triangle(io, 3);
        const radius = io.real();
        const center = io.vec(3);
        const triangleVelocity = io.vec(3);
        const sphereVelocity = io.vec(3);
        const s = Hypersphere.fromCenterRadius(center, radius);
        const r = new IntrSphere3Triangle3FI().find(s, sphereVelocity, t,
            triangleVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }, { exact: true });

    // ======================== IntrAlignedBox2Circle2 =====================

    family.case('IntrAlignedBox2Circle2.test', (io) => {
        const b = box(io, 2);
        const c = sphere(io, 2);
        const r = new IntrAlignedBox2Circle2TI().test(b, c);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrAlignedBox2Circle2.find', (io) => {
        const b = box(io, 2);
        const boxVelocity = io.vec(2);
        const c = sphere(io, 2);
        const circleVelocity = io.vec(2);
        const r = new IntrAlignedBox2Circle2FI().find(b, boxVelocity, c,
            circleVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }, { exact: true });

    family.case('IntrAlignedBox2Circle2.find.aimed', (io) => {
        const b = box(io, 2);
        const radius = io.real();
        const center = io.vec(2);
        const boxVelocity = io.vec(2);
        const circleVelocity = io.vec(2);
        const c = Hypersphere.fromCenterRadius(center, radius);
        const r = new IntrAlignedBox2Circle2FI().find(b, boxVelocity, c,
            circleVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }, { exact: true });

    // ======================= IntrAlignedBox3Cylinder3 ====================

    family.case('IntrAlignedBox3Cylinder3.test', (io) => {
        const b = box(io, 3);
        const cyl = cylinder(io);
        const r = new IntrAlignedBox3Cylinder3TI().test(b, cyl);
        io.outBool(r.intersect);
    }, { exact: true });

    // The aligned-box query is a thin wrapper around the canonical-box query,
    // so it inherits the (U1,-D) sign typo of DoQueryNoZeros and the port's
    // correction of it.
    family.case('IntrAlignedBox3Cylinder3.test.edgeTypoDeviation', (io) => {
        const b = box(io, 3);
        const cyl = cylinder(io);
        const r = new IntrAlignedBox3Cylinder3TI().test(b, cyl);
        io.outBool(r.intersect);
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrCanonicalBox3Cylinder3.h '
            + 'DoQueryNoZeros (U1,-D) sign typo; issue #197'
    });

    // Upstream asserts that the cylinder is finite, so both sides throw.
    family.case('IntrAlignedBox3Cylinder3.test.infinite', (io) => {
        const b = box(io, 3);
        const cyl = cylinder(io);
        const r = new IntrAlignedBox3Cylinder3TI().test(b, cyl);
        io.outBool(r.intersect);
    }, { exact: true });

    // ===================== IntrConvexPolygonHyperplane ===================

    function polygon(io: OracleIO, n: number, dim: number): Vector[] {
        const vertices: Vector[] = [];
        for (let i = 0; i < n; ++i) { vertices.push(io.vec(dim)); }
        return vertices;
    }

    function hyperplane(io: OracleIO, dim: number): Hyperplane {
        const normal = io.vec(dim);
        const constant = io.real();
        return Hyperplane.fromNormalConstant(normal, constant);
    }

    function polygonHyperplaneTest(io: OracleIO, dim: number): void {
        const n = io.integer();
        const poly = polygon(io, n, dim);
        io.integer();  // the recorded 'kind' of the hyperplane construction
        const plane = hyperplane(io, dim);
        const r = new IntrConvexPolygonHyperplaneTI().test(poly, plane);
        io.outBool(r.intersect);
        io.outInt(r.configuration);
    }

    function polygonHyperplaneFind(io: OracleIO, dim: number): void {
        const n = io.integer();
        const poly = polygon(io, n, dim);
        io.integer();  // the recorded 'kind' of the hyperplane construction
        const plane = hyperplane(io, dim);
        const r = new IntrConvexPolygonHyperplaneFI().find(poly, plane);
        io.outInt(r.configuration);
        io.outInt(r.intersection.length);
        for (const p of r.intersection) { io.outVec(p); }
        io.outInt(r.positivePolygon.length);
        for (const p of r.positivePolygon) { io.outVec(p); }
        io.outInt(r.negativePolygon.length);
        for (const p of r.negativePolygon) { io.outVec(p); }
    }

    family.case('IntrConvexPolygonHyperplane.test.2d',
        (io) => { polygonHyperplaneTest(io, 2); }, { exact: true });

    family.case('IntrConvexPolygonHyperplane.find.2d',
        (io) => { polygonHyperplaneFind(io, 2); }, { exact: true });

    family.case('IntrConvexPolygonHyperplane.test.3d',
        (io) => { polygonHyperplaneTest(io, 3); }, { exact: true });

    family.case('IntrConvexPolygonHyperplane.find.3d',
        (io) => { polygonHyperplaneFind(io, 3); }, { exact: true });

    // ======================== IntrAlignedBox3Sphere3 =====================

    family.case('IntrAlignedBox3Sphere3.test', (io) => {
        const b = box(io, 3);
        const s = sphere(io, 3);
        const r = new IntrAlignedBox3Sphere3TI().test(b, s);
        io.outBool(r.intersect);
    }, { exact: true });

    function box3Sphere3Find(io: OracleIO): void {
        const b = box(io, 3);
        const boxVelocity = io.vec(3);
        const center = io.vec(3);
        const radius = io.real();
        const sphereVelocity = io.vec(3);
        const s = Hypersphere.fromCenterRadius(center, radius);
        const r = new IntrAlignedBox3Sphere3FI().find(b, boxVelocity, s,
            sphereVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }

    family.case('IntrAlignedBox3Sphere3.find', box3Sphere3Find, { exact: true });

    family.case('IntrAlignedBox3Sphere3.find.aimed', box3Sphere3Find, { exact: true });

    family.case('IntrAlignedBox3Sphere3.find.probeDeviation', box3Sphere3Find, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrAlignedBox3Sphere3.h '
            + 'DoQueryRayRoundedFace / DoQuery (first-wins probe selection over the '
            + 'pieces of the Minkowski sum); issues #458 and #465'
    });

    // ======================== IntrTriangle3Cylinder3 =====================

    family.case('IntrTriangle3Cylinder3.test', (io) => {
        const t = triangle(io, 3);
        const cyl = cylinder(io);
        const r = new IntrTriangle3Cylinder3TI().test(t, cyl);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrTriangle3Cylinder3.test.degeneratePolygonDeviation', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const radius = io.real();
        const height = io.real();
        const t = triangle(io, 3);
        const cyl = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin, direction), radius, height);
        const r = new IntrTriangle3Cylinder3TI().test(t, cyl);
        io.outBool(r.intersect);
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrTriangle3Cylinder3.h '
            + 'DiskOverlapsPolygon (a degenerate projected polygon is reported as '
            + 'containing the origin); issues #206 and #458'
    });

    family.case('IntrTriangle3Cylinder3.test.infiniteDeviation', (io) => {
        const t = triangle(io, 3);
        const cyl = cylinder(io);
        const r = new IntrTriangle3Cylinder3TI().test(t, cyl);
        io.outBool(r.intersect);
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md cylinder headers, missing IsFinite '
            + 'guard; issues #197, #206 and #255'
    });

    // ========================= IntrEllipse2Ellipse2 ======================

    function outEllipseFI(io: OracleIO, r: IntrEllipse2Ellipse2FIResult,
        tol?: number): void {
        io.outBool(r.intersect);
        io.outInt(r.numPoints === intrEllipse2Ellipse2InfinitePoints
            ? -1 : r.numPoints);
        for (let i = 0; i < 4; ++i) { io.outBool(r.isTransverse[i]); }
        for (let i = 0; i < 4; ++i) { io.outVec(r.points[i], tol); }
    }

    family.case('IntrEllipse2Ellipse2.test', (io) => {
        const e0 = ellipse(io);
        const e1 = ellipse(io);
        io.outInt(new IntrEllipse2Ellipse2TI().test(e0, e1));
    }, { exact: true });

    family.case('IntrEllipse2Ellipse2.test.poleDeviation', (io) => {
        const e0 = ellipse(io);
        const e1 = ellipse(io);
        io.outInt(new IntrEllipse2Ellipse2TI().test(e0, e1));
    }, {
        exact: true,
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrEllipse2Ellipse2.h TIQuery '
            + '(the c_i = 0 terms of f(s) are dropped, losing the two pole critical '
            + 'points); issue #458'
    });

    // The quartic H is built with exact arithmetic (the software fused
    // multiply-add of Functions.ts), but its roots come from
    // RootsPolynomial::SolveQuartic, whose resolvent cubic calls pow, cos and
    // atan2; a 1 ulp libm difference there is amplified near a double root.
    // 1978 of the 2000 deep-run records are bit-identical on every coordinate
    // and the worst scaled error over the remaining ones is 4.3e-10, so the
    // tolerance is 1e-8. The discrete outputs (numPoints, isTransverse) are
    // compared exactly and agreed on every record.
    const quarticTol = 1e-8;

    family.case('IntrEllipse2Ellipse2.find', (io) => {
        const e0 = ellipse(io);
        const e1 = ellipse(io);
        const earlyExit = io.boolean();
        outEllipseFI(io, new IntrEllipse2Ellipse2FI().find(e0, e1, earlyExit),
            quarticTol);
    }, { tol: quarticTol });

    family.case('IntrEllipse2Ellipse2.findStandardForm', (io) => {
        const C0 = io.vec(2);
        const m0 = io.reals(3);
        const C1 = io.vec(2);
        const m1 = io.reals(3);
        const earlyExit = io.boolean();
        const M0 = Matrix.fromArray(2, 2, [m0[0], m0[1], m0[1], m0[2]]);
        const M1 = Matrix.fromArray(2, 2, [m1[0], m1[1], m1[1], m1[2]]);
        outEllipseFI(io, new IntrEllipse2Ellipse2FI()
            .findStandardForm(C0, M0, C1, M1, earlyExit), quarticTol);
    }, { tol: quarticTol });

    family.case('IntrEllipse2Ellipse2.getStandardForm', (io) => {
        const e = ellipse(io);
        const query = new IntrEllipse2Ellipse2FI();
        const { C, M } = query.getStandardForm(e);
        io.outVec(C);
        io.outReal(M.get(0, 0));
        io.outReal(M.get(0, 1));
        io.outReal(M.get(1, 0));
        io.outReal(M.get(1, 1));
        const box0 = query.computeAlignedBox(e);
        const box1 = query.computeAlignedBox(C, M);
        io.outVec(box0.min);
        io.outVec(box0.max);
        io.outVec(box1.min);
        io.outVec(box1.max);
    }, { exact: true });

    family.case('IntrEllipse2Ellipse2.find.divisorDeviation', (io) => {
        const e0 = ellipse(io);
        const e1 = ellipse(io);
        outEllipseFI(io, new IntrEllipse2Ellipse2FI().find(e0, e1, false));
    }, {
        deviation: 'docs/UPSTREAM-FINDINGS.md IntrEllipse2Ellipse2.h CaseE4NotZero '
            + '(divisor == 0 writes both symmetric points to the same slot); issue #250'
    });

    family.finish();
});
