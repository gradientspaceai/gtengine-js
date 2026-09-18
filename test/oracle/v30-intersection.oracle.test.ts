// Replays oracle/cpp/cases/v30-intersection.cpp (verify group 30,
// intersection). Keep the two files in the same order.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { Halfspace } from '../../src/Halfspace.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Ray } from '../../src/Ray.js';
import { Rectangle } from '../../src/Rectangle.js';
import { Segment } from '../../src/Segment.js';
import { Triangle } from '../../src/Triangle.js';
import { Vector } from '../../src/Vector.js';
import {
    IntrAlignedBox2AlignedBox2FI, IntrAlignedBox2AlignedBox2TI
} from '../../src/IntrAlignedBox2AlignedBox2.js';
import { IntrAlignedBox2OrientedBox2TI } from '../../src/IntrAlignedBox2OrientedBox2.js';
import {
    IntrAlignedBox3AlignedBox3FI, IntrAlignedBox3AlignedBox3TI
} from '../../src/IntrAlignedBox3AlignedBox3.js';
import { IntrAlignedBox3OrientedBox3TI } from '../../src/IntrAlignedBox3OrientedBox3.js';
import { IntrCircle2Circle2FI, IntrCircle2Circle2TI } from '../../src/IntrCircle2Circle2.js';
import { IntrHalfspace2Polygon2FI } from '../../src/IntrHalfspace2Polygon2.js';
import { IntrHalfspace3OrientedBox3TI } from '../../src/IntrHalfspace3OrientedBox3.js';
import {
    IntrHalfspace3Segment3FI, IntrHalfspace3Segment3TI
} from '../../src/IntrHalfspace3Segment3.js';
import { IntrHalfspace3Sphere3TI } from '../../src/IntrHalfspace3Sphere3.js';
import {
    IntrHalfspace3Triangle3FI, IntrHalfspace3Triangle3TI
} from '../../src/IntrHalfspace3Triangle3.js';
import { IntrIntervalsFI, IntrIntervalsTI } from '../../src/IntrIntervals.js';
import {
    IntrLine2AlignedBox2FI, IntrLine2AlignedBox2TI,
    defaultIntrLine2AlignedBox2FIResult, defaultIntrLine2AlignedBox2TIResult,
    intrLine2AlignedBox2FIDoQuery, intrLine2AlignedBox2TIDoQuery
} from '../../src/IntrLine2AlignedBox2.js';
import { IntrLine2Line2FI, IntrLine2Line2TI } from '../../src/IntrLine2Line2.js';
import {
    IntrLine2Triangle2FI, IntrLine2Triangle2TI,
    defaultIntrLine2Triangle2FIResult, intrLine2Triangle2FIDoQuery
} from '../../src/IntrLine2Triangle2.js';
import {
    IntrLine3AlignedBox3FI, IntrLine3AlignedBox3TI,
    defaultIntrLine3AlignedBox3FIResult, defaultIntrLine3AlignedBox3TIResult,
    intrLine3AlignedBox3FIDoQuery, intrLine3AlignedBox3TIDoQuery
} from '../../src/IntrLine3AlignedBox3.js';
import { IntrLine3Rectangle3FI, IntrLine3Rectangle3TI } from '../../src/IntrLine3Rectangle3.js';
import {
    IntrLine3Sphere3FI, IntrLine3Sphere3TI,
    defaultIntrLine3Sphere3FIResult, intrLine3Sphere3FIDoQuery
} from '../../src/IntrLine3Sphere3.js';
import { IntrLine3Triangle3FI, IntrLine3Triangle3TI } from '../../src/IntrLine3Triangle3.js';
import {
    IntrOrientedBox2OrientedBox2FI, IntrOrientedBox2OrientedBox2TI
} from '../../src/IntrOrientedBox2OrientedBox2.js';
import { IntrOrientedBox3OrientedBox3TI } from '../../src/IntrOrientedBox3OrientedBox3.js';
import { IntrRay3Triangle3FI, IntrRay3Triangle3TI } from '../../src/IntrRay3Triangle3.js';
import { IntrSegment3Triangle3FI, IntrSegment3Triangle3TI } from '../../src/IntrSegment3Triangle3.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The C++ sentinel std::numeric_limits<int32_t>::max().
const INT32_MAX = 2147483647;

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

function halfspace(io: OracleIO, n: number): Halfspace {
    const normal = io.vec(n);
    const constant = io.real();
    return Halfspace.fromNormalConstant(normal, constant);
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

function circle(io: OracleIO): Hypersphere {
    const center = io.vec(2);
    const radius = io.real();
    return Hypersphere.fromCenterRadius(center, radius);
}

function sphere(io: OracleIO): Hypersphere {
    const center = io.vec(3);
    const radius = io.real();
    return Hypersphere.fromCenterRadius(center, radius);
}

function rectangle3(io: OracleIO): Rectangle {
    const center = io.vec(3);
    const axis0 = io.vec(3);
    const axis1 = io.vec(3);
    // The C++ generator records a full orthonormal frame; the rectangle uses
    // only the first two axes.
    io.vec(3);
    const extent = io.vec(2);
    return Rectangle.fromCenterAxisExtent(center, [axis0, axis1], extent);
}

function polygon2(io: OracleIO): Vector[] {
    const n = io.integer();
    const vertices: Vector[] = [];
    for (let i = 0; i < n; ++i) { vertices.push(io.vec(2)); }
    return vertices;
}

describe('oracle: v30-intersection', () => {
    const family = new OracleFamily('v30-intersection');

    // ============================ IntrIntervals ==========================

    family.case('IntrIntervals.test.finiteFinite', (io) => {
        const interval0 = io.reals(2);
        const interval1 = io.reals(2);
        const r = new IntrIntervalsTI().test(interval0, interval1);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrIntervals.test.finiteSemiInfinite', (io) => {
        const finite = io.reals(2);
        const a = io.real();
        const isPositiveInfinite = io.boolean();
        const r = new IntrIntervalsTI().testFiniteSemiInfinite(finite, a, isPositiveInfinite);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrIntervals.test.semiInfiniteSemiInfinite', (io) => {
        const a0 = io.real();
        const p0 = io.boolean();
        const a1 = io.real();
        const p1 = io.boolean();
        const r = new IntrIntervalsTI().testSemiInfiniteSemiInfinite(a0, p0, a1, p1);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrIntervals.test.dynamic', (io) => {
        const maxTime = io.real();
        const interval0 = io.reals(2);
        const speed0 = io.real();
        const interval1 = io.reals(2);
        const speed1 = io.real();
        const r = new IntrIntervalsTI().testDynamic(maxTime, interval0, speed0, interval1, speed1);
        io.outBool(r.intersect);
        io.outReal(r.firstTime);
        io.outReal(r.lastTime);
    }, { exact: true });

    family.case('IntrIntervals.find.finiteFinite', (io) => {
        const interval0 = io.reals(2);
        const interval1 = io.reals(2);
        const r = new IntrIntervalsFI().find(interval0, interval1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outInt(r.type);
        io.outReal(r.overlap[0]);
        io.outReal(r.overlap[1]);
    }, { exact: true });

    family.case('IntrIntervals.find.finiteSemiInfinite', (io) => {
        const finite = io.reals(2);
        const a = io.real();
        const isPositiveInfinite = io.boolean();
        const r = new IntrIntervalsFI().findFiniteSemiInfinite(finite, a, isPositiveInfinite);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outInt(r.type);
        io.outReal(r.overlap[0]);
        io.outReal(r.overlap[1]);
    }, { exact: true });

    family.case('IntrIntervals.find.semiInfiniteSemiInfinite', (io) => {
        const a0 = io.real();
        const p0 = io.boolean();
        const a1 = io.real();
        const p1 = io.boolean();
        const r = new IntrIntervalsFI().findSemiInfiniteSemiInfinite(a0, p0, a1, p1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outInt(r.type);
        io.outReal(r.overlap[0]);
        io.outReal(r.overlap[1]);
    }, { exact: true });

    family.case('IntrIntervals.findDynamic.soundBranches', (io) => {
        const maxTime = io.real();
        const interval0 = io.reals(2);
        const speed0 = io.real();
        const interval1 = io.reals(2);
        const speed1 = io.real();
        const r = new IntrIntervalsFI().findDynamic(maxTime, interval0, speed0, interval1, speed1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outInt(r.type);
        io.outReal(r.firstTime);
        io.outReal(r.lastTime);
        io.outReal(r.overlap[0]);
        io.outReal(r.overlap[1]);
    }, { exact: true });

    // The port reports the contact point interval0[1] + firstTime * speed0
    // where upstream reports interval0[0] + firstTime * speed0. See
    // docs/UPSTREAM-FINDINGS.md (IntrIntervals.h, dynamic FIQuery) and
    // gradientspaceai/gtengine-js#62.
    family.case('IntrIntervals.findDynamic.leftApproachDeviation', (io) => {
        const maxTime = io.real();
        const interval0 = io.reals(2);
        const speed0 = io.real();
        const interval1 = io.reals(2);
        const speed1 = io.real();
        const r = new IntrIntervalsFI().findDynamic(maxTime, interval0, speed0, interval1, speed1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outInt(r.type);
        io.outReal(r.firstTime);
        io.outReal(r.lastTime);
        io.outReal(r.overlap[0]);
        io.outReal(r.overlap[1]);
    }, { exact: true, deviation: 'UPSTREAM-FINDINGS IntrIntervals.h dynamic FIQuery (#62)' });

    // ===================== aligned box / aligned box =====================

    family.case('IntrAlignedBox2AlignedBox2.test', (io) => {
        const box0 = box(io, 2);
        const box1 = box(io, 2);
        const r = new IntrAlignedBox2AlignedBox2TI().test(box0, box1);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrAlignedBox2AlignedBox2.find', (io) => {
        const box0 = box(io, 2);
        const box1 = box(io, 2);
        const r = new IntrAlignedBox2AlignedBox2FI().find(box0, box1);
        io.outBool(r.intersect);
        io.outVec(r.box.min);
        io.outVec(r.box.max);
    }, { exact: true });

    family.case('IntrAlignedBox3AlignedBox3.test', (io) => {
        const box0 = box(io, 3);
        const box1 = box(io, 3);
        const r = new IntrAlignedBox3AlignedBox3TI().test(box0, box1);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrAlignedBox3AlignedBox3.find', (io) => {
        const box0 = box(io, 3);
        const box1 = box(io, 3);
        const r = new IntrAlignedBox3AlignedBox3FI().find(box0, box1);
        io.outBool(r.intersect);
        io.outVec(r.box.min);
        io.outVec(r.box.max);
    }, { exact: true });

    // ==================== aligned box / oriented box =====================

    family.case('IntrAlignedBox2OrientedBox2.test', (io) => {
        const box0 = box(io, 2);
        const box1 = obox(io, 2);
        const r = new IntrAlignedBox2OrientedBox2TI().test(box0, box1);
        io.outBool(r.intersect);
        io.outInt(r.separating);
    }, { exact: true });

    family.case('IntrAlignedBox3OrientedBox3.test', (io) => {
        const box0 = box(io, 3);
        const box1 = obox(io, 3);
        const epsilon = io.real();
        const r = new IntrAlignedBox3OrientedBox3TI().test(box0, box1, epsilon);
        io.outBool(r.intersect);
        io.outInt(r.separating[0]);
        io.outInt(r.separating[1]);
    }, { exact: true });

    // =================== oriented box / oriented box =====================

    family.case('IntrOrientedBox2OrientedBox2.test', (io) => {
        const box0 = obox(io, 2);
        const box1 = obox(io, 2);
        const r = new IntrOrientedBox2OrientedBox2TI().test(box0, box1);
        io.outBool(r.intersect);
        io.outInt(r.separating);
    }, { exact: true });

    family.case('IntrOrientedBox2OrientedBox2.find', (io) => {
        const box0 = obox(io, 2);
        const box1 = obox(io, 2);
        const r = new IntrOrientedBox2OrientedBox2FI().find(box0, box1);
        io.outBool(r.intersect);
        io.outInt(r.polygon.length);
        for (const p of r.polygon) { io.outVec(p); }
    }, { exact: true });

    family.case('IntrOrientedBox3OrientedBox3.test', (io) => {
        const box0 = obox(io, 3);
        const box1 = obox(io, 3);
        const epsilon = io.real();
        const r = new IntrOrientedBox3OrientedBox3TI().test(box0, box1, epsilon);
        io.outBool(r.intersect);
        io.outInt(r.separating[0]);
        io.outInt(r.separating[1]);
    }, { exact: true });

    // ========================== circle / circle ==========================

    family.case('IntrCircle2Circle2.test', (io) => {
        const circle0 = circle(io);
        const circle1 = circle(io);
        const r = new IntrCircle2Circle2TI().test(circle0, circle1);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrCircle2Circle2.find', (io) => {
        const circle0 = circle(io);
        const circle1 = circle(io);
        const r = new IntrCircle2Circle2FI().find(circle0, circle1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.numIntersections === INT32_MAX) {
            io.outVec(r.circle.center);
            io.outReal(r.circle.radius);
        } else if (r.numIntersections > 0) {
            io.outVec(r.point[0]);
            if (r.numIntersections === 2) { io.outVec(r.point[1]); }
        }
    }, { exact: true });

    // ========================= halfspace queries =========================

    family.case('IntrHalfspace2Polygon2.find', (io) => {
        const polygon = polygon2(io);
        const space = halfspace(io, 2);
        const r = new IntrHalfspace2Polygon2FI().find(space, polygon);
        io.outBool(r.intersect);
        io.outInt(r.polygon.length);
        for (const p of r.polygon) { io.outVec(p); }
    }, { exact: true });

    family.case('IntrHalfspace3OrientedBox3.test', (io) => {
        const box1 = obox(io, 3);
        const space = halfspace(io, 3);
        const r = new IntrHalfspace3OrientedBox3TI().test(space, box1);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrHalfspace3Sphere3.test', (io) => {
        const s = sphere(io);
        const space = halfspace(io, 3);
        const r = new IntrHalfspace3Sphere3TI().test(space, s);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrHalfspace3Segment3.test', (io) => {
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const space = halfspace(io, 3);
        const r = new IntrHalfspace3Segment3TI().test(space, Segment.fromEndpoints(p0, p1));
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrHalfspace3Segment3.find', (io) => {
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const space = halfspace(io, 3);
        const r = new IntrHalfspace3Segment3FI().find(space, Segment.fromEndpoints(p0, p1));
        io.outBool(r.intersect);
        io.outInt(r.numPoints);
        for (let i = 0; i < r.numPoints; ++i) { io.outVec(r.point[i]); }
    }, { exact: true });

    family.case('IntrHalfspace3Triangle3.test', (io) => {
        const tri = triangle(io, 3);
        const space = halfspace(io, 3);
        const r = new IntrHalfspace3Triangle3TI().test(space, tri);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrHalfspace3Triangle3.find', (io) => {
        const tri = triangle(io, 3);
        const space = halfspace(io, 3);
        const r = new IntrHalfspace3Triangle3FI().find(space, tri);
        io.outBool(r.intersect);
        io.outInt(r.numPoints);
        for (let i = 0; i < r.numPoints; ++i) { io.outVec(r.point[i]); }
    }, { exact: true });

    // ============================ line queries ===========================

    family.case('IntrLine2AlignedBox2.test', (io) => {
        const l = line(io, 2);
        const b = box(io, 2);
        const r = new IntrLine2AlignedBox2TI().test(l, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine2AlignedBox2.find', (io) => {
        const l = line(io, 2);
        const b = box(io, 2);
        const r = new IntrLine2AlignedBox2FI().find(l, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    family.case('IntrLine2AlignedBox2.doQuery.ti', (io) => {
        const lineOrigin = io.vec(2);
        const lineDirection = io.vec(2);
        const boxExtent = io.vec(2);
        const r = defaultIntrLine2AlignedBox2TIResult();
        intrLine2AlignedBox2TIDoQuery(lineOrigin, lineDirection, boxExtent, r);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine2AlignedBox2.doQuery.fi', (io) => {
        const lineOrigin = io.vec(2);
        const lineDirection = io.vec(2);
        const boxExtent = io.vec(2);
        const r = defaultIntrLine2AlignedBox2FIResult();
        intrLine2AlignedBox2FIDoQuery(lineOrigin, lineDirection, boxExtent, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrLine3AlignedBox3.test', (io) => {
        const l = line(io, 3);
        const b = box(io, 3);
        const r = new IntrLine3AlignedBox3TI().test(l, b);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3AlignedBox3.find', (io) => {
        const l = line(io, 3);
        const b = box(io, 3);
        const r = new IntrLine3AlignedBox3FI().find(l, b);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrLine3AlignedBox3.doQuery.ti', (io) => {
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const boxExtent = io.vec(3);
        const r = defaultIntrLine3AlignedBox3TIResult();
        intrLine3AlignedBox3TIDoQuery(lineOrigin, lineDirection, boxExtent, r);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3AlignedBox3.doQuery.fi', (io) => {
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const boxExtent = io.vec(3);
        const r = defaultIntrLine3AlignedBox3FIResult();
        intrLine3AlignedBox3FIDoQuery(lineOrigin, lineDirection, boxExtent, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrLine2Line2.test', (io) => {
        const line0 = line(io, 2);
        const line1 = line(io, 2);
        const r = new IntrLine2Line2TI().test(line0, line1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
    }, { exact: true });

    family.case('IntrLine2Line2.find', (io) => {
        const line0 = line(io, 2);
        const line1 = line(io, 2);
        const r = new IntrLine2Line2FI().find(line0, line1);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.line0Parameter[0]);
        io.outReal(r.line0Parameter[1]);
        io.outReal(r.line1Parameter[0]);
        io.outReal(r.line1Parameter[1]);
        io.outVec(r.point);
    }, { exact: true });

    family.case('IntrLine2Triangle2.test', (io) => {
        const l = line(io, 2);
        const tri = triangle(io, 2);
        const r = new IntrLine2Triangle2TI().test(l, tri);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine2Triangle2.find', (io) => {
        const l = line(io, 2);
        const tri = triangle(io, 2);
        const r = new IntrLine2Triangle2FI().find(l, tri);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrLine2Triangle2.doQuery.fi', (io) => {
        const origin = io.vec(2);
        const direction = io.vec(2);
        const tri = triangle(io, 2);
        const r = defaultIntrLine2Triangle2FIResult();
        intrLine2Triangle2FIDoQuery(origin, direction, tri, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrLine3Rectangle3.test', (io) => {
        const l = line(io, 3);
        const rect = rectangle3(io);
        const r = new IntrLine3Rectangle3TI().test(l, rect);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3Rectangle3.find', (io) => {
        const l = line(io, 3);
        const rect = rectangle3(io);
        const r = new IntrLine3Rectangle3FI().find(l, rect);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.rectCoord[0]);
            io.outReal(r.rectCoord[1]);
            io.outReal(r.rectCoord[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    family.case('IntrLine3Sphere3.test', (io) => {
        const l = line(io, 3);
        const s = sphere(io);
        const r = new IntrLine3Sphere3TI().test(l, s);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3Sphere3.doQuery.fi', (io) => {
        const lineOrigin = io.vec(3);
        const lineDirection = io.vec(3);
        const s = sphere(io);
        const r = defaultIntrLine3Sphere3FIResult();
        intrLine3Sphere3FIDoQuery(lineOrigin, lineDirection, s, r);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
    }, { exact: true });

    family.case('IntrLine3Sphere3.find', (io) => {
        const l = line(io, 3);
        const s = sphere(io);
        const r = new IntrLine3Sphere3FI().find(l, s);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.case('IntrLine3Triangle3.test', (io) => {
        const l = line(io, 3);
        const tri = triangle(io, 3);
        const r = new IntrLine3Triangle3TI().test(l, tri);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3Triangle3.find', (io) => {
        const l = line(io, 3);
        const tri = triangle(io, 3);
        const r = new IntrLine3Triangle3FI().find(l, tri);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.triangleBary[0]);
            io.outReal(r.triangleBary[1]);
            io.outReal(r.triangleBary[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    family.case('IntrRay3Triangle3.test', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const tri = triangle(io, 3);
        const r = new IntrRay3Triangle3TI().test(
            Ray.fromOriginDirection(origin, direction), tri);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay3Triangle3.find', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const tri = triangle(io, 3);
        const r = new IntrRay3Triangle3FI().find(
            Ray.fromOriginDirection(origin, direction), tri);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.triangleBary[0]);
            io.outReal(r.triangleBary[1]);
            io.outReal(r.triangleBary[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    family.case('IntrSegment3Triangle3.test', (io) => {
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const tri = triangle(io, 3);
        const r = new IntrSegment3Triangle3TI().test(Segment.fromEndpoints(p0, p1), tri);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment3Triangle3.find', (io) => {
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const tri = triangle(io, 3);
        const r = new IntrSegment3Triangle3FI().find(Segment.fromEndpoints(p0, p1), tri);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.triangleBary[0]);
            io.outReal(r.triangleBary[1]);
            io.outReal(r.triangleBary[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    // ================== constructed hits and tangencies ==================
    // The C++ side aims the line, ray or segment at a point computed from the
    // primitive on a small lattice, so the boundary comparisons are reached
    // with exact equality. See the matching comment in the .cpp.

    // triangle, direction, origin (the C++ helper records them in that order).
    function lineAtTriangle(io: OracleIO): { tri: Triangle; origin: Vector; direction: Vector } {
        const tri = triangle(io, 3);
        const direction = io.vec(3);
        const origin = io.vec(3);
        return { tri, origin, direction };
    }

    family.case('IntrLine3Triangle3.test.throughPoint', (io) => {
        const { tri, origin, direction } = lineAtTriangle(io);
        const r = new IntrLine3Triangle3TI().test(
            Line.fromOriginDirection(origin, direction), tri);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3Triangle3.find.throughPoint', (io) => {
        const { tri, origin, direction } = lineAtTriangle(io);
        const r = new IntrLine3Triangle3FI().find(
            Line.fromOriginDirection(origin, direction), tri);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.triangleBary[0]);
            io.outReal(r.triangleBary[1]);
            io.outReal(r.triangleBary[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    family.case('IntrRay3Triangle3.test.throughPoint', (io) => {
        const { tri, origin, direction } = lineAtTriangle(io);
        const r = new IntrRay3Triangle3TI().test(
            Ray.fromOriginDirection(origin, direction), tri);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrRay3Triangle3.find.throughPoint', (io) => {
        const { tri, origin, direction } = lineAtTriangle(io);
        const r = new IntrRay3Triangle3FI().find(
            Ray.fromOriginDirection(origin, direction), tri);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.triangleBary[0]);
            io.outReal(r.triangleBary[1]);
            io.outReal(r.triangleBary[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    // triangle, direction (recorded but used only to place the endpoints),
    // p0, p1.
    function segmentAtTriangle(io: OracleIO): { tri: Triangle; p0: Vector; p1: Vector } {
        const tri = triangle(io, 3);
        io.vec(3);
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        return { tri, p0, p1 };
    }

    family.case('IntrSegment3Triangle3.test.throughPoint', (io) => {
        const { tri, p0, p1 } = segmentAtTriangle(io);
        const r = new IntrSegment3Triangle3TI().test(Segment.fromEndpoints(p0, p1), tri);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrSegment3Triangle3.find.throughPoint', (io) => {
        const { tri, p0, p1 } = segmentAtTriangle(io);
        const r = new IntrSegment3Triangle3FI().find(Segment.fromEndpoints(p0, p1), tri);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.triangleBary[0]);
            io.outReal(r.triangleBary[1]);
            io.outReal(r.triangleBary[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    function lineAtRectangle(io: OracleIO): { rect: Rectangle; origin: Vector; direction: Vector } {
        const rect = rectangle3(io);
        const direction = io.vec(3);
        const origin = io.vec(3);
        return { rect, origin, direction };
    }

    family.case('IntrLine3Rectangle3.test.throughPoint', (io) => {
        const { rect, origin, direction } = lineAtRectangle(io);
        const r = new IntrLine3Rectangle3TI().test(
            Line.fromOriginDirection(origin, direction), rect);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3Rectangle3.find.throughPoint', (io) => {
        const { rect, origin, direction } = lineAtRectangle(io);
        const r = new IntrLine3Rectangle3FI().find(
            Line.fromOriginDirection(origin, direction), rect);
        io.outBool(r.intersect);
        if (r.intersect) {
            io.outReal(r.parameter);
            io.outReal(r.rectCoord[0]);
            io.outReal(r.rectCoord[1]);
            io.outReal(r.rectCoord[2]);
            io.outVec(r.point);
        }
    }, { exact: true });

    // sphere center, radius, direction, origin.
    function tangentLine(io: OracleIO): { s: Hypersphere; origin: Vector; direction: Vector } {
        const s = sphere(io);
        const direction = io.vec(3);
        const origin = io.vec(3);
        return { s, origin, direction };
    }

    family.case('IntrLine3Sphere3.test.tangent', (io) => {
        const { s, origin, direction } = tangentLine(io);
        const r = new IntrLine3Sphere3TI().test(
            Line.fromOriginDirection(origin, direction), s);
        io.outBool(r.intersect);
    }, { exact: true });

    family.case('IntrLine3Sphere3.find.tangent', (io) => {
        const { s, origin, direction } = tangentLine(io);
        const r = new IntrLine3Sphere3FI().find(
            Line.fromOriginDirection(origin, direction), s);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        if (r.intersect) {
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }, { exact: true });

    family.finish();
});
