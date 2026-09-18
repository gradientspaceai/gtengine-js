// Replays oracle/cpp/cases/v19-distance.cpp (verify group 19, distance).
// Keep the two files in the same order.
//
// Every generator in the C++ file records the same number of doubles in every
// mode, so these replays read the inputs without knowing which mode produced
// them.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { CanonicalBox } from '../../src/CanonicalBox.js';
import { DistAlignedBoxAlignedBox } from '../../src/DistAlignedBoxAlignedBox.js';
import { DistCircle2Circle2 } from '../../src/DistCircle2Circle2.js';
import { DistLine2AlignedBox2, distLine2AlignedBox2DoQuery,
    type DistLine2AlignedBox2Result } from '../../src/DistLine2AlignedBox2.js';
import { DistLine2Circle2 } from '../../src/DistLine2Circle2.js';
import { DistLine2Triangle2 } from '../../src/DistLine2Triangle2.js';
import { DistLine3CanonicalBox3 } from '../../src/DistLine3CanonicalBox3.js';
import { DistLineLine } from '../../src/DistLineLine.js';
import { DistLineRay } from '../../src/DistLineRay.js';
import { DistLineSegment } from '../../src/DistLineSegment.js';
import { DistPoint2Circle2 } from '../../src/DistPoint2Circle2.js';
import { DistPointCanonicalBox } from '../../src/DistPointCanonicalBox.js';
import { DistPointLine } from '../../src/DistPointLine.js';
import { DistPointRay } from '../../src/DistPointRay.js';
import { DistPointRectangle } from '../../src/DistPointRectangle.js';
import { DistPointSegment } from '../../src/DistPointSegment.js';
import { DistPointTriangle } from '../../src/DistPointTriangle.js';
import { DistRayRay } from '../../src/DistRayRay.js';
import { DistRaySegment } from '../../src/DistRaySegment.js';
import { DistSegmentSegment } from '../../src/DistSegmentSegment.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { Ray } from '../../src/Ray.js';
import { Rectangle } from '../../src/Rectangle.js';
import { Segment } from '../../src/Segment.js';
import { Triangle } from '../../src/Triangle.js';
import { Vector, add } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The C++ AlignedBoxInputs records the minimum corner and the size; the
// maximum corner is derived.
function alignedBox(io: OracleIO, n: number): AlignedBox {
    const lo = io.vec(n);
    const size = io.vec(n);
    return AlignedBox.fromMinMax(lo, add(lo, size));
}

describe('oracle: v19-distance', () => {
    const family = new OracleFamily('v19-distance');

    // ---- DistPointLine.h ----
    for (const n of [2, 3]) {
        family.case(`DistPointLine.compute.${n}d`, (io) => {
            const point = io.vec(n);
            const line = Line.fromOriginDirection(io.vec(n), io.vec(n));
            const r = new DistPointLine().compute(point, line);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outReal(r.parameter);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistPointRay.h ----
    for (const n of [2, 3]) {
        family.case(`DistPointRay.compute.${n}d`, (io) => {
            const point = io.vec(n);
            const ray = Ray.fromOriginDirection(io.vec(n), io.vec(n));
            const r = new DistPointRay().compute(point, ray);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outReal(r.parameter);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistPointSegment.h ----
    for (const n of [2, 3]) {
        family.case(`DistPointSegment.compute.${n}d`, (io) => {
            const point = io.vec(n);
            const segment = Segment.fromEndpoints(io.vec(n), io.vec(n));
            const r = new DistPointSegment().compute(point, segment);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outReal(r.parameter);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistPointCanonicalBox.h ----
    for (const n of [2, 3]) {
        family.case(`DistPointCanonicalBox.compute.${n}d`, (io) => {
            const point = io.vec(n);
            const box = CanonicalBox.fromExtent(io.vec(n));
            const r = new DistPointCanonicalBox().compute(point, box);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistPointRectangle.h ----
    for (const n of [2, 3]) {
        family.case(`DistPointRectangle.compute.${n}d`, (io) => {
            const point = io.vec(n);
            const center = io.vec(n);
            const axis = [io.vec(n), io.vec(n)];
            const rectangle =
                Rectangle.fromCenterAxisExtent(center, axis, io.vec(2));
            const r = new DistPointRectangle().compute(point, rectangle);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outReal(r.cartesian[0]);
            io.outReal(r.cartesian[1]);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistPoint2Circle2.h ----
    family.case('DistPoint2Circle2.compute', (io) => {
        const center = io.vec(2);
        const point = io.vec(2);
        const circle = Hypersphere.fromCenterRadius(center, io.real());
        const r = new DistPoint2Circle2().compute(point, circle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
        io.outBool(r.equidistant);
    }, { exact: true });

    // ---- DistLineLine.h ----
    for (const n of [2, 3]) {
        family.case(`DistLineLine.compute.${n}d`, (io) => {
            const line0 = Line.fromOriginDirection(io.vec(n), io.vec(n));
            const line1 = Line.fromOriginDirection(io.vec(n), io.vec(n));
            const r = new DistLineLine().compute(line0, line1);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistLineRay.h ----
    for (const n of [2, 3]) {
        family.case(`DistLineRay.compute.${n}d`, (io) => {
            const line = Line.fromOriginDirection(io.vec(n), io.vec(n));
            const ray = Ray.fromOriginDirection(io.vec(n), io.vec(n));
            const r = new DistLineRay().compute(line, ray);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistLineSegment.h ----
    for (const n of [2, 3]) {
        family.case(`DistLineSegment.compute.${n}d`, (io) => {
            const line = Line.fromOriginDirection(io.vec(n), io.vec(n));
            const segment = Segment.fromEndpoints(io.vec(n), io.vec(n));
            const r = new DistLineSegment().compute(line, segment);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistRayRay.h ----
    for (const n of [2, 3]) {
        family.case(`DistRayRay.compute.${n}d`, (io) => {
            const ray0 = Ray.fromOriginDirection(io.vec(n), io.vec(n));
            const ray1 = Ray.fromOriginDirection(io.vec(n), io.vec(n));
            const r = new DistRayRay().compute(ray0, ray1);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outReal(r.parameter[0]);
            io.outReal(r.parameter[1]);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistRaySegment.h ----
    function raySegment(io: OracleIO, n: number): void {
        const ray = Ray.fromOriginDirection(io.vec(n), io.vec(n));
        const segment = Segment.fromEndpoints(io.vec(n), io.vec(n));
        const r = new DistRaySegment().compute(ray, segment);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    for (const n of [2, 3]) {
        family.case(`DistRaySegment.compute.${n}d`,
            (io) => { raySegment(io, n); }, { exact: true });
    }

    // The port clamps the ray parameter to s0 >= 0; upstream does not. Every
    // record of this case is a configuration where upstream reports s0 < 0.
    family.case('DistRaySegment.compute.3d.deviation',
        (io) => { raySegment(io, 3); },
        { deviation: 'UPSTREAM-FINDINGS DistRaySegment.h, issue #126' });

    // ---- DistSegmentSegment.h ----
    function segmentSegment(io: OracleIO, n: number,
        which: 'compute' | 'computeEndpoints' | 'computeRobust'
            | 'computeRobustEndpoints'): void {
        const p0 = io.vec(n);
        const p1 = io.vec(n);
        const q0 = io.vec(n);
        const q1 = io.vec(n);
        const query = new DistSegmentSegment();
        const r = which === 'compute'
            ? query.compute(Segment.fromEndpoints(p0, p1),
                Segment.fromEndpoints(q0, q1))
            : which === 'computeEndpoints'
                ? query.computeEndpoints(p0, p1, q0, q1)
                : which === 'computeRobust'
                    ? query.computeRobust(Segment.fromEndpoints(p0, p1),
                        Segment.fromEndpoints(q0, q1))
                    : query.computeRobustEndpoints(p0, p1, q0, q1);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    for (const n of [2, 3]) {
        family.case(`DistSegmentSegment.compute.${n}d`,
            (io) => { segmentSegment(io, n, 'compute'); }, { exact: true });
        family.case(`DistSegmentSegment.computeRobust.${n}d`,
            (io) => { segmentSegment(io, n, 'computeRobust'); },
            { exact: true });
    }

    family.case('DistSegmentSegment.computeEndpoints.3d',
        (io) => { segmentSegment(io, 3, 'computeEndpoints'); },
        { exact: true });
    family.case('DistSegmentSegment.computeRobustEndpoints.3d',
        (io) => { segmentSegment(io, 3, 'computeRobustEndpoints'); },
        { exact: true });

    // Upstream replaces an out-of-range endpoint ratio by 1/2; the port clamps
    // it to the nearest endpoint of [0,1].
    family.case('DistSegmentSegment.computeRobust.3d.deviation',
        (io) => { segmentSegment(io, 3, 'computeRobust'); },
        { deviation: 'UPSTREAM-FINDINGS DistSegmentSegment.h, issue #418' });

    // ---- DistPointTriangle.h ----
    function pointTriangle(io: OracleIO, n: number,
        conjugateGradient: boolean): void {
        const point = io.vec(n);
        const triangle = Triangle.fromVertices(io.vec(n), io.vec(n), io.vec(n));
        const query = new DistPointTriangle();
        const r = conjugateGradient
            ? query.useConjugateGradient(point, triangle)
            : query.compute(point, triangle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReals(r.barycentric);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    for (const n of [2, 3]) {
        family.case(`DistPointTriangle.compute.${n}d`,
            (io) => { pointTriangle(io, n, false); }, { exact: true });
        family.case(`DistPointTriangle.useConjugateGradient.${n}d`,
            (io) => { pointTriangle(io, n, true); }, { exact: true });
    }

    // ---- DistAlignedBoxAlignedBox.h ----
    for (const n of [2, 3]) {
        family.case(`DistAlignedBoxAlignedBox.compute.${n}d`, (io) => {
            const box0 = alignedBox(io, n);
            const box1 = alignedBox(io, n);
            const r = new DistAlignedBoxAlignedBox().compute(box0, box1);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outVec(r.closest[0].min);
            io.outVec(r.closest[0].max);
            io.outVec(r.closest[1].min);
            io.outVec(r.closest[1].max);
        }, { exact: true });
    }

    // ---- DistCircle2Circle2.h ----
    family.case('DistCircle2Circle2.compute', (io) => {
        const center0 = io.vec(2);
        const radius0 = io.real();
        const center1 = io.vec(2);
        const radius1 = io.real();
        const r = new DistCircle2Circle2().compute(
            Hypersphere.fromCenterRadius(center0, radius0),
            Hypersphere.fromCenterRadius(center1, radius1));
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outInt(r.numClosestPairs);
        io.outBool(r.concentric);
        io.outBool(r.cocircular);
        for (let j = 0; j < r.numClosestPairs; ++j) {
            io.outVec(r.closest[j][0]);
            io.outVec(r.closest[j][1]);
        }
    }, { exact: true });

    // ---- DistLine2AlignedBox2.h ----
    family.case('DistLine2AlignedBox2.compute', (io) => {
        const box = alignedBox(io, 2);
        const line = Line.fromOriginDirection(io.vec(2), io.vec(2));
        const r = new DistLine2AlignedBox2().compute(line, box);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    // The exported DoQuery (upstream's protected static, granted to the
    // line/oriented-box query by friendship). It works in the box frame, sets
    // only parameter and closest, and reflects the origin and direction in
    // place without undoing the reflection, so those are outputs too.
    family.case('DistLine2AlignedBox2.doQuery', (io) => {
        const extent = io.vec(2);
        const origin = io.vec(2);
        const direction = io.vec(2);
        const result: DistLine2AlignedBox2Result = {
            distance: 0, sqrDistance: 0, parameter: 0,
            closest: [new Vector(2), new Vector(2)]
        };
        distLine2AlignedBox2DoQuery(origin, direction, extent, result);
        io.outReal(result.parameter);
        io.outVec(result.closest[0]);
        io.outVec(result.closest[1]);
        io.outVec(origin);
        io.outVec(direction);
    }, { exact: true });

    // ---- DistLine2Circle2.h ----
    family.case('DistLine2Circle2.compute', (io) => {
        const center = io.vec(2);
        const circle = Hypersphere.fromCenterRadius(center, io.real());
        const line = Line.fromOriginDirection(io.vec(2), io.vec(2));
        const r = new DistLine2Circle2().compute(line, circle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outInt(r.numClosestPairs);
        for (let j = 0; j < r.numClosestPairs; ++j) {
            io.outReal(r.parameter[j]);
            io.outVec(r.closest[j][0]);
            io.outVec(r.closest[j][1]);
        }
    }, { exact: true });

    // ---- DistLine2Triangle2.h ----
    function line2Triangle2(io: OracleIO): void {
        const line = Line.fromOriginDirection(io.vec(2), io.vec(2));
        const triangle =
            Triangle.fromVertices(io.vec(2), io.vec(2), io.vec(2));
        const r = new DistLine2Triangle2().compute(line, triangle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outReals(r.barycentric);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    family.case('DistLine2Triangle2.compute', line2Triangle2, { exact: true });

    // Upstream divides by an edge denominator that has rounded to exactly
    // zero and returns NaN; the port uses the normal-component difference.
    family.case('DistLine2Triangle2.compute.deviation', line2Triangle2,
        { deviation: 'UPSTREAM-FINDINGS DistLine2Triangle2.h, issue #441' });

    // ---- DistLine3CanonicalBox3.h ----
    function line3CanonicalBox3(io: OracleIO): void {
        const box = CanonicalBox.fromExtent(io.vec(3));
        const line = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const r = new DistLine3CanonicalBox3().compute(line, box);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    family.case('DistLine3CanonicalBox3.compute', line3CanonicalBox3,
        { exact: true });

    // Upstream lets the accumulated squared distance go slightly negative and
    // returns NaN; the port clamps it at zero.
    family.case('DistLine3CanonicalBox3.compute.deviation', line3CanonicalBox3,
        { deviation: 'UPSTREAM-FINDINGS DistLine3CanonicalBox3.h, issue #421' });

    family.finish();
});
