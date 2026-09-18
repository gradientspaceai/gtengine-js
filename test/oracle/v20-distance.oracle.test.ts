// Replays oracle/cpp/cases/v20-distance.cpp (verify group 20, distance).
// Keep the two files in the same order.
//
// Every generator in the C++ file records the same number of doubles in every
// mode, so these replays read the inputs without knowing which mode produced
// them.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { Arc2 } from '../../src/Arc2.js';
import { CanonicalBox } from '../../src/CanonicalBox.js';
import { Circle3 } from '../../src/Circle3.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { DistLine2Arc2 } from '../../src/DistLine2Arc2.js';
import { DistLine2OrientedBox2 } from '../../src/DistLine2OrientedBox2.js';
import { DistLine3AlignedBox3 } from '../../src/DistLine3AlignedBox3.js';
import { DistLine3Circle3 } from '../../src/DistLine3Circle3.js';
import { DistLine3OrientedBox3 } from '../../src/DistLine3OrientedBox3.js';
import { DistLine3Rectangle3 } from '../../src/DistLine3Rectangle3.js';
import { DistLine3Triangle3 } from '../../src/DistLine3Triangle3.js';
import { DistPoint2Arc2 } from '../../src/DistPoint2Arc2.js';
import { DistPoint3Circle3 } from '../../src/DistPoint3Circle3.js';
import { DistPoint3Cylinder3 } from '../../src/DistPoint3Cylinder3.js';
import { DistPoint3Frustum3 } from '../../src/DistPoint3Frustum3.js';
import { DistPointAlignedBox } from '../../src/DistPointAlignedBox.js';
import { DistPointOrientedBox } from '../../src/DistPointOrientedBox.js';
import { DistRay2Circle2 } from '../../src/DistRay2Circle2.js';
import { DistRay2Triangle2 } from '../../src/DistRay2Triangle2.js';
import { DistRay3CanonicalBox3 } from '../../src/DistRay3CanonicalBox3.js';
import { DistSegment2Circle2 } from '../../src/DistSegment2Circle2.js';
import { DistSegment2Triangle2 } from '../../src/DistSegment2Triangle2.js';
import { DistSegment3CanonicalBox3 } from '../../src/DistSegment3CanonicalBox3.js';
import { Frustum3 } from '../../src/Frustum3.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { OrientedBox } from '../../src/OrientedBox.js';
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

// The C++ OrientedBoxInputs records center, then the N axes, then the extent.
function orientedBox(io: OracleIO, n: number): OrientedBox {
    const center = io.vec(n);
    const axis: Vector[] = [];
    for (let i = 0; i < n; ++i) { axis.push(io.vec(n)); }
    return OrientedBox.fromCenterAxisExtent(center, axis, io.vec(n));
}

// The C++ ArcInputs records center, radius, end[0], end[1].
function arc2(io: OracleIO): Arc2 {
    const center = io.vec(2);
    const radius = io.real();
    const end0 = io.vec(2);
    return Arc2.fromCenterRadiusEnds(center, radius, end0, io.vec(2));
}

// distance, sqrDistance, numClosestPairs and then, for each pair, the
// parameter and the two closest points. The shared line/ray/segment versus
// circle result in 2D.
function emitCircle2Pairs(io: OracleIO, r: {
    distance: number, sqrDistance: number, numClosestPairs: number,
    parameter: readonly number[], closest: readonly (readonly Vector[])[]
}): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outInt(r.numClosestPairs);
    for (let j = 0; j < r.numClosestPairs; ++j) {
        io.outReal(r.parameter[j]);
        io.outVec(r.closest[j][0]);
        io.outVec(r.closest[j][1]);
    }
}

function emitTriangle2(io: OracleIO, r: {
    distance: number, sqrDistance: number, parameter: number,
    barycentric: readonly number[], closest: readonly Vector[]
}): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter);
    io.outReal(r.barycentric[0]);
    io.outReal(r.barycentric[1]);
    io.outReal(r.barycentric[2]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

function emitOneParameter(io: OracleIO, r: {
    distance: number, sqrDistance: number, parameter: number,
    closest: readonly Vector[]
}): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

describe('oracle: v20-distance', () => {
    const family = new OracleFamily('v20-distance');

    // ---- DistPointAlignedBox.h ----
    for (const n of [2, 3]) {
        family.case(`DistPointAlignedBox.compute.${n}d`, (io) => {
            const point = io.vec(n);
            const box = alignedBox(io, n);
            const r = new DistPointAlignedBox().compute(point, box);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistPointOrientedBox.h ----
    for (const n of [2, 3]) {
        family.case(`DistPointOrientedBox.compute.${n}d`, (io) => {
            const point = io.vec(n);
            const box = orientedBox(io, n);
            const r = new DistPointOrientedBox().compute(point, box);
            io.outReal(r.distance);
            io.outReal(r.sqrDistance);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistPoint2Arc2.h ----
    family.case('DistPoint2Arc2.compute', (io) => {
        const a = arc2(io);
        const point = io.vec(2);
        const r = new DistPoint2Arc2().compute(point, a);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
        io.outBool(r.equidistant);
    }, { exact: true });

    // ---- DistLine2Arc2.h ----
    family.case('DistLine2Arc2.compute', (io) => {
        const a = arc2(io);
        const line = Line.fromOriginDirection(io.vec(2), io.vec(2));
        const r = new DistLine2Arc2().compute(line, a);
        emitCircle2Pairs(io, r);
    }, { exact: true });

    // ---- DistLine2OrientedBox2.h ----
    family.case('DistLine2OrientedBox2.compute', (io) => {
        const box = orientedBox(io, 2);
        const line = Line.fromOriginDirection(io.vec(2), io.vec(2));
        const r = new DistLine2OrientedBox2().compute(line, box);
        emitOneParameter(io, r);
    }, { exact: true });

    // ---- DistLine3AlignedBox3.h ----
    family.case('DistLine3AlignedBox3.compute', (io) => {
        const box = alignedBox(io, 3);
        const line = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const r = new DistLine3AlignedBox3().compute(line, box);
        emitOneParameter(io, r);
    }, { exact: true });

    // ---- DistLine3OrientedBox3.h ----
    // closest[0] is covered by the '.deviation' case below: upstream
    // transforms the world-space line point a second time as if it were in
    // the box frame.
    family.case('DistLine3OrientedBox3.compute', (io) => {
        const box = orientedBox(io, 3);
        const line = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const r = new DistLine3OrientedBox3().compute(line, box);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('DistLine3OrientedBox3.compute.deviation', (io) => {
        const box = orientedBox(io, 3);
        const line = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const r = new DistLine3OrientedBox3().compute(line, box);
        io.outVec(r.closest[0]);
    }, { exact: true, deviation: 'UPSTREAM-FINDINGS DistLine3OrientedBox3.h: '
        + 'closest[0] is written in world space and then transformed again '
        + 'as if it were in the box frame' });

    // ---- DistLine3Rectangle3.h ----
    family.case('DistLine3Rectangle3.compute', (io) => {
        const center = io.vec(3);
        const a0 = io.vec(3);
        const a1 = io.vec(3);
        // The C++ Axes3 records all three frame vectors; the rectangle uses
        // the first two.
        io.vec(3);
        const extent = io.vec(2);
        const rectangle = Rectangle.fromCenterAxisExtent(center, [a0, a1],
            extent);
        const line = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const r = new DistLine3Rectangle3().compute(line, rectangle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outReal(r.cartesian[0]);
        io.outReal(r.cartesian[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    // ---- DistLine3Triangle3.h ----
    family.case('DistLine3Triangle3.compute', (io) => {
        const v0 = io.vec(3);
        const v1 = io.vec(3);
        const triangle = Triangle.fromVertices(v0, v1, io.vec(3));
        const line = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const r = new DistLine3Triangle3().compute(line, triangle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outReal(r.barycentric[0]);
        io.outReal(r.barycentric[1]);
        io.outReal(r.barycentric[2]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    // ---- DistPoint3Circle3.h ----
    family.case('DistPoint3Circle3.compute', (io) => {
        const center = io.vec(3);
        const normal = io.vec(3);
        const circle = Circle3.fromCenterNormalRadius(center, normal,
            io.real());
        const point = io.vec(3);
        const r = new DistPoint3Circle3().compute(point, circle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
        io.outBool(r.equidistant);
    }, { exact: true });

    // ---- DistPoint3Cylinder3.h ----
    family.case('DistPoint3Cylinder3.compute', (io) => {
        const axis = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const radius = io.real();
        const cylinder = Cylinder3.fromAxisRadiusHeight(axis, radius,
            io.real());
        const point = io.vec(3);
        const r = new DistPoint3Cylinder3().compute(point, cylinder);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('DistPoint3Cylinder3.compute.deviation', (io) => {
        const axis = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const radius = io.real();
        // The finite height is recorded but unused; the cylinder is built
        // with the MakeInfiniteCylinder sentinel that follows it.
        io.real();
        const cylinder = Cylinder3.fromAxisRadiusHeight(axis, radius,
            io.real());
        const point = io.vec(3);
        const r = new DistPoint3Cylinder3().compute(point, cylinder);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true, deviation: 'UPSTREAM-FINDINGS DistPoint3Cylinder3.h: '
        + 'the infinite-cylinder test is height == DBL_MAX, but the Cylinder3 '
        + 'sentinel is height = -1, so upstream throws on it' });

    // ---- DistPoint3Frustum3.h ----
    function frustum3(io: OracleIO): Frustum3 {
        const origin = io.vec(3);
        const rVector = io.vec(3);
        const uVector = io.vec(3);
        const dVector = io.vec(3);
        const dMin = io.real();
        const dMax = io.real();
        const uBound = io.real();
        return Frustum3.fromParameters(origin, dVector, uVector, rVector,
            dMin, dMax, uBound, io.real());
    }

    family.case('DistPoint3Frustum3.compute', (io) => {
        const frustum = frustum3(io);
        const point = io.vec(3);
        const r = new DistPoint3Frustum3().compute(point, frustum);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('DistPoint3Frustum3.compute.deviation', (io) => {
        const frustum = frustum3(io);
        const point = io.vec(3);
        const r = new DistPoint3Frustum3().compute(point, frustum);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true, deviation: 'UPSTREAM-FINDINGS DistPoint3Frustum3.h: '
        + 'two of the ten far-edge assignments do not clamp the free '
        + 'coordinate, so the closest point lies outside the frustum' });

    // ---- DistRay2Circle2.h ----
    family.case('DistRay2Circle2.compute', (io) => {
        const center = io.vec(2);
        const circle = Hypersphere.fromCenterRadius(center, io.real());
        const ray = Ray.fromOriginDirection(io.vec(2), io.vec(2));
        const r = new DistRay2Circle2().compute(ray, circle);
        emitCircle2Pairs(io, r);
    }, { exact: true });

    // ---- DistSegment2Circle2.h ----
    family.case('DistSegment2Circle2.compute', (io) => {
        const center = io.vec(2);
        const circle = Hypersphere.fromCenterRadius(center, io.real());
        const segment = Segment.fromEndpoints(io.vec(2), io.vec(2));
        const r = new DistSegment2Circle2().compute(segment, circle);
        emitCircle2Pairs(io, r);
    }, { exact: true });

    // ---- DistRay2Triangle2.h ----
    family.case('DistRay2Triangle2.compute', (io) => {
        const v0 = io.vec(2);
        const v1 = io.vec(2);
        const triangle = Triangle.fromVertices(v0, v1, io.vec(2));
        const ray = Ray.fromOriginDirection(io.vec(2), io.vec(2));
        const r = new DistRay2Triangle2().compute(ray, triangle);
        emitTriangle2(io, r);
    }, { exact: true });

    // ---- DistSegment2Triangle2.h ----
    family.case('DistSegment2Triangle2.compute', (io) => {
        const v0 = io.vec(2);
        const v1 = io.vec(2);
        const triangle = Triangle.fromVertices(v0, v1, io.vec(2));
        const segment = Segment.fromEndpoints(io.vec(2), io.vec(2));
        const r = new DistSegment2Triangle2().compute(segment, triangle);
        emitTriangle2(io, r);
    }, { exact: true });

    // ---- DistRay3CanonicalBox3.h ----
    family.case('DistRay3CanonicalBox3.compute', (io) => {
        const box = CanonicalBox.fromExtent(io.vec(3));
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        const r = new DistRay3CanonicalBox3().compute(ray, box);
        emitOneParameter(io, r);
    }, { exact: true });

    // ---- DistSegment3CanonicalBox3.h ----
    family.case('DistSegment3CanonicalBox3.compute', (io) => {
        const box = CanonicalBox.fromExtent(io.vec(3));
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        const r = new DistSegment3CanonicalBox3().compute(segment, box);
        emitOneParameter(io, r);
    }, { exact: true });

    // ---- DistLine3Circle3.h ----
    function line3Circle3(io: OracleIO): void {
        const center = io.vec(3);
        const normal = io.vec(3);
        const circle = Circle3.fromCenterNormalRadius(center, normal,
            io.real());
        const line = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const r = new DistLine3Circle3().compute(line, circle);
        io.outInt(r.numClosestPairs);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outBool(r.equidistant);
        for (let j = 0; j < r.numClosestPairs; ++j) {
            io.outVec(r.linearClosest[j]);
            io.outVec(r.circularClosest[j]);
        }
    }

    // The closed-form branches (PDFSection411, 412 and 421) use only
    // + - * / and sqrt.
    family.case('DistLine3Circle3.compute.axis', line3Circle3, { exact: true });

    family.case('DistLine3Circle3.compute.axis.deviation', line3Circle3,
        { exact: true, deviation: 'UPSTREAM-FINDINGS DistLine3Circle3.h: '
            + 'Finalize normalizes a possibly-zero projection and then '
            + 'reports the circle center as the closest circle point '
            + '(issue #421)' });

    // PDFSection422. This is the one case of the family that is not compared
    // bit for bit, for two reasons that are both intrinsic to the path:
    //
    // - it is iterative (RootsBisection1 on Phi, up to 4096 halvings), and
    // - the port recovers the critical line parameter as
    //   t = G(tau) - Dot(M,D)/Dot(M,M) instead of upstream's t = tau + s.
    //   That substitution is the documented fix for issue #421 item 3: for a
    //   line near the axis of the circle, s grows like 1/|NxM|^2 while every
    //   bisection bracket has width r*|NxM|/Dot(M,M), so tau approaches -s to
    //   far below ulp(s) and upstream's sum cancels every significant digit.
    //   The defect is one of conditioning, not of a branch, so there is no
    //   exact predicate that separates the inputs upstream gets right from
    //   the ones it gets wrong, and the substitution cannot be guarded. Where
    //   upstream is sound the two expressions agree to a few ulps, which is
    //   what this case measures; the '.nearPerpendicular' case below covers
    //   the regime in which upstream's answer is actually wrong.
    //
    // The generator keeps only configurations on which upstream's own
    // reported distance agrees with an independent reference (the minimum of
    // the point-to-line distance over the circle) to 1e-9 relative.
    // Substituting upstream's t = tau + s into the port makes every record of
    // this case bit-identical, which is how the residual was root-caused.
    family.case('DistLine3Circle3.compute', line3Circle3, { tol: 1e-12 });

    family.case('DistLine3Circle3.compute.tauHat', line3Circle3,
        { exact: true, deviation: 'UPSTREAM-FINDINGS DistLine3Circle3.h: '
            + 'tauHat is missing the division by a2 (issue #247)' });

    family.case('DistLine3Circle3.compute.nearPerpendicular', line3Circle3,
        { exact: true, deviation: 'UPSTREAM-FINDINGS DistLine3Circle3.h: '
            + 't = tau + s cancels away the significant digits of the '
            + 'critical parameter for a near-perpendicular line '
            + '(issue #421 item 3)' });

    family.finish();
});
