// Replays oracle/cpp/cases/v21-distance.cpp (verify group 21, distance).
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
import { ConvexPolyhedron3 } from '../../src/ConvexPolyhedron3.js';
import { DistPlane3CanonicalBox3 } from '../../src/DistPlane3CanonicalBox3.js';
import { DistPoint2Parallelogram2 } from '../../src/DistPoint2Parallelogram2.js';
import { DistPoint3ConvexPolyhedron3 }
    from '../../src/DistPoint3ConvexPolyhedron3.js';
import { DistPointHyperellipsoid }
    from '../../src/DistPointHyperellipsoid.js';
import { DistPointHyperplane } from '../../src/DistPointHyperplane.js';
import { DistRay2AlignedBox2 } from '../../src/DistRay2AlignedBox2.js';
import { DistRay2Arc2 } from '../../src/DistRay2Arc2.js';
import { DistRay2OrientedBox2 } from '../../src/DistRay2OrientedBox2.js';
import { DistRay3AlignedBox3 } from '../../src/DistRay3AlignedBox3.js';
import { DistRay3Circle3 } from '../../src/DistRay3Circle3.js';
import { DistRay3OrientedBox3 } from '../../src/DistRay3OrientedBox3.js';
import { DistRay3Rectangle3 } from '../../src/DistRay3Rectangle3.js';
import { DistRay3Triangle3 } from '../../src/DistRay3Triangle3.js';
import { DistSegment2AlignedBox2 } from '../../src/DistSegment2AlignedBox2.js';
import { DistSegment2Arc2 } from '../../src/DistSegment2Arc2.js';
import { DistSegment2OrientedBox2 }
    from '../../src/DistSegment2OrientedBox2.js';
import { DistSegment3AlignedBox3 } from '../../src/DistSegment3AlignedBox3.js';
import { DistSegment3Circle3 } from '../../src/DistSegment3Circle3.js';
import { DistSegment3OrientedBox3 }
    from '../../src/DistSegment3OrientedBox3.js';
import { DistSegment3Rectangle3 } from '../../src/DistSegment3Rectangle3.js';
import { DistSegment3Triangle3 } from '../../src/DistSegment3Triangle3.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { Matrix, multiplyATB } from '../../src/Matrix.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Parallelogram2 } from '../../src/Parallelogram2.js';
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

// The C++ oriented-box inputs: center, the N axes, the extents.
function orientedBox(io: OracleIO, n: number): OrientedBox {
    const center = io.vec(n);
    const axis: Vector[] = [];
    for (let i = 0; i < n; ++i) { axis.push(io.vec(n)); }
    return OrientedBox.fromCenterAxisExtent(center, axis, io.vec(n));
}

function emitOneParameter(io: OracleIO, r: { distance: number;
    sqrDistance: number; parameter: number; closest: [Vector, Vector] }): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

// The 3D oriented-box queries inherit the closest[0] defect of
// DistLine3OrientedBox3 (UPSTREAM-FINDINGS, issue #421), which the port fixes;
// the ordinary cases therefore do not compare closest[0]. See the '.deviation'
// cases and oracle/reports/v21-distance.md.
function emitNoLinearClosest(io: OracleIO, r: { distance: number;
    sqrDistance: number; parameter: number; closest: [Vector, Vector] }): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter);
    io.outVec(r.closest[1]);
}

function emitTwoClosest(io: OracleIO, r: { distance: number;
    sqrDistance: number; closest: [Vector, Vector] }): void {
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

describe('oracle: v21-distance', () => {
    const family = new OracleFamily('v21-distance');

    // ---- DistRay2AlignedBox2.h, DistRay3AlignedBox3.h ----
    for (const n of [2, 3]) {
        const name = n === 2 ? 'DistRay2AlignedBox2' : 'DistRay3AlignedBox3';
        family.case(`${name}.compute`, (io) => {
            const box = alignedBox(io, n);
            const ray = Ray.fromOriginDirection(io.vec(n), io.vec(n));
            const r = n === 2
                ? new DistRay2AlignedBox2().compute(ray, box)
                : new DistRay3AlignedBox3().compute(ray, box);
            emitOneParameter(io, r);
        }, { exact: true });
    }

    // ---- DistSegment2AlignedBox2.h, DistSegment3AlignedBox3.h ----
    for (const n of [2, 3]) {
        const name = n === 2
            ? 'DistSegment2AlignedBox2' : 'DistSegment3AlignedBox3';
        family.case(`${name}.compute`, (io) => {
            const box = alignedBox(io, n);
            const segment = Segment.fromEndpoints(io.vec(n), io.vec(n));
            const r = n === 2
                ? new DistSegment2AlignedBox2().compute(segment, box)
                : new DistSegment3AlignedBox3().compute(segment, box);
            emitOneParameter(io, r);
        }, { exact: true });
    }

    // ---- DistRay2OrientedBox2.h, DistSegment2OrientedBox2.h ----
    family.case('DistRay2OrientedBox2.compute', (io) => {
        const box = orientedBox(io, 2);
        const ray = Ray.fromOriginDirection(io.vec(2), io.vec(2));
        emitOneParameter(io, new DistRay2OrientedBox2().compute(ray, box));
    }, { exact: true });

    family.case('DistSegment2OrientedBox2.compute', (io) => {
        const box = orientedBox(io, 2);
        const segment = Segment.fromEndpoints(io.vec(2), io.vec(2));
        emitOneParameter(io,
            new DistSegment2OrientedBox2().compute(segment, box));
    }, { exact: true });

    // ---- DistRay3OrientedBox3.h, DistSegment3OrientedBox3.h ----
    family.case('DistRay3OrientedBox3.compute', (io) => {
        const box = orientedBox(io, 3);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        emitNoLinearClosest(io, new DistRay3OrientedBox3().compute(ray, box));
    }, { exact: true });

    family.case('DistRay3OrientedBox3.compute.deviation', (io) => {
        const box = orientedBox(io, 3);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        io.outVec(new DistRay3OrientedBox3().compute(ray, box).closest[0]);
    }, { exact: true, deviation: 'UPSTREAM-FINDINGS DistLine3OrientedBox3.h, '
        + 'issue #421: the world-space line point is written into closest[0] '
        + 'before the loop that maps the box-frame closest points back to the '
        + 'world, so upstream transforms it a second time.' });

    family.case('DistSegment3OrientedBox3.compute', (io) => {
        const box = orientedBox(io, 3);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitNoLinearClosest(io,
            new DistSegment3OrientedBox3().compute(segment, box));
    }, { exact: true });

    family.case('DistSegment3OrientedBox3.compute.deviation', (io) => {
        const box = orientedBox(io, 3);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        io.outVec(new DistSegment3OrientedBox3()
            .compute(segment, box).closest[0]);
    }, { exact: true, deviation: 'UPSTREAM-FINDINGS DistLine3OrientedBox3.h, '
        + 'issue #421: the world-space line point is written into closest[0] '
        + 'before the loop that maps the box-frame closest points back to the '
        + 'world, so upstream transforms it a second time.' });

    // ---- DistRay3Rectangle3.h, DistSegment3Rectangle3.h ----
    // The C++ RectangleInputs records the center, three orthonormal axes (the
    // third is unused) and the two extents.
    function rectangle(io: OracleIO): Rectangle {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3)];
        io.vec(3);
        return Rectangle.fromCenterAxisExtent(center, axis, io.vec(2));
    }

    function emitCartesian(io: OracleIO, r: { distance: number;
        sqrDistance: number; parameter: number; cartesian: number[];
        closest: [Vector, Vector] }): void {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outReal(r.cartesian[0]);
        io.outReal(r.cartesian[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    function emitBarycentric(io: OracleIO, r: { distance: number;
        sqrDistance: number; parameter: number; barycentric: number[];
        closest: [Vector, Vector] }): void {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outReal(r.barycentric[0]);
        io.outReal(r.barycentric[1]);
        io.outReal(r.barycentric[2]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    family.case('DistRay3Rectangle3.compute', (io) => {
        const rect = rectangle(io);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        emitCartesian(io, new DistRay3Rectangle3().compute(ray, rect));
    }, { exact: true });

    family.case('DistSegment3Rectangle3.compute', (io) => {
        const rect = rectangle(io);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitCartesian(io,
            new DistSegment3Rectangle3().compute(segment, rect));
    }, { exact: true });

    family.case('DistRay3Triangle3.compute', (io) => {
        const triangle = Triangle.fromVertices(io.vec(3), io.vec(3),
            io.vec(3));
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        emitBarycentric(io, new DistRay3Triangle3().compute(ray, triangle));
    }, { exact: true });

    family.case('DistSegment3Triangle3.compute', (io) => {
        const triangle = Triangle.fromVertices(io.vec(3), io.vec(3),
            io.vec(3));
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitBarycentric(io,
            new DistSegment3Triangle3().compute(segment, triangle));
    }, { exact: true });

    // ---- DistPointHyperplane.h ----
    for (const n of [2, 3, 4]) {
        family.case(`DistPointHyperplane.compute.${n}d`, (io) => {
            const normal = io.vec(n);
            const plane = Hyperplane.fromNormalOrigin(normal, io.vec(n));
            const r = new DistPointHyperplane().compute(io.vec(n), plane);
            io.outReal(r.distance);
            io.outReal(r.signedDistance);
            io.outVec(r.closest[0]);
            io.outVec(r.closest[1]);
        }, { exact: true });
    }

    // ---- DistPlane3CanonicalBox3.h ----
    family.case('DistPlane3CanonicalBox3.compute', (io) => {
        const normal = io.vec(3);
        const plane = Hyperplane.fromNormalOrigin(normal, io.vec(3));
        const box = CanonicalBox.fromExtent(io.vec(3));
        emitTwoClosest(io, new DistPlane3CanonicalBox3().compute(plane, box));
    }, { exact: true });

    // ---- DistPoint2Parallelogram2.h ----
    family.case('DistPoint2Parallelogram2.compute', (io) => {
        const center = io.vec(2);
        const axis = [io.vec(2), io.vec(2)];
        const pgm = Parallelogram2.fromCenterAxis(center, axis);
        emitTwoClosest(io,
            new DistPoint2Parallelogram2().compute(io.vec(2), pgm));
    }, { exact: true });

    family.case('DistPoint2Parallelogram2.compute.degenerate', (io) => {
        const center = io.vec(2);
        const axis = [io.vec(2), io.vec(2)];
        const pgm = Parallelogram2.fromCenterAxis(center, axis);
        emitTwoClosest(io,
            new DistPoint2Parallelogram2().compute(io.vec(2), pgm));
    }, { exact: true });

    family.case('DistPoint2Parallelogram2.getMinimizer', (io) => {
        const B = new Matrix(2, 2);
        B.setCol(0, io.vec(2));
        B.setCol(1, io.vec(2));
        const A = multiplyATB(B, B);
        io.outVec(new DistPoint2Parallelogram2().getMinimizer(A, io.vec(2)));
    }, { exact: true });

    // ---- DistPointHyperellipsoid.h ----
    for (const n of [2, 3]) {
        family.case(`DistPointHyperellipsoid.compute.${n}d`, (io) => {
            const center = io.vec(n);
            const axis: Vector[] = [];
            for (let i = 0; i < n; ++i) { axis.push(io.vec(n)); }
            const hyperellipsoid =
                Hyperellipsoid.fromCenterAxisExtent(center, axis, io.vec(n));
            emitTwoClosest(io, new DistPointHyperellipsoid()
                .compute(io.vec(n), hyperellipsoid));
        }, { exact: true });

        family.case(`DistPointHyperellipsoid.computeAxisAligned.${n}d`, (io) => {
            const extent = io.vec(n);
            emitTwoClosest(io, new DistPointHyperellipsoid()
                .computeAxisAligned(io.vec(n), extent));
        }, { exact: true });
    }

    // ---- DistPoint3ConvexPolyhedron3.h ----
    function indices(io: OracleIO): number[] {
        const count = io.integer();
        const out: number[] = [];
        for (let i = 0; i < count; ++i) { out.push(io.integer()); }
        return out;
    }

    function vertices(io: OracleIO, count: number): Vector[] {
        const out: Vector[] = [];
        for (let i = 0; i < count; ++i) { out.push(io.vec(3)); }
        return out;
    }

    function emitConvexPolyhedron(io: OracleIO, r: {
        queryIsSuccessful: boolean; numLCPIterations: number; distance: number;
        sqrDistance: number; closest: [Vector, Vector] }): void {
        io.outBool(r.queryIsSuccessful);
        io.outInt(r.numLCPIterations);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    family.case('DistPoint3ConvexPolyhedron3.compute.tetrahedron', (io) => {
        const polyhedron = new ConvexPolyhedron3(vertices(io, 4), indices(io),
            true, true);
        emitConvexPolyhedron(io,
            new DistPoint3ConvexPolyhedron3().compute(io.vec(3), polyhedron));
    }, { exact: true });

    family.case('DistPoint3ConvexPolyhedron3.compute.box', (io) => {
        const polyhedron = new ConvexPolyhedron3(vertices(io, 8), indices(io),
            true, true);
        emitConvexPolyhedron(io, new DistPoint3ConvexPolyhedron3(12)
            .compute(io.vec(3), polyhedron));
    }, { exact: true });

    family.case('DistPoint3ConvexPolyhedron3.compute.noPlanes', (io) => {
        const point = io.vec(3);
        const polyhedron = new ConvexPolyhedron3(vertices(io, 3), indices(io),
            true, true);
        emitConvexPolyhedron(io,
            new DistPoint3ConvexPolyhedron3().compute(point, polyhedron));
    }, { exact: true });

    family.case('DistPoint3ConvexPolyhedron3.setMaxLCPIterations', (io) => {
        const polyhedron = new ConvexPolyhedron3(vertices(io, 4), indices(io),
            true, true);
        const point = io.vec(3);
        const query = new DistPoint3ConvexPolyhedron3();
        query.setMaxLCPIterations(io.integer());
        emitConvexPolyhedron(io, query.compute(point, polyhedron));
    }, { exact: true });

    // ---- DistRay2Arc2.h, DistSegment2Arc2.h ----
    function arc(io: OracleIO): Arc2 {
        const center = io.vec(2);
        const radius = io.real();
        return Arc2.fromCenterRadiusEnds(center, radius, io.vec(2), io.vec(2));
    }

    function emitArc(io: OracleIO, r: { distance: number; sqrDistance: number;
        numClosestPairs: number; parameter: [number, number];
        closest: [[Vector, Vector], [Vector, Vector]] }): void {
        io.outInt(r.numClosestPairs);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        for (let j = 0; j < r.numClosestPairs; ++j) {
            io.outReal(r.parameter[j]);
            io.outVec(r.closest[j][0]);
            io.outVec(r.closest[j][1]);
        }
    }

    family.case('DistRay2Arc2.compute', (io) => {
        const a = arc(io);
        const ray = Ray.fromOriginDirection(io.vec(2), io.vec(2));
        emitArc(io, new DistRay2Arc2().compute(ray, a));
    }, { exact: true });

    family.case('DistSegment2Arc2.compute', (io) => {
        const a = arc(io);
        const segment = Segment.fromEndpoints(io.vec(2), io.vec(2));
        emitArc(io, new DistSegment2Arc2().compute(segment, a));
    }, { exact: true });

    // ---- DistRay3Circle3.h, DistSegment3Circle3.h ----
    function circle3(io: OracleIO): Circle3 {
        const center = io.vec(3);
        const normal = io.vec(3);
        return Circle3.fromCenterNormalRadius(center, normal, io.real());
    }

    function emitCircle3(io: OracleIO, r: { numClosestPairs: number;
        distance: number; sqrDistance: number; equidistant: boolean;
        linearClosest: [Vector, Vector];
        circularClosest: [Vector, Vector] }): void {
        io.outInt(r.numClosestPairs);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outBool(r.equidistant);
        for (let j = 0; j < r.numClosestPairs; ++j) {
            io.outVec(r.linearClosest[j]);
            io.outVec(r.circularClosest[j]);
        }
    }

    // The closed-form branches (PDFSection411, 412, 421) are pure arithmetic
    // and are compared bit for bit.
    family.case('DistRay3Circle3.compute.axis', (io) => {
        const circle = circle3(io);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        emitCircle3(io, new DistRay3Circle3().compute(ray, circle));
    }, { exact: true });

    family.case('DistSegment3Circle3.compute.axis', (io) => {
        const circle = circle3(io);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitCircle3(io, new DistSegment3Circle3().compute(segment, circle));
    }, { exact: true });

    const finalizeDeviation = 'UPSTREAM-FINDINGS DistLine3Circle3.h, issue '
        + '#421: Finalize normalizes the in-plane component of a critical '
        + 'line point without checking that it is nonzero, so upstream '
        + 'reports the circle center as a closest circle point and the '
        + 'distance to the center instead of to the circle.';

    family.case('DistRay3Circle3.compute.axis.deviation', (io) => {
        const circle = circle3(io);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        emitCircle3(io, new DistRay3Circle3().compute(ray, circle));
    }, { exact: true, deviation: finalizeDeviation });

    family.case('DistSegment3Circle3.compute.axis.deviation', (io) => {
        const circle = circle3(io);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitCircle3(io, new DistSegment3Circle3().compute(segment, circle));
    }, { exact: true, deviation: finalizeDeviation });

    // The generic (bisection) path of DistLine3Circle3, reached through the
    // ray and segment wrappers. It is compared with the default scaled
    // tolerance 1e-12 because the port recovers the critical line parameter as
    // t = G(tau) - Dot(M,D)/Dot(M,M) instead of upstream's t = tau + s, the
    // documented conditioning fix of issue #421 item 3 (group 20 root-caused
    // the same residual for DistLine3Circle3.compute). The generator keeps
    // only configurations on which upstream's own answer is verifiably
    // correct.
    family.case('DistRay3Circle3.compute', (io) => {
        const circle = circle3(io);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        emitCircle3(io, new DistRay3Circle3().compute(ray, circle));
    });

    family.case('DistSegment3Circle3.compute', (io) => {
        const circle = circle3(io);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitCircle3(io, new DistSegment3Circle3().compute(segment, circle));
    });

    // numClosestPairs == 2: the tie arm of SelectClosestPoint, where the two
    // candidates are the segment endpoints, so the output comes from two
    // point-circle queries alone and is exact.
    family.case('DistSegment3Circle3.compute.twoPairs', (io) => {
        const circle = circle3(io);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitCircle3(io, new DistSegment3Circle3().compute(segment, circle));
    }, { exact: true });

    family.finish();
});
