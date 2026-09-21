// Replays oracle/cpp/cases/v14-containment.cpp. Keep the two files in the
// same order; each case here reads the inputs the C++ case recorded, in the
// same order, and emits the outputs the C++ case recorded, in the same order.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { Arc2 } from '../../src/Arc2.js';
import { Capsule, type Capsule3 } from '../../src/Capsule.js';
import { Circle3 } from '../../src/Circle3.js';
import { Cone } from '../../src/Cone.js';
import {
    getContainerAlignedBox, inContainerAlignedBox, mergeContainersAlignedBox
} from '../../src/ContAlignedBox.js';
import { getContainerAlignedBox2Arc2 } from '../../src/ContAlignedBox2Arc2.js';
import {
    getContainerCapsule3, inContainerCapsule3, inContainerCapsuleCapsule3,
    inContainerSphereCapsule3, mergeContainersCapsule3
} from '../../src/ContCapsule3.js';
import {
    getContainerCircle2, inContainerCircle2, mergeContainersCircle2
} from '../../src/ContCircle2.js';
import { inContainerCone } from '../../src/ContCone.js';
import { getContainerCylinder3, inContainerCylinder3 } from '../../src/ContCylinder3.js';
import {
    getContainerEllipse2, inContainerEllipse2, mergeContainersEllipse2
} from '../../src/ContEllipse2.js';
import { getContainerEllipse2MinCR } from '../../src/ContEllipse2MinCR.js';
import {
    getContainerEllipsoid3, inContainerEllipsoid3, mergeContainersEllipsoid3
} from '../../src/ContEllipsoid3.js';
import { getContainerEllipsoid3MinCR } from '../../src/ContEllipsoid3MinCR.js';
import { getContainerLozenge3, inContainerLozenge3 } from '../../src/ContLozenge3.js';
import {
    getContainerOrientedBox2, inContainerOrientedBox2, mergeContainersOrientedBox2
} from '../../src/ContOrientedBox2.js';
import {
    getContainerOrientedBox3, inContainerOrientedBox3, mergeContainersOrientedBox3
} from '../../src/ContOrientedBox3.js';
import { PointInPolygon2 } from '../../src/ContPointInPolygon2.js';
import {
    PointInPolyhedron3, PointInPolyhedron3Face, PointInPolyhedron3FaceType
} from '../../src/ContPointInPolyhedron3.js';
import { circumscribeCircle2, inscribeCircle2 } from '../../src/ContScribeCircle2.js';
import {
    circumscribeCircle3, circumscribeSphere3, inscribeCircle3, inscribeSphere3
} from '../../src/ContScribeCircle3Sphere3.js';
import {
    getContainerSphere3, inContainerSphere3, mergeContainersSphere3
} from '../../src/ContSphere3.js';
import { inContainerTetrahedron3 } from '../../src/ContTetrahedron3.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Hyperplane, type Plane3 } from '../../src/Hyperplane.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { Lozenge3 } from '../../src/Lozenge3.js';
import { Matrix } from '../../src/Matrix.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Ray } from '../../src/Ray.js';
import { Rectangle } from '../../src/Rectangle.js';
import { Segment } from '../../src/Segment.js';
import { Tetrahedron3 } from '../../src/Tetrahedron3.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- shared readers, mirroring the C++ generators ----

function cloud(io: OracleIO, n: number, count: number): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) { points.push(io.vec(n)); }
    return points;
}

function frame(io: OracleIO, n: number): Vector[] {
    const axis: Vector[] = [];
    for (let i = 0; i < n; ++i) { axis.push(io.vec(n)); }
    return axis;
}

function abox(io: OracleIO, n: number): AlignedBox {
    const min = io.vec(n);
    const max = io.vec(n);
    return AlignedBox.fromMinMax(min, max);
}

function obox(io: OracleIO, n: number): OrientedBox {
    const center = io.vec(n);
    const axis = frame(io, n);
    const extent = io.vec(n);
    return OrientedBox.fromCenterAxisExtent(center, axis, extent);
}

function ellip(io: OracleIO, n: number): Hyperellipsoid {
    const center = io.vec(n);
    const axis = frame(io, n);
    const extent = io.vec(n);
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

function sph(io: OracleIO, n: number): Hypersphere {
    const center = io.vec(n);
    const radius = io.real();
    return Hypersphere.fromCenterRadius(center, radius);
}

function cap3(io: OracleIO): Capsule3 {
    const p0 = io.vec(3);
    const p1 = io.vec(3);
    const radius = io.real();
    return Capsule.fromSegmentRadius(Segment.fromEndpoints(p0, p1), radius);
}

// The C++ Cn<N> records the ray, then the angle and the six values Cone's
// SetAngle derives from it with cos/sin/tan, then the height pair with -1
// standing for an infinite maximum height. The replay assigns the recorded
// trigonometry instead of recomputing it, so the C math library never enters
// the compared computation.
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

function outABox(io: OracleIO, box: AlignedBox): void {
    io.outVec(box.min);
    io.outVec(box.max);
}

function outSphere(io: OracleIO, s: Hypersphere): void {
    io.outVec(s.center);
    io.outReal(s.radius);
}

function outOBox(io: OracleIO, box: OrientedBox): void {
    io.outVec(box.center);
    for (const a of box.axis) { io.outVec(a); }
    io.outVec(box.extent);
}

function outEllip(io: OracleIO, e: Hyperellipsoid): void {
    io.outVec(e.center);
    for (const a of e.axis) { io.outVec(a); }
    io.outVec(e.extent);
}

function outCapsule(io: OracleIO, c: Capsule3): void {
    io.outVec(c.segment.p[0]);
    io.outVec(c.segment.p[1]);
    io.outReal(c.radius);
}

function matrixFromCols(axis: readonly Vector[]): Matrix {
    const m = Matrix.zero(axis.length, axis.length);
    for (let j = 0; j < axis.length; ++j) { m.setCol(j, axis[j]); }
    return m;
}

describe('oracle: v14-containment', () => {
    const family = new OracleFamily('v14-containment');

    // ============================ ContAlignedBox =========================

    const alignedBoxGetContainer = (n: number) => (io: OracleIO): void => {
        const count = io.integer();
        const points = cloud(io, n, count);
        const box = getContainerAlignedBox(points);
        io.outBool(box !== null);
        outABox(io, box!);
    };

    family.case('ContAlignedBox.getContainer.2d',
        alignedBoxGetContainer(2), { exact: true });
    family.case('ContAlignedBox.getContainer.3d',
        alignedBoxGetContainer(3), { exact: true });

    const alignedBoxInContainer = (n: number) => (io: OracleIO): void => {
        const box = abox(io, n);
        const p = io.vec(n);
        io.outBool(inContainerAlignedBox(p, box));
    };

    family.case('ContAlignedBox.inContainer.2d',
        alignedBoxInContainer(2), { exact: true });
    family.case('ContAlignedBox.inContainer.3d',
        alignedBoxInContainer(3), { exact: true });

    family.case('ContAlignedBox.mergeContainers.3d', (io) => {
        const box0 = abox(io, 3);
        const box1 = abox(io, 3);
        const merge = mergeContainersAlignedBox(box0, box1);
        io.outBool(true);
        outABox(io, merge);
    }, { exact: true });

    // The C++ case records the six min/max values of each box as 0/1 flags
    // (0 means +0, 1 means -0), so that stdMin / stdMax are evaluated on
    // signed zeros; Math.min / Math.max would order -0 below +0.
    family.case('ContAlignedBox.mergeContainers.signedZero', (io) => {
        const read = (): AlignedBox => {
            const min = new Vector(3);
            const max = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                min.set(i, io.integer() === 0 ? 0 : -0);
                max.set(i, io.integer() === 0 ? 0 : -0);
            }
            return AlignedBox.fromMinMax(min, max);
        };
        const box0 = read();
        const box1 = read();
        const merge = mergeContainersAlignedBox(box0, box1);
        io.outBool(true);
        outABox(io, merge);
    }, { exact: true });

    // ====================== ContCircle2 / ContSphere3 ====================

    family.case('ContCircle2.getContainer', (io) => {
        const count = io.integer();
        const points = cloud(io, 2, count);
        io.outBool(true);
        outSphere(io, getContainerCircle2(points));
    }, { exact: true });

    family.case('ContCircle2.inContainer', (io) => {
        const circle = sph(io, 2);
        const p = io.vec(2);
        io.outBool(inContainerCircle2(p, circle));
    }, { exact: true });

    family.case('ContCircle2.mergeContainers', (io) => {
        const c0 = sph(io, 2);
        const c1 = sph(io, 2);
        io.outBool(true);
        outSphere(io, mergeContainersCircle2(c0, c1));
    }, { exact: true });

    family.case('ContSphere3.getContainer', (io) => {
        const count = io.integer();
        const points = cloud(io, 3, count);
        io.outBool(true);
        outSphere(io, getContainerSphere3(points));
    }, { exact: true });

    family.case('ContSphere3.inContainer', (io) => {
        const sphere = sph(io, 3);
        const p = io.vec(3);
        io.outBool(inContainerSphere3(p, sphere));
    }, { exact: true });

    family.case('ContSphere3.mergeContainers', (io) => {
        const s0 = sph(io, 3);
        const s1 = sph(io, 3);
        io.outBool(true);
        outSphere(io, mergeContainersSphere3(s0, s1));
    }, { exact: true });

    // ========================= ContPointInPolygon2 =======================

    const polygonCase = (query: (pip: PointInPolygon2, p: Vector) => boolean) =>
        (io: OracleIO): void => {
            const m = io.integer();
            const poly = cloud(io, 2, 2 * m);
            const p = io.vec(2);
            io.outBool(query(new PointInPolygon2(poly), p));
        };

    family.case('ContPointInPolygon2.contains.convex',
        polygonCase((pip, p) => pip.contains(p)), { exact: true });

    family.case('ContPointInPolygon2.contains.nonconvex', (io) => {
        const poly = cloud(io, 2, 6);
        const p = io.vec(2);
        io.outBool(new PointInPolygon2(poly).contains(p));
    }, { exact: true });

    family.case('ContPointInPolygon2.containsConvexOrderN',
        polygonCase((pip, p) => pip.containsConvexOrderN(p)), { exact: true });

    family.case('ContPointInPolygon2.containsConvexOrderLogN',
        polygonCase((pip, p) => pip.containsConvexOrderLogN(p)), { exact: true });

    family.case('ContPointInPolygon2.containsQuadrilateral',
        polygonCase((pip, p) => pip.containsQuadrilateral(p)), { exact: true });

    // ========================= ContAlignedBox2Arc2 =======================

    family.case('ContAlignedBox2Arc2.getContainer', (io) => {
        const center = io.vec(2);
        const radius = io.real();
        const end0 = io.vec(2);
        const end1 = io.vec(2);
        const arc = Arc2.fromCenterRadiusEnds(center, radius, end0, end1);
        io.outBool(true);
        outABox(io, getContainerAlignedBox2Arc2(arc));
    }, { exact: true });

    // ============================ ContCapsule3 ===========================

    family.case('ContCapsule3.getContainer', (io) => {
        const count = io.integer();
        const points = cloud(io, 3, count);
        io.outBool(true);
        outCapsule(io, getContainerCapsule3(points));
    }, { exact: true });

    family.case('ContCapsule3.inContainer.point', (io) => {
        const capsule = cap3(io);
        const p = io.vec(3);
        io.outBool(inContainerCapsule3(p, capsule));
    }, { exact: true });

    family.case('ContCapsule3.inContainer.sphere', (io) => {
        const capsule = cap3(io);
        const sphere = sph(io, 3);
        io.outBool(inContainerSphereCapsule3(sphere, capsule));
    }, { exact: true });

    // The C++ case records the CONTAINER capsule first, then the test capsule.
    family.case('ContCapsule3.inContainer.capsule', (io) => {
        const capsule = cap3(io);
        const testCapsule = cap3(io);
        io.outBool(inContainerCapsuleCapsule3(testCapsule, capsule));
    }, { exact: true });

    family.case('ContCapsule3.mergeContainers', (io) => {
        const c0 = cap3(io);
        const c1 = cap3(io);
        io.outBool(true);
        outCapsule(io, mergeContainersCapsule3(c0, c1));
    }, { exact: true });

    // ============================ ContCylinder3 ==========================

    family.case('ContCylinder3.getContainer', (io) => {
        const count = io.integer();
        const points = cloud(io, 3, count);
        const cylinder = getContainerCylinder3(points);
        io.outBool(true);
        io.outVec(cylinder.axis.origin);
        io.outVec(cylinder.axis.direction);
        io.outReal(cylinder.radius);
        io.outReal(cylinder.height);
    }, { exact: true });

    family.case('ContCylinder3.inContainer', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const radius = io.real();
        const height = io.real();
        const cylinder = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin, direction), radius, height);
        const p = io.vec(3);
        io.outBool(inContainerCylinder3(p, cylinder));
    }, { exact: true });

    // ============================= ContLozenge3 ==========================

    const outLozengeNoCenter = (io: OracleIO, l: Lozenge3): void => {
        io.outReal(l.radius);
        io.outVec(l.rectangle.axis[0]);
        io.outVec(l.rectangle.axis[1]);
        io.outReal(l.rectangle.extent.get(0));
        io.outReal(l.rectangle.extent.get(1));
    };

    // The rectangle centre is compared in the two cases below, not here: the
    // port centres the rectangle on the midpoint of the fitted parameter
    // interval where upstream uses a corner (issue #174).
    family.case('ContLozenge3.getContainer', (io) => {
        const count = io.integer();
        const points = cloud(io, 3, count);
        const lozenge = getContainerLozenge3(points);
        io.outBool(true);
        outLozengeNoCenter(io, lozenge);
    }, { exact: true });

    family.case('ContLozenge3.getContainer.cornerDeviation', (io) => {
        const count = io.integer();
        const points = cloud(io, 3, count);
        io.outVec(getContainerLozenge3(points).rectangle.center);
    }, { deviation: 'docs/UPSTREAM-FINDINGS.md ContLozenge3.h, issue #174' });

    // Coincident points are the only configuration that reaches upstream's
    // "container is a sphere" branch, where its own centre formula is the
    // midpoint expression the port uses everywhere, so the centre is compared
    // here. See the C++ case for why the offsets are necessarily zero.
    family.case('ContLozenge3.getContainer.sphereBranch', (io) => {
        const count = io.integer();
        const points = cloud(io, 3, count);
        const lozenge = getContainerLozenge3(points);
        io.outBool(true);
        io.outVec(lozenge.rectangle.center);
        outLozengeNoCenter(io, lozenge);
    }, { exact: true });

    family.case('ContLozenge3.inContainer', (io) => {
        const center = io.vec(3);
        const axis = frame(io, 3);
        const extent = io.vec(2);
        const radius = io.real();
        const lozenge = Lozenge3.fromRectangleRadius(
            Rectangle.fromCenterAxisExtent(center, [axis[0], axis[1]], extent),
            radius);
        const p = io.vec(3);
        io.outBool(inContainerLozenge3(p, lozenge));
    }, { exact: true });

    // ================= ContOrientedBox2 / ContOrientedBox3 ===============

    const oboxGetContainer = (n: number) => (io: OracleIO): void => {
        const count = io.integer();
        const points = cloud(io, n, count);
        const box = n === 2 ? getContainerOrientedBox2(points)
            : getContainerOrientedBox3(points);
        io.outBool(box !== null);
        outOBox(io, box!);
    };

    family.case('ContOrientedBox2.getContainer', oboxGetContainer(2), { exact: true });
    family.case('ContOrientedBox3.getContainer', oboxGetContainer(3), { exact: true });

    const oboxInContainer = (n: number) => (io: OracleIO): void => {
        const box = obox(io, n);
        const p = io.vec(n);
        io.outBool(n === 2 ? inContainerOrientedBox2(p, box)
            : inContainerOrientedBox3(p, box));
    };

    family.case('ContOrientedBox2.inContainer', oboxInContainer(2), { exact: true });
    family.case('ContOrientedBox2.inContainer.leftHanded',
        oboxInContainer(2), { exact: true });
    family.case('ContOrientedBox3.inContainer', oboxInContainer(3), { exact: true });
    family.case('ContOrientedBox3.inContainer.leftHanded',
        oboxInContainer(3), { exact: true });

    const oboxMerge = (n: number) => (io: OracleIO): void => {
        const box0 = obox(io, n);
        const box1 = obox(io, n);
        io.outBool(true);
        outOBox(io, n === 2 ? mergeContainersOrientedBox2(box0, box1)
            : mergeContainersOrientedBox3(box0, box1));
    };

    family.case('ContOrientedBox2.mergeContainers', oboxMerge(2), { exact: true });
    family.case('ContOrientedBox3.mergeContainers', oboxMerge(3), { exact: true });
    family.case('ContOrientedBox3.mergeContainers.leftHanded',
        oboxMerge(3), { exact: true });

    // =================== ContEllipse2 / ContEllipsoid3 ===================

    const ellipGetContainer = (n: number) => (io: OracleIO): void => {
        const count = io.integer();
        const points = cloud(io, n, count);
        const e = n === 2 ? getContainerEllipse2(points)
            : getContainerEllipsoid3(points);
        io.outBool(e !== null);
        outEllip(io, e!);
    };

    family.case('ContEllipse2.getContainer', ellipGetContainer(2), { exact: true });
    family.case('ContEllipsoid3.getContainer', ellipGetContainer(3), { exact: true });

    const ellipInContainer = (n: number) => (io: OracleIO): void => {
        const e = ellip(io, n);
        const p = io.vec(n);
        io.outBool(n === 2 ? inContainerEllipse2(p, e) : inContainerEllipsoid3(p, e));
    };

    family.case('ContEllipse2.inContainer', ellipInContainer(2), { exact: true });
    family.case('ContEllipse2.inContainer.leftHanded',
        ellipInContainer(2), { exact: true });
    family.case('ContEllipsoid3.inContainer', ellipInContainer(3), { exact: true });
    family.case('ContEllipsoid3.inContainer.leftHanded',
        ellipInContainer(3), { exact: true });

    const ellipMerge = (n: number) => (io: OracleIO): void => {
        const e0 = ellip(io, n);
        const e1 = ellip(io, n);
        io.outBool(true);
        outEllip(io, n === 2 ? mergeContainersEllipse2(e0, e1)
            : mergeContainersEllipsoid3(e0, e1));
    };

    family.case('ContEllipse2.mergeContainers', ellipMerge(2), { exact: true });
    family.case('ContEllipsoid3.mergeContainers', ellipMerge(3), { exact: true });
    family.case('ContEllipsoid3.mergeContainers.leftHanded',
        ellipMerge(3), { exact: true });

    // =============================== ContCone ============================

    const coneCase = (n: number) => (io: OracleIO): void => {
        const c = cone(io, n);
        const p = io.vec(n);
        io.outBool(inContainerCone(p, c));
    };

    family.case('ContCone.inContainer.3d', coneCase(3), { exact: true });
    family.case('ContCone.inContainer.2d', coneCase(2), { exact: true });

    // =========================== ContEllipse2MinCR =======================

    const minCR2 = (io: OracleIO): void => {
        const C = io.vec(2);
        const axis = frame(io, 2);
        const count = io.integer();
        const points = cloud(io, 2, count);
        const D = getContainerEllipse2MinCR(points, C, matrixFromCols(axis));
        io.outReal(D[0]);
        io.outReal(D[1]);
    };

    family.case('ContEllipse2MinCR.compute', minCR2, { exact: true });
    family.case('ContEllipse2MinCR.compute.verticalLineDeviation', minCR2,
        { deviation: 'docs/UPSTREAM-FINDINGS.md ContEllipse2MinCR.h, issue #234' });

    // ========================== ContEllipsoid3MinCR ======================

    const minCR3 = (io: OracleIO): void => {
        const C = io.vec(3);
        const axis = frame(io, 3);
        const count = io.integer();
        const points = cloud(io, 3, count);
        const D = getContainerEllipsoid3MinCR(points, C, matrixFromCols(axis));
        io.outReal(D[0]);
        io.outReal(D[1]);
        io.outReal(D[2]);
    };

    family.case('ContEllipsoid3MinCR.compute', minCR3, { exact: true });
    // Upstream's LogAssert(numer >= 0) fires on every record of this case; the
    // port clamps the slack to zero as the adjacent upstream comment
    // prescribes and returns a feasible result instead of throwing.
    family.case('ContEllipsoid3MinCR.compute.assertDeviation', minCR3,
        { deviation: 'docs/UPSTREAM-FINDINGS.md ContEllipsoid3MinCR.h, issue #409' });

    // ============ ContScribeCircle2 / ContScribeCircle3Sphere3 ===========

    family.case('ContScribeCircle2.circumscribe', (io) => {
        const v = cloud(io, 2, 3);
        const circle = circumscribeCircle2(v[0], v[1], v[2]);
        io.outBool(circle !== null);
        if (circle !== null) { outSphere(io, circle); }
    }, { exact: true });

    family.case('ContScribeCircle2.inscribe', (io) => {
        const v = cloud(io, 2, 3);
        const circle = inscribeCircle2(v[0], v[1], v[2]);
        io.outBool(circle !== null);
        if (circle !== null) { outSphere(io, circle); }
    }, { exact: true });

    const outCircle3 = (io: OracleIO, c: Circle3): void => {
        io.outVec(c.center);
        io.outVec(c.normal);
        io.outReal(c.radius);
    };

    family.case('ContScribeCircle3Sphere3.circumscribeCircle3', (io) => {
        const v = cloud(io, 3, 3);
        const circle = circumscribeCircle3(v[0], v[1], v[2]);
        io.outBool(circle !== null);
        if (circle !== null) { outCircle3(io, circle); }
    }, { exact: true });

    family.case('ContScribeCircle3Sphere3.circumscribeSphere3', (io) => {
        const v = cloud(io, 3, 4);
        const sphere = circumscribeSphere3(v[0], v[1], v[2], v[3]);
        io.outBool(sphere !== null);
        if (sphere !== null) { outSphere(io, sphere); }
    }, { exact: true });

    family.case('ContScribeCircle3Sphere3.inscribeCircle3', (io) => {
        const v = cloud(io, 3, 3);
        const circle = inscribeCircle3(v[0], v[1], v[2]);
        io.outBool(circle !== null);
        if (circle !== null) { outCircle3(io, circle); }
    }, { exact: true });

    family.case('ContScribeCircle3Sphere3.inscribeSphere3', (io) => {
        const v = cloud(io, 3, 4);
        const sphere = inscribeSphere3(v[0], v[1], v[2], v[3]);
        io.outBool(sphere !== null);
        if (sphere !== null) { outSphere(io, sphere); }
    }, { exact: true });

    // =========================== ContTetrahedron3 ========================

    family.case('ContTetrahedron3.inContainer', (io) => {
        const v = cloud(io, 3, 4);
        const tetra = Tetrahedron3.fromArray(v);
        const p = io.vec(3);
        io.outBool(inContainerTetrahedron3(p, tetra));
    }, { exact: true });

    // ======================= ContPointInPolyhedron3 ======================

    // The C++ case records the vertex list, then each face as
    // (index count, indices, plane normal, plane constant, triangle count,
    // triangles), then the ray count and directions, then the test point.
    const polyhedron = (io: OracleIO): {
        vertices: Vector[], faces: PointInPolyhedron3Face[]
    } => {
        const numVertices = io.integer();
        const vertices = cloud(io, 3, numVertices);
        const numFaces = io.integer();
        const faces: PointInPolyhedron3Face[] = [];
        for (let i = 0; i < numFaces; ++i) {
            const numIndices = io.integer();
            const indices: number[] = [];
            for (let k = 0; k < numIndices; ++k) { indices.push(io.integer()); }
            const normal = io.vec(3);
            const constant = io.real();
            const plane: Plane3 = Hyperplane.fromNormalConstant(normal, constant);
            const numTriangles = io.integer();
            const triangles: number[] = [];
            for (let k = 0; k < numTriangles; ++k) { triangles.push(io.integer()); }
            faces.push(PointInPolyhedron3Face.fromIndicesPlane(
                indices, plane, triangles));
        }
        return { vertices, faces };
    };

    const polyhedronCase = (type: PointInPolyhedron3FaceType,
        method: number | 'index') => (io: OracleIO): void => {
        const { vertices, faces } = polyhedron(io);
        const numRays = 1 + 2 * io.integer();
        const directions: Vector[] = [];
        for (let j = 0; j < numRays; ++j) { directions.push(io.vec(3)); }
        const p = io.vec(3);
        const m = method === 'index' ? 1 + (io.index % 2) : method;
        const query = new PointInPolyhedron3(type, vertices, faces, directions, m);
        io.outBool(query.contains(p));
    };

    const TRIANGLE = PointInPolyhedron3FaceType.TRIANGLE;
    const CONVEX = PointInPolyhedron3FaceType.CONVEX;
    const SIMPLE = PointInPolyhedron3FaceType.SIMPLE;

    family.case('ContPointInPolyhedron3.contains.triangle',
        polyhedronCase(TRIANGLE, 0), { exact: true });
    family.case('ContPointInPolyhedron3.contains.convex0',
        polyhedronCase(CONVEX, 0), { exact: true });
    family.case('ContPointInPolyhedron3.contains.convex1',
        polyhedronCase(CONVEX, 1), { exact: true });
    family.case('ContPointInPolyhedron3.contains.convex2',
        polyhedronCase(CONVEX, 2), { exact: true });
    family.case('ContPointInPolyhedron3.contains.simple0',
        polyhedronCase(SIMPLE, 0), { exact: true });
    family.case('ContPointInPolyhedron3.contains.simple1',
        polyhedronCase(SIMPLE, 1), { exact: true });
    family.case('ContPointInPolyhedron3.contains.unsupported',
        polyhedronCase(SIMPLE, 2), { exact: true });

    const quadDeviation = 'docs/UPSTREAM-FINDINGS.md ContPointInPolyhedron3.h, issue #343';
    family.case('ContPointInPolyhedron3.contains.convex0QuadDeviation',
        polyhedronCase(CONVEX, 0), { deviation: quadDeviation });
    family.case('ContPointInPolyhedron3.contains.convex12QuadDeviation',
        polyhedronCase(CONVEX, 'index'), { deviation: quadDeviation });
    family.case('ContPointInPolyhedron3.contains.simple1QuadDeviation',
        polyhedronCase(SIMPLE, 1), { deviation: quadDeviation });

    family.finish();
});
