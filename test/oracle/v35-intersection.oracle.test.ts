// Replays oracle/cpp/cases/v35-intersection.cpp. Keep the two files in the
// same order: each case body reads the inputs in the order the C++ generator
// recorded them and pushes the outputs in the order the C++ case recorded
// them.
//
// Every query of group 35 uses only + - * / sqrt fabs min max and
// comparisons, so every ordinary case is declared { exact: true }. The
// sin/cos/tan that build frames, unit directions and cone angles are applied
// on the C++ side to unrecorded draws; only the resulting vectors and
// trigonometric values are recorded, so the C math library never enters the
// compared computation.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { CanonicalBox } from '../../src/CanonicalBox.js';
import { Cone } from '../../src/Cone.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import {
    IntrOrientedBox3Cone3TI
} from '../../src/IntrOrientedBox3Cone3.js';
import { IntrRay3Cone3FI } from '../../src/IntrRay3Cone3.js';
import {
    IntrRay3Plane3FI, IntrRay3Plane3TI, defaultIntrRay3Plane3FIResult,
    intrRay3Plane3FIDoQuery
} from '../../src/IntrRay3Plane3.js';
import { IntrSegment3Cone3FI } from '../../src/IntrSegment3Cone3.js';
import {
    IntrSegment3Plane3FI, IntrSegment3Plane3TI,
    defaultIntrSegment3Plane3FIResult, intrSegment3Plane3FIDoQuery
} from '../../src/IntrSegment3Plane3.js';
import {
    IntrSphere3Cone3FI, IntrSphere3Cone3TI
} from '../../src/IntrSphere3Cone3.js';
import {
    IntrTetrahedron3Tetrahedron3TI, intrTetrahedron3Tetrahedron3InvalidIndex
} from '../../src/IntrTetrahedron3Tetrahedron3.js';
import {
    IntrTriangle2Triangle2FI, IntrTriangle2Triangle2TI
} from '../../src/IntrTriangle2Triangle2.js';
import {
    IntrTriangle3AlignedBox3FI, IntrTriangle3AlignedBox3TI,
    intrTriangle3BoxFacePlanes
} from '../../src/IntrTriangle3AlignedBox3.js';
import {
    IntrTriangle3CanonicalBox3FI, IntrTriangle3CanonicalBox3TI
} from '../../src/IntrTriangle3CanonicalBox3.js';
import {
    IntrTriangle3OrientedBox3FI, IntrTriangle3OrientedBox3TI
} from '../../src/IntrTriangle3OrientedBox3.js';
import {
    IntrTriangle3Triangle3FI, IntrTriangle3Triangle3TI
} from '../../src/IntrTriangle3Triangle3.js';
import type { IntrLine3Cone3FIResult } from '../../src/IntrLine3Cone3.js';
import type { IntrLine3Plane3FIResult } from '../../src/IntrLine3Plane3.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Ray } from '../../src/Ray.js';
import { Segment } from '../../src/Segment.js';
import { Tetrahedron3 } from '../../src/Tetrahedron3.js';
import { Triangle } from '../../src/Triangle.js';
import { Vector, add, mul, sub } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- readers that mirror the C++ generator helpers ----

function plane3(io: OracleIO): Hyperplane {
    const normal = io.vec(3);
    const constant = io.real();
    return Hyperplane.fromNormalConstant(normal, constant);
}

function triangle3(io: OracleIO): Triangle {
    const v0 = io.vec(3);
    const v1 = io.vec(3);
    const v2 = io.vec(3);
    return Triangle.fromVertices(v0, v1, v2);
}

function triangle2(io: OracleIO): Triangle {
    const v0 = io.vec(2);
    const v1 = io.vec(2);
    const v2 = io.vec(2);
    return Triangle.fromVertices(v0, v1, v2);
}

function alignedBox3(io: OracleIO): AlignedBox {
    // The C++ side records the centre and the extent and forms min and max.
    const center = io.vec(3);
    const extent = io.vec(3);
    return AlignedBox.fromMinMax(sub(center, extent), add(center, extent));
}

function orientedBox3(io: OracleIO): OrientedBox {
    const center = io.vec(3);
    const axis0 = io.vec(3);
    const axis1 = io.vec(3);
    const axis2 = io.vec(3);
    const extent = io.vec(3);
    return OrientedBox.fromCenterAxisExtent(center, [axis0, axis1, axis2],
        extent);
}

// Cone3::SetAngle's derived members, recomputed from the recorded cosine,
// sine and tangent exactly as SetAngle computes them. Cone3::angle is not
// read by any query of this group.
function setConeTrig(cone: Cone, cosA: number, sinA: number,
    tanA: number): void {
    cone.angle = 0;
    cone.cosAngle = cosA;
    cone.sinAngle = sinA;
    cone.tanAngle = tanA;
    cone.cosAngleSqr = cosA * cosA;
    cone.sinAngleSqr = sinA * sinA;
    cone.invSinAngle = 1 / sinA;
}

function cone3(io: OracleIO): Cone {
    const origin = io.vec(3);
    const direction = io.vec(3);
    const cosA = io.real();
    const sinA = io.real();
    const tanA = io.real();
    const hmin = io.real();
    const hmax = io.real();
    const cone = new Cone(3);
    cone.ray.origin = origin;
    cone.ray.direction = direction;
    setConeTrig(cone, cosA, sinA, tanA);
    if (hmax === -1) {
        cone.makeInfiniteTruncatedCone(hmin);
    } else {
        cone.makeConeFrustum(hmin, hmax);
    }
    return cone;
}

// ---- output helpers ----

function emitLinePlaneFI(io: OracleIO, r: IntrLine3Plane3FIResult): void {
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter);
    io.outVec(r.point);
}

function emitPolygon(io: OracleIO, poly: readonly Vector[]): void {
    io.outInt(poly.length);
    for (const p of poly) { io.outVec(p); }
}

function emitTriangleBoxFI(io: OracleIO,
    r: { insidePolygon: Vector[], outsidePolygons: Vector[][] }): void {
    emitPolygon(io, r.insidePolygon);
    io.outInt(r.outsidePolygons.length);
    for (const poly of r.outsidePolygons) { emitPolygon(io, poly); }
}

function emitConeFI(io: OracleIO, r: IntrLine3Cone3FIResult): void {
    io.outBool(r.intersect);
    io.outInt(r.type);
    for (let i = 0; i < 2; ++i) {
        io.outReal(r.t[i].x[0] as number);
        io.outReal(r.t[i].x[1] as number);
        io.outReal(r.t[i].d);
    }
    for (let i = 0; i < 2; ++i) {
        for (let j = 0; j < 3; ++j) {
            io.outReal(r.P[i][j].x[0] as number);
            io.outReal(r.P[i][j].x[1] as number);
            io.outReal(r.P[i][j].d);
        }
    }
}

describe('oracle: v35-intersection', () => {
    const family = new OracleFamily('v35-intersection');

    // ---------------------------- IntrRay3Plane3 --------------------------

    family.case('IntrRay3Plane3.test', (io) => {
        const plane = plane3(io);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        io.outBool(new IntrRay3Plane3TI().test(ray, plane).intersect);
    }, { exact: true });

    family.case('IntrRay3Plane3.find', (io) => {
        const plane = plane3(io);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        emitLinePlaneFI(io, new IntrRay3Plane3FI().find(ray, plane));
    }, { exact: true });

    family.case('IntrRay3Plane3.doQuery.fi', (io) => {
        const plane = plane3(io);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = defaultIntrRay3Plane3FIResult();
        intrRay3Plane3FIDoQuery(origin, direction, plane, r);
        emitLinePlaneFI(io, r);
    }, { exact: true });

    family.case('IntrRay3Plane3.find.throughLatticePoint', (io) => {
        const plane = plane3(io);
        const direction = io.vec(3);
        const origin = io.vec(3);
        const ray = Ray.fromOriginDirection(origin, direction);
        emitLinePlaneFI(io, new IntrRay3Plane3FI().find(ray, plane));
    }, { exact: true });

    family.case('IntrRay3Plane3.test.throughLatticePoint', (io) => {
        const plane = plane3(io);
        const direction = io.vec(3);
        const origin = io.vec(3);
        const ray = Ray.fromOriginDirection(origin, direction);
        io.outBool(new IntrRay3Plane3TI().test(ray, plane).intersect);
    }, { exact: true });

    // -------------------------- IntrSegment3Plane3 ------------------------

    family.case('IntrSegment3Plane3.test', (io) => {
        const plane = plane3(io);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        io.outBool(new IntrSegment3Plane3TI().test(segment, plane).intersect);
    }, { exact: true });

    family.case('IntrSegment3Plane3.find', (io) => {
        const plane = plane3(io);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitLinePlaneFI(io, new IntrSegment3Plane3FI().find(segment, plane));
    }, { exact: true });

    family.case('IntrSegment3Plane3.doQuery.fi', (io) => {
        const plane = plane3(io);
        const segOrigin = io.vec(3);
        const segDirection = io.vec(3);
        const segExtent = io.real();
        const r = defaultIntrSegment3Plane3FIResult();
        intrSegment3Plane3FIDoQuery(segOrigin, segDirection, segExtent, plane, r);
        emitLinePlaneFI(io, r);
    }, { exact: true });

    family.case('IntrSegment3Plane3.find.throughLatticePoint', (io) => {
        const plane = plane3(io);
        io.vec(3);  // the direction the endpoints were built from
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitLinePlaneFI(io, new IntrSegment3Plane3FI().find(segment, plane));
    }, { exact: true });

    family.case('IntrSegment3Plane3.test.throughLatticePoint', (io) => {
        const plane = plane3(io);
        io.vec(3);  // the direction the endpoints were built from
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        io.outBool(new IntrSegment3Plane3TI().test(segment, plane).intersect);
    }, { exact: true });

    // ------------------------ IntrTriangle2Triangle2 ----------------------

    family.case('IntrTriangle2Triangle2.test', (io) => {
        const triangle0 = triangle2(io);
        const triangle1 = triangle2(io);
        io.outBool(new IntrTriangle2Triangle2TI()
            .test(triangle0, triangle1).intersect);
    }, { exact: true });

    family.case('IntrTriangle2Triangle2.find', (io) => {
        const triangle0 = triangle2(io);
        const triangle1 = triangle2(io);
        emitPolygon(io, new IntrTriangle2Triangle2FI()
            .find(triangle0, triangle1).intersection);
    }, { exact: true });

    family.case('IntrTriangle2Triangle2.test.sharedVertex', (io) => {
        const triangle0 = triangle2(io);
        const triangle1 = triangle2(io);
        io.outBool(new IntrTriangle2Triangle2TI()
            .test(triangle0, triangle1).intersect);
    }, { exact: true });

    family.case('IntrTriangle2Triangle2.find.sharedVertex', (io) => {
        const triangle0 = triangle2(io);
        const triangle1 = triangle2(io);
        emitPolygon(io, new IntrTriangle2Triangle2FI()
            .find(triangle0, triangle1).intersection);
    }, { exact: true });

    // --------------------- IntrTetrahedron3Tetrahedron3 -------------------
    //
    // 'separating' carries std::numeric_limits<size_t>::max() on the C++ side
    // and Number.MAX_SAFE_INTEGER here; the sentinel is not representable as
    // a double, so both sides emit -1 for it (the value itself is a
    // documented port convention, not a computed quantity).

    const tetraCase = (io: OracleIO): void => {
        const epsilon = io.real();
        const tetra0 = Tetrahedron3.fromVertices(io.vec(3), io.vec(3),
            io.vec(3), io.vec(3));
        const tetra1 = Tetrahedron3.fromVertices(io.vec(3), io.vec(3),
            io.vec(3), io.vec(3));
        const r = new IntrTetrahedron3Tetrahedron3TI()
            .test(tetra0, tetra1, epsilon);
        io.outBool(r.intersect);
        for (let i = 0; i < 2; ++i) {
            const s = r.separating[i];
            io.outInt(s === intrTetrahedron3Tetrahedron3InvalidIndex ? -1 : s);
        }
    };

    family.case('IntrTetrahedron3Tetrahedron3.test', tetraCase,
        { exact: true });
    family.case('IntrTetrahedron3Tetrahedron3.test.smallEdges', tetraCase,
        { exact: true });

    // The port normalises the edge vectors before comparing against the
    // cosine cutoff; upstream does not, so the whole edge-edge phase is
    // skipped for long edges. docs/UPSTREAM-FINDINGS.md
    // IntrTetrahedron3Tetrahedron3.h, issue #307.
    family.case('IntrTetrahedron3Tetrahedron3.test.edgeCutoffDeviation',
        tetraCase, { deviation: 'UPSTREAM-FINDINGS IntrTetrahedron3Tetrahedron3.h edge-edge cutoff; issue #307' });

    // The port tests projection-interval disjointness on the cross-product
    // axis; upstream tests which side of a plane through one edge endpoint the
    // other tetrahedron lies on. docs/UPSTREAM-FINDINGS.md
    // IntrTetrahedron3Tetrahedron3.h, issue #307.
    family.case('IntrTetrahedron3Tetrahedron3.test.edgeSeparationDeviation',
        tetraCase, { deviation: 'UPSTREAM-FINDINGS IntrTetrahedron3Tetrahedron3.h edge-edge separation test; issue #307' });

    // ------------------------- IntrTriangle3*Box3 -------------------------

    const canonicalTest = (io: OracleIO): void => {
        const box = CanonicalBox.fromExtent(io.vec(3));
        const triangle = triangle3(io);
        io.outBool(new IntrTriangle3CanonicalBox3TI()
            .test(triangle, box).intersect);
    };
    const canonicalFind = (io: OracleIO): void => {
        const box = CanonicalBox.fromExtent(io.vec(3));
        const triangle = triangle3(io);
        emitTriangleBoxFI(io,
            new IntrTriangle3CanonicalBox3FI().find(triangle, box));
    };

    family.case('IntrTriangle3CanonicalBox3.test', canonicalTest,
        { exact: true });
    family.case('IntrTriangle3CanonicalBox3.find', canonicalFind,
        { exact: true });
    family.case('IntrTriangle3CanonicalBox3.test.onFaces', canonicalTest,
        { exact: true });
    family.case('IntrTriangle3CanonicalBox3.find.onFaces', canonicalFind,
        { exact: true });

    const alignedTest = (io: OracleIO): void => {
        const box = alignedBox3(io);
        const triangle = triangle3(io);
        io.outBool(new IntrTriangle3AlignedBox3TI()
            .test(triangle, box).intersect);
    };
    const alignedFind = (io: OracleIO): void => {
        const box = alignedBox3(io);
        const triangle = triangle3(io);
        emitTriangleBoxFI(io,
            new IntrTriangle3AlignedBox3FI().find(triangle, box));
    };

    family.case('IntrTriangle3AlignedBox3.test', alignedTest, { exact: true });
    family.case('IntrTriangle3AlignedBox3.find', alignedFind, { exact: true });
    family.case('IntrTriangle3AlignedBox3.test.onFaces', alignedTest,
        { exact: true });
    family.case('IntrTriangle3AlignedBox3.find.onFaces', alignedFind,
        { exact: true });

    const orientedTest = (io: OracleIO): void => {
        const box = orientedBox3(io);
        const triangle = triangle3(io);
        io.outBool(new IntrTriangle3OrientedBox3TI()
            .test(triangle, box).intersect);
    };
    const orientedFind = (io: OracleIO): void => {
        const box = orientedBox3(io);
        const triangle = triangle3(io);
        emitTriangleBoxFI(io,
            new IntrTriangle3OrientedBox3FI().find(triangle, box));
    };

    // The port's exported 'intrTriangle3BoxFacePlanes' against upstream's
    // inline expressions in the three FI queries.
    family.case('IntrTriangle3OrientedBox3.facePlanes', (io) => {
        const box = orientedBox3(io);
        const planes = intrTriangle3BoxFacePlanes(box.center, box.axis,
            box.extent);
        for (const plane of planes) {
            io.outVec(plane.normal);
            io.outReal(plane.constant);
        }
    }, { exact: true });

    family.case('IntrTriangle3OrientedBox3.test', orientedTest,
        { exact: true });
    family.case('IntrTriangle3OrientedBox3.find', orientedFind,
        { exact: true });
    family.case('IntrTriangle3OrientedBox3.test.onFaces', orientedTest,
        { exact: true });
    family.case('IntrTriangle3OrientedBox3.find.onFaces', orientedFind,
        { exact: true });

    // ------------------------ IntrTriangle3Triangle3 ----------------------

    const stationaryTest = (io: OracleIO): void => {
        const triangle0 = triangle3(io);
        const triangle1 = triangle3(io);
        const r = new IntrTriangle3Triangle3TI().test(triangle0, triangle1);
        io.outBool(r.intersect);
        io.outReal(r.contactTime);
    };
    const stationaryFind = (io: OracleIO): void => {
        const triangle0 = triangle3(io);
        const triangle1 = triangle3(io);
        const r = new IntrTriangle3Triangle3FI().find(triangle0, triangle1);
        io.outBool(r.intersect);
        io.outReal(r.contactTime);
        emitPolygon(io, r.intersection);
    };

    family.case('IntrTriangle3Triangle3.test', stationaryTest,
        { exact: true });
    family.case('IntrTriangle3Triangle3.find', stationaryFind,
        { exact: true });
    family.case('IntrTriangle3Triangle3.test.onPlane', stationaryTest,
        { exact: true });
    family.case('IntrTriangle3Triangle3.find.onPlane', stationaryFind,
        { exact: true });
    family.case('IntrTriangle3Triangle3.test.coplanar', stationaryTest,
        { exact: true });

    const movingTest = (io: OracleIO): void => {
        const triangle0 = triangle3(io);
        const triangle1 = triangle3(io);
        const velocity0 = io.vec(3);
        const velocity1 = io.vec(3);
        const tMax = io.real();
        const r = new IntrTriangle3Triangle3TI()
            .testDynamic(tMax, triangle0, velocity0, triangle1, velocity1);
        io.outBool(r.intersect);
        io.outReal(r.contactTime);
    };
    const movingFind = (io: OracleIO): void => {
        const triangle0 = triangle3(io);
        const triangle1 = triangle3(io);
        const velocity0 = io.vec(3);
        const velocity1 = io.vec(3);
        const tMax = io.real();
        const r = new IntrTriangle3Triangle3FI()
            .findDynamic(tMax, triangle0, velocity0, triangle1, velocity1);
        io.outBool(r.intersect);
        io.outReal(r.contactTime);
        emitPolygon(io, r.intersection);
    };

    family.case('IntrTriangle3Triangle3.test.moving', movingTest,
        { exact: true });
    family.case('IntrTriangle3Triangle3.find.moving', movingFind,
        { exact: true });
    family.case('IntrTriangle3Triangle3.test.moving.coplanar', movingTest,
        { exact: true });
    family.case('IntrTriangle3Triangle3.find.moving.coplanar', movingFind,
        { exact: true });

    // Upstream decides "parallel" with fabs(Dot(N0,N1)) < 1 on unnormalised
    // edge cross products; the port uses |Cross(N0,N1)|^2 > 0, the criterion
    // upstream's own stationary query uses. docs/UPSTREAM-FINDINGS.md
    // IntrTriangle3Triangle3.h moving-triangle overloads, issue #334.
    family.case('IntrTriangle3Triangle3.test.moving.parallelDeviation',
        movingTest, { deviation: 'UPSTREAM-FINDINGS IntrTriangle3Triangle3.h moving parallel test; issue #334' });
    family.case('IntrTriangle3Triangle3.find.moving.parallelDeviation',
        movingFind, { deviation: 'UPSTREAM-FINDINGS IntrTriangle3Triangle3.h moving parallel test; issue #334' });

    // --------------------------- IntrSphere3Cone3 -------------------------

    const sphereConeTest = (io: OracleIO): void => {
        const cone = cone3(io);
        const sphere = Hypersphere.fromCenterRadius(io.vec(3), io.real());
        io.outBool(new IntrSphere3Cone3TI().test(sphere, cone).intersect);
    };

    family.case('IntrSphere3Cone3.test.infiniteCone', sphereConeTest,
        { exact: true });
    family.case('IntrSphere3Cone3.test.infiniteTruncatedCone', sphereConeTest,
        { exact: true });
    family.case('IntrSphere3Cone3.test.finiteCone', sphereConeTest,
        { exact: true });
    family.case('IntrSphere3Cone3.test.coneFrustum', sphereConeTest,
        { exact: true });

    // The FI 'point' is the field the port corrects (upstream omits the cone
    // vertex), so the general case compares 'intersect' only.
    family.case('IntrSphere3Cone3.find', (io) => {
        const cone = cone3(io);
        const sphere = Hypersphere.fromCenterRadius(io.vec(3), io.real());
        io.outBool(new IntrSphere3Cone3FI().find(sphere, cone).intersect);
    }, { exact: true });

    // With the cone vertex exactly at the origin, upstream's 't * D' and the
    // port's 'V + t * D' are bit-identical, so the point is compared too.
    family.case('IntrSphere3Cone3.find.originVertex', (io) => {
        const cone = cone3(io);
        const sphere = Hypersphere.fromCenterRadius(io.vec(3), io.real());
        const r = new IntrSphere3Cone3FI().find(sphere, cone);
        io.outBool(r.intersect);
        io.outVec(r.point);
    }, { exact: true });

    // docs/UPSTREAM-FINDINGS.md IntrSphere3Cone3.h FI query, issue #307.
    family.case('IntrSphere3Cone3.find.vertexDeviation', (io) => {
        const cone = new Cone(3);
        cone.ray.origin = io.vec(3);
        cone.ray.direction = io.vec(3);
        const cosA = io.real();
        const sinA = io.real();
        const tanA = io.real();
        setConeTrig(cone, cosA, sinA, tanA);
        cone.makeInfiniteCone();
        const sphere = Hypersphere.fromCenterRadius(io.vec(3), io.real());
        const r = new IntrSphere3Cone3FI().find(sphere, cone);
        io.outBool(r.intersect);
        io.outVec(r.point);
    }, { deviation: 'UPSTREAM-FINDINGS IntrSphere3Cone3.h FI point omits the cone vertex; issue #307' });

    // ------------------ IntrRay3Cone3 / IntrSegment3Cone3 -----------------

    const rayConeFind = (io: OracleIO): void => {
        const cone = cone3(io);
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        emitConeFI(io, new IntrRay3Cone3FI().find(ray, cone));
    };
    const segmentConeFind = (io: OracleIO): void => {
        const cone = cone3(io);
        const segment = Segment.fromEndpoints(io.vec(3), io.vec(3));
        emitConeFI(io, new IntrSegment3Cone3FI().find(segment, cone));
    };

    family.case('IntrRay3Cone3.find', rayConeFind, { exact: true });
    family.case('IntrSegment3Cone3.find', segmentConeFind, { exact: true });
    family.case('IntrSegment3Cone3.find.throughInterior', segmentConeFind,
        { exact: true });

    // Both inherit the port's corrections of IntrLine3Cone3 for a line
    // through the cone vertex. docs/UPSTREAM-FINDINGS.md IntrLine3Cone3.h
    // CaseC2NotZeroDiscrZero and CaseC2NotZeroDiscrPos Block 3, issues #304
    // and #465.
    family.case('IntrRay3Cone3.find.throughVertexDeviation', rayConeFind,
        { deviation: 'UPSTREAM-FINDINGS IntrLine3Cone3.h through-vertex analysis; issues #304, #465' });
    family.case('IntrSegment3Cone3.find.throughVertexDeviation',
        segmentConeFind,
        { deviation: 'UPSTREAM-FINDINGS IntrLine3Cone3.h through-vertex analysis; issues #304, #465' });

    // ------------------------ IntrOrientedBox3Cone3 -----------------------

    family.case('IntrOrientedBox3Cone3.test', (io) => {
        const box = orientedBox3(io);
        const cone = cone3(io);
        io.outBool(new IntrOrientedBox3Cone3TI().test(box, cone).intersect);
    }, { exact: true });

    // Three calls on one query object: a clipping configuration, a
    // fully-inside-the-slab configuration and the first one again. Upstream
    // leaves stale adjacency bits behind, so its third answer drops the
    // clipped candidate edges. docs/UPSTREAM-FINDINGS.md
    // IntrAlignedBox3Cone3.h BoxFullyInConeSlab, issue #301.
    family.case('IntrOrientedBox3Cone3.test.staleAdjacencyDeviation', (io) => {
        const box = orientedBox3(io);
        const coneA = cone3(io);
        const coneB = cone3(io);
        const query = new IntrOrientedBox3Cone3TI();
        io.outBool(query.test(box, coneA).intersect);
        io.outBool(query.test(box, coneB).intersect);
        io.outBool(query.test(box, coneA).intersect);
    }, { deviation: 'UPSTREAM-FINDINGS IntrAlignedBox3Cone3.h BoxFullyInConeSlab stale adjacency; issue #301' });

    family.finish();
});
