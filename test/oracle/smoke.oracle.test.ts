// Replays oracle/cpp/cases/smoke.cpp. Keep the two files in the same order.
import { describe } from 'vitest';
import { AxisAngle } from '../../src/AxisAngle.js';
import { DistPointTriangle } from '../../src/DistPointTriangle.js';
import { DistSegmentSegment } from '../../src/DistSegmentSegment.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { IntrLine3Sphere3FI } from '../../src/IntrLine3Sphere3.js';
import { IntrRay3Sphere3FI } from '../../src/IntrRay3Sphere3.js';
import { Line } from '../../src/Line.js';
import { Ray } from '../../src/Ray.js';
import { Rotation } from '../../src/Rotation.js';
import { Segment } from '../../src/Segment.js';
import { SymmetricEigensolver3x3 } from '../../src/SymmetricEigensolver3x3.js';
import { Triangle } from '../../src/Triangle.js';
import { OracleFamily, type OracleIO } from './harness.js';

function segment3(io: OracleIO): Segment {
    return Segment.fromEndpoints(io.vec(3), io.vec(3));
}

describe('oracle: smoke', () => {
    const family = new OracleFamily('smoke');

    family.case('DistSegmentSegment.compute.3d', (io) => {
        const r = new DistSegmentSegment().compute(segment3(io), segment3(io));
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('DistSegmentSegment.computeRobust.3d', (io) => {
        const r = new DistSegmentSegment().computeRobust(segment3(io), segment3(io));
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('DistPointTriangle.compute.3d', (io) => {
        const point = io.vec(3);
        const triangle = Triangle.fromVertices(io.vec(3), io.vec(3), io.vec(3));
        const r = new DistPointTriangle().compute(point, triangle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReals(r.barycentric);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }, { exact: true });

    family.case('IntrLine3Sphere3.find', (io) => {
        const line = Line.fromOriginDirection(io.vec(3), io.vec(3));
        const sphere = Hypersphere.fromCenterRadius(io.vec(3), io.real());
        const r = new IntrLine3Sphere3FI().find(line, sphere);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    family.case('IntrRay3Sphere3.find', (io) => {
        const ray = Ray.fromOriginDirection(io.vec(3), io.vec(3));
        const sphere = Hypersphere.fromCenterRadius(io.vec(3), io.real());
        const r = new IntrRay3Sphere3FI().find(ray, sphere);
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        for (let i = 0; i < r.numIntersections; ++i) {
            io.outReal(r.parameter[i]);
            io.outVec(r.point[i]);
        }
    }, { exact: true });

    family.case('SymmetricEigensolver3x3.solve', (io) => {
        const a = io.reals(6);
        const aggressive = io.boolean();
        const sortType = io.integer();
        const r = new SymmetricEigensolver3x3().solve(
            a[0], a[1], a[2], a[3], a[4], a[5], aggressive, sortType);
        io.outInt(r.iterations);
        io.outReals(r.evals);
        for (const row of r.evecs) { io.outReals(row); }
    }, { exact: true });

    // sin/cos come from the MSVC runtime in C++ and from V8 here, so this
    // case is compared with a tolerance.
    family.case('Rotation.axisAngleToMatrixAndQuaternion', (io) => {
        const rotation = Rotation.fromAxisAngle(new AxisAngle(io.vec(3), io.real()));
        io.outMat(rotation.toMatrix());
        const q = rotation.toQuaternion();
        for (let i = 0; i < 4; ++i) { io.outReal(q.get(i)); }
    });

    family.finish();
});
