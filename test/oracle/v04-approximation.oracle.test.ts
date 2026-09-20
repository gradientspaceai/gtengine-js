// Replays oracle/cpp/cases/v04-approximation.cpp (verify group 4,
// approximation). Keep the two files in the same order.
import { describe } from 'vitest';
import { ApprCylinder3 } from '../../src/ApprCylinder3.js';
import { ApprEllipse2 } from '../../src/ApprEllipse2.js';
import { approximateEllipseByArcs } from '../../src/ApprEllipseByArcs.js';
import { ApprEllipsoid3 } from '../../src/ApprEllipsoid3.js';
import {
    ApprGreatArc3, ApprGreatCircle3
} from '../../src/ApprGreatCircle3.js';
import { ApprParabola2 } from '../../src/ApprParabola2.js';
import { ApprParaboloid3 } from '../../src/ApprParaboloid3.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// Mirrors Points2 / Points3 of the C++ case file: the count is recorded
// first, then that many points, in every generator mode.
function points(io: OracleIO, dimension: number): Vector[] {
    const n = io.integer();
    const list: Vector[] = [];
    for (let i = 0; i < n; ++i) { list.push(io.vec(dimension)); }
    return list;
}

describe('oracle: v04-approximation', () => {
    const family = new OracleFamily('v04-approximation');

    // ---- ApprParabola2 -----------------------------------------------------
    // A 3x3 closed-form inverse and a sum of products: arithmetic only.

    family.case('ApprParabola2.fit', (io) => {
        const r = ApprParabola2.fit(points(io, 2));
        io.outBool(r.success);
        io.outReals(r.u);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParabola2.fitRobust', (io) => {
        const r = ApprParabola2.fitRobust(points(io, 2));
        io.outBool(r.success);
        io.outVec(r.average);
        io.outReals(r.v);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParabola2.fit.throw', (io) => {
        const r = ApprParabola2.fit(points(io, 2));
        io.outBool(r.success);
        io.outReals(r.u);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParabola2.fitRobust.throw', (io) => {
        const r = ApprParabola2.fitRobust(points(io, 2));
        io.outBool(r.success);
        io.outVec(r.average);
        io.outReals(r.v);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    // ---- ApprParaboloid3 ---------------------------------------------------
    // An LDL^T decomposition of a 6x6 system: arithmetic only.

    family.case('ApprParaboloid3.fit', (io) => {
        const r = ApprParaboloid3.fit(points(io, 3));
        io.outBool(r.success);
        io.outReals(r.u);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParaboloid3.fitRobust', (io) => {
        const r = ApprParaboloid3.fitRobust(points(io, 3));
        io.outBool(r.success);
        io.outVec(r.average);
        io.outReals(r.v);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParaboloid3.fit.throw', (io) => {
        const r = ApprParaboloid3.fit(points(io, 3));
        io.outBool(r.success);
        io.outReals(r.u);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParaboloid3.fitRobust.throw', (io) => {
        const r = ApprParaboloid3.fitRobust(points(io, 3));
        io.outBool(r.success);
        io.outVec(r.average);
        io.outReals(r.v);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    // ---- ApprGreatCircle3 --------------------------------------------------

    family.case('ApprGreatCircle3.compute', (io) => {
        io.outVec(new ApprGreatCircle3().compute(points(io, 3)));
    }, { exact: true });

    // atan2 on every sample; the generator keeps the largest-gap search
    // libm-stable (see the C++ comment).
    family.case('ApprGreatArc3.compute', (io) => {
        const r = new ApprGreatArc3().compute(points(io, 3));
        io.outVec(r.normal);
        io.outVec(r.arcEnd0);
        io.outVec(r.arcEnd1);
    });

    // ---- ApprEllipseByArcs -------------------------------------------------

    function arcs(io: OracleIO): void {
        const a = io.real();
        const b = io.real();
        const numArcs = io.integer();
        const r = approximateEllipseByArcs(a, b, numArcs);
        io.outBool(r !== null);
        io.outInt(r === null ? 0 : r.points.length);
        if (r !== null) { for (const p of r.points) { io.outVec(p); } }
        io.outInt(r === null ? 0 : r.centers.length);
        if (r !== null) {
            for (const c of r.centers) { io.outVec(c); }
            io.outReals(r.radii);
        }
    }

    // std::pow is on the path for every intermediate point.
    family.case('ApprEllipseByArcs.approximate', arcs);

    family.case('ApprEllipseByArcs.approximate.deviation', (io) => {
        io.integer();  // the adjacent-double offset k, recorded by C++
        arcs(io);
    }, { deviation: 'issue #322: Circumscribe failure discarded in the '
        + 'intermediate-arc loop, leaving the previous arc stored' });

    // ---- ApprEllipse2 / ApprEllipsoid3 -------------------------------------
    // RootsPolynomial.solveCubic (pow, atan2, sin, cos) is on the centre
    // update, so these are tolerance cases. The iteration count is fixed by
    // the caller, so it cannot drift.

    function hyperellipsoid(io: OracleIO, n: number): Hyperellipsoid {
        const center = io.vec(n);
        const axis: Vector[] = [];
        for (let i = 0; i < n; ++i) { axis.push(io.vec(n)); }
        const extent = io.vec(n);
        return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
    }

    function emitHyperellipsoid(io: OracleIO, error: number,
        e: Hyperellipsoid): void {
        io.outReal(error);
        io.outVec(e.center);
        for (const a of e.axis) { io.outVec(a); }
        io.outVec(e.extent);
    }

    family.case('ApprEllipse2.compute.box', (io) => {
        const pts = points(io, 2);
        const numIterations = io.integer();
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(pts, numIterations, false, ellipse);
        emitHyperellipsoid(io, error, ellipse);
    });

    family.case('ApprEllipse2.compute.ellipse', (io) => {
        const pts = points(io, 2);
        const numIterations = io.integer();
        const ellipse = hyperellipsoid(io, 2);
        const error = new ApprEllipse2().compute(pts, numIterations, true, ellipse);
        emitHyperellipsoid(io, error, ellipse);
    });

    family.case('ApprEllipsoid3.compute.box', (io) => {
        const pts = points(io, 3);
        const numIterations = io.integer();
        const ellipsoid = new Hyperellipsoid(3);
        const error = new ApprEllipsoid3().compute(pts, numIterations, false, ellipsoid);
        emitHyperellipsoid(io, error, ellipsoid);
    });

    family.case('ApprEllipsoid3.compute.ellipsoid', (io) => {
        const pts = points(io, 3);
        const numIterations = io.integer();
        const ellipsoid = hyperellipsoid(io, 3);
        const error = new ApprEllipsoid3().compute(pts, numIterations, true, ellipsoid);
        emitHyperellipsoid(io, error, ellipsoid);
    });

    // ---- ApprCylinder3 -----------------------------------------------------

    function emitCylinder(io: OracleIO, cylinder: Cylinder3): void {
        io.outVec(cylinder.axis.origin);
        io.outVec(cylinder.axis.direction);
        io.outReal(cylinder.radius);
        io.outReal(cylinder.height);
    }

    family.case('ApprCylinder3.compute.eigenIndex', (io) => {
        const pts = points(io, 3);
        const fitter = ApprCylinder3.fromEigenIndex(io.integer());
        const cylinder = new Cylinder3();
        io.outReal(fitter.compute(pts, cylinder));
        emitCylinder(io, cylinder);
    }, { exact: true });

    family.case('ApprCylinder3.compute.specifiedAxis', (io) => {
        const pts = points(io, 3);
        const fitter = ApprCylinder3.fromCylinderAxis(io.vec(3));
        const cylinder = new Cylinder3();
        io.outReal(fitter.compute(pts, cylinder));
        emitCylinder(io, cylinder);
    }, { exact: true });

    family.case('ApprCylinder3.compute.throw', (io) => {
        const pts = points(io, 3);
        const fitter = ApprCylinder3.fromCylinderAxis(io.vec(3));
        const cylinder = new Cylinder3();
        io.outReal(fitter.compute(pts, cylinder));
        emitCylinder(io, cylinder);
    }, { exact: true });

    // cos and sin generate the candidate directions of the hemisphere search.
    family.case('ApprCylinder3.compute.hemisphere', (io) => {
        const pts = points(io, 3);
        const numThetaSamples = io.integer();
        const numPhiSamples = io.integer();
        const fitter = ApprCylinder3.fromHemisphereSearch(0, numThetaSamples,
            numPhiSamples, true);
        const cylinder = new Cylinder3();
        io.outReal(fitter.compute(pts, cylinder));
        emitCylinder(io, cylinder);
    });

    family.case('ApprCylinder3.computeMesh', (io) => {
        const pts = points(io, 3);
        const numThetaSamples = io.integer();
        const numPhiSamples = io.integer();
        // The same triangle fan the C++ case builds from the point count.
        const indices: number[] = [];
        for (let t = 0; t < pts.length - 2; ++t) { indices.push(0, t + 1, t + 2); }
        const fitter = ApprCylinder3.fromHemisphereSearch(0, numThetaSamples,
            numPhiSamples, false);
        const cylinder = new Cylinder3();
        fitter.computeMesh(pts, indices, cylinder);
        emitCylinder(io, cylinder);
    });

    family.finish();
});
