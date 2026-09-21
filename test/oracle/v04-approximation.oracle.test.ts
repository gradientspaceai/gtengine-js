// Replays oracle/cpp/cases/v04-approximation.cpp (verify group 4,
// approximation). Keep the two files in the same order.
import { describe } from 'vitest';
import { ApprCone3 } from '../../src/ApprCone3.js';
import type { ApprCone3Parameters } from '../../src/ApprCone3.js';
import {
    ApprCone3EllipseAndPoints, ApprCone3EllipseAndPointsControl,
    ApprCone3ExtractEllipses
} from '../../src/ApprCone3EllipseAndPoints.js';
import { ApprCylinder3 } from '../../src/ApprCylinder3.js';
import { ApprEllipse2 } from '../../src/ApprEllipse2.js';
import { approximateEllipseByArcs } from '../../src/ApprEllipseByArcs.js';
import { ApprEllipsoid3 } from '../../src/ApprEllipsoid3.js';
import {
    ApprGreatArc3, ApprGreatCircle3
} from '../../src/ApprGreatCircle3.js';
import { ApprParabola2 } from '../../src/ApprParabola2.js';
import { ApprParaboloid3 } from '../../src/ApprParaboloid3.js';
import { Cone } from '../../src/Cone.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { Ellipse3 } from '../../src/Ellipse3.js';
import { ApprTorus3 } from '../../src/ApprTorus3.js';
import type { ApprTorus3Parameters } from '../../src/ApprTorus3.js';
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

    // One to three points make the initial oriented box degenerate, so M has
    // a repeated eigenvalue whose eigenvectors are not a function of the
    // data; only the error, the centre and the extents are emitted.
    family.case('ApprEllipsoid3.compute.degenerateBox', (io) => {
        const pts = points(io, 3);
        const numIterations = io.integer();
        const ellipsoid = new Hyperellipsoid(3);
        const error = new ApprEllipsoid3().compute(pts, numIterations, false,
            ellipsoid);
        io.outReal(error);
        io.outVec(ellipsoid.center);
        io.outVec(ellipsoid.extent);
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

    // ---- ApprCone3 ---------------------------------------------------------
    // The generator fixes the initial cone angle to a dyadic value whose
    // cosine is bit-identical in the MSVC runtime and in V8, so the only libm
    // call left in the iterated cases is the final acos. Everything else is
    // emitted with outRealExact and held to bit-identity.

    // The initial cone: vertex, unit axis and the agreed angle.
    function initialCone(io: OracleIO): ApprCone3Parameters {
        const vertex = io.vec(3);
        const axis = io.vec(3);
        io.integer();  // the angle index j, recorded by C++
        const angle = io.real();
        return { vertex, axis, angle };
    }

    function emitCone(io: OracleIO, cone: ApprCone3Parameters): void {
        io.outVecExact(cone.vertex);
        io.outVecExact(cone.axis);
        io.outReal(cone.angle);  // acos
    }

    family.case('ApprCone3.gaussNewton.initialGuess', (io) => {
        const pts = points(io, 3);
        const maxIterations = io.integer();
        const cone = initialCone(io);
        const r = new ApprCone3().computeGaussNewton(pts, maxIterations,
            0, 0, true, cone);
        io.outInt(r.numIterations);
        io.outBool(r.converged);
        io.outRealExact(r.minError);
        io.outRealExact(r.minErrorDifference);
        io.outRealExact(r.minUpdateLength);
        io.outVecExact(r.minLocation);
        emitCone(io, cone);
    });

    // The C++ case records maxIterations after the initial cone, because it
    // is chosen by probing upstream's own control flow (see the C++ comment).
    function levenbergMarquardtCone(io: OracleIO): void {
        const pts = points(io, 3);
        const maxAdjustments = io.integer();
        const lambdaFactor = io.real();
        const lambdaAdjust = io.real();
        const cone = initialCone(io);
        const maxIterations = io.integer();
        const r = new ApprCone3().computeLevenbergMarquardt(pts, maxIterations,
            0, 0, lambdaFactor, lambdaAdjust, maxAdjustments, true, cone);
        io.outInt(r.numIterations);
        io.outInt(r.numAdjustments);
        io.outBool(r.converged);
        io.outRealExact(r.minError);
        io.outRealExact(r.minErrorDifference);
        io.outRealExact(r.minUpdateLength);
        io.outVecExact(r.minLocation);
        emitCone(io, cone);
    }

    family.case('ApprCone3.levenbergMarquardt.initialGuess',
        levenbergMarquardtCone);

    family.case('ApprCone3.levenbergMarquardt.staleResidual.deviation',
        levenbergMarquardtCone,
        { deviation: 'issue #261: LevenbergMarquardtMinimizer::DoIteration '
            + 'builds -J^T*F from the residual at the previously rejected '
            + 'candidate' });

    // atan2 in ComputeInitialCone, then the cos/acos round trip.
    family.case('ApprCone3.gaussNewton.computeInitialCone', (io) => {
        const pts = points(io, 3);
        const cone: ApprCone3Parameters =
            { vertex: new Vector(3), axis: new Vector(3), angle: 0 };
        const r = new ApprCone3().computeGaussNewton(pts, 0, 0, 0, false, cone);
        io.outInt(r.numIterations);
        io.outBool(r.converged);
        io.outReal(r.minError);
        io.outVec(r.minLocation);
        io.outVec(cone.vertex);
        io.outVec(cone.axis);
        io.outReal(cone.angle);
    });

    family.case('ApprCone3.levenbergMarquardt.computeInitialCone', (io) => {
        const pts = points(io, 3);
        const cone: ApprCone3Parameters =
            { vertex: new Vector(3), axis: new Vector(3), angle: 0 };
        const r = new ApprCone3().computeLevenbergMarquardt(pts, 0, 0, 0,
            1e-3, 10, 2, false, cone);
        io.outInt(r.numIterations);
        io.outInt(r.numAdjustments);
        io.outBool(r.converged);
        io.outReal(r.minError);
        io.outVec(r.minLocation);
        io.outVec(cone.vertex);
        io.outVec(cone.axis);
        io.outReal(cone.angle);
    });

    // ---- ApprTorus3 --------------------------------------------------------
    // Every F and J evaluation calls sin and cos, the non-iterative fit goes
    // through RootsPolynomial.solveCubic, and the spherical angles of the
    // initial guess come from atan2 and acos, so all of these are tolerance
    // cases.

    function torusInput(io: OracleIO): ApprTorus3Parameters {
        const C = io.vec(3);
        const N = io.vec(3);
        const r0 = io.real();
        const r1 = io.real();
        return { C, N, r0, r1 };
    }

    function emitTorus(io: OracleIO, t: ApprTorus3Parameters): void {
        io.outVec(t.C);
        io.outVec(t.N);
        io.outReal(t.r0);
        io.outReal(t.r1);
    }

    function newTorus(): ApprTorus3Parameters {
        return { C: new Vector(3), N: new Vector(3), r0: 0, r1: 0 };
    }

    family.case('ApprTorus3.compute', (io) => {
        const pts = points(io, 3);
        const torus = newTorus();
        const r = new ApprTorus3().compute(pts, torus);
        io.outBool(r.success);
        io.outReal(r.error);
        emitTorus(io, torus);
    });

    family.case('ApprTorus3.gaussNewton.initialGuess', (io) => {
        const pts = points(io, 3);
        const maxIterations = io.integer();
        const torus = torusInput(io);
        const r = new ApprTorus3().computeGaussNewton(pts, maxIterations,
            0, 0, true, torus);
        io.outInt(r.numIterations);
        io.outBool(r.converged);
        io.outReal(r.minError);
        io.outReal(r.minErrorDifference);
        io.outReal(r.minUpdateLength);
        io.outVec(r.minLocation);
        emitTorus(io, torus);
        // Every F and J evaluation calls sin and cos on the current
        // parameters, so the residual and the Jacobian differ in the last bit
        // and the Cholesky solve of the 7x7 normal equations multiplies that
        // by its condition number. The same minimizer, driven from a
        // libm-free initial guess by ApprCone3, is bit-identical over the
        // whole deep run, which is what places the cause in sin/cos rather
        // than in the port. Measured maximum over 2000 records: 1.25e-8.
    }, { tol: 1e-7 });

    family.case('ApprTorus3.levenbergMarquardt.initialGuess', (io) => {
        const pts = points(io, 3);
        const maxAdjustments = io.integer();
        const lambdaFactor = io.real();
        const lambdaAdjust = io.real();
        const torus = torusInput(io);
        const maxIterations = io.integer();
        const r = new ApprTorus3().computeLevenbergMarquardt(pts, maxIterations,
            0, 0, lambdaFactor, lambdaAdjust, maxAdjustments, true, torus);
        io.outInt(r.numIterations);
        io.outInt(r.numAdjustments);
        io.outBool(r.converged);
        io.outReal(r.minError);
        io.outReal(r.minErrorDifference);
        io.outReal(r.minUpdateLength);
        io.outVec(r.minLocation);
        emitTorus(io, torus);
        // The same sin/cos amplification as the Gauss-Newton case above; the
        // damped normal equations are better conditioned. Measured maximum
        // over 2000 records: 9.11e-12.
    }, { tol: 1e-10 });

    family.case('ApprTorus3.gaussNewton.computeInitialTorus', (io) => {
        const pts = points(io, 3);
        const torus = newTorus();
        const r = new ApprTorus3().computeGaussNewton(pts, 0, 0, 0, false, torus);
        io.outInt(r.numIterations);
        io.outBool(r.converged);
        io.outReal(r.minError);
        io.outVec(r.minLocation);
        emitTorus(io, torus);
    });

    family.case('ApprTorus3.levenbergMarquardt.computeInitialTorus', (io) => {
        const pts = points(io, 3);
        const torus = newTorus();
        const r = new ApprTorus3().computeLevenbergMarquardt(pts, 0, 0, 0,
            1e-3, 10, 2, false, torus);
        io.outInt(r.numIterations);
        io.outInt(r.numAdjustments);
        io.outBool(r.converged);
        io.outReal(r.minError);
        io.outVec(r.minLocation);
        emitTorus(io, torus);
    });

    // ---- ApprCone3EllipseAndPoints -----------------------------------------
    // Minimize1's bracket search compares sin/cos-derived error values, so
    // this is a tolerance case; mode 1 puts the samples exactly on the cone
    // that (ellipse, theta) determines so the minimum is deep.

    function ellipse3(io: OracleIO): Ellipse3 {
        const center = io.vec(3);
        const normal = io.vec(3);
        const axis0 = io.vec(3);
        const axis1 = io.vec(3);
        const extent = io.vec(2);
        return Ellipse3.fromCenterNormalAxisExtent(center, normal,
            [axis0, axis1], extent);
    }

    function control(io: OracleIO): ApprCone3EllipseAndPointsControl {
        const c = new ApprCone3EllipseAndPointsControl();
        c.maxSubdivisions = io.integer();
        c.maxBisections = io.integer();
        c.epsilon = io.real();
        c.tolerance = io.real();
        c.padding = io.real();
        c.penalty = io.real();
        return c;
    }

    function emitCone3(io: OracleIO, cone: Cone): void {
        io.outVec(cone.ray.origin);
        io.outVec(cone.ray.direction);
        io.outReal(cone.angle);
        io.outReal(cone.cosAngle);
        io.outReal(cone.sinAngle);
        io.outReal(cone.tanAngle);
        io.outReal(cone.cosAngleSqr);
        io.outReal(cone.sinAngleSqr);
        io.outReal(cone.invSinAngle);
    }

    // Minimize1's parabolic-interpolation loop compares F values that sin and
    // cos make 1 ulp apart, so the bracket can close on a slightly different
    // theta; near the minimum F is flat, so a 1 ulp difference in F moves the
    // argmin far more than it moves F. maxBisections is kept at 1 or 2 so the
    // bracket never closes on that flat region; the deep run then agrees to
    // 7.0e-14, inside the default tolerance.
    family.case('ApprCone3EllipseAndPoints.fit', (io) => {
        const e = ellipse3(io);
        const pts = points(io, 3);
        emitCone3(io, ApprCone3EllipseAndPoints.fit(e, pts, control(io)));
    });

    family.case('ApprCone3EllipseAndPoints.fit.deviation', (io) => {
        const e = ellipse3(io);
        const pts = points(io, 3);
        emitCone3(io, ApprCone3EllipseAndPoints.fit(e, pts, control(io)));
    }, { deviation: 'issue #349: ComputeCone divides by the ellipse extent a '
        + 'without validating it' });

    family.case('ApprCone3ExtractEllipses.extract', (io) => {
        const n = io.integer();
        const pts: Vector[] = [];
        for (let i = 0; i < 2 * n; ++i) { pts.push(io.vec(3)); }
        const cosAngleEpsilon = io.real();
        const boxExtentEpsilon = io.real();
        const extractor = new ApprCone3ExtractEllipses();
        const ellipses = extractor.extract(pts, boxExtentEpsilon,
            cosAngleEpsilon);
        io.outInt(extractor.getPlanes().length);
        for (const plane of extractor.getPlanes()) {
            io.outVec(plane.normal);
            io.outReal(plane.constant);
        }
        io.outInt(extractor.getIndices().length);
        for (const list of extractor.getIndices()) {
            io.outInt(list.length);
            for (const index of list) { io.outInt(index); }
        }
        // Only the plane normal of each ellipse is compared; see the C++
        // comment for why the ApprEllipse2-derived centre, axes and extents
        // are not.
        io.outInt(ellipses.length);
        for (const e of ellipses) { io.outVec(e.normal); }
        // 2000 records x 2 planes x 1024 ApprEllipse2 iterations.
    }, { timeout: 180000 });

    family.case('ApprCone3EllipseAndPoints.fit.throw', (io) => {
        const e = ellipse3(io);
        const pts = points(io, 3);
        emitCone3(io, ApprCone3EllipseAndPoints.fit(e, pts, control(io)));
    });

    family.finish();
});
