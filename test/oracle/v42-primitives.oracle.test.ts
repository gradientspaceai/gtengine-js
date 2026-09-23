// Replays oracle/cpp/cases/v42-primitives.cpp. Keep the two files in the same
// order: every case body reads the inputs in the order the C++ generator
// recorded them and emits the outputs in the order the C++ case recorded them.
import { describe } from 'vitest';
import { AlignedBox } from '../../src/AlignedBox.js';
import { Arc2 } from '../../src/Arc2.js';
import { CanonicalBox } from '../../src/CanonicalBox.js';
import { Capsule } from '../../src/Capsule.js';
import { Circle3 } from '../../src/Circle3.js';
import { Cylinder3 } from '../../src/Cylinder3.js';
import { Ellipse3 } from '../../src/Ellipse3.js';
import { Frustum3 } from '../../src/Frustum3.js';
import { Halfspace } from '../../src/Halfspace.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Line } from '../../src/Line.js';
import { Lozenge3 } from '../../src/Lozenge3.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Parallelepiped3 } from '../../src/Parallelepiped3.js';
import { Parallelogram2 } from '../../src/Parallelogram2.js';
import { Polyhedron3 } from '../../src/Polyhedron3.js';
import { Ray } from '../../src/Ray.js';
import { Rectangle } from '../../src/Rectangle.js';
import { Sector2 } from '../../src/Sector2.js';
import { Segment } from '../../src/Segment.js';
import { SegmentMesh } from '../../src/SegmentMesh.js';
import { Torus3 } from '../../src/Torus3.js';
import { Triangle } from '../../src/Triangle.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The six comparison methods every value type of this group defines, mirroring
// upstream's operator==, !=, <, <=, > and >=.
interface Comparable<T> {
    equals(other: T): boolean;
    notEquals(other: T): boolean;
    lessThan(other: T): boolean;
    lessThanOrEqual(other: T): boolean;
    greaterThan(other: T): boolean;
    greaterThanOrEqual(other: T): boolean;
}

// Mirrors the C++ EmitComparisons.
function emitComparisons<T extends Comparable<T>>(io: OracleIO, a: T, b: T): void {
    io.outBool(a.equals(b));
    io.outBool(a.notEquals(b));
    io.outBool(a.lessThan(b));
    io.outBool(a.lessThanOrEqual(b));
    io.outBool(a.greaterThan(b));
    io.outBool(a.greaterThanOrEqual(b));
    io.outBool(b.lessThan(a));
    io.outBool(b.lessThanOrEqual(a));
    io.outBool(a.equals(a));
    io.outBool(a.lessThan(a));
    io.outBool(a.lessThanOrEqual(a));
    io.outBool(a.greaterThanOrEqual(a));
}

// Mirrors the C++ ComparePair: the K members of the first object then the K
// members of the second.
function comparePair(io: OracleIO, k: number): [number[], number[]] {
    const a = io.reals(k);
    const b = io.reals(k);
    return [a, b];
}

// Mirrors the C++ Unpack.
function unpack(src: readonly number[], offset: number, n: number): Vector {
    const v = new Vector(n);
    for (let i = 0; i < n; ++i) {
        v.values[i] = src[offset + i];
    }
    return v;
}

function emitVectors(io: OracleIO, vectors: readonly Vector[]): void {
    for (const v of vectors) {
        io.outVec(v);
    }
}

// Mirrors the C++ EmitTemplatedDefaults: the default state of every type of
// the group that is templated on the dimension.
function emitTemplatedDefaults(io: OracleIO, n: number): void {
    const box = new AlignedBox(n);
    io.outVec(box.min);
    io.outVec(box.max);
    const canonicalBox = new CanonicalBox(n);
    io.outVec(canonicalBox.extent);
    const halfspace = new Halfspace(n);
    io.outVec(halfspace.normal);
    io.outReal(halfspace.constant);
    const hypersphere = new Hypersphere(n);
    io.outVec(hypersphere.center);
    io.outReal(hypersphere.radius);
    const line = new Line(n);
    io.outVec(line.origin);
    io.outVec(line.direction);
    const ray = new Ray(n);
    io.outVec(ray.origin);
    io.outVec(ray.direction);
    const orientedBox = new OrientedBox(n);
    io.outVec(orientedBox.center);
    emitVectors(io, orientedBox.axis);
    io.outVec(orientedBox.extent);
    const rectangle = new Rectangle(n);
    io.outVec(rectangle.center);
    io.outVec(rectangle.axis[0]);
    io.outVec(rectangle.axis[1]);
    io.outVec(rectangle.extent);
    const segment = new Segment(n);
    io.outVec(segment.p[0]);
    io.outVec(segment.p[1]);
    const triangle = new Triangle(n);
    emitVectors(io, triangle.v);
    const capsule = new Capsule(n);
    io.outVec(capsule.segment.p[0]);
    io.outVec(capsule.segment.p[1]);
    io.outReal(capsule.radius);
}

describe('oracle: v42-primitives', () => {
    const family = new OracleFamily('v42-primitives');

    // -----------------------------------------------------------------
    // Comparison operators.
    // -----------------------------------------------------------------

    family.case('AlignedBox.comparisons.3d', (io) => {
        const [a, b] = comparePair(io, 6);
        const boxA = AlignedBox.fromMinMax(unpack(a, 0, 3), unpack(a, 3, 3));
        const boxB = AlignedBox.fromMinMax(unpack(b, 0, 3), unpack(b, 3, 3));
        emitComparisons(io, boxA, boxB);
    }, { exact: true });

    family.case('CanonicalBox.comparisons.2d', (io) => {
        const [a, b] = comparePair(io, 2);
        const boxA = CanonicalBox.fromExtent(unpack(a, 0, 2));
        const boxB = CanonicalBox.fromExtent(unpack(b, 0, 2));
        emitComparisons(io, boxA, boxB);
    }, { exact: true });

    family.case('Halfspace.comparisons.3d', (io) => {
        const [a, b] = comparePair(io, 4);
        const hsA = Halfspace.fromNormalConstant(unpack(a, 0, 3), a[3]);
        const hsB = Halfspace.fromNormalConstant(unpack(b, 0, 3), b[3]);
        emitComparisons(io, hsA, hsB);
    }, { exact: true });

    family.case('Hypersphere.comparisons.3d', (io) => {
        const [a, b] = comparePair(io, 4);
        const sphereA = Hypersphere.fromCenterRadius(unpack(a, 0, 3), a[3]);
        const sphereB = Hypersphere.fromCenterRadius(unpack(b, 0, 3), b[3]);
        emitComparisons(io, sphereA, sphereB);
    }, { exact: true });

    family.case('Line.comparisons.3d', (io) => {
        const [a, b] = comparePair(io, 6);
        const lineA = Line.fromOriginDirection(unpack(a, 0, 3), unpack(a, 3, 3));
        const lineB = Line.fromOriginDirection(unpack(b, 0, 3), unpack(b, 3, 3));
        emitComparisons(io, lineA, lineB);
    }, { exact: true });

    family.case('Ray.comparisons.2d', (io) => {
        const [a, b] = comparePair(io, 4);
        const rayA = Ray.fromOriginDirection(unpack(a, 0, 2), unpack(a, 2, 2));
        const rayB = Ray.fromOriginDirection(unpack(b, 0, 2), unpack(b, 2, 2));
        emitComparisons(io, rayA, rayB);
    }, { exact: true });

    family.case('OrientedBox.comparisons.3d', (io) => {
        const [a, b] = comparePair(io, 15);
        const build = (s: number[]): OrientedBox => OrientedBox.fromCenterAxisExtent(
            unpack(s, 0, 3),
            [unpack(s, 3, 3), unpack(s, 6, 3), unpack(s, 9, 3)],
            unpack(s, 12, 3));
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Rectangle.comparisons.3d', (io) => {
        const [a, b] = comparePair(io, 11);
        const build = (s: number[]): Rectangle => Rectangle.fromCenterAxisExtent(
            unpack(s, 0, 3), [unpack(s, 3, 3), unpack(s, 6, 3)], unpack(s, 9, 2));
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Segment.comparisons.3d', (io) => {
        const [a, b] = comparePair(io, 6);
        const segA = Segment.fromEndpoints(unpack(a, 0, 3), unpack(a, 3, 3));
        const segB = Segment.fromEndpoints(unpack(b, 0, 3), unpack(b, 3, 3));
        emitComparisons(io, segA, segB);
    }, { exact: true });

    family.case('Triangle.comparisons.2d', (io) => {
        // B goes through the vertex-array constructor and A through the
        // three-vertex one, so both are covered.
        const [a, b] = comparePair(io, 6);
        const triangleA = Triangle.fromVertices(
            unpack(a, 0, 2), unpack(a, 2, 2), unpack(a, 4, 2));
        const triangleB = Triangle.fromVertexArray(
            [unpack(b, 0, 2), unpack(b, 2, 2), unpack(b, 4, 2)]);
        emitComparisons(io, triangleA, triangleB);
    }, { exact: true });

    family.case('Arc2.comparisons', (io) => {
        const [a, b] = comparePair(io, 7);
        const build = (s: number[]): Arc2 => Arc2.fromCenterRadiusEnds(
            unpack(s, 0, 2), s[2], unpack(s, 3, 2), unpack(s, 5, 2));
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Capsule.comparisons.3d', (io) => {
        const [a, b] = comparePair(io, 7);
        const build = (s: number[]): Capsule => Capsule.fromSegmentRadius(
            Segment.fromEndpoints(unpack(s, 0, 3), unpack(s, 3, 3)), s[6]);
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Circle3.comparisons', (io) => {
        const [a, b] = comparePair(io, 7);
        const build = (s: number[]): Circle3 => Circle3.fromCenterNormalRadius(
            unpack(s, 0, 3), unpack(s, 3, 3), s[6]);
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Cylinder3.comparisons', (io) => {
        const [a, b] = comparePair(io, 8);
        const build = (s: number[]): Cylinder3 => Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(unpack(s, 0, 3), unpack(s, 3, 3)), s[6], s[7]);
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Ellipse3.comparisons', (io) => {
        const [a, b] = comparePair(io, 14);
        const build = (s: number[]): Ellipse3 => Ellipse3.fromCenterNormalAxisExtent(
            unpack(s, 0, 3), unpack(s, 3, 3),
            [unpack(s, 6, 3), unpack(s, 9, 3)], unpack(s, 12, 2));
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Frustum3.comparisons', (io) => {
        const [a, b] = comparePair(io, 16);
        const build = (s: number[]): Frustum3 => Frustum3.fromParameters(
            unpack(s, 0, 3), unpack(s, 3, 3), unpack(s, 6, 3), unpack(s, 9, 3),
            s[12], s[13], s[14], s[15]);
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Lozenge3.comparisons', (io) => {
        const [a, b] = comparePair(io, 12);
        const build = (s: number[]): Lozenge3 => Lozenge3.fromRectangleRadius(
            Rectangle.fromCenterAxisExtent(unpack(s, 0, 3),
                [unpack(s, 3, 3), unpack(s, 6, 3)], unpack(s, 9, 2)), s[11]);
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Parallelepiped3.comparisons', (io) => {
        // The axes are unconstrained here, so the objects are built by member
        // assignment: the C++ case does the same to avoid the constructor's
        // right-handedness assert, which the comparisons do not involve.
        const [a, b] = comparePair(io, 12);
        const build = (s: number[]): Parallelepiped3 => {
            const p = new Parallelepiped3();
            p.center = unpack(s, 0, 3);
            p.axis = [unpack(s, 3, 3), unpack(s, 6, 3), unpack(s, 9, 3)];
            return p;
        };
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Parallelogram2.comparisons', (io) => {
        const [a, b] = comparePair(io, 6);
        const build = (s: number[]): Parallelogram2 => {
            const p = new Parallelogram2();
            p.center = unpack(s, 0, 2);
            p.axis = [unpack(s, 2, 2), unpack(s, 4, 2)];
            return p;
        };
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Sector2.comparisons', (io) => {
        // Upstream compares vertex, radius, direction and angle but not the
        // derived cosAngle/sinAngle, which stay at their defaults on both
        // sides (the C++ case assigns the members directly rather than calling
        // SetAngle, which would call std::cos).
        const [a, b] = comparePair(io, 6);
        const build = (s: number[]): Sector2 => {
            const sector = new Sector2();
            sector.vertex = unpack(s, 0, 2);
            sector.radius = s[2];
            sector.direction = unpack(s, 3, 2);
            sector.angle = s[5];
            return sector;
        };
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    family.case('Torus3.comparisons', (io) => {
        const [a, b] = comparePair(io, 14);
        const build = (s: number[]): Torus3 => Torus3.fromCenterFrameRadii(
            unpack(s, 0, 3), unpack(s, 3, 3), unpack(s, 6, 3), unpack(s, 9, 3),
            s[12], s[13]);
        emitComparisons(io, build(a), build(b));
    }, { exact: true });

    // -----------------------------------------------------------------
    // Boxes, rectangles, segments.
    // -----------------------------------------------------------------

    family.case('AlignedBox.centeredFormAndVertices.3d', (io) => {
        const min = io.vec(3);
        const max = io.vec(3);
        const box = AlignedBox.fromMinMax(min, max);
        const form = box.getCenteredForm();
        io.outVec(form.center);
        io.outVec(form.extent);
        emitVectors(io, box.getVertices());
    }, { exact: true });

    family.case('AlignedBox.centeredFormAndVertices.2d', (io) => {
        const min = io.vec(2);
        const max = io.vec(2);
        const box = AlignedBox.fromMinMax(min, max);
        const form = box.getCenteredForm();
        io.outVec(form.center);
        io.outVec(form.extent);
        emitVectors(io, box.getVertices());
    }, { exact: true });

    family.case('CanonicalBox.getVertices.3d', (io) => {
        const box = CanonicalBox.fromExtent(io.vec(3));
        emitVectors(io, box.getVertices());
    }, { exact: true });

    family.case('CanonicalBox.getVertices.2d', (io) => {
        const box = CanonicalBox.fromExtent(io.vec(2));
        emitVectors(io, box.getVertices());
    }, { exact: true });

    family.case('OrientedBox.getVertices.3d', (io) => {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3), io.vec(3)];
        const extent = io.vec(3);
        const box = OrientedBox.fromCenterAxisExtent(center, axis, extent);
        emitVectors(io, box.getVertices());
    }, { exact: true });

    family.case('OrientedBox.getVertices.2d', (io) => {
        const center = io.vec(2);
        const axis = [io.vec(2), io.vec(2)];
        const extent = io.vec(2);
        const box = OrientedBox.fromCenterAxisExtent(center, axis, extent);
        emitVectors(io, box.getVertices());
    }, { exact: true });

    family.case('Rectangle.getVertices.3d', (io) => {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3)];
        const extent = io.vec(2);
        const rectangle = Rectangle.fromCenterAxisExtent(center, axis, extent);
        emitVectors(io, rectangle.getVertices());
    }, { exact: true });

    family.case('Rectangle.getVertices.2d', (io) => {
        const center = io.vec(2);
        const axis = [io.vec(2), io.vec(2)];
        const extent = io.vec(2);
        const rectangle = Rectangle.fromCenterAxisExtent(center, axis, extent);
        emitVectors(io, rectangle.getVertices());
    }, { exact: true });

    family.case('Segment.getCenteredForm.3d', (io) => {
        const p0 = io.vec(3);
        const p1 = io.vec(3);
        const segment = Segment.fromEndpoints(p0, p1);
        const form = segment.getCenteredForm();
        io.outVec(form.center);
        io.outVec(form.direction);
        io.outReal(form.extent);
        const round = Segment.fromCenteredForm(form.center, form.direction,
            form.extent);
        io.outVec(round.p[0]);
        io.outVec(round.p[1]);
    }, { exact: true });

    family.case('Segment.getCenteredForm.2d', (io) => {
        const p0 = io.vec(2);
        const p1 = io.vec(2);
        const segment = Segment.fromEndpoints(p0, p1);
        const form = segment.getCenteredForm();
        io.outVec(form.center);
        io.outVec(form.direction);
        io.outReal(form.extent);
        const round = Segment.fromCenteredForm(form.center, form.direction,
            form.extent);
        io.outVec(round.p[0]);
        io.outVec(round.p[1]);
    }, { exact: true });

    family.case('Segment.setCenteredForm.3d', (io) => {
        const center = io.vec(3);
        const direction = io.vec(3);
        const extent = io.real();
        const segment = new Segment(3);
        segment.setCenteredForm(center, direction, extent);
        io.outVec(segment.p[0]);
        io.outVec(segment.p[1]);
        const copy = Segment.fromPointArray([segment.p[0], segment.p[1]]);
        io.outVec(copy.p[0]);
        io.outVec(copy.p[1]);
    }, { exact: true });

    // -----------------------------------------------------------------
    // Parallelepiped3 and Parallelogram2.
    // -----------------------------------------------------------------

    family.case('Parallelepiped3.getVertices', (io) => {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3), io.vec(3)];
        const parallelepiped = Parallelepiped3.fromCenterAxis(center, axis);
        emitVectors(io, parallelepiped.getVertices());
    }, { exact: true });

    family.case('Parallelepiped3.rightHandedAssert', (io) => {
        const center = io.vec(3);
        const axis = [io.vec(3), io.vec(3), io.vec(3)];
        const parallelepiped = Parallelepiped3.fromCenterAxis(center, axis);
        emitVectors(io, parallelepiped.getVertices());
    }, { exact: true });

    family.case('Parallelogram2.getVertices', (io) => {
        const center = io.vec(2);
        const axis = [io.vec(2), io.vec(2)];
        const parallelogram = Parallelogram2.fromCenterAxis(center, axis);
        emitVectors(io, parallelogram.getVertices());
    }, { exact: true });

    family.case('Parallelogram2.rightHandedAssert', (io) => {
        const center = io.vec(2);
        const axis = [io.vec(2), io.vec(2)];
        const parallelogram = Parallelogram2.fromCenterAxis(center, axis);
        emitVectors(io, parallelogram.getVertices());
    }, { exact: true });

    // -----------------------------------------------------------------
    // Frustum3.
    // -----------------------------------------------------------------

    family.case('Frustum3.updateAndComputeVertices', (io) => {
        const origin = io.vec(3);
        const rVector = io.vec(3);
        const uVector = io.vec(3);
        const dVector = io.vec(3);
        const dMin = io.real();
        const dMax = io.real();
        const uBound = io.real();
        const rBound = io.real();
        const frustum = Frustum3.fromParameters(origin, dVector, uVector,
            rVector, dMin, dMax, uBound, rBound);
        io.outReal(frustum.getDRatio());
        io.outReal(frustum.getMTwoUF());
        io.outReal(frustum.getMTwoRF());
        emitVectors(io, frustum.computeVertices());
        frustum.uBound = io.real();
        frustum.update();
        io.outReal(frustum.getDRatio());
        io.outReal(frustum.getMTwoUF());
        io.outReal(frustum.getMTwoRF());
    }, { exact: true });

    // -----------------------------------------------------------------
    // Arc2 and Sector2.
    // -----------------------------------------------------------------

    family.case('Arc2.contains', (io) => {
        const center = io.vec(2);
        const radius = io.real();
        const e0 = io.vec(2);
        const e1 = io.vec(2);
        const p = io.vec(2);
        const epsKind = io.integer();
        const epsilon = epsKind === 0 ? -1
            : (epsKind === 1 ? 0 : (epsKind === 2 ? 1e-9 : 0.5));
        const arc = Arc2.fromCenterRadiusEnds(center, radius, e0, e1);
        io.outBool(arc.containsOnCircle(p));
        io.outBool(arc.contains(p, epsilon));
    }, { exact: true });

    family.case('Sector2.contains', (io) => {
        const sector = new Sector2();
        sector.vertex = io.vec(2);
        sector.radius = io.real();
        sector.direction = io.vec(2);
        sector.angle = io.real();
        sector.cosAngle = io.real();
        sector.sinAngle = io.real();
        const p = [io.vec(2), io.vec(2), io.vec(2)];
        for (const point of p) {
            io.outBool(sector.contains(point));
        }
    }, { exact: true });

    // std::cos and std::sin in the MSVC runtime against V8's: compared with
    // the default scaled tolerance of 1e-12.
    family.case('Sector2.setAngle', (io) => {
        const sector = new Sector2();
        sector.setAngle(io.real());
        io.outRealExact(sector.angle);
        io.outReal(sector.cosAngle);
        io.outReal(sector.sinAngle);
    });

    // -----------------------------------------------------------------
    // Torus3: std::cos/std::sin in Evaluate and std::atan2 in GetParameters,
    // so both carry the default scaled tolerance of 1e-12.
    // -----------------------------------------------------------------

    family.case('Torus3.evaluate', (io) => {
        const center = io.vec(3);
        const direction0 = io.vec(3);
        const direction1 = io.vec(3);
        const normal = io.vec(3);
        const radius0 = io.real();
        const radius1 = io.real();
        const u = io.real();
        const v = io.real();
        const maxOrder = io.integer();
        const torus = Torus3.fromCenterFrameRadii(center, direction0,
            direction1, normal, radius0, radius1);
        const jet = torus.evaluate(u, v, maxOrder);
        io.outInt(jet.length);
        emitVectors(io, jet);
    });

    family.case('Torus3.getParameters', (io) => {
        const center = io.vec(3);
        const direction0 = io.vec(3);
        const direction1 = io.vec(3);
        const normal = io.vec(3);
        const radius0 = io.real();
        const radius1 = io.real();
        const torus = Torus3.fromCenterFrameRadii(center, direction0,
            direction1, normal, radius0, radius1);
        const x = io.vec(3);
        const parameters = torus.getParameters(x);
        io.outReal(parameters.u);
        io.outReal(parameters.v);
    });

    // -----------------------------------------------------------------
    // Cylinder3.
    // -----------------------------------------------------------------

    family.case('Cylinder3.finiteInfinite', (io) => {
        const origin = io.vec(3);
        const direction = io.vec(3);
        const radius = io.real();
        const height = io.real();
        const axis = Line.fromOriginDirection(origin, direction);
        const cylinder = Cylinder3.fromAxisRadiusHeight(axis, radius, height);
        io.outReal(cylinder.height);
        io.outBool(cylinder.isFinite());
        io.outBool(cylinder.isInfinite());
        cylinder.makeFiniteCylinder(io.real());
        io.outReal(cylinder.height);
        io.outBool(cylinder.isFinite());
        cylinder.makeInfiniteCylinder();
        io.outReal(cylinder.height);
        io.outBool(cylinder.isFinite());
        io.outBool(cylinder.isInfinite());
    }, { exact: true });

    // -----------------------------------------------------------------
    // Polyhedron3.
    // -----------------------------------------------------------------

    family.case('Polyhedron3.queries', (io) => {
        const shape = io.integer();
        const extra = io.integer();
        const poolSize = (shape === 0 ? 4 : (shape === 1 ? 8 : 4)) + extra;
        const pool: Vector[] = [];
        for (let k = 0; k < poolSize; ++k) {
            pool.push(io.vec(3));
        }
        const numIndices = io.integer();
        const indices: number[] = [];
        for (let i = 0; i < numIndices; ++i) {
            indices.push(io.integer());
        }
        const counterClockwise = io.boolean();
        const polyhedron = new Polyhedron3(pool, numIndices, indices,
            counterClockwise);
        io.outBool(polyhedron.isValid());
        io.outBool(polyhedron.counterClockwise());
        io.outInt(polyhedron.getVertices().length);
        io.outInt(polyhedron.getUniqueIndices().length);
        for (const index of polyhedron.getUniqueIndices()) {
            io.outInt(index);
        }
        io.outInt(polyhedron.getIndices().length);
        io.outVec(polyhedron.computeVertexAverage());
        io.outReal(polyhedron.computeSurfaceArea());
        io.outReal(polyhedron.computeVolume());
    }, { exact: true });

    family.case('Polyhedron3.invalidInput', (io) => {
        const poolSize = io.integer();
        const pool: Vector[] = [];
        for (let k = 0; k < poolSize; ++k) {
            pool.push(io.vec(3));
        }
        const numIndices = io.integer();
        const indices: number[] = [];
        for (let i = 0; i < numIndices; ++i) {
            indices.push(io.integer());
        }
        const counterClockwise = io.boolean();
        const polyhedron = new Polyhedron3(pool, numIndices, indices,
            counterClockwise);
        io.outBool(polyhedron.isValid());
        io.outBool(polyhedron.counterClockwise());
        io.outInt(polyhedron.getUniqueIndices().length);
        io.outInt(polyhedron.getIndices().length);
        io.outVec(polyhedron.computeVertexAverage());
        io.outReal(polyhedron.computeSurfaceArea());
        io.outReal(polyhedron.computeVolume());
    }, { exact: true });

    // -----------------------------------------------------------------
    // SegmentMesh.
    // -----------------------------------------------------------------

    const segmentMeshCase = (io: OracleIO, n: number): void => {
        const kind = io.integer();
        const numVertices = io.integer();
        const vertices: Vector[] = [];
        for (let i = 0; i < numVertices; ++i) {
            vertices.push(io.vec(n));
        }
        const validate = io.boolean();
        const numPairs = io.integer();
        const indices: [number, number][] = [];
        for (let i = 0; i < numPairs; ++i) {
            const i0 = io.integer();
            const i1 = io.integer();
            indices.push([i0, i1]);
        }

        let mesh = new SegmentMesh();
        if (kind === 1) {
            mesh = SegmentMesh.fromDisjoint(vertices);
        } else if (kind === 2) {
            mesh = SegmentMesh.fromContiguous(vertices, true);
        } else if (kind === 3) {
            mesh = SegmentMesh.fromContiguous(vertices, false);
        } else if (kind === 4) {
            mesh = SegmentMesh.fromIndexed(vertices, indices, validate);
        }

        io.outInt(mesh.getTopology());
        io.outInt(mesh.getVertices().length);
        emitVectors(io, mesh.getVertices());
        io.outInt(mesh.getIndices().length);
        for (const pair of mesh.getIndices()) {
            io.outInt(pair[0]);
            io.outInt(pair[1]);
        }
    };

    family.case('SegmentMesh.constructors.2d', (io) => {
        segmentMeshCase(io, 2);
    }, { exact: true });

    family.case('SegmentMesh.constructors.3d', (io) => {
        segmentMeshCase(io, 3);
    }, { exact: true });

    // -----------------------------------------------------------------
    // The default constructors of every type of the group.
    // -----------------------------------------------------------------

    family.case('Primitives.defaultConstructors', (io) => {
        const n = io.integer();
        emitTemplatedDefaults(io, n);

        const arc = new Arc2();
        io.outVec(arc.center);
        io.outReal(arc.radius);
        io.outVec(arc.end[0]);
        io.outVec(arc.end[1]);
        const circle = new Circle3();
        io.outVec(circle.center);
        io.outVec(circle.normal);
        io.outReal(circle.radius);
        const cylinder = new Cylinder3();
        io.outVec(cylinder.axis.origin);
        io.outVec(cylinder.axis.direction);
        io.outReal(cylinder.radius);
        io.outReal(cylinder.height);
        io.outBool(cylinder.isFinite());
        const ellipse = new Ellipse3();
        io.outVec(ellipse.center);
        io.outVec(ellipse.normal);
        io.outVec(ellipse.axis[0]);
        io.outVec(ellipse.axis[1]);
        io.outVec(ellipse.extent);
        const frustum = new Frustum3();
        io.outVec(frustum.origin);
        io.outVec(frustum.dVector);
        io.outVec(frustum.uVector);
        io.outVec(frustum.rVector);
        io.outReal(frustum.dMin);
        io.outReal(frustum.dMax);
        io.outReal(frustum.uBound);
        io.outReal(frustum.rBound);
        io.outReal(frustum.getDRatio());
        io.outReal(frustum.getMTwoUF());
        io.outReal(frustum.getMTwoRF());
        emitVectors(io, frustum.computeVertices());
        const lozenge = new Lozenge3();
        io.outVec(lozenge.rectangle.center);
        io.outVec(lozenge.rectangle.axis[0]);
        io.outVec(lozenge.rectangle.axis[1]);
        io.outVec(lozenge.rectangle.extent);
        io.outReal(lozenge.radius);
        const parallelepiped = new Parallelepiped3();
        io.outVec(parallelepiped.center);
        emitVectors(io, parallelepiped.axis);
        emitVectors(io, parallelepiped.getVertices());
        const parallelogram = new Parallelogram2();
        io.outVec(parallelogram.center);
        emitVectors(io, parallelogram.axis);
        emitVectors(io, parallelogram.getVertices());
        const sector = new Sector2();
        io.outVec(sector.vertex);
        io.outReal(sector.radius);
        io.outVec(sector.direction);
        io.outReal(sector.angle);
        io.outReal(sector.cosAngle);
        io.outReal(sector.sinAngle);
        const torus = new Torus3();
        io.outVec(torus.center);
        io.outVec(torus.direction0);
        io.outVec(torus.direction1);
        io.outVec(torus.normal);
        io.outReal(torus.radius0);
        io.outReal(torus.radius1);
        const mesh = new SegmentMesh();
        io.outInt(mesh.getTopology());
        io.outInt(mesh.getVertices().length);
        io.outInt(mesh.getIndices().length);
    }, { exact: true });

    family.finish();
});
