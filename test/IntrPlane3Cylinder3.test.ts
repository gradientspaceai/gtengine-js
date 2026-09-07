import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3Cylinder3TI,
    IntrPlane3Cylinder3FI,
    IntrPlane3Cylinder3FIResultType,
    defaultIntrPlane3Cylinder3FIResult
} from '../src/IntrPlane3Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, plane as arbPlane, positive,
    rotationFrame, unitVector, wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

// A cylinder with the given axis line, radius and height; a negative height
// means an infinite cylinder.
function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const d = Vector.fromArray(direction);
    normalize(d);
    const C = Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(Vector.fromArray(origin), d), radius, 1);
    if (height < 0) {
        C.makeInfiniteCylinder();
    }
    else {
        C.makeFiniteCylinder(height);
    }
    return C;
}

// The squared distance from X to the cylinder axis line.
function sqrRadialDistance(C: Cylinder3, X: Vector): number {
    const d = sub(X, C.axis.origin);
    const along = dot(d, C.axis.direction);
    const radial = sub(d, mul(along, C.axis.direction));
    return dot(radial, radial);
}

const ti = new IntrPlane3Cylinder3TI();
const fi = new IntrPlane3Cylinder3FI();

describe('IntrPlane3Cylinder3', () => {
    it('has a default result with no intersection', () => {
        const result = defaultIntrPlane3Cylinder3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.type)
            .toBe(IntrPlane3Cylinder3FIResultType.noIntersection);
        expect(result.line[0].origin.values).toEqual([0, 0, 0]);
        expect(result.trimLine[1].origin.values).toEqual([0, 0, 0]);
    });

    it('tests an infinite cylinder against parallel and oblique planes', () => {
        const C = cylinder([0, 0, 0], [0, 0, 1], 1, -1);
        // A plane that is not parallel to the axis always intersects.
        expect(ti.test(plane([0, 0, 1], [0, 0, 100]), C).intersect).toBe(true);
        // Planes parallel to the axis intersect only within the radius.
        expect(ti.test(plane([1, 0, 0], [0.5, 0, 0]), C).intersect).toBe(true);
        expect(ti.test(plane([1, 0, 0], [1, 0, 0]), C).intersect).toBe(true);
        expect(ti.test(plane([1, 0, 0], [1.0001, 0, 0]), C).intersect)
            .toBe(false);
    });

    it('tests a finite cylinder with the projection-radius formula', () => {
        // A cylinder of radius 1 and height 4 along z, centered at the origin,
        // spans z in [-2,2].
        const C = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        expect(ti.test(plane([0, 0, 1], [0, 0, 2]), C).intersect).toBe(true);
        expect(ti.test(plane([0, 0, 1], [0, 0, 2.0001]), C).intersect)
            .toBe(false);
        expect(ti.test(plane([1, 0, 0], [1, 0, 0]), C).intersect).toBe(true);
        expect(ti.test(plane([1, 0, 0], [1.0001, 0, 0]), C).intersect)
            .toBe(false);
    });

    it('reports a circle for a plane perpendicular to the axis', () => {
        const C = cylinder([0, 0, 0], [0, 0, 1], 2, -1);
        const result = fi.find(plane([0, 0, 1], [0, 0, 3]), C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrPlane3Cylinder3FIResultType.circle);
        expect(result.ellipse.extent.values[0]).toBeCloseTo(2, 10);
        expect(result.ellipse.extent.values[1]).toBeCloseTo(2, 10);
        expect(result.ellipse.center.values[2]).toBeCloseTo(3, 12);
        expect(result.ellipse.normal.values).toEqual([0, 0, 1]);
    });

    it('reports an ellipse with the expected semi-axes for a tilted plane', () => {
        // A plane tilted by 'a' from the axis-perpendicular plane cuts the
        // radius-r cylinder in an ellipse with semi-axes r and r/cos(a).
        const r = 1.5;
        const a = Math.PI / 6;
        const C = cylinder([0, 0, 0], [0, 0, 1], r, -1);
        const P = plane([0, Math.sin(a), Math.cos(a)], [0, 0, 0]);
        const result = fi.find(P, C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrPlane3Cylinder3FIResultType.ellipse);
        const extents = [result.ellipse.extent.values[0],
            result.ellipse.extent.values[1]].sort((x, y) => x - y);
        expect(extents[0]).toBeCloseTo(r, 8);
        expect(extents[1]).toBeCloseTo(r / Math.cos(a), 8);

        // Points on the ellipse are on the plane and on the cylinder.
        for (let k = 0; k < 8; ++k) {
            const t = 2 * Math.PI * k / 8;
            const X = add(result.ellipse.center,
                add(mul(result.ellipse.extent.values[0] * Math.cos(t),
                    result.ellipse.axis[0]),
                    mul(result.ellipse.extent.values[1] * Math.sin(t),
                        result.ellipse.axis[1])));
            expect(dot(P.normal, X) - P.constant).toBeCloseTo(0, 8);
            expect(sqrRadialDistance(C, X)).toBeCloseTo(r * r, 7);
        }
    });

    it('reports two parallel lines for a plane containing the axis direction', () => {
        const C = cylinder([0, 0, 0], [0, 0, 1], 1, -1);
        const result = fi.find(plane([1, 0, 0], [0.5, 0, 0]), C);
        expect(result.intersect).toBe(true);
        expect(result.type)
            .toBe(IntrPlane3Cylinder3FIResultType.parallelLines);
        const y = Math.sqrt(1 - 0.25);
        for (const L of result.line) {
            expect(L.origin.values[0]).toBeCloseTo(0.5, 12);
            expect(Math.abs(L.origin.values[1])).toBeCloseTo(y, 12);
            expect(L.direction.values).toEqual([0, 0, 1]);
            expect(sqrRadialDistance(C, L.origin)).toBeCloseTo(1, 10);
        }
        expect(result.line[0].origin.values[1])
            .not.toBeCloseTo(result.line[1].origin.values[1], 6);
    });

    it('reports a single line for a tangent plane', () => {
        const C = cylinder([0, 0, 0], [0, 0, 1], 1, -1);
        const result = fi.find(plane([1, 0, 0], [1, 0, 0]), C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrPlane3Cylinder3FIResultType.singleLine);
        expect(result.line[0].origin.values[0]).toBeCloseTo(1, 12);
        expect(result.line[0].origin.values[1]).toBeCloseTo(0, 12);
        expect(result.line[0].direction.values).toEqual([0, 0, 1]);
        // The second line is untouched.
        expect(result.line[1].origin.values).toEqual([0, 0, 0]);
    });

    it('reports no intersection for a separated parallel plane', () => {
        const C = cylinder([0, 0, 0], [0, 0, 1], 1, -1);
        const result = fi.find(plane([1, 0, 0], [2, 0, 0]), C);
        expect(result.intersect).toBe(false);
        expect(result.type)
            .toBe(IntrPlane3Cylinder3FIResultType.noIntersection);
    });

    it('computes trim lines on the end planes of a finite cylinder', () => {
        // A cylinder of height 4 along z centered at the origin has end
        // planes z = -2 and z = 2. The plane y = z meets them in the lines
        // y = -2 and y = 2.
        const C = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const a = Math.SQRT1_2;
        const P = plane([0, -a, a], [0, 0, 0]);
        const result = fi.find(P, C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrPlane3Cylinder3FIResultType.ellipse);
        // The trim lines lie in both the query plane and the end planes.
        for (let i = 0; i < 2; ++i) {
            const L = result.trimLine[i];
            const expectedZ = (i === 0 ? -2 : 2);
            expect(dot(P.normal, L.origin) - P.constant).toBeCloseTo(0, 10);
            expect(L.origin.values[2]).toBeCloseTo(expectedZ, 10);
            expect(L.direction.values[2]).toBeCloseTo(0, 12);
            expect(length(L.direction)).toBeCloseTo(1, 10);
            // Walking along the trim line stays in both planes.
            const X = add(L.origin, mul(3, L.direction));
            expect(dot(P.normal, X) - P.constant).toBeCloseTo(0, 10);
            expect(X.values[2]).toBeCloseTo(expectedZ, 10);
        }
    });

    it('leaves the trim lines at their defaults when the axis is parallel', () => {
        const C = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(plane([1, 0, 0], [0.5, 0, 0]), C);
        expect(result.type)
            .toBe(IntrPlane3Cylinder3FIResultType.parallelLines);
        expect(result.trimLine[0].origin.values).toEqual([0, 0, 0]);
        expect(result.trimLine[1].origin.values).toEqual([0, 0, 0]);
    });

    it('agrees with the TI query and the cylinder equation on random inputs', () => {
        let state = 616161;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numEllipses = 0;
        let numLines = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const radius = 0.5 + Math.abs(rand()) * 2;
            const C = cylinder([rand(), rand(), rand()],
                [rand(), rand(), rand() + 0.001], radius, -1);
            const P = plane([rand(), rand(), rand()],
                [rand() * 2, rand() * 2, rand() * 2]);
            const result = fi.find(P, C);
            expect(ti.test(P, C).intersect).toBe(result.intersect);

            if (result.type === IntrPlane3Cylinder3FIResultType.ellipse
                || result.type === IntrPlane3Cylinder3FIResultType.circle) {
                ++numEllipses;
                for (let k = 0; k < 5; ++k) {
                    const t = 2 * Math.PI * k / 5;
                    const X = add(result.ellipse.center,
                        add(mul(result.ellipse.extent.values[0] * Math.cos(t),
                            result.ellipse.axis[0]),
                            mul(result.ellipse.extent.values[1] * Math.sin(t),
                                result.ellipse.axis[1])));
                    expect(dot(P.normal, X) - P.constant).toBeCloseTo(0, 7);
                    expect(sqrRadialDistance(C, X))
                        .toBeCloseTo(radius * radius, 6);
                }
            }
            else if (result.type
                === IntrPlane3Cylinder3FIResultType.parallelLines) {
                ++numLines;
                for (const L of result.line) {
                    expect(dot(P.normal, L.origin) - P.constant)
                        .toBeCloseTo(0, 8);
                    expect(sqrRadialDistance(C, L.origin))
                        .toBeCloseTo(radius * radius, 7);
                    expect(dot(P.normal, L.direction)).toBeCloseTo(0, 10);
                }
            }
        }
        expect(numEllipses).toBeGreaterThan(200);
        expect(numLines).toBeGreaterThanOrEqual(0);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrPlane3Cylinder3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrPlane3Cylinder3 verification', () => {
    const tiQ = new IntrPlane3Cylinder3TI();
    const fiQ = new IntrPlane3Cylinder3FI();

    const arbFiniteCylinder = fc.tuple(wellScaledVector(3), unitVector(3),
        positive(4), positive(6))
        .map(([c, w, r, h]) => Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(c, w), r, h));

    const arbInfiniteCylinder = fc.tuple(wellScaledVector(3), unitVector(3),
        positive(4))
        .map(([c, w, r]) => Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(c, w), r, -1));

    // An orthonormal pair spanning the plane orthogonal to w.
    function basis(w: Vector): [Vector, Vector] {
        const helper = Math.abs(w.values[0]) < 0.9
            ? Vector.fromArray([1, 0, 0]) : Vector.fromArray([0, 1, 0]);
        const u = cross(w, helper);
        normalize(u);
        return [u, cross(w, u)];
    }

    // Distance from X to the infinite cylinder axis.
    function axisDistance(cyl: Cylinder3, X: Vector): number {
        const d = sub(X, cyl.axis.origin);
        const h = dot(d, cyl.axis.direction);
        return length(sub(d, mul(h, cyl.axis.direction)));
    }

    it('TI agrees with the extreme signed distances of the finite cylinder', () => {
        // The maximum of Dot(N, X - C) over the solid finite cylinder is
        //   (h/2)*|Dot(N,W)| + r*|N - Dot(N,W)*W| = (h/2)*|Dot(N,W)|
        //     + r*|Cross(N,W)|
        // for a unit-length N and W. Sampling the boundary reproduces it.
        check(fc.tuple(arbPlane(3), arbFiniteCylinder), ([P, C]) => {
            const [u, v] = basis(C.axis.direction);
            let lo = Number.POSITIVE_INFINITY;
            let hi = Number.NEGATIVE_INFINITY;
            for (const s of [-0.5, 0.5]) {
                for (let i = 0; i < 180; ++i) {
                    const a = (2 * Math.PI * i) / 180;
                    const X = add(C.axis.origin, add(
                        mul(s * C.height, C.axis.direction),
                        add(mul(C.radius * Math.cos(a), u),
                            mul(C.radius * Math.sin(a), v))));
                    const sd = dot(P.normal, X) - P.constant;
                    lo = Math.min(lo, sd);
                    hi = Math.max(hi, sd);
                }
            }
            const got = tiQ.test(P, C).intersect;
            const scale = 1 + C.radius + C.height;
            if (lo < -1e-6 * scale && hi > 1e-6 * scale) {
                expect(got).toBe(true);
            }
            if (lo > 1e-6 * scale || hi < -1e-6 * scale) {
                expect(got).toBe(false);
            }
        }, 60);
    });

    it('TI treats an infinite cylinder as unbounded along its axis', () => {
        check(fc.tuple(arbPlane(3), arbInfiniteCylinder), ([P, C]) => {
            const dotNW = dot(P.normal, C.axis.direction);
            const d = Math.abs(dot(P.normal,
                sub(C.axis.origin, P.origin)));
            expect(tiQ.test(P, C).intersect)
                .toBe(dotNW !== 0 ? true : d <= C.radius);
        });
    });

    it('the FI ellipse lies on the plane and on the infinite cylinder', () => {
        check(fc.tuple(arbPlane(3), arbInfiniteCylinder), ([P, C]) => {
            const r = fiQ.find(P, C);
            if (r.type !== IntrPlane3Cylinder3FIResultType.ellipse
                && r.type !== IntrPlane3Cylinder3FIResultType.circle) {
                return;
            }
            expect(r.intersect).toBe(true);
            const dotNW = Math.abs(dot(P.normal, C.axis.direction));
            // The ellipse is ill-conditioned when the plane is nearly
            // parallel to the axis (its major extent goes to infinity).
            if (dotNW < 1e-2) { return; }
            const scale = 1 + C.radius / dotNW + length(C.axis.origin)
                + length(P.origin);
            for (let i = 0; i < 24; ++i) {
                const a = (2 * Math.PI * i) / 24;
                const X = add(r.ellipse.center, add(
                    mul(r.ellipse.extent.values[0] * Math.cos(a),
                        r.ellipse.axis[0]),
                    mul(r.ellipse.extent.values[1] * Math.sin(a),
                        r.ellipse.axis[1])));
                expectClose(dot(P.normal, X) - P.constant, 0, 1e-8 * scale, 0);
                expectClose(axisDistance(C, X), C.radius, 1e-8 * scale, 1e-8);
            }
            // The ellipse normal is the plane normal and the minor extent of
            // the section equals the cylinder radius.
            expectVectorClose(r.ellipse.normal, P.normal, 1e-15, 1e-15);
            const minor = Math.min(r.ellipse.extent.values[0],
                r.ellipse.extent.values[1]);
            expectClose(minor, C.radius, 1e-8 * scale, 1e-8);
            // The major extent is radius / |Dot(N,W)|.
            const major = Math.max(r.ellipse.extent.values[0],
                r.ellipse.extent.values[1]);
            expectClose(major, C.radius / dotNW, 1e-7 * scale, 1e-7);
        }, 60);
    });

    it('reports a near-circular section when the plane is orthogonal to the axis', () => {
        // Upstream chooses CIRCLE over ELLIPSE with the exact comparison
        // 'ellipse2.extent[0] != ellipse2.extent[1]'; the two extents of an
        // orthogonal section agree only to rounding, so the reported type is
        // CIRCLE exactly when they come out bit-equal. The port preserves the
        // exact comparison; the geometry is checked instead of the label.
        check(arbInfiniteCylinder, C => {
            const P = Hyperplane.fromNormalOrigin(C.axis.direction,
                C.axis.origin);
            const r = fiQ.find(P, C);
            expect(r.intersect).toBe(true);
            const e0 = r.ellipse.extent.values[0];
            const e1 = r.ellipse.extent.values[1];
            expect(r.type).toBe(e0 === e1
                ? IntrPlane3Cylinder3FIResultType.circle
                : IntrPlane3Cylinder3FIResultType.ellipse);
            expectClose(e0, e1, 1e-12, 1e-10);
            expectClose(e0, C.radius, 1e-10, 1e-9);
        });
    });

    it('labels an exactly orthogonal axis-aligned section a circle', () => {
        const C = Cylinder3.fromAxisRadiusHeight(Line.fromOriginDirection(
            Vector.zero(3), Vector.fromArray([0, 0, 1])), 2, -1);
        const P = Hyperplane.fromNormalOrigin(Vector.fromArray([0, 0, 1]),
            Vector.zero(3));
        const r = fiQ.find(P, C);
        expect(r.type).toBe(IntrPlane3Cylinder3FIResultType.circle);
        expect(r.ellipse.extent.values[0]).toBe(r.ellipse.extent.values[1]);
        expectClose(r.ellipse.extent.values[0], 2, 1e-12, 1e-12);
    });

    it('the FI parallel lines lie on the plane and on the cylinder', () => {
        // The parallel branch is selected by the exact test Dot(N,W) == 0, so
        // the axis direction and the plane normal must be an exactly
        // orthonormal pair; these are (their components are exact binary
        // fractions whose sums of squares are exactly 1).
        const exactPairs: [Vector, Vector][] = [
            [Vector.fromArray([0, 0, 1]), Vector.fromArray([1, 0, 0])],
            [Vector.fromArray([0, 1, 0]), Vector.fromArray([0, 0, 1])],
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 0.6, 0.8])],
            [Vector.fromArray([0.6, 0.8, 0]), Vector.fromArray([-0.8, 0.6, 0])],
            [Vector.fromArray([0, 0.8, 0.6]), Vector.fromArray([0, -0.6, 0.8])],
            [Vector.fromArray([0.8, 0, -0.6]), Vector.fromArray([0.6, 0, 0.8])]
        ];
        check(fc.tuple(fc.constantFrom(...exactPairs), wellScaledVector(3),
            positive(4), wellScaled(-3, 3)),
            ([[W, u], origin, radius, offset]) => {
                expect(dot(W, u)).toBe(0);
                const C = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(origin, W), radius, -1);
                const P = Hyperplane.fromNormalOrigin(u,
                    add(origin, mul(offset, u)));
                const r = fiQ.find(P, C);
                const d = Math.abs(offset);
                if (d > radius) {
                    expect(r.intersect).toBe(false);
                    expect(r.type)
                        .toBe(IntrPlane3Cylinder3FIResultType.noIntersection);
                    return;
                }
                if (d === radius) {
                    expect(r.type)
                        .toBe(IntrPlane3Cylinder3FIResultType.singleLine);
                } else {
                    expect(r.type)
                        .toBe(IntrPlane3Cylinder3FIResultType.parallelLines);
                }
                const scale = 1 + radius + length(origin);
                const numLines = r.type
                    === IntrPlane3Cylinder3FIResultType.singleLine ? 1 : 2;
                for (let i = 0; i < numLines; ++i) {
                    expectVectorClose(r.line[i].direction, W, 1e-15, 1e-15);
                    for (const t of [-2, 0, 3]) {
                        const X = add(r.line[i].origin,
                            mul(t, r.line[i].direction));
                        expectClose(dot(P.normal, X) - P.constant, 0,
                            1e-10 * scale, 0);
                        expectClose(axisDistance(C, X), radius,
                            1e-9 * scale, 1e-9);
                    }
                }
                // The trim lines are untouched in the parallel branch.
                expect(r.trimLine[0].direction.values).toEqual([1, 0, 0]);
                expect(r.trimLine[1].direction.values).toEqual([1, 0, 0]);
            });
    });

    it('the trim lines lie on the plane and on the cylinder end planes', () => {
        check(fc.tuple(arbPlane(3), arbFiniteCylinder), ([P, C]) => {
            const r = fiQ.find(P, C);
            if (r.type !== IntrPlane3Cylinder3FIResultType.ellipse
                && r.type !== IntrPlane3Cylinder3FIResultType.circle) {
                return;
            }
            const W = C.axis.direction;
            const dotNW = Math.abs(dot(P.normal, W));
            if (dotNW < 1e-2 || dotNW > 1 - 1e-6) { return; }
            const offset = mul(0.5 * C.height, W);
            const ends = [sub(C.axis.origin, offset),
                add(C.axis.origin, offset)];
            const scale = 1 + C.radius + C.height + length(C.axis.origin)
                + length(P.origin);
            for (let i = 0; i < 2; ++i) {
                const L = r.trimLine[i];
                for (const t of [-1, 0, 2]) {
                    const X = add(L.origin, mul(t, L.direction));
                    // On the query plane ...
                    expectClose(dot(P.normal, X) - P.constant, 0,
                        1e-9 * scale, 0);
                    // ... and on the cylinder end plane.
                    expectClose(dot(W, sub(X, ends[i])), 0, 1e-9 * scale, 0);
                }
            }
        }, 60);
    });

    it('reports no intersection for a separated finite cylinder', () => {
        check(arbFiniteCylinder, C => {
            const W = C.axis.direction;
            const [u] = basis(W);
            // A plane orthogonal to the axis, above the cylinder top.
            const P = Hyperplane.fromNormalOrigin(W,
                add(C.axis.origin, mul(0.5 * C.height + 1, W)));
            const r = fiQ.find(P, C);
            expect(tiQ.test(P, C).intersect).toBe(false);
            expect(r.intersect).toBe(false);
            expect(r.type).toBe(IntrPlane3Cylinder3FIResultType.noIntersection);
            // A plane parallel to the axis, beyond the cylinder wall.
            const P2 = Hyperplane.fromNormalOrigin(u,
                add(C.axis.origin, mul(C.radius + 1, u)));
            expect(tiQ.test(P2, C).intersect).toBe(false);
            expect(fiQ.find(P2, C).intersect).toBe(false);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbPlane(3), arbInfiniteCylinder, rotationFrame(3),
            wellScaledVector(3)),
            ([P, C, R, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(R[0], v), dot(R[1], v), dot(R[2], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const P2 = Hyperplane.fromNormalOrigin(rot(P.normal),
                    map(P.origin));
                const C2 = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(map(C.axis.origin),
                        rot(C.axis.direction)), C.radius, C.height);
                const dotNW = Math.abs(dot(P.normal, C.axis.direction));
                // The parallel branch is an exact test on Dot(N,W); a
                // rotation perturbs it off zero.
                if (dotNW < 1e-6) { return; }
                const r0 = fiQ.find(P, C);
                const r1 = fiQ.find(P2, C2);
                expect(r1.intersect).toBe(r0.intersect);
                if (r0.type === IntrPlane3Cylinder3FIResultType.ellipse
                    || r0.type === IntrPlane3Cylinder3FIResultType.circle) {
                    expectVectorClose(r1.ellipse.center, map(r0.ellipse.center),
                        1e-7, 1e-7);
                    expectClose(r1.ellipse.extent.values[0],
                        r0.ellipse.extent.values[0], 1e-7, 1e-7);
                    expectClose(r1.ellipse.extent.values[1],
                        r0.ellipse.extent.values[1], 1e-7, 1e-7);
                }
            });
    });

    it('leaves the unused result members at their default values', () => {
        // Upstream value-initializes the Line3 and Ellipse3 members, which
        // calls their default constructors (origin 0 / direction (1,0,0);
        // center 0, normal (0,0,1), axes (1,0,0) and (0,1,0), extents 1) --
        // not "all zero members" as the upstream Result comment claims.
        const d = defaultIntrPlane3Cylinder3FIResult();
        expect(d.line[0].origin.values).toEqual([0, 0, 0]);
        expect(d.line[0].direction.values).toEqual([1, 0, 0]);
        expect(d.ellipse.center.values).toEqual([0, 0, 0]);
        expect(d.ellipse.normal.values).toEqual([0, 0, 1]);
        expect(d.ellipse.extent.values).toEqual([1, 1]);
        expect(d.trimLine[1].direction.values).toEqual([1, 0, 0]);

        // A separated configuration leaves them untouched.
        const C = Cylinder3.fromAxisRadiusHeight(Line.fromOriginDirection(
            Vector.zero(3), Vector.fromArray([0, 0, 1])), 1, 2);
        const P = Hyperplane.fromNormalOrigin(Vector.fromArray([0, 0, 1]),
            Vector.fromArray([0, 0, 5]));
        const r = fiQ.find(P, C);
        expect(r.type).toBe(IntrPlane3Cylinder3FIResultType.noIntersection);
        expect(r.ellipse.extent.values).toEqual([1, 1]);
        expect(r.trimLine[0].direction.values).toEqual([1, 0, 0]);
    });
});
