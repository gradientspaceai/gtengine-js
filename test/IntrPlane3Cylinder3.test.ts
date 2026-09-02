import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3';
import { Hyperplane } from '../src/Hyperplane';
import {
    IntrPlane3Cylinder3TI,
    IntrPlane3Cylinder3FI,
    IntrPlane3Cylinder3FIResultType,
    defaultIntrPlane3Cylinder3FIResult
} from '../src/IntrPlane3Cylinder3';
import { Line } from '../src/Line';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector';

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
