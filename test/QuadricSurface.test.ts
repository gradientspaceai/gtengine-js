import { describe, expect, it } from 'vitest';
import { QuadricSurface, QuadricSurfaceClassification } from '../src/QuadricSurface.js';
import { Matrix } from '../src/Matrix.js';
import { Vector } from '../src/Vector.js';

const C = QuadricSurfaceClassification;

// A quadric from its 10 implicit coefficients
// q0 + q1*x + q2*y + q3*z + q4*x^2 + q5*x*y + q6*x*z + q7*y^2 + q8*y*z + q9*z^2.
function fromQ(q: readonly number[]): QuadricSurface {
    return QuadricSurface.fromCoefficients(q);
}

// The known quadrics used by several tests, given as (name, q, expected).
const KNOWN: Array<{ name: string, q: number[], expected: QuadricSurfaceClassification }> = [
    // A has three nonzero eigenvalues.
    { name: 'unit sphere', q: [-1, 0, 0, 0, 1, 0, 0, 1, 0, 1], expected: C.ELLIPSOID },
    { name: 'ellipsoid', q: [-36, 0, 0, 0, 9, 0, 0, 4, 0, 36], expected: C.ELLIPSOID },
    { name: 'imaginary ellipsoid', q: [1, 0, 0, 0, 1, 0, 0, 1, 0, 1], expected: C.NO_SOLUTION },
    { name: 'point', q: [0, 0, 0, 0, 1, 0, 0, 1, 0, 1], expected: C.POINT },
    { name: 'elliptic cone', q: [0, 0, 0, 0, 1, 0, 0, 1, 0, -1], expected: C.ELLIPTIC_CONE },
    { name: 'hyperboloid of one sheet', q: [-1, 0, 0, 0, 1, 0, 0, 1, 0, -1], expected: C.HYPERBOLOID_ONE_SHEET },
    { name: 'hyperboloid of two sheets', q: [1, 0, 0, 0, 1, 0, 0, 1, 0, -1], expected: C.HYPERBOLOID_TWO_SHEETS },
    // A has two nonzero eigenvalues.
    { name: 'elliptic cylinder', q: [-1, 0, 0, 0, 1, 0, 0, 4, 0, 0], expected: C.ELLIPTIC_CYLINDER },
    { name: 'imaginary cylinder', q: [1, 0, 0, 0, 1, 0, 0, 4, 0, 0], expected: C.NO_SOLUTION },
    { name: 'line', q: [0, 0, 0, 0, 1, 0, 0, 4, 0, 0], expected: C.LINE },
    { name: 'hyperbolic cylinder', q: [-1, 0, 0, 0, 1, 0, 0, -1, 0, 0], expected: C.HYPERBOLIC_CYLINDER },
    { name: 'two intersecting planes', q: [0, 0, 0, 0, 1, 0, 0, -1, 0, 0], expected: C.TWO_PLANES },
    { name: 'elliptic paraboloid', q: [0, 0, 0, -1, 1, 0, 0, 1, 0, 0], expected: C.ELLIPTIC_PARABOLOID },
    { name: 'hyperbolic paraboloid', q: [0, 0, 0, -1, 1, 0, 0, -1, 0, 0], expected: C.HYPERBOLIC_PARABOLOID },
    // A has one nonzero eigenvalue.
    { name: 'parabolic cylinder', q: [0, 0, -1, 0, 1, 0, 0, 0, 0, 0], expected: C.PARABOLIC_CYLINDER },
    { name: 'two parallel planes', q: [-1, 0, 0, 0, 1, 0, 0, 0, 0, 0], expected: C.TWO_PLANES },
    { name: 'two imaginary parallel planes', q: [1, 0, 0, 0, 1, 0, 0, 0, 0, 0], expected: C.NO_SOLUTION },
    { name: 'double plane', q: [0, 0, 0, 0, 1, 0, 0, 0, 0, 0], expected: C.PLANE },
    { name: 'negative double plane', q: [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], expected: C.PLANE },
    { name: 'negative two parallel planes', q: [1, 0, 0, 0, -1, 0, 0, 0, 0, 0], expected: C.TWO_PLANES },
    { name: 'negative imaginary planes', q: [-1, 0, 0, 0, -1, 0, 0, 0, 0, 0], expected: C.NO_SOLUTION },
    // A is zero.
    { name: 'plane', q: [-1, 1, 2, 3, 0, 0, 0, 0, 0, 0], expected: C.PLANE },
    { name: 'entire space', q: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], expected: C.ENTIRE_SPACE },
    { name: 'empty space', q: [5, 0, 0, 0, 0, 0, 0, 0, 0, 0], expected: C.NO_SOLUTION }
];

// Substitute X = M*Y + t into q and return the coefficients of the quadric in
// Y. M must be invertible; then the two quadrics are affinely equivalent and
// have the same classification.
function changeOfVariables(q: readonly number[], M: number[][],
    t: number[]): number[] {
    const A = [
        [q[4], 0.5 * q[5], 0.5 * q[6]],
        [0.5 * q[5], q[7], 0.5 * q[8]],
        [0.5 * q[6], 0.5 * q[8], q[9]]
    ];
    const b = [q[1], q[2], q[3]];
    const c = q[0];

    // First the translation: A stays, b -> 2*A*t + b, c -> t^T*A*t + b^T*t + c.
    const At = [0, 1, 2].map(i => A[i][0] * t[0] + A[i][1] * t[1] + A[i][2] * t[2]);
    const bT = [0, 1, 2].map(i => 2 * At[i] + b[i]);
    const cT = t[0] * At[0] + t[1] * At[1] + t[2] * At[2] +
        b[0] * t[0] + b[1] * t[1] + b[2] * t[2] + c;

    // Then the linear map: A -> M^T*A*M, b -> M^T*b, c unchanged.
    const AM = [0, 1, 2].map(i => [0, 1, 2].map(j =>
        A[i][0] * M[0][j] + A[i][1] * M[1][j] + A[i][2] * M[2][j]));
    const A2 = [0, 1, 2].map(i => [0, 1, 2].map(j =>
        M[0][i] * AM[0][j] + M[1][i] * AM[1][j] + M[2][i] * AM[2][j]));
    const b2 = [0, 1, 2].map(j => M[0][j] * bT[0] + M[1][j] * bT[1] + M[2][j] * bT[2]);

    return [cT, b2[0], b2[1], b2[2], A2[0][0], 2 * A2[0][1], 2 * A2[0][2],
        A2[1][1], 2 * A2[1][2], A2[2][2]];
}

describe('QuadricSurface', () => {
    it('default construction is the zero quadric', () => {
        const surface = new QuadricSurface();
        expect(surface.getC()).toBe(0);
        expect(surface.getA().values).toEqual(new Array<number>(9).fill(0));
        expect(surface.getB().values).toEqual([0, 0, 0]);
        expect(surface.getQ()).toEqual(new Array<number>(10).fill(0));
        expect(surface.getClassification()).toBe(C.ENTIRE_SPACE);
    });

    it('fromCoefficients builds (A, b, c) and getQ inverts it', () => {
        const q = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        const surface = fromQ(q);
        const A = surface.getA();
        expect(A.get(0, 0)).toBe(5);
        expect(A.get(0, 1)).toBe(3);
        expect(A.get(1, 0)).toBe(3);
        expect(A.get(0, 2)).toBe(3.5);
        expect(A.get(2, 0)).toBe(3.5);
        expect(A.get(1, 1)).toBe(8);
        expect(A.get(1, 2)).toBe(4.5);
        expect(A.get(2, 1)).toBe(4.5);
        expect(A.get(2, 2)).toBe(10);
        expect(surface.getB().values).toEqual([2, 3, 4]);
        expect(surface.getC()).toBe(1);
        expect(surface.getQ()).toEqual(q);
    });

    it('fromMatrix stores copies of the inputs', () => {
        const A = Matrix.fromArray(3, 3, [1, 0, 0, 0, 2, 0, 0, 0, 3]);
        const b = Vector.fromArray([4, 5, 6]);
        const surface = QuadricSurface.fromMatrix(A, b, 7);
        A.set(0, 0, 100);
        b.values[0] = 100;
        expect(surface.getA().get(0, 0)).toBe(1);
        expect(surface.getB().values[0]).toBe(4);
        expect(surface.getC()).toBe(7);
        expect(surface.getQ()).toEqual([7, 4, 5, 6, 1, 0, 0, 2, 0, 3]);
    });

    it('rejects inputs of the wrong size', () => {
        expect(() => QuadricSurface.fromCoefficients([1, 2, 3])).toThrow('Invalid size.');
        expect(() => QuadricSurface.fromMatrix(new Matrix(2, 2),
            Vector.fromArray([1, 2]), 0)).toThrow('Invalid size.');
        expect(() => QuadricSurface.fromMatrix(new Matrix(3, 3),
            Vector.fromArray([1, 2]), 0)).toThrow('Invalid size.');
    });

    it('evaluates the function and its derivatives at known values', () => {
        // 0 = -14 + x + 2*y + 3*z + 4*x^2 + 5*x*y + 6*x*z + 7*y^2 + 8*y*z + 9*z^2
        const q = [-14, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        const surface = fromQ(q);
        const p = Vector.fromArray([1, -2, 3]);
        const x = 1, y = -2, z = 3;
        const expectedF = q[0] + q[1] * x + q[2] * y + q[3] * z + q[4] * x * x +
            q[5] * x * y + q[6] * x * z + q[7] * y * y + q[8] * y * z + q[9] * z * z;
        expect(surface.f(p)).toBeCloseTo(expectedF, 12);

        // Gradient by hand.
        expect(surface.fx(p)).toBeCloseTo(q[1] + 2 * q[4] * x + q[5] * y + q[6] * z, 12);
        expect(surface.fy(p)).toBeCloseTo(q[2] + q[5] * x + 2 * q[7] * y + q[8] * z, 12);
        expect(surface.fz(p)).toBeCloseTo(q[3] + q[6] * x + q[8] * y + 2 * q[9] * z, 12);

        // Hessian is constant.
        expect(surface.fxx()).toBe(2 * q[4]);
        expect(surface.fxy()).toBe(q[5]);
        expect(surface.fxz()).toBe(q[6]);
        expect(surface.fyy()).toBe(2 * q[7]);
        expect(surface.fyz()).toBe(q[8]);
        expect(surface.fzz()).toBe(2 * q[9]);
    });

    it('has a gradient matching central differences at random points', () => {
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const h = 1e-5;
        for (let trial = 0; trial < 200; ++trial) {
            const q = new Array<number>(10).fill(0).map(() => 4 * rand() - 2);
            const surface = fromQ(q);
            const p = Vector.fromArray([4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2]);
            const grad = [surface.fx(p), surface.fy(p), surface.fz(p)];
            for (let k = 0; k < 3; ++k) {
                const pp = p.clone(), pm = p.clone();
                pp.values[k] += h;
                pm.values[k] -= h;
                const fd = (surface.f(pp) - surface.f(pm)) / (2 * h);
                expect(fd).toBeCloseTo(grad[k], 6);
            }
        }
    });

    it('f vanishes on constructed points of a sphere, cylinder and cone', () => {
        // Sphere of radius 3 centered at (1,2,3).
        const sphere = fromQ([1 + 4 + 9 - 9, -2, -4, -6, 1, 0, 0, 1, 0, 1]);
        for (const p of [[4, 2, 3], [1, 5, 3], [1, 2, 6], [-2, 2, 3]]) {
            expect(sphere.f(Vector.fromArray(p))).toBeCloseTo(0, 12);
        }
        expect(sphere.getClassification()).toBe(C.ELLIPSOID);

        // Circular cylinder x^2 + y^2 = 4 about the z axis.
        const cylinder = fromQ([-4, 0, 0, 0, 1, 0, 0, 1, 0, 0]);
        for (const p of [[2, 0, 0], [0, 2, 5], [-2, 0, -3]]) {
            expect(cylinder.f(Vector.fromArray(p))).toBeCloseTo(0, 12);
        }
        expect(cylinder.getClassification()).toBe(C.ELLIPTIC_CYLINDER);

        // Cone x^2 + y^2 = z^2.
        const cone = fromQ([0, 0, 0, 0, 1, 0, 0, 1, 0, -1]);
        for (const p of [[1, 0, 1], [0, -3, 3], [3, 4, -5]]) {
            expect(cone.f(Vector.fromArray(p))).toBeCloseTo(0, 12);
        }
        expect(cone.getClassification()).toBe(C.ELLIPTIC_CONE);
    });

    it('classifies the known quadrics', () => {
        for (const testCase of KNOWN) {
            const actual = fromQ(testCase.q).getClassification();
            expect(`${testCase.name}: ${C[actual]}`)
                .toBe(`${testCase.name}: ${C[testCase.expected]}`);
        }
    });

    it('caches the classification', () => {
        const surface = fromQ([-1, 0, 0, 0, 1, 0, 0, 1, 0, 1]);
        expect(surface.getClassification()).toBe(C.ELLIPSOID);
        expect(surface.getClassification()).toBe(C.ELLIPSOID);
    });

    it('classifies quadrics in rotated frames the same way', () => {
        // An exact orthogonal frame scaled by 5, built from the (3,4,5)
        // Pythagorean triple so that all entries are integers.
        const rotation = [
            [3, -4, 0],
            [4, 3, 0],
            [0, 0, 5]
        ];
        const translation = [2, -3, 1];
        for (const testCase of KNOWN) {
            const q = changeOfVariables(testCase.q, rotation, translation);
            const actual = fromQ(q).getClassification();
            expect(`${testCase.name}: ${C[actual]}`)
                .toBe(`${testCase.name}: ${C[testCase.expected]}`);
        }
    });

    it('classifies quadrics after invertible affine changes of variables', () => {
        // Affine equivalence preserves every classification category, so the
        // classification must be invariant under these substitutions. The
        // matrices are integer valued with determinant 1 or -1.
        const maps: Array<{ M: number[][], t: number[] }> = [
            { M: [[1, 1, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0] },
            { M: [[1, 0, 0], [2, 1, 0], [-1, 3, 1]], t: [1, 1, 1] },
            { M: [[0, 0, 1], [1, 0, 0], [0, 1, 0]], t: [-2, 5, 7] },
            { M: [[2, 1, 1], [1, 1, 0], [1, 0, 0]], t: [0, -1, 4] },
            { M: [[1, 0, 0], [0, 1, 0], [0, 0, -1]], t: [3, 0, -3] }
        ];
        for (const testCase of KNOWN) {
            for (const { M, t } of maps) {
                const q = changeOfVariables(testCase.q, M, t);
                const actual = fromQ(q).getClassification();
                expect(`${testCase.name}: ${C[actual]}`)
                    .toBe(`${testCase.name}: ${C[testCase.expected]}`);
            }
        }
    });

    it('classifies a paraboloid, cylinder and cone in a general position', () => {
        // z = x^2 + y^2 rotated so the axis is (1,1,1)/sqrt(3) is not exactly
        // representable, so instead verify the classification of quadrics
        // whose axes are the integer directions of an orthogonal frame.
        // u = (1,1,0), v = (1,-1,0), w = (0,0,1) with |u|^2 = |v|^2 = 2.
        // f = (u.X)^2 + (v.X)^2 - 2*(w.X) = 2*x^2 + 2*y^2 - 2*z.
        const paraboloid = fromQ([0, 0, 0, -2, 2, 0, 0, 2, 0, 0]);
        expect(paraboloid.getClassification()).toBe(C.ELLIPTIC_PARABOLOID);
        // f = (u.X)^2 - (v.X)^2 - 2*(w.X) = 4*x*y - 2*z.
        const saddle = fromQ([0, 0, 0, -2, 0, 4, 0, 0, 0, 0]);
        expect(saddle.getClassification()).toBe(C.HYPERBOLIC_PARABOLOID);
        // f = (u.X)^2 - (w.X)^2 = x^2 + 2*x*y + y^2 - z^2, rank 2 with signs
        // (+,-), b = 0, so the solution set is two intersecting planes.
        const twoPlanes = fromQ([0, 0, 0, 0, 1, 2, 0, 1, 0, -1]);
        expect(twoPlanes.getClassification()).toBe(C.TWO_PLANES);
        // The same with a nonzero constant is a hyperbolic cylinder.
        const hypCylinder = fromQ([-1, 0, 0, 0, 1, 2, 0, 1, 0, -1]);
        expect(hypCylinder.getClassification()).toBe(C.HYPERBOLIC_CYLINDER);
    });

    it('handles a rank-1 A whose first two rows are zero', () => {
        // Only the (2,2) entry of A is nonzero: z^2 - 1 = 0 is two planes,
        // z^2 = 0 is one plane, z^2 + y = 0 is a parabolic cylinder.
        expect(fromQ([-1, 0, 0, 0, 0, 0, 0, 0, 0, 1]).getClassification())
            .toBe(C.TWO_PLANES);
        expect(fromQ([0, 0, 0, 0, 0, 0, 0, 0, 0, 1]).getClassification())
            .toBe(C.PLANE);
        expect(fromQ([0, 0, 1, 0, 0, 0, 0, 0, 0, 1]).getClassification())
            .toBe(C.PARABOLIC_CYLINDER);
        // Only the (1,1) entry of A is nonzero.
        expect(fromQ([-4, 0, 0, 0, 0, 0, 0, 1, 0, 0]).getClassification())
            .toBe(C.TWO_PLANES);
        expect(fromQ([0, 0, 0, 1, 0, 0, 0, 1, 0, 0]).getClassification())
            .toBe(C.PARABOLIC_CYLINDER);
    });

    it('handles a rank-2 A whose first row is zero', () => {
        // y^2 + z^2 - 1 = 0 is a circular cylinder about the x axis.
        expect(fromQ([-1, 0, 0, 0, 0, 0, 0, 1, 0, 1]).getClassification())
            .toBe(C.ELLIPTIC_CYLINDER);
        // y^2 - z^2 - 1 = 0 is a hyperbolic cylinder about the x axis.
        expect(fromQ([-1, 0, 0, 0, 0, 0, 0, 1, 0, -1]).getClassification())
            .toBe(C.HYPERBOLIC_CYLINDER);
        // y^2 + z^2 - x = 0 is an elliptic paraboloid about the x axis.
        expect(fromQ([0, -1, 0, 0, 0, 0, 0, 1, 0, 1]).getClassification())
            .toBe(C.ELLIPTIC_PARABOLOID);
    });

    it('classifies exactly where a floating-point test would be ambiguous', () => {
        // The cone x^2 + y^2 - z^2 = 0 and the nearby hyperboloids differ
        // only in the last bits of the constant term.
        const eps = Number.EPSILON;
        expect(fromQ([0, 0, 0, 0, 1, 0, 0, 1, 0, -1]).getClassification())
            .toBe(C.ELLIPTIC_CONE);
        expect(fromQ([-eps, 0, 0, 0, 1, 0, 0, 1, 0, -1]).getClassification())
            .toBe(C.HYPERBOLOID_ONE_SHEET);
        expect(fromQ([eps, 0, 0, 0, 1, 0, 0, 1, 0, -1]).getClassification())
            .toBe(C.HYPERBOLOID_TWO_SHEETS);

        // A sphere of radius sqrt(eps) is still an ellipsoid, not a point.
        expect(fromQ([-eps, 0, 0, 0, 1, 0, 0, 1, 0, 1]).getClassification())
            .toBe(C.ELLIPSOID);
        expect(fromQ([0, 0, 0, 0, 1, 0, 0, 1, 0, 1]).getClassification())
            .toBe(C.POINT);
    });
});
