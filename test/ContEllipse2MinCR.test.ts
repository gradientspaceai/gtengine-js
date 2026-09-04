import { describe, it, expect } from 'vitest';
import { getContainerEllipse2MinCR } from '../src/ContEllipse2MinCR.js';
import { Matrix } from '../src/Matrix.js';
import { Vector, add, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, fc, latticeVector, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';

function v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// The rotation matrix whose columns are the ellipse axes.
function rot(angle: number): Matrix {
    const c = Math.cos(angle), s = Math.sin(angle);
    return Matrix.fromArray(2, 2, [c, -s, s, c]);
}

// The quadratic form (X-C)^T R D R^T (X-C).
function form(X: Vector, C: Vector, R: Matrix, D: readonly number[]): number {
    const d = sub(X, C);
    let sum = 0;
    for (let j = 0; j < 2; ++j) {
        let u = 0;
        for (let i = 0; i < 2; ++i) {
            u += R.get(i, j) * d.values[i];
        }
        sum += D[j] * u * u;
    }
    return sum;
}

// The constraint coefficients A[i] = (u^2, v^2) with (u,v) = R^T*(P[i]-C).
function constraints(points: readonly Vector[], C: Vector,
    R: Matrix): number[][] {
    return points.map(P => {
        const d = sub(P, C);
        const a: number[] = [];
        for (let j = 0; j < 2; ++j) {
            let u = 0;
            for (let i = 0; i < 2; ++i) {
                u += R.get(i, j) * d.values[i];
            }
            a.push(u * u);
        }
        return a;
    });
}

// An independent solver for the same problem: maximize D[0]*D[1] subject to
// A[i][0]*D[0] + A[i][1]*D[1] <= 1 and D >= 0. The objective is strictly
// increasing in each variable, so the maximizer lies on the boundary of the
// polytope with either one or two active constraints. This enumerates every
// such candidate (tangency point of a single line, intersection of a pair of
// lines), discards the infeasible ones and keeps the largest product. It
// shares no algorithmic idea with the ported hull walk.
function maxProductByEnumeration(A: readonly number[][]): number[] {
    const feasible = (D: readonly number[]): boolean => {
        if (D[0] < 0 || D[1] < 0 || !isFinite(D[0]) || !isFinite(D[1])) {
            return false;
        }
        for (const a of A) {
            if (a[0] * D[0] + a[1] * D[1] > 1 + 1e-9) {
                return false;
            }
        }
        return true;
    };

    let best = [0, 0];
    const consider = (D: number[]): void => {
        if (feasible(D) && D[0] * D[1] > best[0] * best[1]) {
            best = D;
        }
    };

    for (let i = 0; i < A.length; ++i) {
        // Maximum of x*y on the single line A[i][0]*x + A[i][1]*y = 1.
        if (A[i][0] > 0 && A[i][1] > 0) {
            consider([0.5 / A[i][0], 0.5 / A[i][1]]);
        }
        for (let j = i + 1; j < A.length; ++j) {
            // Intersection of the two constraint lines.
            const det = A[i][0] * A[j][1] - A[j][0] * A[i][1];
            if (det !== 0) {
                consider([(A[j][1] - A[i][1]) / det,
                    (A[i][0] - A[j][0]) / det]);
            }
        }
    }
    return best;
}

function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1103515245) + 12345) >>> 0;
        return s / 4294967296;
    };
}

describe('getContainerEllipse2MinCR', () => {
    it('recovers the extents of an axis-aligned ellipse sample', () => {
        // Points sampled on x^2/9 + y^2/4 = 1, including the four vertices.
        // The constraints from (3,0) and (0,2) bound D[0] <= 1/9 and
        // D[1] <= 1/4, and D = (1/9, 1/4) satisfies every other constraint
        // with equality, so it is the unique maximizer of D[0]*D[1].
        const a = 3, b = 2;
        const points: Vector[] = [];
        for (let i = 0; i < 32; ++i) {
            const t = 2 * Math.PI * i / 32;
            points.push(v(a * Math.cos(t), b * Math.sin(t)));
        }
        const C = v(0, 0);
        const R = rot(0);
        const D = getContainerEllipse2MinCR(points, C, R);
        expect(D[0]).toBeCloseTo(1 / (a * a), 12);
        expect(D[1]).toBeCloseTo(1 / (b * b), 12);
    });

    it('recovers the extents in a rotated, translated frame', () => {
        const a = 5, b = 0.75, angle = 0.7, C = v(-2, 3);
        const R = rot(angle);
        const points: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            const t = 2 * Math.PI * i / 24;
            const u = a * Math.cos(t), w = b * Math.sin(t);
            // X = C + u*axis0 + w*axis1.
            points.push(v(
                C.values[0] + u * R.get(0, 0) + w * R.get(0, 1),
                C.values[1] + u * R.get(1, 0) + w * R.get(1, 1)));
        }
        const D = getContainerEllipse2MinCR(points, C, R);
        expect(D[0]).toBeCloseTo(1 / (a * a), 12);
        expect(D[1]).toBeCloseTo(1 / (b * b), 12);
        for (const P of points) {
            expect(form(P, C, R, D)).toBeLessThanOrEqual(1 + 1e-10);
        }
    });

    it('handles a single point (tangency at the midpoint of the line)', () => {
        // One constraint a0*x+a1*y = 1; the maximum of x*y on it is at
        // x = 1/(2*a0), y = 1/(2*a1).
        const C = v(1, 1);
        const R = rot(0);
        const P = v(4, 3);  // (u,v) = (3,2), A = (9,4)
        const D = getContainerEllipse2MinCR([P], C, R);
        expect(D[0]).toBeCloseTo(0.5 / 9, 15);
        expect(D[1]).toBeCloseTo(0.5 / 4, 15);
        // The point is on the boundary of the ellipse.
        expect(form(P, C, R, D)).toBeCloseTo(1, 12);
    });

    it('is bounded by the two axis-extreme points of a box sample', () => {
        // Points (+-3, +-2): a single distinct constraint (9,4) survives the
        // redundant-constraint elimination.
        const points = [v(3, 2), v(-3, 2), v(3, -2), v(-3, -2)];
        const D = getContainerEllipse2MinCR(points, v(0, 0), rot(0));
        expect(D[0]).toBeCloseTo(0.5 / 9, 15);
        expect(D[1]).toBeCloseTo(0.5 / 4, 15);
        for (const P of points) {
            expect(form(P, v(0, 0), rot(0), D)).toBeCloseTo(1, 12);
        }
    });

    it('eliminates redundant constraints (interior points do not matter)',
        () => {
            const outer = [v(3, 0), v(0, 2), v(2, 1.4), v(-2.5, -1)];
            const D0 = getContainerEllipse2MinCR(outer, v(0, 0), rot(0));
            const withInterior = outer.concat(
                [v(0.1, 0.2), v(-0.5, 0.3), v(1, 0.1), v(0, 0)]);
            const D1 = getContainerEllipse2MinCR(withInterior, v(0, 0),
                rot(0));
            expect(D1[0]).toBeCloseTo(D0[0], 14);
            expect(D1[1]).toBeCloseTo(D0[1], 14);
        });

    it('matches an independent exhaustive-candidate solver on random clouds',
        () => {
            const rnd = makeRandom(12345);
            for (let trial = 0; trial < 40; ++trial) {
                const angle = 2 * Math.PI * rnd();
                const R = rot(angle);
                const C = v(4 * rnd() - 2, 4 * rnd() - 2);
                const numPoints = 3 + Math.floor(12 * rnd());
                const points: Vector[] = [];
                for (let i = 0; i < numPoints; ++i) {
                    points.push(v(C.values[0] + 6 * rnd() - 3,
                        C.values[1] + 6 * rnd() - 3));
                }
                // Guarantee both constraint components are positive so the
                // problem is bounded.
                points.push(v(C.values[0] + 1.5, C.values[1] + 1.1));
                points.push(v(C.values[0] - 1.3, C.values[1] - 1.7));

                const D = getContainerEllipse2MinCR(points, C, R);
                expect(D[0]).toBeGreaterThan(0);
                expect(D[1]).toBeGreaterThan(0);

                // Every point is inside or on the ellipse.
                for (const P of points) {
                    expect(form(P, C, R, D)).toBeLessThanOrEqual(1 + 1e-10);
                }

                // At least one constraint is active (the fit is tight).
                let maxForm = 0;
                for (const P of points) {
                    maxForm = Math.max(maxForm, form(P, C, R, D));
                }
                expect(maxForm).toBeGreaterThan(1 - 1e-8);

                // The product matches the independent maximizer.
                const A = constraints(points, C, R);
                const E = maxProductByEnumeration(A);
                expect(D[0] * D[1]).toBeGreaterThan(E[0] * E[1] * (1 - 1e-9));
                expect(D[0] * D[1]).toBeLessThan(E[0] * E[1] * (1 + 1e-9));
            }
        });

    it('handles a point on the second ellipse axis (vertical constraint)',
        () => {
            // The point (3,0) yields the constraint 9*D[0] <= 1, whose line
            // is vertical (its D[1] coefficient is 0). Upstream divides by
            // that zero coefficient and returns NaN; see the port note in
            // ContEllipse2MinCR.ts.
            const points = [v(3, 0), v(0, 2), v(2, 1.4), v(-2.5, -1)];
            const C = v(0, 0), R = rot(0);
            const D = getContainerEllipse2MinCR(points, C, R);
            expect(Number.isNaN(D[0])).toBe(false);
            expect(Number.isNaN(D[1])).toBe(false);
            expect(D[0]).toBeCloseTo(1 / 9, 14);
            expect(D[1]).toBeCloseTo(1 / 4, 14);
            const E = maxProductByEnumeration(constraints(points, C, R));
            expect(D[0] * D[1]).toBeCloseTo(E[0] * E[1], 14);
            for (const P of points) {
                expect(form(P, C, R, D)).toBeLessThanOrEqual(1 + 1e-10);
            }
        });

    it('throws when the problem is unbounded along an axis', () => {
        // All points on the first axis: every constraint has A[i][1] = 0, so
        // D[1] is unbounded and upstream's LogAssert on iYMin fires.
        const points = [v(1, 0), v(-2, 0), v(3, 0)];
        expect(() => getContainerEllipse2MinCR(points, v(0, 0), rot(0)))
            .toThrow('Unexpected condition.');
    });

    it('throws when every point is the center', () => {
        const points = [v(1, 1), v(1, 1)];
        expect(() => getContainerEllipse2MinCR(points, v(1, 1), rot(0)))
            .toThrow('Unexpected condition.');
    });

    it('throws for an empty point set', () => {
        expect(() => getContainerEllipse2MinCR([], v(0, 0), rot(0)))
            .toThrow('no points');
    });

    it('validates the dimensions of its inputs', () => {
        expect(() => getContainerEllipse2MinCR([v(1, 1)],
            Vector.fromArray([0, 0, 0]), rot(0))).toThrow('2D');
        expect(() => getContainerEllipse2MinCR([Vector.fromArray([1, 1, 1])],
            v(0, 0), rot(0))).toThrow('2D');
        expect(() => getContainerEllipse2MinCR([v(1, 1)], v(0, 0),
            new Matrix(3, 3))).toThrow('2x2');
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContEllipse2MinCR.h semantics.
// ---------------------------------------------------------------------------

describe('ContEllipse2MinCR verification', () => {
    // Lattice clouds keep the constraint coefficients A[i] = (u^2, v^2)
    // exact for the identity frame, which keeps the comparison against the
    // enumeration solver meaningful. The problem must be bounded in both
    // variables, so at least one point must be off each ellipse axis.
    // Lattice data in the identity frame: the constraint coefficients
    // A[i] = ((x-cx)^2, (y-cy)^2) are then small exact integers, so no
    // constraint line is nearly parallel to another (their determinants are
    // nonzero integers) and no coefficient is a tiny nonzero number. That
    // conditioning matters: the upstream hull walk evaluates
    // y0 = (1 - a0*x0)/b0 on the final line, which cancels catastrophically
    // when b0 is a tiny nonzero number, and the products compared below then
    // lose most of their significant digits. (That is a property of the
    // upstream algorithm, not of the port; the containment property below
    // uses arbitrary frames and still holds.)
    const problem = fc.tuple(
        fc.array(latticeVector(2, -6, 6), { minLength: 1, maxLength: 8 }),
        latticeVector(2, -3, 3))
        .map(([points, C]) => ({
            points, C,
            R: Matrix.fromArray(2, 2, [1, 0, 0, 1])
        }))
        .filter(({ points, C }) =>
            // Bounded in both variables: some point off each ellipse axis.
            points.some(p => p.get(0) !== C.get(0))
            && points.some(p => p.get(1) !== C.get(1)));

    // Arbitrary center and frame; used only for properties that do not
    // compare two floating-point computations of the same optimum.
    const rotatedProblem = fc.tuple(
        fc.array(latticeVector(2, -6, 6), { minLength: 1, maxLength: 8 }),
        rotationFrame(2), wellScaledVector(2, -3, 3))
        .map(([points, frame, C]) => ({
            points, C,
            R: Matrix.fromArray(2, 2, [
                frame[0].get(0), frame[1].get(0),
                frame[0].get(1), frame[1].get(1)])
        }))
        .filter(({ points, C, R }) => {
            // Bounded in both variables, and no coefficient that is tiny but
            // nonzero. A point that lies (nearly) on an ellipse axis gives a
            // constraint line that is (nearly) vertical or horizontal, and
            // the upstream walk then evaluates y0 = (1 - a0*x0)/b0 with
            // a0*x0 within rounding of 1 and b0 within rounding of 0: the
            // answer loses every significant digit. That is a property of
            // the algorithm, not of the port (the exactly-degenerate case is
            // the one upstream divides by zero on, covered by the #234 test
            // below).
            const A = constraints(points, C, R);
            const maxAll = Math.max(...A.map(a => Math.max(a[0], a[1])));
            return maxAll > 0
                && A.every(a => a.every(x => x === 0 || x >= 1e-4 * maxAll))
                && A.some(a => a[0] > 0) && A.some(a => a[1] > 0);
        });

    // The design claim: every input point satisfies the quadratic form. The
    // hull walk ends on an active constraint, so the tightest point evaluates
    // to 1 up to rounding of the constraint coefficients.
    it('every input point is inside the resulting ellipse', () => {
        check(rotatedProblem, ({ points, C, R }:
            { points: Vector[], C: Vector, R: Matrix }) => {
            const D = getContainerEllipse2MinCR(points, C, R);
            expect(D[0]).toBeGreaterThanOrEqual(0);
            expect(D[1]).toBeGreaterThanOrEqual(0);
            for (const P of points) {
                expect(form(P, C, R, D)).toBeLessThanOrEqual(1 + 1e-9);
            }
            // Tight: at least one point is on the ellipse.
            expect(points.some(P => Math.abs(form(P, C, R, D) - 1) <= 1e-9))
                .toBe(true);
        });
    });

    // The product D[0]*D[1] is the maximum of the linear program, which the
    // independent enumeration solver finds by a completely different route.
    it('attains the maximum product found by enumeration', () => {
        check(problem, ({ points, C, R }:
            { points: Vector[], C: Vector, R: Matrix }) => {
            const D = getContainerEllipse2MinCR(points, C, R);
            const A = constraints(points, C, R);
            const best = maxProductByEnumeration(A);
            const got = D[0] * D[1];
            const want = best[0] * best[1];
            // The enumeration solver accepts candidates that violate a
            // constraint by up to 1e-9, so its optimum can differ from the
            // exact one by a comparable relative amount; both routines also
            // form the same intersection points by different expressions.
            expect(got).toBeGreaterThanOrEqual(want * (1 - 1e-9) - 1e-12);
            expect(got).toBeLessThanOrEqual(want * (1 + 1e-9) + 1e-12);
        });
    });

    // Transformations that leave the constraint coefficients
    // A[i] = ((u,v) of P[i]-C, squared) bit-identical must leave D
    // bit-identical: translating the points and the center together by an
    // integer vector, reflecting a coordinate of both, and negating a column
    // of R (an ellipse axis). A general rotation is *not* such a
    // transformation: it perturbs the coefficients in the last bits, and a
    // point that lands near an ellipse axis then produces a nearly vertical
    // constraint line, where the final y0 = (1 - a0*x0)/b0 of the walk
    // cancels catastrophically. Rotational equivariance is therefore not a
    // property of this algorithm.
    it('is invariant under exact translations, reflections and axis flips',
        () => {
            check(fc.tuple(problem, latticeVector(2, -20, 20),
                fc.integer({ min: 0, max: 1 })),
                ([{ points, C, R }, t, axis]:
                    [{ points: Vector[], C: Vector, R: Matrix }, Vector,
                        number]) => {
                    const D0 = getContainerEllipse2MinCR(points, C, R);

                    const shifted = getContainerEllipse2MinCR(
                        points.map(P => add(P, t)), add(C, t), R);
                    expect(shifted[0]).toBe(D0[0]);
                    expect(shifted[1]).toBe(D0[1]);

                    const flip = (p: Vector): Vector => {
                        const q = p.clone();
                        q.set(axis, -q.get(axis));
                        return q;
                    };
                    const mirrored = getContainerEllipse2MinCR(
                        points.map(flip), flip(C), R);
                    expect(mirrored[0]).toBe(D0[0]);
                    expect(mirrored[1]).toBe(D0[1]);

                    // Negating one column of R negates that component of
                    // R^T*(P-C), which the squaring undoes.
                    const R2 = R.clone();
                    for (let r = 0; r < 2; ++r) {
                        R2.set(r, axis, -R.get(r, axis));
                    }
                    const flippedAxis = getContainerEllipse2MinCR(points, C, R2);
                    expect(flippedAxis[0]).toBe(D0[0]);
                    expect(flippedAxis[1]).toBe(D0[1]);
                });
        });

    // Scaling the data about the center by a power of two scales the
    // constraint coefficients by s^2 exactly, so D scales by 1/s^2 exactly.
    it('scales as 1/s^2 when the data is scaled by s about the center', () => {
        check(fc.tuple(problem, fc.integer({ min: -4, max: 4 })),
            ([{ points, C, R }, e]:
                [{ points: Vector[], C: Vector, R: Matrix }, number]) => {
                const s = Math.pow(2, e);
                const scaled = points.map(P => add(C, mul(s, sub(P, C))));
                const D0 = getContainerEllipse2MinCR(points, C, R);
                const D1 = getContainerEllipse2MinCR(scaled, C, R);
                expectClose(D1[0] * s * s, D0[0], 1e-12, 1e-9);
                expectClose(D1[1] * s * s, D0[1], 1e-12, 1e-9);
            });
    });

    // Regression pin for upstream issue #234: a point on the second ellipse
    // axis produces a vertical constraint (b = 0) that the hull walk can step
    // onto, where upstream evaluates 0/0. The concrete case from the issue.
    it('does not produce NaN on a vertical constraint (#234)', () => {
        const C = v(0, 0);
        const R = Matrix.fromArray(2, 2, [1, 0, 0, 1]);
        const points = [v(3, 0), v(0, 2), v(2, 1.4), v(-2.5, -1)];
        const D = getContainerEllipse2MinCR(points, C, R);
        expect(Number.isFinite(D[0])).toBe(true);
        expect(Number.isFinite(D[1])).toBe(true);
        expectClose(D[0], 1 / 9, 1e-12, 1e-12);
        expectClose(D[1], 1 / 4, 1e-12, 1e-12);
        for (const P of points) {
            expect(form(P, C, R, D)).toBeLessThanOrEqual(1 + 1e-12);
        }
    });
});
