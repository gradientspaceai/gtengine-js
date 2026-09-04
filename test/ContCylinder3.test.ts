import { describe, it, expect } from 'vitest';
import {
    getContainerCylinder3,
    inContainerCylinder3
} from '../src/ContCylinder3.js';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, dot, length, mul, sub } from '../src/Vector.js';
import { add } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import {
    check, expectClose, fc, rotationFrame, seededRandom, unitVector,
    wellScaledVector
} from './helpers/arbitraries.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Independent containment check with a tolerance, since points on the
// cylinder boundary satisfy the upstream strict tests only up to rounding.
function containedWithin(point: Vector, cylinder: Cylinder3,
    tol: number): boolean {
    const diff = sub(point, cylinder.axis.origin);
    const zProj = dot(diff, cylinder.axis.direction);
    if (Math.abs(zProj) * 2 > cylinder.height + tol) {
        return false;
    }
    const xyProj = sub(diff, mul(zProj, cylinder.axis.direction));
    return length(xyProj) <= cylinder.radius + tol;
}

function makeRandom(seed: number): () => number {
    let s = seed;
    return () => {
        s = (s * 1103515245 + 12345) % 2147483648;
        return s / 2147483648;
    };
}

describe('getContainerCylinder3', () => {
    it('recovers a cylinder from points on its lateral surface', () => {
        const points: Vector[] = [];
        for (const z of [-2, 2]) {
            for (let k = 0; k < 16; ++k) {
                const a = (2 * Math.PI * k) / 16;
                points.push(v(2 * Math.cos(a), 2 * Math.sin(a), z));
            }
        }

        const cylinder = getContainerCylinder3(points);
        // The dominant variance is along z, so the fitted axis is +/- e2.
        expect(Math.abs(cylinder.axis.direction.values[2])).toBeCloseTo(1, 10);
        expect(cylinder.axis.direction.values[0]).toBeCloseTo(0, 10);
        expect(cylinder.axis.direction.values[1]).toBeCloseTo(0, 10);
        expect(cylinder.axis.origin.values[0]).toBeCloseTo(0, 10);
        expect(cylinder.axis.origin.values[1]).toBeCloseTo(0, 10);
        expect(cylinder.axis.origin.values[2]).toBeCloseTo(0, 10);
        expect(cylinder.radius).toBeCloseTo(2, 10);
        expect(cylinder.height).toBeCloseTo(4, 10);
    });

    it('contains every input point (random clouds)', () => {
        const rand = makeRandom(9871);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 40; ++i) {
                // An elongated cloud so the least-squares axis is stable.
                points.push(v(
                    10 * (rand() - 0.5),
                    2 * (rand() - 0.5),
                    2 * (rand() - 0.5)));
            }
            const cylinder = getContainerCylinder3(points);
            expect(cylinder.radius).toBeGreaterThan(0);
            expect(cylinder.height).toBeGreaterThan(0);
            for (const p of points) {
                expect(containedWithin(p, cylinder, 1e-9)).toBe(true);
            }
        }
    });

    it('handles collinear points (zero radius)', () => {
        const points = [v(0, 0, 0), v(1, 1, 1), v(2, 2, 2), v(-1, -1, -1)];
        const cylinder = getContainerCylinder3(points);
        expect(cylinder.radius).toBeCloseTo(0, 10);
        expect(cylinder.height).toBeCloseTo(3 * Math.sqrt(3), 10);
        expect(cylinder.axis.origin.values[0]).toBeCloseTo(0.5, 10);
        expect(cylinder.axis.origin.values[1]).toBeCloseTo(0.5, 10);
        expect(cylinder.axis.origin.values[2]).toBeCloseTo(0.5, 10);
    });

    it('handles coincident points (degenerate cylinder)', () => {
        const points = [v(3, -1, 2), v(3, -1, 2), v(3, -1, 2)];
        const cylinder = getContainerCylinder3(points);
        expect(cylinder.radius).toBeCloseTo(0, 12);
        expect(cylinder.height).toBeCloseTo(0, 12);
        expect(cylinder.axis.origin.values[0]).toBeCloseTo(3, 12);
        expect(cylinder.axis.origin.values[1]).toBeCloseTo(-1, 12);
        expect(cylinder.axis.origin.values[2]).toBeCloseTo(2, 12);
    });

    it('throws for an empty point set', () => {
        expect(() => getContainerCylinder3([])).toThrow();
    });

    it('throws for non-3D points', () => {
        expect(() => getContainerCylinder3([Vector.fromArray([1, 2])]))
            .toThrow();
    });
});

describe('inContainerCylinder3', () => {
    const cylinder = Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(v(0, 0, 0), v(0, 0, 1)), 2, 4);

    it('accepts interior points', () => {
        expect(inContainerCylinder3(v(0, 0, 0), cylinder)).toBe(true);
        expect(inContainerCylinder3(v(1, 1, 1), cylinder)).toBe(true);
    });

    it('accepts boundary points (the boundary is part of the cylinder)', () => {
        expect(inContainerCylinder3(v(2, 0, 0), cylinder)).toBe(true);
        expect(inContainerCylinder3(v(0, 0, 2), cylinder)).toBe(true);
        expect(inContainerCylinder3(v(0, 2, -2), cylinder)).toBe(true);
    });

    it('rejects points outside the radius or the height', () => {
        expect(inContainerCylinder3(v(2.0001, 0, 0), cylinder)).toBe(false);
        expect(inContainerCylinder3(v(0, 0, 2.0001), cylinder)).toBe(false);
        expect(inContainerCylinder3(v(0, 0, -2.0001), cylinder)).toBe(false);
    });

    it('works for an oblique axis', () => {
        const d = v(1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3));
        const oblique = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(v(1, 2, 3), d), 1, 2);
        // A point one unit along the axis from the origin, on the axis.
        expect(inContainerCylinder3(v(1 + d.values[0], 2 + d.values[1],
            3 + d.values[2]), oblique)).toBe(true);
        // Two units along the axis is outside the half-height of 1.
        expect(inContainerCylinder3(v(1 + 2 * d.values[0], 2 + 2 * d.values[1],
            3 + 2 * d.values[2]), oblique)).toBe(false);
    });

    it('throws for non-3D points', () => {
        expect(() => inContainerCylinder3(Vector.fromArray([1, 2]), cylinder))
            .toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContCylinder3.h semantics.
// ---------------------------------------------------------------------------

describe('ContCylinder3 verification', () => {
    // Elongated clouds so the least-squares axis (a principal eigenvector) is
    // well conditioned; see the ContCapsule3 verification block.
    // The parameters t are spread evenly over [-5, 5] rather than drawn at
    // random, so the variance along the axis (about 8) always dominates the
    // transverse variance (at most 0.09) and the principal direction is never
    // ambiguous. Random t values can all collapse to zero, which makes the
    // eigenvector arbitrary and every equivariance property meaningless.
    const elongatedCloud = fc.tuple(unitVector(3), wellScaledVector(3, -4, 4),
        fc.array(wellScaledVector(3, -0.3, 0.3),
            { minLength: 4, maxLength: 12 }))
        .map(([axis, origin, noises]) =>
            noises.map((noise, i) => {
                const t = -5 + (10 * i) / (noises.length - 1);
                return add(add(origin, mul(t, axis)), noise);
            }));

    // The design claim of GetContainer: the cylinder contains every input
    // point. The radius and the height are both maxima over the same
    // projections that inContainer recomputes, so only rounding separates
    // the two evaluations.
    it('the fitted cylinder contains every input point', () => {
        check(elongatedCloud, (points: Vector[]) => {
            const cyl = getContainerCylinder3(points);
            for (const p of points) {
                expect(containedWithin(p, cyl, 1e-9)).toBe(true);
            }
        });
    });

    // Cross-check the fitted radius and height against independent maxima
    // computed from the returned axis.
    it('radius and height are the extremes of the projections', () => {
        check(elongatedCloud, (points: Vector[]) => {
            const cyl = getContainerCylinder3(points);
            const dir = cyl.axis.direction;
            if (length(dir) < 0.5) {
                return;    // degenerate cloud, no axis direction
            }
            let maxRadial = 0;
            let wMin = Infinity, wMax = -Infinity;
            for (const p of points) {
                const diff = sub(p, cyl.axis.origin);
                const w = dot(diff, dir);
                wMin = Math.min(wMin, w);
                wMax = Math.max(wMax, w);
                maxRadial = Math.max(maxRadial,
                    length(sub(diff, mul(w, dir))));
            }
            expectClose(cyl.radius, maxRadial, 1e-9, 1e-9);
            expectClose(cyl.height, wMax - wMin, 1e-9, 1e-9);
            // The origin returned is the midpoint of the projection interval,
            // so the interval is symmetric about it.
            expectClose(wMin + wMax, 0, 1e-9, 1e-9);
        });
    });

    // Rigid motions: same radius and height, transformed axis.
    it('is equivariant under rigid motions', () => {
        check(fc.tuple(elongatedCloud, rotationFrame(3), wellScaledVector(3)),
            ([points, frame, t]: [Vector[], Vector[], Vector]) => {
                const xform = (p: Vector): Vector =>
                    add(add(add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])),
                        mul(p.get(2), frame[2])), t);
                const c0 = getContainerCylinder3(points);
                const c1 = getContainerCylinder3(points.map(xform));
                expectClose(c1.radius, c0.radius, 1e-8, 1e-8);
                expectClose(c1.height, c0.height, 1e-8, 1e-8);
                const wantOrigin = xform(c0.axis.origin);
                expect(length(sub(c1.axis.origin, wantOrigin)))
                    .toBeLessThanOrEqual(1e-7);
                // The direction comes from an eigen decomposition, so it may
                // be negated.
                const wantDir = add(add(
                    mul(c0.axis.direction.get(0), frame[0]),
                    mul(c0.axis.direction.get(1), frame[1])),
                    mul(c0.axis.direction.get(2), frame[2]));
                expect(Math.min(length(sub(c1.axis.direction, wantDir)),
                    length(add(c1.axis.direction, wantDir))))
                    .toBeLessThanOrEqual(1e-7);
            });
    });

    // inContainer against an independent brute-force evaluation. Points near
    // the boundary are skipped: the two evaluations differ only by rounding
    // there, and the upstream test uses strict comparisons.
    it('inContainer agrees with a brute-force cylinder test', () => {
        check(fc.tuple(wellScaledVector(3, -3, 3), unitVector(3),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0.2, max: 6, noNaN: true }),
            wellScaledVector(3, -8, 8)),
            ([origin, dir, radius, height, p]:
                [Vector, Vector, number, number, Vector]) => {
                const cyl = new Cylinder3();
                cyl.axis = Line.fromOriginDirection(origin, dir);
                cyl.radius = radius;
                cyl.height = height;
                const diff = sub(p, origin);
                const w = dot(diff, dir);
                const radial = length(sub(diff, mul(w, dir)));
                const slackAxial = 0.5 * height - Math.abs(w);
                const slackRadial = radius - radial;
                if (Math.abs(slackAxial) < 1e-9 || Math.abs(slackRadial) < 1e-9) {
                    return;
                }
                expect(inContainerCylinder3(p, cyl))
                    .toBe(slackAxial > 0 && slackRadial > 0);
            });
    });

    // Points sampled inside the fitted cylinder of a cloud are reported
    // inside, which exercises inContainer against getContainer output.
    it('interior samples of the fitted cylinder are reported inside', () => {
        const rand = seededRandom(0x5eedc0);
        check(elongatedCloud, (points: Vector[]) => {
            const cyl = getContainerCylinder3(points);
            if (length(cyl.axis.direction) < 0.5 || cyl.radius < 1e-6
                || cyl.height < 1e-6) {
                return;
            }
            const basis = [cyl.axis.direction.clone(), new Vector(3),
                new Vector(3)];
            computeOrthogonalComplement3(1, basis);
            for (let k = 0; k < 20; ++k) {
                const w = (rand() - 0.5) * cyl.height * 0.9;
                const r = rand() * cyl.radius * 0.9;
                const a = 2 * Math.PI * rand();
                const q = add(cyl.axis.origin,
                    add(mul(w, basis[0]),
                        add(mul(r * Math.cos(a), basis[1]),
                            mul(r * Math.sin(a), basis[2]))));
                expect(inContainerCylinder3(q, cyl)).toBe(true);
            }
        }, 50);
    });
});
