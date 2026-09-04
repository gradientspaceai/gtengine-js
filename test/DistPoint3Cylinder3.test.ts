import { describe, expect, it } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { DistPoint3Cylinder3 } from '../src/DistPoint3Cylinder3.js';
import { Line } from '../src/Line.js';
import {
    Vector, add, dot, getOrthogonal, length, mul, normalize, sub
} from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';
import { cross } from '../src/Vector3.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const d = v(...direction);
    normalize(d);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(v(...origin), d), radius, height);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistPoint3Cylinder3', () => {
    const query = new DistPoint3Cylinder3();

    it('reports zero distance for a point inside a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 2, 4);
        const result = query.compute(v(0.5, 0.5, 1), c);
        expect(result.distance).toBe(0);
        expect(result.closest[1].values[0]).toBeCloseTo(0.5, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0.5, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports zero distance for a point on the cylinder axis', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 2, 4);
        const result = query.compute(v(0, 0, 1), c);
        expect(result.distance).toBe(0);
    });

    it('measures a point outside the wall of an infinite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 1);
        c.makeInfiniteCylinder();
        const result = query.compute(v(5, 0, 100), c);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(100, 12);
    });

    it('reports zero distance inside an infinite cylinder at any height',
        () => {
            const c = cylinder([0, 0, 0], [0, 0, 1], 1, 1);
            c.makeInfiniteCylinder();
            const result = query.compute(v(0.5, 0, -1000), c);
            expect(result.distance).toBe(0);
        });

    it('clamps to the cap of a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 2, 4);
        // Directly above the axis, beyond the +z cap at z = 2.
        const result = query.compute(v(0, 0, 6), c);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(2, 12);
    });

    it('measures the rim of a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 2);
        // Radially 3 units outside the wall and 4 units above the +z cap.
        const result = query.compute(v(4, 0, 5), c);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('measures the -z cap of a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 2, 4);
        const result = query.compute(v(0.5, 0.5, -7), c);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(-2, 12);
    });

    it('rejects a nonpositive radius', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 0, 2);
        expect(() => query.compute(v(1, 1, 1), c)).toThrow(
            /positive radius/);
    });

    it('rejects a zero height for a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 0);
        expect(() => query.compute(v(1, 1, 1), c)).toThrow(
            /positive height/);
    });

    it('agrees with a dense sampling of a tilted finite cylinder', () => {
        const rnd = makeRandom(2718);
        const origin = v(0.5, -1, 0.25);
        const dir = v(1, 2, 3);
        normalize(dir);
        const radius = 1.25;
        const height = 3;
        const c = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin, dir), radius, height);

        // An orthonormal basis of the plane perpendicular to the axis.
        const U = v(-dir.values[1], dir.values[0], 0);
        normalize(U);
        const W = cross(dir, U);

        for (let trial = 0; trial < 30; ++trial) {
            const p = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const result = query.compute(p, c);

            // The reported closest point is in the solid cylinder.
            const delta = sub(result.closest[1], origin);
            const h = dot(delta, dir);
            const radial = sub(delta, mul(h, dir));
            expect(Math.abs(h)).toBeLessThanOrEqual(0.5 * height + 1e-9);
            expect(Math.sqrt(dot(radial, radial))).toBeLessThanOrEqual(
                radius + 1e-9);

            // The reported closest point realizes the reported distance.
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 9);

            // No sampled solid-cylinder point is closer.
            const nR = 10, nA = 40, nH = 12;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= nR; ++i) {
                const r = radius * i / nR;
                for (let j = 0; j < nA; ++j) {
                    const a = 2 * Math.PI * j / nA;
                    for (let k = 0; k <= nH; ++k) {
                        const hh = height * (k / nH - 0.5);
                        const q = add(origin, add(mul(hh, dir),
                            add(mul(r * Math.cos(a), U),
                                mul(r * Math.sin(a), W))));
                        const f = sub(p, q);
                        best = Math.min(best, dot(f, f));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistPoint3Cylinder3.h.
// ---------------------------------------------------------------------------

describe('DistPoint3Cylinder3 verification', () => {
    const query = new DistPoint3Cylinder3();

    const finiteArb = fc.tuple(wellScaledVector(3, -5, 5), unitVector(3),
        finite(0.1, 4), finite(0.1, 8))
        .map(([o, d, radius, height]) => Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(o, d), radius, height));

    const infiniteArb = fc.tuple(wellScaledVector(3, -5, 5), unitVector(3),
        finite(0.1, 4))
        .map(([o, d, radius]) => {
            const c = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(o, d), radius, 1);
            c.makeInfiniteCylinder();
            return c;
        });

    // Cylinder coordinates of a point: the signed distance along the axis and
    // the distance from the axis.
    function cylinderCoords(p: Vector, cyl: Cylinder3):
        { axial: number, radial: number } {
        const delta = sub(p, cyl.axis.origin);
        const axial = dot(cyl.axis.direction, delta);
        const radial = length(sub(delta, mul(axial, cyl.axis.direction)));
        return { axial, radial };
    }

    // A solid finite cylinder is the product of a disk and an interval, and
    // the two factors are orthogonal, so the distance is the Euclidean
    // combination of the per-factor distances. This is an independent closed
    // form for the query.
    function closedForm(p: Vector, cyl: Cylinder3): number {
        const { axial, radial } = cylinderCoords(p, cyl);
        const dr = Math.max(radial - cyl.radius, 0);
        if (cyl.isInfinite()) {
            return dr;
        }
        const dz = Math.max(Math.abs(axial) - 0.5 * cyl.height, 0);
        return Math.sqrt(dr * dr + dz * dz);
    }

    it('matches the closed form for finite cylinders', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), finiteArb), ([p, cyl]) => {
            const r = query.compute(p, cyl);
            expectClose(r.distance, closedForm(p, cyl), 1e-9, 1e-9);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-9, 1e-9);
            expectVectorClose(r.closest[0], p, 0, 0);
            expect(r.closest[0]).not.toBe(p);
            // closest[1] is a cylinder point.
            const c = cylinderCoords(r.closest[1], cyl);
            expect(c.radial).toBeLessThanOrEqual(cyl.radius + 1e-9);
            expect(Math.abs(c.axial))
                .toBeLessThanOrEqual(0.5 * cyl.height + 1e-9);
        });
    });

    it('matches the closed form for infinite cylinders', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), infiniteArb),
            ([p, cyl]) => {
                const r = query.compute(p, cyl);
                expectClose(r.distance, closedForm(p, cyl), 1e-9, 1e-9);
                const c = cylinderCoords(r.closest[1], cyl);
                expect(c.radial).toBeLessThanOrEqual(cyl.radius + 1e-9);
            });
    });

    it('accepts the Cylinder3 infinite sentinel (upstream issue #187)', () => {
        // Upstream tests height == numeric_limits<T>::max(), but
        // MakeInfiniteCylinder sets height = -1, so upstream falls into the
        // finite branch and trips its "positive height" assertion.
        const cyl = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(v(0, 0, 0), v(0, 0, 1)), 2, 1);
        cyl.makeInfiniteCylinder();
        expect(cyl.height).toBe(-1);
        const r = query.compute(v(5, 0, 1000), cyl);
        expect(r.distance).toBeCloseTo(3, 12);
        expectVectorClose(r.closest[1], v(2, 0, 1000), 1e-9, 1e-9);
    });

    it('reports zero distance for points inside the cylinder', () => {
        check(fc.tuple(finiteArb, finite(0, 0.99), finite(-0.99, 0.99),
            finite(-Math.PI, Math.PI)), ([cyl, ru, zu, angle]) => {
            const u = getOrthogonal(cyl.axis.direction, true);
            const w = cross(cyl.axis.direction, u);
            const radial = ru * cyl.radius;
            const p = add(cyl.axis.origin,
                add(mul(zu * 0.5 * cyl.height, cyl.axis.direction),
                    add(mul(radial * Math.cos(angle), u),
                        mul(radial * Math.sin(angle), w))));
            const r = query.compute(p, cyl);
            expectClose(r.distance, 0, 1e-9, 1e-9);
            expectVectorClose(r.closest[1], p, 1e-9, 1e-9);
        });
    });

    it('handles points on the cylinder axis', () => {
        check(fc.tuple(finiteArb, finite(-8, 8)), ([cyl, t]) => {
            const p = add(cyl.axis.origin, mul(t, cyl.axis.direction));
            const r = query.compute(p, cyl);
            const expected = Math.max(Math.abs(t) - 0.5 * cyl.height, 0);
            expectClose(r.distance, expected, 1e-9, 1e-9);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), finiteArb,
            rotationFrame(3), wellScaledVector(3, -6, 6)),
            ([p, cyl, frame, shift]) => {
                const rot = (q: Vector): Vector =>
                    add(add(mul(q.values[0], frame[0]),
                        mul(q.values[1], frame[1])),
                        mul(q.values[2], frame[2]));
                const moved = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(add(shift, rot(cyl.axis.origin)),
                        rot(cyl.axis.direction)), cyl.radius, cyl.height);
                const r0 = query.compute(p, cyl);
                const r1 = query.compute(add(shift, rot(p)), moved);
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            });
    });

    it('rejects a non-positive radius or height', () => {
        const axis = Line.fromOriginDirection(v(0, 0, 0), v(0, 0, 1));
        expect(() => query.compute(v(1, 0, 0),
            Cylinder3.fromAxisRadiusHeight(axis, 0, 1)))
            .toThrow('The cylinder must have a positive radius.');
        expect(() => query.compute(v(1, 0, 0),
            Cylinder3.fromAxisRadiusHeight(axis, 1, 0)))
            .toThrow('The cylinder must have a positive height.');
    });
});
