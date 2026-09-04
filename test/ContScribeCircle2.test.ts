import { describe, it, expect } from 'vitest';
import { circumscribeCircle2, inscribeCircle2 } from '../src/ContScribeCircle2.js';
import { Vector, add, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, fc, latticeVector, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';
import { dotPerp } from '../src/Vector2.js';

function v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Distance from a point to the line through p0 and p1.
function lineDistance(c: Vector, p0: Vector, p1: Vector): number {
    const e = sub(p1, p0);
    return Math.abs(dotPerp(e, sub(c, p0))) / length(e);
}

describe('ContScribeCircle2', () => {
    describe('circumscribeCircle2', () => {
        it('computes the circumcircle of a right triangle', () => {
            // Legs 3 and 4, hypotenuse 5; the circumcenter is the midpoint of
            // the hypotenuse and the circumradius is 2.5.
            const circle = circumscribeCircle2(v(0, 0), v(3, 0), v(0, 4));
            expect(circle).not.toBeNull();
            expect(circle!.center.values[0]).toBeCloseTo(1.5, 12);
            expect(circle!.center.values[1]).toBeCloseTo(2, 12);
            expect(circle!.radius).toBeCloseTo(2.5, 12);
        });

        it('computes the circumcircle of an equilateral triangle', () => {
            // Side 1 centered at the origin; circumradius = 1/sqrt(3).
            const a = v(1, 0);
            const b = v(Math.cos(2 * Math.PI / 3), Math.sin(2 * Math.PI / 3));
            const c = v(Math.cos(4 * Math.PI / 3), Math.sin(4 * Math.PI / 3));
            const circle = circumscribeCircle2(a, b, c);
            expect(circle).not.toBeNull();
            expect(length(circle!.center)).toBeCloseTo(0, 12);
            expect(circle!.radius).toBeCloseTo(1, 12);
        });

        it('returns null for collinear points', () => {
            expect(circumscribeCircle2(v(0, 0), v(1, 1), v(2, 2))).toBeNull();
            expect(circumscribeCircle2(v(0, 0), v(0, 0), v(1, 3))).toBeNull();
        });

        it('rejects non-2D points', () => {
            expect(() => circumscribeCircle2(Vector.fromArray([0, 0, 0]),
                v(1, 0), v(0, 1))).toThrow();
        });
    });

    describe('inscribeCircle2', () => {
        it('computes the incircle of a 3-4-5 right triangle', () => {
            // The inradius of a right triangle is (a+b-c)/2 = (3+4-5)/2 = 1
            // and the incenter is (r,r) for this placement.
            const circle = inscribeCircle2(v(0, 0), v(3, 0), v(0, 4));
            expect(circle).not.toBeNull();
            expect(circle!.radius).toBeCloseTo(1, 12);
            expect(circle!.center.values[0]).toBeCloseTo(1, 12);
            expect(circle!.center.values[1]).toBeCloseTo(1, 12);
        });

        it('returns null for degenerate triangles', () => {
            expect(inscribeCircle2(v(0, 0), v(0, 0), v(0, 0))).toBeNull();
            expect(inscribeCircle2(v(0, 0), v(1, 1), v(2, 2))).toBeNull();
        });
    });

    it('the circumcircle passes through all three vertices and the incircle '
        + 'is tangent to all three edges (randomized)', () => {
        let seed = 424242;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rv = () => v(20 * rand() - 10, 20 * rand() - 10);

        let tested = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const a = rv(), b = rv(), c = rv();
            const twiceArea = Math.abs(dotPerp(sub(b, a), sub(c, a)));
            if (twiceArea < 1) {
                continue; // too close to degenerate for tight tolerances
            }
            ++tested;

            const circum = circumscribeCircle2(a, b, c);
            expect(circum).not.toBeNull();
            for (const p of [a, b, c]) {
                expect(length(sub(p, circum!.center)))
                    .toBeCloseTo(circum!.radius, 8);
            }

            const inc = inscribeCircle2(a, b, c);
            expect(inc).not.toBeNull();
            expect(lineDistance(inc!.center, a, b)).toBeCloseTo(inc!.radius, 8);
            expect(lineDistance(inc!.center, b, c)).toBeCloseTo(inc!.radius, 8);
            expect(lineDistance(inc!.center, c, a)).toBeCloseTo(inc!.radius, 8);

            // The incircle radius equals Area/s.
            const s = 0.5 * (length(sub(b, a)) + length(sub(c, b))
                + length(sub(a, c)));
            expect(inc!.radius).toBeCloseTo(0.5 * twiceArea / s, 8);
        }
        expect(tested).toBeGreaterThan(400);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContScribeCircle2.h semantics.
// ---------------------------------------------------------------------------

describe('ContScribeCircle2 verification', () => {
    // Lattice triangles with a nonzero exact signed area, so the circumcircle
    // and incircle are well defined and the degenerate branches are reached
    // only by the explicitly collinear cases below.
    const triangle2 = fc.tuple(latticeVector(2, -6, 6), latticeVector(2, -6, 6),
        latticeVector(2, -6, 6))
        .filter(([a, b, c]) => Math.abs(dotPerp(sub(b, a), sub(c, a))) >= 2);

    // The circumcircle is equidistant from all three vertices; that distance
    // is the radius.
    it('the circumcircle passes through all three vertices', () => {
        check(triangle2, ([a, b, c]: Vector[]) => {
            const circle = circumscribeCircle2(a, b, c)!;
            expect(circle).not.toBeNull();
            for (const p of [a, b, c]) {
                expectClose(length(sub(p, circle.center)), circle.radius,
                    1e-9, 1e-9);
            }
        });
    });

    // The incircle is tangent to all three edge lines: the distance from its
    // center to each line equals the radius, and the center is inside the
    // triangle.
    it('the incircle is tangent to all three edges from inside', () => {
        check(triangle2, ([a, b, c]: Vector[]) => {
            const circle = inscribeCircle2(a, b, c)!;
            expect(circle).not.toBeNull();
            for (const [p, q] of [[a, b], [b, c], [c, a]] as const) {
                expectClose(lineDistance(circle.center, p, q), circle.radius,
                    1e-9, 1e-9);
            }
            // Inside: the center is on the same side of each directed edge as
            // the opposite vertex.
            for (const [p, q, r] of [[a, b, c], [b, c, a], [c, a, b]] as const) {
                const e = sub(q, p);
                const s1 = dotPerp(e, sub(circle.center, p));
                const s2 = dotPerp(e, sub(r, p));
                expect(s1 * s2).toBeGreaterThan(0);
            }
        });
    });

    // Rigid motions: the radii are invariant and the centers follow.
    it('both circles are equivariant under rigid motions', () => {
        check(fc.tuple(triangle2, rotationFrame(2), wellScaledVector(2)),
            ([[a, b, c], frame, t]: [Vector[], Vector[], Vector]) => {
                const xform = (p: Vector): Vector =>
                    add(add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])), t);
                for (const f of [circumscribeCircle2, inscribeCircle2]) {
                    const c0 = f(a, b, c)!;
                    const c1 = f(xform(a), xform(b), xform(c))!;
                    // The circumradius of a lattice triangle can be large when
                    // the triangle is thin, and the rotation costs a few ulps
                    // of that magnitude.
                    expectClose(c1.radius, c0.radius, 1e-9, 1e-8);
                    expect(length(sub(c1.center, xform(c0.center))))
                        .toBeLessThanOrEqual(1e-8 * (1 + c0.radius));
                }
            });
    });

    // Degenerate input: collinear lattice points make the 2x2 system singular
    // (circumscribe) and give a zero inradius (inscribe). Both return null.
    it('returns null for collinear triples', () => {
        check(fc.tuple(latticeVector(2, -6, 6), latticeVector(2, -6, 6),
            fc.integer({ min: -4, max: 4 }), fc.integer({ min: -4, max: 4 })),
            ([a, d, s, t]: [Vector, Vector, number, number]) => {
                if (d.get(0) === 0 && d.get(1) === 0) {
                    return;
                }
                const b = add(a, mul(s, d));
                const c = add(a, mul(t, d));
                expect(circumscribeCircle2(a, b, c)).toBeNull();
                expect(inscribeCircle2(a, b, c)).toBeNull();
            });
    });

    // Upstream issue #292 item 6: Inscribe writes a degenerate circle into the
    // output before returning false. The port returns null instead, so the
    // caller cannot mistake the degenerate circle for a result.
    it('returns null rather than a degenerate circle (#292)', () => {
        expect(inscribeCircle2(v(0, 0), v(0, 0), v(0, 0))).toBeNull();
        expect(inscribeCircle2(v(0, 0), v(1, 0), v(2, 0))).toBeNull();
    });
});
