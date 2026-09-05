import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, sub, dot, length, normalize } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { IntrCylinder3Cylinder3TI } from '../src/IntrCylinder3Cylinder3.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function cylinder(c: Vector, d: Vector, r: number, h: number): Cylinder3 {
    const dir = d.clone();
    normalize(dir);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(c, dir), r, h);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The radius of the projection of the cylinder onto the unit direction D.
function projectionRadius(cyl: Cylinder3, d: Vector): number {
    return cyl.radius * length(cross(cyl.axis.direction, d))
        + 0.5 * cyl.height * Math.abs(dot(cyl.axis.direction, d));
}

// Verify that the reported direction really separates the two cylinders.
function separates(c0: Cylinder3, c1: Cylinder3, d: Vector): boolean {
    const unit = d.clone();
    normalize(unit);
    const delta = sub(c1.axis.origin, c0.axis.origin);
    return Math.abs(dot(delta, unit))
        > projectionRadius(c0, unit) + projectionRadius(c1, unit);
}

// Brute-force containment in the solid finite cylinder.
function inCylinder(p: Vector, cyl: Cylinder3): boolean {
    const diff = sub(p, cyl.axis.origin);
    const t = dot(diff, cyl.axis.direction);
    if (Math.abs(t) > 0.5 * cyl.height) {
        return false;
    }
    const radial = v3(0, 0, 0);
    for (let i = 0; i < 3; ++i) {
        radial.values[i] = diff.values[i]
            - t * cyl.axis.direction.values[i];
    }
    return length(radial) <= cyl.radius;
}

describe('IntrCylinder3Cylinder3', () => {
    const ti = new IntrCylinder3Cylinder3TI(1, 16, 8);

    it('rejects invalid angle counts', () => {
        expect(() => new IntrCylinder3Cylinder3TI(1, 0, 4)).toThrow();
        expect(() => new IntrCylinder3Cylinder3TI(1, 4, 0)).toThrow();
    });

    it('keeps the requested thread count for API compatibility', () => {
        expect(new IntrCylinder3Cylinder3TI(4, 8, 8).numThreads).toBe(4);
    });

    it('reports no separation for coincident axes origins', () => {
        const c0 = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
        const c1 = cylinder(v3(0, 0, 0), v3(1, 0, 0), 1, 2);
        const result = ti.test(c0, c1);
        expect(result.separated).toBe(false);
        expect(result.separatingDirection.values).toEqual([0, 0, 0]);
    });

    it('separates parallel cylinders stacked along the axis', () => {
        const c0 = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
        const c1 = cylinder(v3(0, 0, 5), v3(0, 0, 1), 1, 2);
        const result = ti.test(c0, c1);
        expect(result.separated).toBe(true);
        expect(result.separatingDirection.values).toEqual([0, 0, 1]);
        expect(separates(c0, c1, result.separatingDirection)).toBe(true);
    });

    it('separates parallel cylinders offset radially', () => {
        const c0 = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
        const c1 = cylinder(v3(5, 0, 0), v3(0, 0, 1), 1, 2);
        const result = ti.test(c0, c1);
        expect(result.separated).toBe(true);
        expect(separates(c0, c1, result.separatingDirection)).toBe(true);
    });

    it('does not separate overlapping parallel cylinders', () => {
        const c0 = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
        const c1 = cylinder(v3(1, 0, 0.5), v3(0, 0, 1), 1, 2);
        expect(ti.test(c0, c1).separated).toBe(false);
    });

    it('separates perpendicular cylinders that are far apart', () => {
        const c0 = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
        const c1 = cylinder(v3(0, 8, 0), v3(1, 0, 0), 1, 2);
        const result = ti.test(c0, c1);
        expect(result.separated).toBe(true);
        expect(separates(c0, c1, result.separatingDirection)).toBe(true);
    });

    it('does not separate crossing perpendicular cylinders', () => {
        const c0 = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 4);
        const c1 = cylinder(v3(0, 0, 0.5), v3(1, 0, 0), 1, 4);
        expect(ti.test(c0, c1).separated).toBe(false);
    });

    it('reports a genuine separating direction whenever it reports one',
        () => {
            const rnd = makeRandom(31337);
            let numSeparated = 0;
            for (let k = 0; k < 250; ++k) {
                const d0 = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
                const d1 = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
                if (dot(d0, d0) < 1e-6 || dot(d1, d1) < 1e-6) {
                    continue;
                }
                const c0 = cylinder(v3(0, 0, 0), d0, 0.2 + rnd(),
                    0.4 + rnd() * 2);
                const c1 = cylinder(
                    v3(rnd() * 8 - 4, rnd() * 8 - 4, rnd() * 8 - 4),
                    d1, 0.2 + rnd(), 0.4 + rnd() * 2);
                const result = ti.test(c0, c1);
                if (result.separated) {
                    // Allow a tiny tolerance for the sampled directions.
                    const unit = result.separatingDirection.clone();
                    normalize(unit);
                    const delta = sub(c1.axis.origin, c0.axis.origin);
                    const gap = Math.abs(dot(delta, unit))
                        - projectionRadius(c0, unit)
                        - projectionRadius(c1, unit);
                    expect(gap).toBeGreaterThan(-1e-9);
                    ++numSeparated;
                }
            }
            expect(numSeparated).toBeGreaterThan(0);
        });

    it('never reports separation when the cylinders actually overlap', () => {
        const rnd = makeRandom(60649);
        let numOverlap = 0;
        for (let k = 0; k < 200; ++k) {
            const d0 = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            const d1 = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d0, d0) < 1e-6 || dot(d1, d1) < 1e-6) {
                continue;
            }
            const c0 = cylinder(v3(0, 0, 0), d0, 0.3 + rnd(),
                0.5 + rnd() * 2);
            const c1 = cylinder(
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                d1, 0.3 + rnd(), 0.5 + rnd() * 2);
            if (length(sub(c1.axis.origin, c0.axis.origin)) === 0) {
                continue;
            }

            // Sample points of cylinder1 and test them against cylinder0.
            const a = Math.abs(c1.axis.direction.values[0]) < 0.9
                ? v3(1, 0, 0) : v3(0, 1, 0);
            const u = cross(a, c1.axis.direction);
            normalize(u);
            const w = cross(c1.axis.direction, u);

            let overlaps = false;
            for (let s = 0; s <= 8 && !overlaps; ++s) {
                const t = c1.height * (-0.5 + s / 8);
                for (let ri = 0; ri <= 3 && !overlaps; ++ri) {
                    const rr = (c1.radius * ri) / 3;
                    for (let ai = 0; ai < 16 && !overlaps; ++ai) {
                        const th = (2 * Math.PI * ai) / 16;
                        const p = v3(0, 0, 0);
                        for (let c = 0; c < 3; ++c) {
                            p.values[c] = c1.axis.origin.values[c]
                                + t * c1.axis.direction.values[c]
                                + rr * Math.cos(th) * u.values[c]
                                + rr * Math.sin(th) * w.values[c];
                        }
                        if (inCylinder(p, c0)) {
                            overlaps = true;
                        }
                    }
                }
            }

            if (overlaps) {
                expect(ti.test(c0, c1).separated).toBe(false);
                ++numOverlap;
            }
        }
        expect(numOverlap).toBeGreaterThan(0);
    });

    it('finds more separations with a finer hemisphere sampling', () => {
        // A coarse sampling can miss a separating direction; a finer one
        // must not lose a separation that the coarse one found.
        const coarse = new IntrCylinder3Cylinder3TI(1, 8, 4);
        const fine = new IntrCylinder3Cylinder3TI(1, 64, 32);
        const rnd = makeRandom(4096);
        for (let k = 0; k < 60; ++k) {
            const c0 = cylinder(v3(0, 0, 0),
                v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1),
                0.2 + rnd(), 0.4 + rnd() * 2);
            const c1 = cylinder(
                v3(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3),
                v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1),
                0.2 + rnd(), 0.4 + rnd() * 2);
            if (coarse.test(c0, c1).separated) {
                expect(fine.test(c0, c1).separated).toBe(true);
            }
        }
    });

    it('throws for an infinite cylinder', () => {
        const c0 = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
        const c1 = cylinder(v3(5, 0, 0), v3(0, 0, 1), 1, 2);
        c1.makeInfiniteCylinder();
        expect(() => ti.test(c0, c1)).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, seededRandom
} from './helpers/arbitraries.js';
import { add, mul } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

const cylArb = fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
    fc.double({ min: 0.2, max: 2, noNaN: true }),
    fc.double({ min: 0.3, max: 4, noNaN: true }))
    .map(([o, w, r, h]) => Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(o, w), r, h));

// The upstream separation functional for a candidate direction D: the
// cylinders are separated by D when this is negative.
function separationTest(c0: Cylinder3, c1: Cylinder3, d: Vector): number {
    const delta = sub(c1.axis.origin, c0.axis.origin);
    return c0.radius * length(cross(c0.axis.direction, d))
        + c1.radius * length(cross(c1.axis.direction, d))
        + 0.5 * c0.height * Math.abs(dot(c0.axis.direction, d))
        + 0.5 * c1.height * Math.abs(dot(c1.axis.direction, d))
        - Math.abs(dot(delta, d));
}

// Sample the surface of a cylinder (which contains the extreme points along
// any direction).
function sampleCylinder(c: Cylinder3, nz: number, na: number): Vector[] {
    const basis = [c.axis.direction.clone(), new Vector(3), new Vector(3)];
    computeOrthogonalComplement3(1, basis);
    const half = 0.5 * c.height;
    const pts: Vector[] = [];
    for (let i = 0; i <= nz; ++i) {
        const z = -half + (2 * half * i) / nz;
        for (let k = 0; k < na; ++k) {
            const a = (2 * Math.PI * k) / na;
            pts.push(add(c.axis.origin, add(mul(z, basis[0]),
                add(mul(c.radius * Math.cos(a), basis[1]),
                    mul(c.radius * Math.sin(a), basis[2])))));
        }
    }
    return pts;
}

describe('IntrCylinder3Cylinder3 verification', () => {
    const tiq = new IntrCylinder3Cylinder3TI(1, 16, 8);

    it('a reported separating direction really separates', () => {
        check(fc.tuple(cylArb, cylArb), ([c0, c1]) => {
            const r = tiq.test(c0, c1);
            if (!r.separated) {
                return;
            }
            const d = r.separatingDirection;
            expect(length(d)).toBeGreaterThan(0);
            expect(separationTest(c0, c1, d)).toBeLessThan(0);
            // The sampled point sets are on opposite sides of some plane
            // with normal d.
            const a = sampleCylinder(c0, 6, 24).map(p => dot(d, p));
            const b = sampleCylinder(c1, 6, 24).map(p => dot(d, p));
            const gap = Math.min(...b) - Math.max(...a);
            const gap2 = Math.min(...a) - Math.max(...b);
            expect(Math.max(gap, gap2)).toBeGreaterThan(0);
        });
    });

    it('is symmetric under argument swap', () => {
        check(fc.tuple(cylArb, cylArb), ([c0, c1]) => {
            const a = tiq.test(c0, c1);
            const b = tiq.test(c1, c0);
            if (a.separated !== b.separated) {
                // The hemisphere sampling is not swap invariant (the basis is
                // built from Delta, which changes sign), so a direction found
                // for one order may be missed in the other. Accept only when
                // the functional is essentially zero for the found direction.
                const found = a.separated ? a : b;
                expect(Math.abs(separationTest(c0, c1,
                    found.separatingDirection))).toBeLessThan(1e-6);
                return;
            }
            expect(a.separated).toBe(b.separated);
        });
    });

    it('coincident axis origins are never separated', () => {
        check(fc.tuple(cylArb, unitVector(3),
            fc.double({ min: 0.2, max: 2, noNaN: true }),
            fc.double({ min: 0.3, max: 4, noNaN: true })),
            ([c0, w, r, h]) => {
                const c1 = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(c0.axis.origin, w), r, h);
                const res = tiq.test(c0, c1);
                expect(res.separated).toBe(false);
                expect(res.separatingDirection.values).toEqual([0, 0, 0]);
            });
    });

    it('cylinders moved far apart are separated', () => {
        check(fc.tuple(cylArb, cylArb, unitVector(3)), ([c0, c1, d]) => {
            const reach = c0.radius + c1.radius + 0.5 * c0.height
                + 0.5 * c1.height + 1;
            const moved = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(
                    add(c0.axis.origin, mul(2 * reach, d)),
                    c1.axis.direction), c1.radius, c1.height);
            const res = tiq.test(c0, moved);
            expect(res.separated).toBe(true);
            expect(separationTest(c0, moved, res.separatingDirection))
                .toBeLessThan(0);
        });
    });

    it('overlapping cylinders are never reported as separated', () => {
        const rnd = seededRandom(0x2f9a13cb);
        for (let trial = 0; trial < 120; ++trial) {
            const mk = (o: Vector): Cylinder3 => {
                const w = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                    rnd() * 2 - 1]);
                normalize(w);
                return Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(o, w), 0.3 + rnd(),
                    0.5 + rnd() * 2);
            };
            const o0 = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            const c0 = mk(o0);
            // Place the second axis origin inside the first cylinder, so the
            // solids certainly overlap.
            const c1 = mk(o0);
            const res = tiq.test(c0, c1);
            expect(res.separated).toBe(false);
        }
    }, 30000);

    it('is invariant under a common rigid motion', () => {
        check(fc.tuple(cylArb, cylArb, rotationFrame(3),
            wellScaledVector(3, -3, 3)), ([c0, c1, R, tr]) => {
            const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
            const xf = (c: Cylinder3): Cylinder3 =>
                Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(add(tr, rot(c.axis.origin)),
                        rot(c.axis.direction)), c.radius, c.height);
            const a = tiq.test(c0, c1);
            const b = tiq.test(xf(c0), xf(c1));
            if (a.separated !== b.separated) {
                // The hemisphere sampling is defined in a frame built from
                // Delta, which rotates with the configuration, but the sample
                // directions are only approximately equivariant; accept a
                // difference only for a marginal separation.
                const found = a.separated ? a.separatingDirection
                    : b.separatingDirection;
                const t = a.separated ? separationTest(c0, c1, found)
                    : separationTest(xf(c0), xf(c1), found);
                expect(Math.abs(t)).toBeLessThan(1e-3);
                return;
            }
            expect(a.separated).toBe(b.separated);
        });
    });

    it('throws for infinite cylinders and invalid angle counts', () => {
        const c = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(Vector.fromArray([0, 0, 0]),
                Vector.fromArray([0, 0, 1])), 1, 2);
        const inf = c.clone();
        inf.makeInfiniteCylinder();
        expect(() => tiq.test(inf, c)).toThrow();
        expect(() => new IntrCylinder3Cylinder3TI(1, 0, 4)).toThrow();
        expect(() => new IntrCylinder3Cylinder3TI(1, 4, 0)).toThrow();
    });
});
