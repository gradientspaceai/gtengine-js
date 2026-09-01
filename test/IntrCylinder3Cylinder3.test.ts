import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3';
import { Line } from '../src/Line';
import { Vector, sub, dot, length, normalize } from '../src/Vector';
import { cross } from '../src/Vector3';
import { IntrCylinder3Cylinder3TI } from '../src/IntrCylinder3Cylinder3';

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
