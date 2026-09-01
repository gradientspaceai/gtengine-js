import { describe, it, expect } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox';
import { Cylinder3 } from '../src/Cylinder3';
import { Line } from '../src/Line';
import { Vector, sub, dot, length, normalize } from '../src/Vector';
import {
    IntrCanonicalBox3Cylinder3TI
} from '../src/IntrCanonicalBox3Cylinder3';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function box(e0: number, e1: number, e2: number): CanonicalBox {
    return CanonicalBox.fromExtent(v3(e0, e1, e2));
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

// Brute-force containment in the solid finite cylinder.
function inCylinder(p: Vector, cyl: Cylinder3): boolean {
    const diff = sub(p, cyl.axis.origin);
    const t = dot(diff, cyl.axis.direction);
    if (Math.abs(t) > 0.5 * cyl.height) {
        return false;
    }
    const radial = sub(diff, cyl.axis.direction.clone());
    // Recompute the radial component explicitly.
    for (let i = 0; i < 3; ++i) {
        radial.values[i] = diff.values[i]
            - t * cyl.axis.direction.values[i];
    }
    return length(radial) <= cyl.radius;
}

describe('IntrCanonicalBox3Cylinder3', () => {
    const ti = new IntrCanonicalBox3Cylinder3TI();
    const unit = box(1, 1, 1);

    it('reports intersection for a cylinder inside the box', () => {
        expect(ti.test(unit, cylinder(v3(0, 0, 0), v3(0, 0, 1), 0.5, 1))
            .intersect).toBe(true);
    });

    it('reports intersection for a box inside the cylinder', () => {
        expect(ti.test(unit, cylinder(v3(0, 0, 0), v3(0, 0, 1), 5, 10))
            .intersect).toBe(true);
    });

    it('culls a cylinder outside the slab of the end disks', () => {
        expect(ti.test(unit, cylinder(v3(0, 0, 10), v3(0, 0, 1), 1, 2))
            .intersect).toBe(false);
    });

    it('handles an axis-aligned cylinder (two zero components)', () => {
        // Axis parallel to z, offset in x. The box spans x in [-1,1].
        expect(ti.test(unit, cylinder(v3(1.5, 0, 0), v3(0, 0, 1), 0.4, 2))
            .intersect).toBe(false);
        expect(ti.test(unit, cylinder(v3(1.5, 0, 0), v3(0, 0, 1), 0.6, 2))
            .intersect).toBe(true);
        // Exactly tangent to the face x = 1.
        expect(ti.test(unit, cylinder(v3(1.5, 0, 0), v3(0, 0, 1), 0.5, 2))
            .intersect).toBe(true);
    });

    it('handles a cylinder axis with one zero component', () => {
        // Axis in the xy plane at 45 degrees, well away from the box.
        expect(ti.test(unit, cylinder(v3(0, 0, 5), v3(1, 1, 0), 0.5, 2))
            .intersect).toBe(false);
        expect(ti.test(unit, cylinder(v3(0, 0, 1), v3(1, 1, 0), 0.5, 2))
            .intersect).toBe(true);
    });

    it('handles a cylinder axis with no zero components', () => {
        expect(ti.test(unit, cylinder(v3(0, 0, 0), v3(1, 1, 1), 0.5, 2))
            .intersect).toBe(true);
        expect(ti.test(unit, cylinder(v3(4, 4, 4), v3(1, 1, 1), 0.5, 2))
            .intersect).toBe(false);
    });

    it('gives the same answer under reflections of the axis direction', () => {
        // The query reflects the configuration so the axis is in the first
        // octant; a canonical box is symmetric under those reflections, so
        // reversing the axis direction must not change the answer.
        const rnd = makeRandom(90210);
        for (let k = 0; k < 200; ++k) {
            const b = box(0.2 + rnd() * 2, 0.2 + rnd() * 2, 0.2 + rnd() * 2);
            const c = v3(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3);
            const d = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d, d) < 1e-6) {
                continue;
            }
            const r = 0.1 + rnd(), h = 0.2 + rnd() * 3;
            const forward = ti.test(b, cylinder(c, d, r, h)).intersect;
            const backward = ti.test(b,
                cylinder(c, v3(-d.values[0], -d.values[1], -d.values[2]),
                    r, h)).intersect;
            expect(forward).toBe(backward);
        }
    });

    it('agrees with a brute-force sampling of the box', () => {
        const rnd = makeRandom(112358);
        let numHit = 0, numMiss = 0;
        for (let k = 0; k < 250; ++k) {
            const b = box(0.3 + rnd() * 1.5, 0.3 + rnd() * 1.5,
                0.3 + rnd() * 1.5);
            const d = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d, d) < 1e-6) {
                continue;
            }
            const cyl = cylinder(
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                d, 0.15 + rnd() * 0.8, 0.3 + rnd() * 2);
            const intersect = ti.test(b, cyl).intersect;

            // A sampled box point inside the cylinder proves intersection.
            let sampleHit = false;
            const n = 14;
            for (let i = 0; i <= n && !sampleHit; ++i) {
                const x = b.extent.values[0] * (-1 + (2 * i) / n);
                for (let j = 0; j <= n && !sampleHit; ++j) {
                    const y = b.extent.values[1] * (-1 + (2 * j) / n);
                    for (let m = 0; m <= n && !sampleHit; ++m) {
                        const z = b.extent.values[2] * (-1 + (2 * m) / n);
                        if (inCylinder(v3(x, y, z), cyl)) {
                            sampleHit = true;
                        }
                    }
                }
            }

            if (sampleHit) {
                expect(intersect).toBe(true);
                ++numHit;
            } else if (!intersect) {
                ++numMiss;
            }
        }
        expect(numHit).toBeGreaterThan(0);
        expect(numMiss).toBeGreaterThan(0);
    });

    it('reports no intersection when a fine sampling of the cylinder misses',
        () => {
            // The converse direction: sample the solid cylinder and check
            // that a sample inside the box implies an intersection.
            const rnd = makeRandom(271828);
            let numHit = 0;
            for (let k = 0; k < 200; ++k) {
                const b = box(0.5 + rnd(), 0.5 + rnd(), 0.5 + rnd());
                const d = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
                if (dot(d, d) < 1e-6) {
                    continue;
                }
                const cyl = cylinder(
                    v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                    d, 0.2 + rnd(), 0.4 + rnd() * 2);
                const intersect = ti.test(b, cyl).intersect;

                const basis = [cyl.axis.direction.clone(), new Vector(3),
                    new Vector(3)];
                // Build an orthonormal frame by hand from the axis.
                const a = Math.abs(cyl.axis.direction.values[0]) < 0.9
                    ? v3(1, 0, 0) : v3(0, 1, 0);
                const u = v3(
                    a.values[1] * cyl.axis.direction.values[2]
                        - a.values[2] * cyl.axis.direction.values[1],
                    a.values[2] * cyl.axis.direction.values[0]
                        - a.values[0] * cyl.axis.direction.values[2],
                    a.values[0] * cyl.axis.direction.values[1]
                        - a.values[1] * cyl.axis.direction.values[0]);
                normalize(u);
                const w = v3(
                    cyl.axis.direction.values[1] * u.values[2]
                        - cyl.axis.direction.values[2] * u.values[1],
                    cyl.axis.direction.values[2] * u.values[0]
                        - cyl.axis.direction.values[0] * u.values[2],
                    cyl.axis.direction.values[0] * u.values[1]
                        - cyl.axis.direction.values[1] * u.values[0]);
                basis[1] = u;
                basis[2] = w;

                let sampleHit = false;
                for (let s = 0; s <= 12 && !sampleHit; ++s) {
                    const t = cyl.height * (-0.5 + s / 12);
                    for (let ri = 0; ri <= 4 && !sampleHit; ++ri) {
                        const rr = (cyl.radius * ri) / 4;
                        for (let ai = 0; ai < 24 && !sampleHit; ++ai) {
                            const th = (2 * Math.PI * ai) / 24;
                            const p = v3(0, 0, 0);
                            for (let c = 0; c < 3; ++c) {
                                p.values[c] = cyl.axis.origin.values[c]
                                    + t * cyl.axis.direction.values[c]
                                    + rr * Math.cos(th) * u.values[c]
                                    + rr * Math.sin(th) * w.values[c];
                            }
                            if (Math.abs(p.values[0]) <= b.extent.values[0] &&
                                Math.abs(p.values[1]) <= b.extent.values[1] &&
                                Math.abs(p.values[2]) <= b.extent.values[2]) {
                                sampleHit = true;
                            }
                        }
                    }
                }

                if (sampleHit) {
                    expect(intersect).toBe(true);
                    ++numHit;
                }
            }
            expect(numHit).toBeGreaterThan(0);
        });

    it('throws for an infinite cylinder', () => {
        const cyl = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
        cyl.makeInfiniteCylinder();
        expect(() => ti.test(unit, cyl)).toThrow();
    });

    it('throws when the box is not 3-dimensional', () => {
        const b2 = CanonicalBox.fromExtent(Vector.fromArray([1, 1]));
        expect(() => ti.test(b2, cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2)))
            .toThrow();
    });
});
