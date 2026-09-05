import { describe, it, expect } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, sub, dot, length, normalize } from '../src/Vector.js';
import {
    IntrCanonicalBox3Cylinder3TI
} from '../src/IntrCanonicalBox3Cylinder3.js';

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

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, unitVector, wellScaledVector, seededRandom
} from './helpers/arbitraries.js';
import { add, mul } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

const canonBox = fc.tuple(fc.double({ min: 0.2, max: 3, noNaN: true }),
    fc.double({ min: 0.2, max: 3, noNaN: true }),
    fc.double({ min: 0.2, max: 3, noNaN: true }))
    .map(([e0, e1, e2]) =>
        CanonicalBox.fromExtent(Vector.fromArray([e0, e1, e2])));

const cylV31 = fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
    fc.double({ min: 0.1, max: 2, noNaN: true }),
    fc.double({ min: 0.2, max: 4, noNaN: true }))
    .map(([o, w, r, h]) => Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(o, w), r, h));

function inBox(box: CanonicalBox, p: Vector, tol: number): boolean {
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(p.values[i]) > box.extent.values[i] + tol) {
            return false;
        }
    }
    return true;
}

function inCyl(c: Cylinder3, p: Vector, tol: number): boolean {
    const diff = sub(p, c.axis.origin);
    const z = dot(diff, c.axis.direction);
    if (Math.abs(z) > 0.5 * c.height + tol) {
        return false;
    }
    const radial = length(sub(diff, mul(z, c.axis.direction)));
    return radial <= c.radius + tol;
}

describe('IntrCanonicalBox3Cylinder3 verification', () => {
    const tiq = new IntrCanonicalBox3Cylinder3TI();

    it('reports an intersection whenever a sampled common point exists', () => {
        const rnd = seededRandom(0x64bd2ae1);
        for (let trial = 0; trial < 150; ++trial) {
            const box = CanonicalBox.fromExtent(Vector.fromArray(
                [0.3 + rnd() * 2, 0.3 + rnd() * 2, 0.3 + rnd() * 2]));
            const w = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (length(w) < 0.3) {
                continue;
            }
            normalize(w);
            const c = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(Vector.fromArray([rnd() * 6 - 3,
                    rnd() * 6 - 3, rnd() * 6 - 3]), w),
                0.2 + rnd() * 1.2, 0.4 + rnd() * 2.5);
            const got = tiq.test(box, c).intersect;
            // Sample the solid cylinder and look for a point inside the box.
            const basis = [w.clone(), new Vector(3), new Vector(3)];
            computeOrthogonalComplement3(1, basis);
            let common = false;
            const half = 0.5 * c.height;
            for (let i = 0; i <= 12 && !common; ++i) {
                const z = -half + (2 * half * i) / 12;
                for (let j = 0; j <= 6 && !common; ++j) {
                    const rr = (c.radius * j) / 6;
                    for (let k = 0; k < 24; ++k) {
                        const a = (2 * Math.PI * k) / 24;
                        const p = add(c.axis.origin, add(mul(z, basis[0]),
                            add(mul(rr * Math.cos(a), basis[1]),
                                mul(rr * Math.sin(a), basis[2]))));
                        if (inBox(box, p, 0)) {
                            common = true;
                            break;
                        }
                    }
                }
            }
            if (!common) {
                // Also sample the box for a point inside the cylinder.
                for (let i = 0; i <= 12 && !common; ++i) {
                    for (let j = 0; j <= 12 && !common; ++j) {
                        for (let k = 0; k <= 12; ++k) {
                            const p = Vector.fromArray([
                                box.extent.values[0] * (2 * i / 12 - 1),
                                box.extent.values[1] * (2 * j / 12 - 1),
                                box.extent.values[2] * (2 * k / 12 - 1)]);
                            if (inCyl(c, p, 0)) {
                                common = true;
                                break;
                            }
                        }
                    }
                }
            }
            if (common) {
                expect(got).toBe(true);
            }
        }
    }, 30000);

    it('is invariant under reflections of the coordinate axes', () => {
        check(fc.tuple(canonBox, cylV31, fc.boolean(), fc.boolean(),
            fc.boolean()), ([box, c, f0, f1, f2]) => {
            const sgn = [f0 ? -1 : 1, f1 ? -1 : 1, f2 ? -1 : 1];
            const refl = (v: Vector): Vector => Vector.fromArray(
                [0, 1, 2].map(i => sgn[i] * v.values[i]));
            const c2 = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(refl(c.axis.origin),
                    refl(c.axis.direction)), c.radius, c.height);
            expect(tiq.test(box, c).intersect)
                .toBe(tiq.test(box, c2).intersect);
        });
    });

    it('is invariant under a coordinate permutation of box and cylinder', () => {
        check(fc.tuple(canonBox, cylV31, fc.constantFrom(
            [0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1], [1, 0, 2], [2, 1, 0])),
            ([box, c, perm]) => {
                const permute = (v: Vector): Vector =>
                    Vector.fromArray(perm.map(i => v.values[i]));
                const box2 = CanonicalBox.fromExtent(permute(box.extent));
                const c2 = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(permute(c.axis.origin),
                        permute(c.axis.direction)), c.radius, c.height);
                expect(tiq.test(box, c).intersect)
                    .toBe(tiq.test(box2, c2).intersect);
            });
    });

    it('a cylinder whose axis passes through the box intersects', () => {
        check(fc.tuple(canonBox, unitVector(3),
            fc.double({ min: 0.1, max: 2, noNaN: true }),
            fc.double({ min: -0.8, max: 0.8, noNaN: true }),
            fc.double({ min: -0.8, max: 0.8, noNaN: true }),
            fc.double({ min: -0.8, max: 0.8, noNaN: true })),
            ([box, w, r, s0, s1, s2]) => {
                const inside = Vector.fromArray([
                    s0 * box.extent.values[0], s1 * box.extent.values[1],
                    s2 * box.extent.values[2]]);
                const c = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(inside, w), r, 1);
                expect(tiq.test(box, c).intersect).toBe(true);
            });
    });

    it('a cylinder translated far away does not intersect', () => {
        check(fc.tuple(canonBox, cylV31, unitVector(3)), ([box, c, d]) => {
            const reach = length(box.extent) + c.radius + 0.5 * c.height + 1;
            const moved = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(mul(2 * reach, d),
                    c.axis.direction), c.radius, c.height);
            expect(tiq.test(box, moved).intersect).toBe(false);
        });
    });

    it('a cylinder that contains the box intersects', () => {
        check(fc.tuple(canonBox, unitVector(3)), ([box, w]) => {
            const reach = length(box.extent);
            const c = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(Vector.fromArray([0, 0, 0]), w),
                reach + 1, 2 * reach + 2);
            expect(tiq.test(box, c).intersect).toBe(true);
        });
    });
});
