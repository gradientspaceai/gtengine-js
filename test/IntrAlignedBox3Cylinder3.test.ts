import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, normalize, sub } from '../src/Vector.js';
import { IntrAlignedBox3Cylinder3TI } from '../src/IntrAlignedBox3Cylinder3.js';
import { IntrCanonicalBox3Cylinder3TI } from '../src/IntrCanonicalBox3Cylinder3.js';
import { add, dot, mul } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import {
    check, fc, positive, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function box(min: Vector, max: Vector): AlignedBox {
    return AlignedBox.fromMinMax(min, max);
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

// The delegation the ported query performs, computed independently here.
function viaCanonical(b: AlignedBox, cyl: Cylinder3): boolean {
    const { center, extent } = b.getCenteredForm();
    const cbox = CanonicalBox.fromExtent(extent);
    const translated = cyl.clone();
    translated.axis.origin = sub(translated.axis.origin, center);
    return new IntrCanonicalBox3Cylinder3TI().test(cbox, translated).intersect;
}

describe('IntrAlignedBox3Cylinder3TI', () => {
    const query = new IntrAlignedBox3Cylinder3TI();

    it('reports intersection for a cylinder inside a box', () => {
        const b = box(v3(-1, -1, -1), v3(1, 1, 1));
        const cyl = cylinder(v3(0, 0, 0), v3(0, 0, 1), 0.5, 1);
        expect(query.test(b, cyl).intersect).toBe(true);
    });

    it('reports no intersection for a distant cylinder', () => {
        const b = box(v3(-1, -1, -1), v3(1, 1, 1));
        const cyl = cylinder(v3(10, 10, 10), v3(0, 0, 1), 0.5, 1);
        expect(query.test(b, cyl).intersect).toBe(false);
    });

    it('is translation invariant (box and cylinder moved together)', () => {
        const b = box(v3(-1, -2, -3), v3(2, 1, 4));
        const cyl = cylinder(v3(1.5, 0.5, 3.5), v3(1, 1, 1), 0.75, 2);
        const shift = v3(100, -50, 25);
        const b2 = box(
            Vector.fromArray([b.min.values[0] + shift.values[0],
                b.min.values[1] + shift.values[1],
                b.min.values[2] + shift.values[2]]),
            Vector.fromArray([b.max.values[0] + shift.values[0],
                b.max.values[1] + shift.values[1],
                b.max.values[2] + shift.values[2]]));
        const cyl2 = cyl.clone();
        cyl2.axis.origin = Vector.fromArray([
            cyl.axis.origin.values[0] + shift.values[0],
            cyl.axis.origin.values[1] + shift.values[1],
            cyl.axis.origin.values[2] + shift.values[2]]);
        expect(query.test(b2, cyl2).intersect).toBe(query.test(b, cyl).intersect);
    });

    it('separates a cylinder that just misses a face', () => {
        // Axis along z, so the cylinder is a disk of radius 1 in xy; its
        // closest approach to the box face x = 1 is at x = 2.5 - 1 = 1.5.
        const b = box(v3(-1, -1, -1), v3(1, 1, 1));
        const cyl = cylinder(v3(2.5, 0, 0), v3(0, 0, 1), 1, 2);
        expect(query.test(b, cyl).intersect).toBe(false);

        // Move it so that the disk boundary touches x = 1.
        const touching = cylinder(v3(2, 0, 0), v3(0, 0, 1), 1, 2);
        expect(query.test(b, touching).intersect).toBe(true);
    });

    it('handles a cylinder outside the box slab', () => {
        const b = box(v3(-1, -1, -1), v3(1, 1, 1));
        const cyl = cylinder(v3(0, 0, 5), v3(0, 0, 1), 0.5, 2);
        expect(query.test(b, cyl).intersect).toBe(false);
    });

    it('throws for an infinite cylinder', () => {
        const b = box(v3(-1, -1, -1), v3(1, 1, 1));
        const cyl = cylinder(v3(0, 0, 0), v3(0, 0, 1), 0.5, 1);
        cyl.makeInfiniteCylinder();
        expect(() => query.test(b, cyl)).toThrow();
    });

    it('throws for a non-3D box', () => {
        const b = AlignedBox.fromMinMax(Vector.fromArray([-1, -1]),
            Vector.fromArray([1, 1]));
        const cyl = cylinder(v3(0, 0, 0), v3(0, 0, 1), 0.5, 1);
        expect(() => query.test(b, cyl)).toThrow();
    });

    it('matches the canonical-box query after the transform (randomized)', () => {
        const rand = makeRandom(20250901);
        let mismatches = 0;
        let numIntersect = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const cx = 4 * rand() - 2, cy = 4 * rand() - 2, cz = 4 * rand() - 2;
            const ex = 0.25 + rand(), ey = 0.25 + rand(), ez = 0.25 + rand();
            const b = box(v3(cx - ex, cy - ey, cz - ez),
                v3(cx + ex, cy + ey, cz + ez));
            const cyl = cylinder(
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2),
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1),
                0.1 + rand(), 0.2 + 2 * rand());
            const actual = query.test(b, cyl).intersect;
            if (actual !== viaCanonical(b, cyl)) {
                ++mismatches;
            }
            if (actual) {
                ++numIntersect;
            }
        }
        expect(mismatches).toBe(0);
        // Sanity: the random configurations exercise both outcomes.
        expect(numIntersect).toBeGreaterThan(20);
        expect(numIntersect).toBeLessThan(380);
    });

    it('does not modify the input cylinder', () => {
        const b = box(v3(-1, -1, -1), v3(3, 3, 3));
        const cyl = cylinder(v3(1, 1, 1), v3(0, 0, 1), 0.5, 1);
        query.test(b, cyl);
        expect(cyl.axis.origin.values).toEqual([1, 1, 1]);
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrAlignedBox3Cylinder3.h.
// ---------------------------------------------------------------------------

describe('IntrAlignedBox3Cylinder3 verification', () => {
    const ti = new IntrAlignedBox3Cylinder3TI();
    const cbTi = new IntrCanonicalBox3Cylinder3TI();

    const arbBox = fc.tuple(wellScaledVector(3, -3, 3),
        fc.array(positive(2, 0.2), { minLength: 3, maxLength: 3 }))
        .map(([c, e]) => AlignedBox.fromMinMax(
            Vector.fromArray([c.get(0) - e[0], c.get(1) - e[1],
                c.get(2) - e[2]]),
            Vector.fromArray([c.get(0) + e[0], c.get(1) + e[1],
                c.get(2) + e[2]])));
    const arbCylinder = fc.tuple(wellScaledVector(3, -3, 3), unitVector(3),
        positive(2, 0.25), positive(3, 0.5))
        .map(([o, d, r, h]) => Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(o, d), r, h));
    const arbPair = fc.tuple(arbBox, arbCylinder);

    // Signed containment functions: negative strictly inside the solid.
    function inBox(p: Vector, b: AlignedBox): number {
        let worst = -Number.MAX_VALUE;
        for (let i = 0; i < 3; ++i) {
            worst = Math.max(worst, b.min.get(i) - p.get(i),
                p.get(i) - b.max.get(i));
        }
        return worst;
    }

    function inCylinder(p: Vector, cyl: Cylinder3): number {
        const diff = sub(p, cyl.axis.origin);
        const z = dot(diff, cyl.axis.direction);
        const radial = sub(diff, mul(cyl.axis.direction, z));
        return Math.max(Math.sqrt(dot(radial, radial)) - cyl.radius,
            Math.abs(z) - 0.5 * cyl.height);
    }

    it('reports an intersection when a sampled box point is well inside the'
        + ' cylinder', () => {
        check(arbPair, ([b, cyl]) => {
            const n = 10;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    for (let k = 0; k <= n; ++k) {
                        const p = Vector.fromArray([
                            b.min.get(0) + (i / n) * (b.max.get(0)
                                - b.min.get(0)),
                            b.min.get(1) + (j / n) * (b.max.get(1)
                                - b.min.get(1)),
                            b.min.get(2) + (k / n) * (b.max.get(2)
                                - b.min.get(2))]);
                        if (inCylinder(p, cyl) < -1e-6) {
                            expect(ti.test(b, cyl).intersect).toBe(true);
                            return;
                        }
                    }
                }
            }
        }, 60);
    });

    it('reports an intersection when a sampled cylinder point is well inside'
        + ' the box', () => {
        check(arbPair, ([b, cyl]) => {
            const basis = [cyl.axis.direction.clone(), Vector.zero(3),
                Vector.zero(3)];
            computeOrthogonalComplement3(1, basis);
            for (let i = 0; i <= 6; ++i) {
                for (let j = 0; j <= 6; ++j) {
                    for (let k = 0; k < 12; ++k) {
                        const z = (i / 6 - 0.5) * cyl.height;
                        const rho = (j / 6) * cyl.radius;
                        const a = (2 * Math.PI * k) / 12;
                        const p = add(cyl.axis.origin,
                            add(mul(z, basis[0]),
                                add(mul(rho * Math.cos(a), basis[1]),
                                    mul(rho * Math.sin(a), basis[2]))));
                        if (inBox(p, b) < -1e-6) {
                            expect(ti.test(b, cyl).intersect).toBe(true);
                            return;
                        }
                    }
                }
            }
        }, 60);
    });

    it('reports no intersection when the box is outside the cylinder slab',
        () => {
            // A sound separation certificate: every box vertex is strictly on
            // the same side of a cylinder end plane, or every box vertex is
            // farther than the radius from the axis and the whole box is on
            // one side of a plane through the axis (so the nearest point is a
            // vertex). The first certificate alone is enough here.
            check(arbPair, ([b, cyl]) => {
                const vertices = b.getVertices();
                let allAbove = true, allBelow = true;
                for (const v of vertices) {
                    const z = dot(sub(v, cyl.axis.origin), cyl.axis.direction);
                    allAbove = allAbove && z > 0.5 * cyl.height + 1e-9;
                    allBelow = allBelow && z < -0.5 * cyl.height - 1e-9;
                }
                if (allAbove || allBelow) {
                    expect(ti.test(b, cyl).intersect).toBe(false);
                }
            });
        });

    it('agrees with the canonical-box query on the translated problem', () => {
        check(arbPair, ([b, cyl]) => {
            const { center, extent } = b.getCenteredForm();
            const cbox = CanonicalBox.fromExtent(extent);
            const moved = cyl.clone();
            moved.axis.origin = sub(cyl.axis.origin, center);
            expect(ti.test(b, cyl).intersect)
                .toBe(cbTi.test(cbox, moved).intersect);
        });
    });

    it('does not modify the caller cylinder', () => {
        check(arbPair, ([b, cyl]) => {
            const before = cyl.axis.origin.clone();
            ti.test(b, cyl);
            expect(cyl.axis.origin.equals(before)).toBe(true);
        });
    });

    it('is equivariant under translation', () => {
        check(fc.tuple(arbPair, wellScaledVector(3, -4, 4)),
            ([[b, cyl], shift]) => {
                const b2 = AlignedBox.fromMinMax(add(b.min, shift),
                    add(b.max, shift));
                const cyl2 = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(add(cyl.axis.origin, shift),
                        cyl.axis.direction.clone()),
                    cyl.radius, cyl.height);
                expect(ti.test(b2, cyl2).intersect)
                    .toBe(ti.test(b, cyl).intersect);
            });
    });

    it('rejects infinite cylinders and non-3D boxes', () => {
        const infinite = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(v3(0, 0, 0), v3(0, 0, 1)), 1, -1);
        const unit = box(v3(-1, -1, -1), v3(1, 1, 1));
        expect(() => ti.test(unit, infinite))
            .toThrow('Infinite cylinders are not yet supported.');
        const flat = AlignedBox.fromMinMax(Vector.zero(2), Vector.zero(2));
        expect(() => ti.test(flat, cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2)))
            .toThrow('mismatched sizes');
    });
});
