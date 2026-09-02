import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { CanonicalBox } from '../src/CanonicalBox';
import { Cylinder3 } from '../src/Cylinder3';
import { Line } from '../src/Line';
import { Vector, normalize, sub } from '../src/Vector';
import { IntrAlignedBox3Cylinder3TI } from '../src/IntrAlignedBox3Cylinder3';
import { IntrCanonicalBox3Cylinder3TI } from '../src/IntrCanonicalBox3Cylinder3';

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
