import { describe, it, expect } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { Cylinder3 } from '../src/Cylinder3.js';
import { IntrCanonicalBox3Cylinder3TI } from '../src/IntrCanonicalBox3Cylinder3.js';
import { IntrOrientedBox3Cylinder3TI } from '../src/IntrOrientedBox3Cylinder3.js';
import { Line } from '../src/Line.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function cylinder(origin: Vector, direction: Vector, radius: number,
    height: number): Cylinder3 {
    const d = direction.clone();
    normalize(d);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(origin, d), radius, height);
}

function orthonormalFrame(w: Vector): Vector[] {
    const v = [w.clone(), Vector.zero(3), Vector.zero(3)];
    normalize(v[0]);
    computeOrthogonalComplement3(1, v, false);
    return [v[0], v[1], v[2]];
}

// Map a point of the canonical (origin-centered, axis-aligned) frame into the
// world frame defined by 'center' and the orthonormal 'axes'.
function toWorld(p: Vector, center: Vector, axes: Vector[]): Vector {
    let q = center.clone();
    for (let d = 0; d < 3; ++d) {
        q = add(q, mul(p.values[d], axes[d]));
    }
    return q;
}

function rotateDirection(d: Vector, axes: Vector[]): Vector {
    let q = Vector.zero(3);
    for (let i = 0; i < 3; ++i) {
        q = add(q, mul(d.values[i], axes[i]));
    }
    return q;
}

describe('IntrOrientedBox3Cylinder3TI', () => {
    const query = new IntrOrientedBox3Cylinder3TI();
    const cbQuery = new IntrCanonicalBox3Cylinder3TI();
    const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

    it('matches the canonical-box query when the box is axis-aligned at the origin', () => {
        const extent = vec(1, 2, 3);
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0), unitAxes,
            extent);
        const cbox = CanonicalBox.fromExtent(extent);

        const cases: Cylinder3[] = [
            cylinder(vec(0, 0, 0), vec(0, 0, 1), 0.5, 2),
            cylinder(vec(10, 0, 0), vec(0, 0, 1), 0.5, 2),
            cylinder(vec(1.4, 0, 0), vec(0, 0, 1), 0.5, 2),
            cylinder(vec(1.6, 0, 0), vec(0, 0, 1), 0.5, 2),
            cylinder(vec(0, 0, 5), vec(1, 1, 1), 1, 4),
            cylinder(vec(3, 4, 5), vec(1, -1, 0), 0.25, 1)
        ];
        for (const c of cases) {
            expect(query.test(box, c).intersect)
                .toBe(cbQuery.test(cbox, c).intersect);
        }
    });

    it('detects the obvious inside and far-away configurations', () => {
        const box = OrientedBox.fromCenterAxisExtent(vec(5, -3, 2), unitAxes,
            vec(1, 1, 1));
        // A cylinder centered inside the box.
        expect(query.test(box,
            cylinder(vec(5, -3, 2), vec(0, 0, 1), 0.2, 0.5)).intersect)
            .toBe(true);
        // A cylinder far away.
        expect(query.test(box,
            cylinder(vec(100, 100, 100), vec(0, 0, 1), 1, 1)).intersect)
            .toBe(false);
    });

    it('is invariant under a rigid motion of both objects', () => {
        const axes = orthonormalFrame(vec(0.4, -0.6, 0.7));
        const center = vec(2, -1, 4);
        const extent = vec(1, 2, 0.5);
        const alignedBox = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0),
            unitAxes, extent);
        const movedBox = OrientedBox.fromCenterAxisExtent(center, axes,
            extent);

        let seed = 314159;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        let numHits = 0, numMisses = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const o = vec(rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4);
            const d = vec(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1);
            if (d.values.every(v => Math.abs(v) < 1e-3)) {
                continue;
            }
            const radius = 0.2 + rand() * 1.5;
            const height = 0.2 + rand() * 3;

            const c0 = cylinder(o, d, radius, height);
            const c1 = cylinder(toWorld(c0.axis.origin, center, axes),
                rotateDirection(c0.axis.direction, axes), radius, height);

            const r0 = query.test(alignedBox, c0).intersect;
            const r1 = query.test(movedBox, c1).intersect;
            expect(r1).toBe(r0);
            if (r0) {
                ++numHits;
            } else {
                ++numMisses;
            }
        }
        expect(numHits).toBeGreaterThan(5);
        expect(numMisses).toBeGreaterThan(5);
    });

    it('rejects infinite cylinders and non-3D boxes', () => {
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0), unitAxes,
            vec(1, 1, 1));
        const infinite = cylinder(vec(0, 0, 0), vec(0, 0, 1), 1, -1);
        expect(() => query.test(box, infinite)).toThrow();

        const box2 = new OrientedBox(2);
        expect(() => query.test(box2,
            cylinder(vec(0, 0, 0), vec(0, 0, 1), 1, 1))).toThrow();
    });
});
