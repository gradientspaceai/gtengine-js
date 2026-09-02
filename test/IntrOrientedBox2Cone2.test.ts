import { describe, it, expect } from 'vitest';
import { Cone } from '../src/Cone';
import { IntrRay2OrientedBox2TI } from '../src/IntrRay2OrientedBox2';
import {
    IntrOrientedBox2Cone2TI,
    defaultIntrOrientedBox2Cone2TIResult
} from '../src/IntrOrientedBox2Cone2';
import { OrientedBox } from '../src/OrientedBox';
import { Ray } from '../src/Ray';
import { Vector, dot, normalize, sub } from '../src/Vector';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function obox(center: number[], angle: number, extent: number[]): OrientedBox {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center),
        [vec(c, s), vec(-s, c)], Vector.fromArray(extent));
}

function cone2(origin: number[], direction: number[], angle: number): Cone {
    const c = new Cone(2);
    const d = Vector.fromArray(direction);
    normalize(d);
    c.ray = Ray.fromOriginDirection(Vector.fromArray(origin), d);
    c.setAngle(angle);
    return c;
}

// An independent evaluation of the documented criterion in world coordinates:
// the cone axis ray meets the box, or some box vertex P satisfies
// Dot(D,P-V) > 0 and (Dot(D,P-V))^2 > |P-V|^2 * cosAngle^2.
function reference(box: OrientedBox, cone: Cone): boolean {
    const rbQuery = new IntrRay2OrientedBox2TI();
    if (rbQuery.test(cone.ray, box).intersect) {
        return true;
    }
    for (const P of box.getVertices()) {
        const diff = sub(P, cone.ray.origin);
        const num = dot(cone.ray.direction, diff);
        if (num > 0) {
            if (num * num > dot(diff, diff) * cone.cosAngle * cone.cosAngle) {
                return true;
            }
        }
    }
    return false;
}

const ti = new IntrOrientedBox2Cone2TI();

describe('IntrOrientedBox2Cone2', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrOrientedBox2Cone2TIResult().intersect).toBe(false);
    });

    it('detects a box straddling the cone axis', () => {
        const C = cone2([0, 0], [1, 0], Math.PI / 6);
        expect(ti.test(obox([5, 0], 0, [1, 1]), C).intersect).toBe(true);
    });

    it('rejects a box behind the cone apex', () => {
        const C = cone2([0, 0], [1, 0], Math.PI / 6);
        expect(ti.test(obox([-5, 0], 0, [1, 1]), C).intersect).toBe(false);
    });

    it('rejects a box outside the cone angle and accepts one inside', () => {
        // A 30-degree half-angle cone along +x. At x = 5 the cone half-width
        // is 5*tan(30 deg) = 2.887.
        const C = cone2([0, 0], [1, 0], Math.PI / 6);
        // A small box centered well above the cone.
        expect(ti.test(obox([5, 6], 0, [0.5, 0.5]), C).intersect).toBe(false);
        // A small box centered inside the cone.
        expect(ti.test(obox([5, 2], 0, [0.2, 0.2]), C).intersect).toBe(true);
    });

    it('does not report a box that only touches the cone boundary', () => {
        // The 45-degree cone boundary through the origin is the line y = x.
        // A box whose only cone-side corner is exactly on y = x touches but
        // does not overlap.
        const C = cone2([0, 0], [1, 0], Math.PI / 4);
        const B = obox([3, 5], 0, [1, 1]);  // corner (4, 4) is on y = x
        expect(ti.test(B, C).intersect).toBe(false);
        // Nudging the box towards the axis produces an overlap.
        expect(ti.test(obox([3, 4.99], 0, [1, 1]), C).intersect).toBe(true);
    });

    it('accounts for the box orientation', () => {
        const C = cone2([0, 0], [1, 0], Math.PI / 12);  // 15 degrees
        // An axis-aligned square well above the narrow cone misses it.
        expect(ti.test(obox([10, 4], 0, [1, 1]), C).intersect).toBe(false);
        // Rotating it 45 degrees pushes a corner down to y = 4 - sqrt(2),
        // which is inside the cone (10*tan(15 deg) = 2.679).
        expect(ti.test(obox([10, 4], Math.PI / 4, [1, 1]), C).intersect)
            .toBe(true);
    });

    it('agrees with a world-coordinate evaluation on random inputs', () => {
        let state = 30313;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const C = cone2([rand() * 2, rand() * 2],
                [rand(), rand() + 0.001],
                0.15 + Math.abs(rand()) * 1.2);
            const B = obox([rand() * 6, rand() * 6], rand() * Math.PI,
                [0.3 + Math.abs(rand()) * 2, 0.3 + Math.abs(rand()) * 2]);
            const actual = ti.test(B, C).intersect;
            expect(actual).toBe(reference(B, C));
            if (actual) {
                ++numHits;
            }
        }
        expect(numHits).toBeGreaterThan(50);
        expect(numHits).toBeLessThan(450);
    });
});
