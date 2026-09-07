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
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0),
            [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)], vec(1, 1, 1));
        const infinite = cylinder(vec(0, 0, 0), vec(0, 0, 1), 1, -1);
        expect(() => query.test(box, infinite)).toThrow();

        const box2 = new OrientedBox(2);
        expect(() => query.test(box2,
            cylinder(vec(0, 0, 0), vec(0, 0, 1), 1, 1))).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrOrientedBox3Cylinder3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { dot, sub } from '../src/Vector.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no frame component is a
// subnormal that would underflow when squared.
const angleB = () => wellScaled(-Math.PI, Math.PI);

// R = Rz(a)*Ry(b)*Rx(c); the columns are the box axes.
function rotFrameB(a: number, b: number, c: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const cc = Math.cos(c), sc = Math.sin(c);
    return [
        vec(ca * cb, sa * cb, -sb),
        vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc, cb * sc),
        vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc, cb * cc)];
}

function unitDirB(th: number, ph: number): Vector {
    const d = vec(Math.cos(th) * Math.cos(ph), Math.sin(th) * Math.cos(ph),
        Math.sin(ph));
    normalize(d);
    return d;
}

const boxCylinder = fc.tuple(
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    angleB(), angleB(), angleB(),
    fc.double({ min: 0.3, max: 2.5, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.3, max: 2.5, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.3, max: 2.5, noNaN: true, noDefaultInfinity: true }),
    wellScaled(-4, 4), wellScaled(-4, 4), wellScaled(-4, 4),
    angleB(), angleB(),
    fc.double({ min: 0.2, max: 2, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.4, max: 4, noNaN: true, noDefaultInfinity: true })
).map(([cx, cy, cz, a, b, c, e0, e1, e2, ox, oy, oz, th, ph, r, h]) => ({
    box: OrientedBox.fromCenterAxisExtent(vec(cx, cy, cz), rotFrameB(a, b, c),
        vec(e0, e1, e2)),
    cyl: cylinder(vec(ox, oy, oz), unitDirB(th, ph), r, h)
}));

function insideBox(box: OrientedBox, p: Vector): boolean {
    const diff = sub(p, box.center);
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(dot(diff, box.axis[i])) > box.extent.values[i]) {
            return false;
        }
    }
    return true;
}

function insideCylinder(c: Cylinder3, p: Vector): boolean {
    const diff = sub(p, c.axis.origin);
    const z = dot(diff, c.axis.direction);
    if (Math.abs(z) > 0.5 * c.height) { return false; }
    const radial = sub(diff, mul(z, c.axis.direction));
    return dot(radial, radial) <= c.radius * c.radius;
}

describe('IntrOrientedBox3Cylinder3 verification', () => {
    const q = new IntrOrientedBox3Cylinder3TI();
    const cbq = new IntrCanonicalBox3Cylinder3TI();

    it('matches the canonical-box query on the transformed configuration',
        () => {
            check(boxCylinder, ({ box, cyl }) => {
                const cbox = CanonicalBox.fromExtent(box.extent);
                const diff = sub(cyl.axis.origin, box.center);
                const o = vec(dot(box.axis[0], diff), dot(box.axis[1], diff),
                    dot(box.axis[2], diff));
                const d = vec(dot(box.axis[0], cyl.axis.direction),
                    dot(box.axis[1], cyl.axis.direction),
                    dot(box.axis[2], cyl.axis.direction));
                const tc = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(o, d), cyl.radius, cyl.height);
                expect(q.test(box, cyl).intersect)
                    .toBe(cbq.test(cbox, tc).intersect);
            });
        });

    it('any sampled common point of the two solids forces intersect', () => {
        const rnd = seededRandom(0x8ac31f7);
        let witnessed = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const box = OrientedBox.fromCenterAxisExtent(
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rotFrameB(rnd() * 6, rnd() * 6, rnd() * 6),
                vec(0.3 + rnd() * 2, 0.3 + rnd() * 2, 0.3 + rnd() * 2));
            const axis = unitDirB(rnd() * 6, rnd() * 6 - 3);
            const cyl = cylinder(
                vec(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3), axis,
                0.2 + rnd() * 1.5, 0.5 + rnd() * 3);

            // Sample the solid cylinder and look for a point inside the box.
            let common = false;
            for (let k = 0; k < 400 && !common; ++k) {
                const z = (rnd() - 0.5) * cyl.height;
                const rad = cyl.radius * Math.sqrt(rnd());
                const ang = rnd() * 2 * Math.PI;
                const basis = [axis.clone(), Vector.zero(3), Vector.zero(3)];
                computeOrthogonalComplement3(1, basis, false);
                const p = add(cyl.axis.origin,
                    add(mul(z, axis),
                        add(mul(rad * Math.cos(ang), basis[1]),
                            mul(rad * Math.sin(ang), basis[2]))));
                if (insideBox(box, p)) { common = true; }
            }
            if (common) {
                expect(q.test(box, cyl).intersect).toBe(true);
                ++witnessed;
            }
        }
        expect(witnessed).toBeGreaterThan(20);
    }, 60000);

    it('a cylinder centred on a sampled box point always intersects', () => {
        check(fc.tuple(boxCylinder,
            fc.double({ min: -0.9, max: 0.9, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true,
                noDefaultInfinity: true })),
        ([{ box, cyl }, u0, u1, u2]) => {
            const p = add(box.center,
                add(mul(u0 * box.extent.values[0], box.axis[0]),
                    add(mul(u1 * box.extent.values[1], box.axis[1]),
                        mul(u2 * box.extent.values[2], box.axis[2]))));
            const moved = cylinder(p, cyl.axis.direction, cyl.radius,
                cyl.height);
            expect(insideBox(box, p)).toBe(true);
            expect(insideCylinder(moved, p)).toBe(true);
            expect(q.test(box, moved).intersect).toBe(true);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(boxCylinder, angleB(), angleB(), angleB(),
            wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3)),
        ([{ box, cyl }, a1, a2, a3, tx, ty, tz]) => {
            const fr = rotFrameB(a1, a2, a3);
            const rot = (v: Vector): Vector => add(mul(v.get(0), fr[0]),
                add(mul(v.get(1), fr[1]), mul(v.get(2), fr[2])));
            const t = vec(tx, ty, tz);
            const xf = (v: Vector): Vector => add(rot(v), t);
            const box2 = OrientedBox.fromCenterAxisExtent(xf(box.center),
                [rot(box.axis[0]), rot(box.axis[1]), rot(box.axis[2])],
                box.extent);
            const cyl2 = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(xf(cyl.axis.origin),
                    rot(cyl.axis.direction)), cyl.radius, cyl.height);
            const a = q.test(box, cyl).intersect;
            const b = q.test(box2, cyl2).intersect;
            if (a === b) { return; }
            // A disagreement is allowed only for a configuration that is
            // grazing: shrinking or growing the cylinder radius by 1e-6 must
            // change the original answer too.
            const shrunk = cylinder(cyl.axis.origin, cyl.axis.direction,
                Math.max(cyl.radius - 1e-6, 1e-9), cyl.height);
            const grown = cylinder(cyl.axis.origin, cyl.axis.direction,
                cyl.radius + 1e-6, cyl.height);
            expect(q.test(box, shrunk).intersect
                !== q.test(box, grown).intersect).toBe(true);
        });
    });

    it('separates a cylinder moved beyond the combined reach', () => {
        check(boxCylinder, ({ box, cyl }) => {
            const far = cylinder(add(cyl.axis.origin, vec(1000, 0, 0)),
                cyl.axis.direction, cyl.radius, cyl.height);
            expect(q.test(box, far).intersect).toBe(false);
        });
    });

    it('rejects an infinite cylinder and a non-3D box', () => {
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0),
            [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)], vec(1, 1, 1));
        const c = cylinder(vec(0, 0, 0), vec(0, 0, 1), 1, 2);
        c.makeInfiniteCylinder();
        expect(() => q.test(box, c)).toThrow(
            'Infinite cylinders are not yet supported.');

        const box2 = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(() => q.test(box2,
            cylinder(vec(0, 0, 0), vec(0, 0, 1), 1, 2))).toThrow();
    });
});
