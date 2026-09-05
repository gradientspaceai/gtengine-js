import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import { IntrHalfspace3OrientedBox3TI } from '../src/IntrHalfspace3OrientedBox3.js';
import { length } from '../src/Vector.js';
import { check, fc, positive, rotationFrame, unitVector, wellScaled } from './helpers/arbitraries.js';

function halfspace(nx: number, ny: number, nz: number, c: number): Halfspace {
    const n = Vector.fromArray([nx, ny, nz]);
    normalize(n);
    return Halfspace.fromNormalConstant(n, c);
}

function rotationAxes(ax: number, ay: number, az: number): Vector[] {
    const cx = Math.cos(ax), sx = Math.sin(ax);
    const cy = Math.cos(ay), sy = Math.sin(ay);
    const cz = Math.cos(az), sz = Math.sin(az);
    const m = [
        [cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx],
        [sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx],
        [-sy, cy * sx, cy * cx]
    ];
    return [
        Vector.fromArray([m[0][0], m[1][0], m[2][0]]),
        Vector.fromArray([m[0][1], m[1][1], m[2][1]]),
        Vector.fromArray([m[0][2], m[1][2], m[2][2]])
    ];
}

function box(center: number[], axes: Vector[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center), axes,
        Vector.fromArray(extent));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrHalfspace3OrientedBox3', () => {
    const ti = new IntrHalfspace3OrientedBox3TI();
    const identity = rotationAxes(0, 0, 0);

    it('detects an axis-aligned box inside, outside and straddling', () => {
        const h = halfspace(0, 0, 1, 0);  // z >= 0
        expect(ti.test(h, box([0, 0, 5], identity, [1, 1, 1])).intersect).toBe(true);
        expect(ti.test(h, box([0, 0, -5], identity, [1, 1, 1])).intersect).toBe(false);
        expect(ti.test(h, box([0, 0, 0], identity, [1, 1, 1])).intersect).toBe(true);
    });

    it('treats face contact as an intersection and finds the exact threshold', () => {
        const h = halfspace(0, 0, 1, 0);
        expect(ti.test(h, box([0, 0, -1], identity, [1, 1, 1])).intersect).toBe(true);
        expect(ti.test(h, box([0, 0, -1.0000001], identity, [1, 1, 1])).intersect)
            .toBe(false);
    });

    it('uses the projected radius of a rotated box', () => {
        // A unit cube rotated 45 degrees about z has half-width sqrt(2) along
        // x, so it reaches the plane x >= 0 when its center is at -sqrt(2).
        const h = halfspace(1, 0, 0, 0);
        const axes = rotationAxes(0, 0, Math.PI / 4);
        expect(ti.test(h, box([-Math.SQRT2 + 1e-12, 0, 0], axes, [1, 1, 1])).intersect)
            .toBe(true);
        expect(ti.test(h, box([-Math.SQRT2 - 1e-6, 0, 0], axes, [1, 1, 1])).intersect)
            .toBe(false);
        // Without the rotation, contact happens at -1 instead.
        expect(ti.test(h, box([-1.2, 0, 0], identity, [1, 1, 1])).intersect)
            .toBe(false);
    });

    it('handles a degenerate box (a point)', () => {
        const h = halfspace(1, 1, 1, 0);
        const axes = rotationAxes(0.2, 0.3, 0.4);
        expect(ti.test(h, box([0, 0, 0], axes, [0, 0, 0])).intersect).toBe(true);
        expect(ti.test(h, box([-1, -1, -1], axes, [0, 0, 0])).intersect).toBe(false);
    });

    it('agrees with a vertex-enumeration oracle on random configurations', () => {
        const rand = makeRandom(112233);
        let numIn = 0, numOut = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const h = halfspace(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1,
                4 * rand() - 2);
            const b = box([4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2],
                rotationAxes(2 * Math.PI * rand(), 2 * Math.PI * rand(),
                    2 * Math.PI * rand()),
                [0.1 + rand(), 0.1 + rand(), 0.1 + rand()]);

            // Oracle: the box intersects the halfspace iff some vertex has a
            // nonnegative signed distance (the box is a convex hull of its
            // vertices).
            let oracle = false;
            for (const v of b.getVertices()) {
                if (dot(h.normal, v) - h.constant >= 0) {
                    oracle = true;
                }
            }

            const intersect = ti.test(h, b).intersect;
            expect(intersect).toBe(oracle);
            if (intersect) {
                ++numIn;
            } else {
                ++numOut;
            }
        }
        expect(numIn).toBeGreaterThan(50);
        expect(numOut).toBeGreaterThan(50);
    });
});

describe('IntrHalfspace3OrientedBox3 verification', () => {
    const ti = new IntrHalfspace3OrientedBox3TI();

    const hsArb = fc.tuple(unitVector(3), wellScaled(-6, 6))
        .map(([n, c]) => Halfspace.fromNormalConstant(n, c));
    const boxArb = fc.tuple(
        fc.array(wellScaled(-6, 6), { minLength: 3, maxLength: 3 }),
        rotationFrame(3),
        fc.array(positive(3), { minLength: 3, maxLength: 3 }))
        .map(([c, axes, ext]) => box(c, axes, ext));

    it('intersect equals "some box vertex is in the halfspace"', () => {
        check(fc.tuple(hsArb, boxArb), ([h, b]) => {
            // A convex body meets a halfspace iff its support point in the
            // normal direction does; for a box that support point is one of
            // the eight vertices. This uses the box vertices instead of the
            // query's projection-radius algebra.
            let best = -Infinity;
            for (const v of b.getVertices()) {
                best = Math.max(best, dot(h.normal, v) - h.constant);
            }
            const scale = 1 + Math.abs(h.constant) + length(b.center)
                + b.extent.values[0] + b.extent.values[1] + b.extent.values[2];
            if (Math.abs(best) < 1e-12 * scale) {
                return;    // the box exactly touches the plane
            }
            expect(ti.test(h, b).intersect).toBe(best >= 0);
        });
    });

    it('a sampled box point inside the halfspace forces intersect = true', () => {
        const rnd = makeRandom(0x2c19bd);
        check(fc.tuple(hsArb, boxArb), ([h, b]) => {
            if (ti.test(h, b).intersect) {
                return;
            }
            const e = b.extent.values;
            for (let k = 0; k < 200; ++k) {
                const p = [b.center.values[0], b.center.values[1],
                    b.center.values[2]];
                for (let d = 0; d < 3; ++d) {
                    const s = (2 * rnd() - 1) * e[d];
                    p[0] += s * b.axis[d].values[0];
                    p[1] += s * b.axis[d].values[1];
                    p[2] += s * b.axis[d].values[2];
                }
                const f = h.normal.values[0] * p[0] + h.normal.values[1] * p[1]
                    + h.normal.values[2] * p[2] - h.constant;
                expect(f).toBeLessThan(1e-9);
            }
        }, 60);
    }, 30000);

    it('is equivariant under a rigid motion of the halfspace and box', () => {
        check(fc.tuple(hsArb, boxArb, rotationFrame(3),
            fc.array(wellScaled(-4, 4), { minLength: 3, maxLength: 3 })),
            ([h, b, R, t]) => {
                const apply = (v: Vector): number[] => {
                    const out = [0, 0, 0];
                    for (let d = 0; d < 3; ++d) {
                        out[d] = R[0].values[d] * v.values[0]
                            + R[1].values[d] * v.values[1]
                            + R[2].values[d] * v.values[2];
                    }
                    return out;
                };
                const n2 = apply(h.normal);
                const c2 = apply(b.center);
                const h2 = Halfspace.fromNormalConstant(Vector.fromArray(n2),
                    h.constant + n2[0] * t[0] + n2[1] * t[1] + n2[2] * t[2]);
                const b2 = box([c2[0] + t[0], c2[1] + t[1], c2[2] + t[2]],
                    [Vector.fromArray(apply(b.axis[0])),
                        Vector.fromArray(apply(b.axis[1])),
                        Vector.fromArray(apply(b.axis[2]))],
                    [...b.extent.values]);
                let best = -Infinity;
                for (const v of b.getVertices()) {
                    best = Math.max(best, dot(h.normal, v) - h.constant);
                }
                const scale = 1 + Math.abs(h.constant) + length(b.center)
                    + b.extent.values[0] + b.extent.values[1]
                    + b.extent.values[2];
                if (Math.abs(best) < 1e-9 * scale) {
                    return;
                }
                expect(ti.test(h2, b2).intersect).toBe(ti.test(h, b).intersect);
            });
    });

    it('a zero-extent box is the point test on its center', () => {
        const identity = [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
            Vector.fromArray([0, 0, 1])];
        const h = halfspace(0, 0, 1, 2);   // z >= 2
        expect(ti.test(h, box([0, 0, 3], identity, [0, 0, 0])).intersect)
            .toBe(true);
        expect(ti.test(h, box([0, 0, 2], identity, [0, 0, 0])).intersect)
            .toBe(true);   // on the plane counts (closed halfspace)
        expect(ti.test(h, box([0, 0, 1], identity, [0, 0, 0])).intersect)
            .toBe(false);
    });
});
