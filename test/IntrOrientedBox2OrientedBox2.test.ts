import { describe, it, expect } from 'vitest';
import {
    IntrOrientedBox2OrientedBox2TI,
    IntrOrientedBox2OrientedBox2FI
} from '../src/IntrOrientedBox2OrientedBox2.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

const ti = new IntrOrientedBox2OrientedBox2TI();
const fi = new IntrOrientedBox2OrientedBox2FI();

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function rotatedAxes(angle: number): Vector[] {
    const c = Math.cos(angle), s = Math.sin(angle);
    return [v2(c, s), v2(-s, c)];
}

function makeBox(center: Vector, angle: number, extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(center, rotatedAxes(angle),
        Vector.fromArray(extent));
}

// The signed area of a polygon listed counterclockwise (shoelace formula).
function signedArea(polygon: Vector[]): number {
    let area = 0;
    const n = polygon.length;
    for (let i = 0; i < n; ++i) {
        const p = polygon[i], q = polygon[(i + 1) % n];
        area += p.values[0] * q.values[1] - q.values[0] * p.values[1];
    }
    return 0.5 * area;
}

// Test whether a point lies in a box, with a tolerance.
function inBox(p: Vector, box: OrientedBox, tolerance: number): boolean {
    const delta = sub(p, box.center);
    for (let k = 0; k < 2; ++k) {
        if (Math.abs(dot(delta, box.axis[k]))
            > box.extent.values[k] + tolerance) {
            return false;
        }
    }
    return true;
}

describe('IntrOrientedBox2OrientedBox2TI', () => {
    const unit = makeBox(v2(0, 0), 0, [1, 1]);

    it('reports coincident boxes as intersecting', () => {
        expect(ti.test(unit, unit).intersect).toBe(true);
    });

    it('reports partially overlapping boxes as intersecting', () => {
        const other = makeBox(v2(1, 0.5), 0, [1, 1]);
        expect(ti.test(unit, other).intersect).toBe(true);
    });

    it('treats edge-touching boxes as intersecting (solid, closed boxes)', () => {
        const touching = makeBox(v2(2, 0), 0, [1, 1]);
        expect(ti.test(unit, touching).intersect).toBe(true);
    });

    it('reports each of the four separating axes', () => {
        // box0.axis[0]
        let result = ti.test(unit, makeBox(v2(2.5, 0), 0, [1, 1]));
        expect(result.intersect).toBe(false);
        expect(result.separating).toBe(0);

        // box0.axis[1]
        result = ti.test(unit, makeBox(v2(0, 2.5), 0, [1, 1]));
        expect(result.intersect).toBe(false);
        expect(result.separating).toBe(1);

        // box1.axis[0]: a diamond (45 degrees) pushed away along its own
        // first axis. The box0 axes do not separate.
        const diamond = makeBox(v2(1.5, 1.5), Math.PI / 4, [0.5, 0.5]);
        result = ti.test(unit, diamond);
        expect(result.intersect).toBe(false);
        expect(result.separating).toBe(2);

        // box1.axis[1]: the same diamond, pushed along its second axis.
        const diamond2 = makeBox(v2(-1.5, 1.5), Math.PI / 4, [0.5, 0.5]);
        result = ti.test(unit, diamond2);
        expect(result.intersect).toBe(false);
        expect(result.separating).toBe(3);
    });

    it('handles zero-extent boxes (points)', () => {
        expect(ti.test(unit, makeBox(v2(0.5, 0.5), 0, [0, 0])).intersect)
            .toBe(true);
        expect(ti.test(unit, makeBox(v2(1, 1), 0, [0, 0])).intersect)
            .toBe(true);
        expect(ti.test(unit, makeBox(v2(1.5, 0), 0, [0, 0])).intersect)
            .toBe(false);
    });

    it('is symmetric in its arguments', () => {
        const rand = makeRandom(31415);
        for (let trial = 0; trial < 300; ++trial) {
            const box1 = makeBox(v2(4 * rand() - 2, 4 * rand() - 2),
                6 * rand(), [0.2 + rand(), 0.2 + rand()]);
            expect(ti.test(unit, box1).intersect)
                .toBe(ti.test(box1, unit).intersect);
        }
    });

    it('agrees with a point-sampling oracle', () => {
        const rand = makeRandom(271828);
        let found = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const e0 = [0.2 + rand(), 0.2 + rand()];
            const box0 = makeBox(v2(2 * rand() - 1, 2 * rand() - 1),
                6 * rand(), e0);
            const box1 = makeBox(v2(2 * rand() - 1, 2 * rand() - 1),
                6 * rand(), [0.2 + rand(), 0.2 + rand()]);

            let common = false;
            for (let s = 0; s < 200 && !common; ++s) {
                let p = box0.center.clone();
                for (let k = 0; k < 2; ++k) {
                    p = add(p, mul(e0[k] * (2 * rand() - 1), box0.axis[k]));
                }
                common = inBox(p, box1, 0);
            }

            if (common) {
                ++found;
                expect(ti.test(box0, box1).intersect).toBe(true);
            }
        }
        expect(found).toBeGreaterThan(20);
    });
});

describe('IntrOrientedBox2OrientedBox2FI', () => {
    const unit = makeBox(v2(0, 0), 0, [1, 1]);

    it('returns box0 when box0 is contained in box1', () => {
        const outer = makeBox(v2(0, 0), 0, [3, 3]);
        const result = fi.find(unit, outer);
        expect(result.intersect).toBe(true);
        expect(result.polygon).toHaveLength(4);
        expect(Math.abs(signedArea(result.polygon))).toBeCloseTo(4, 12);
        // The polygon is box0's vertices in counterclockwise order.
        expect(result.polygon[0].values).toEqual([-1, -1]);
        expect(result.polygon[1].values).toEqual([1, -1]);
        expect(result.polygon[2].values).toEqual([1, 1]);
        expect(result.polygon[3].values).toEqual([-1, 1]);
    });

    it('returns box0 for coincident boxes', () => {
        const result = fi.find(unit, unit);
        expect(result.intersect).toBe(true);
        expect(Math.abs(signedArea(result.polygon))).toBeCloseTo(4, 12);
    });

    it('clips box0 against an overlapping box1', () => {
        // The overlap of [-1,1]^2 and [0,2]x[-1,1] is [0,1]x[-1,1], area 2.
        const other = makeBox(v2(1, 0), 0, [1, 1]);
        const result = fi.find(unit, other);
        expect(result.intersect).toBe(true);
        expect(Math.abs(signedArea(result.polygon))).toBeCloseTo(2, 12);
        for (const p of result.polygon) {
            expect(inBox(p, unit, 1e-12)).toBe(true);
            expect(inBox(p, other, 1e-12)).toBe(true);
        }
    });

    it('clips a rotated box to an octagon-like polygon', () => {
        const diamond = makeBox(v2(0, 0), Math.PI / 4, [1.2, 1.2]);
        const result = fi.find(unit, diamond);
        expect(result.intersect).toBe(true);
        expect(result.polygon.length).toBeGreaterThan(4);
        const area = Math.abs(signedArea(result.polygon));
        // The clipped area is less than box0's area of 4 and more than half.
        expect(area).toBeLessThan(4);
        expect(area).toBeGreaterThan(2);
        for (const p of result.polygon) {
            expect(inBox(p, unit, 1e-9)).toBe(true);
            expect(inBox(p, diamond, 1e-9)).toBe(true);
        }
    });

    it('reports no intersection for separated boxes', () => {
        const apart = makeBox(v2(5, 0), 0, [1, 1]);
        const result = fi.find(unit, apart);
        expect(result.intersect).toBe(false);
        expect(result.polygon).toEqual([]);
    });

    it('reports no intersection for edge-touching boxes, unlike the TI query', () => {
        // The clipping query treats a degenerate (zero-area) overlap as no
        // intersection, because no box0 vertex is strictly inside box1. The
        // TI query, which uses closed-box separation tests, reports true.
        const touching = makeBox(v2(2, 0), 0, [1, 1]);
        expect(ti.test(unit, touching).intersect).toBe(true);
        expect(fi.find(unit, touching).intersect).toBe(false);
    });
});

describe('IntrOrientedBox2OrientedBox2 consistency', () => {
    it('produces a polygon inside both boxes and implies TI intersection', () => {
        const rand = makeRandom(161803);
        let polygons = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const box0 = makeBox(v2(2 * rand() - 1, 2 * rand() - 1),
                6 * rand(), [0.2 + rand(), 0.2 + rand()]);
            const box1 = makeBox(v2(2 * rand() - 1, 2 * rand() - 1),
                6 * rand(), [0.2 + rand(), 0.2 + rand()]);

            const tiResult = ti.test(box0, box1);
            const fiResult = fi.find(box0, box1);

            if (fiResult.intersect) {
                ++polygons;
                // A nonempty clip implies the boxes intersect.
                expect(tiResult.intersect).toBe(true);
                expect(fiResult.polygon.length).toBeGreaterThanOrEqual(3);
                for (const p of fiResult.polygon) {
                    expect(inBox(p, box0, 1e-9)).toBe(true);
                    expect(inBox(p, box1, 1e-9)).toBe(true);
                }
                // The clipped polygon is convex and counterclockwise.
                expect(signedArea(fiResult.polygon)).toBeGreaterThan(-1e-12);
            }
            else if (tiResult.intersect) {
                // The queries differ only for degenerate (zero-area) overlap.
                const delta = sub(box1.center, box0.center);
                let touching = false;
                for (let k = 0; k < 2; ++k) {
                    const gap0 = Math.abs(dot(delta, box0.axis[k]))
                        - (box0.extent.values[k]
                            + box1.extent.values[0] * Math.abs(dot(box0.axis[k], box1.axis[0]))
                            + box1.extent.values[1] * Math.abs(dot(box0.axis[k], box1.axis[1])));
                    const gap1 = Math.abs(dot(delta, box1.axis[k]))
                        - (box1.extent.values[k]
                            + box0.extent.values[0] * Math.abs(dot(box1.axis[k], box0.axis[0]))
                            + box0.extent.values[1] * Math.abs(dot(box1.axis[k], box0.axis[1])));
                    if (Math.abs(gap0) < 1e-9 || Math.abs(gap1) < 1e-9) {
                        touching = true;
                    }
                }
                expect(touching).toBe(true);
            }
        }
        expect(polygons).toBeGreaterThan(20);
    });
});
