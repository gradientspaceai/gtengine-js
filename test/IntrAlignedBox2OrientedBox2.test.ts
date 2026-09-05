import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector } from '../src/Vector.js';
import { IntrAlignedBox2OrientedBox2TI } from '../src/IntrAlignedBox2OrientedBox2.js';
import { IntrOrientedBox2OrientedBox2TI } from '../src/IntrOrientedBox2OrientedBox2.js';
import { check, fc, positive, seededRandom, wellScaled } from './helpers/arbitraries.js';

function alignedBox(minX: number, minY: number, maxX: number, maxY: number): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray([minX, minY]),
        Vector.fromArray([maxX, maxY]));
}

function orientedBox(cx: number, cy: number, angle: number,
    e0: number, e1: number): OrientedBox {
    const c = Math.cos(angle), s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(
        Vector.fromArray([cx, cy]),
        [Vector.fromArray([c, s]), Vector.fromArray([-s, c])],
        Vector.fromArray([e0, e1]));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// An independent oracle: clip the aligned box's rectangle against the four
// half planes of the oriented box (Sutherland-Hodgman) and return the area of
// the clipped convex polygon.
function overlapArea(box0: AlignedBox, box1: OrientedBox): number {
    let poly: number[][] = [
        [box0.min.values[0], box0.min.values[1]],
        [box0.max.values[0], box0.min.values[1]],
        [box0.max.values[0], box0.max.values[1]],
        [box0.min.values[0], box0.max.values[1]]
    ];
    const c = box1.center.values;
    for (let i = 0; i < 2; ++i) {
        const u = box1.axis[i].values;
        const e = box1.extent.values[i];
        for (const sign of [1, -1]) {
            // Keep the points with sign * dot(u, p - c) <= e.
            const f = (p: number[]) =>
                e - sign * (u[0] * (p[0] - c[0]) + u[1] * (p[1] - c[1]));
            const next: number[][] = [];
            for (let k = 0; k < poly.length; ++k) {
                const p = poly[k], q = poly[(k + 1) % poly.length];
                const fp = f(p), fq = f(q);
                if (fp >= 0) {
                    next.push(p);
                }
                if ((fp > 0 && fq < 0) || (fp < 0 && fq > 0)) {
                    const t = fp / (fp - fq);
                    next.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
                }
            }
            poly = next;
            if (poly.length === 0) {
                return 0;
            }
        }
    }
    let area = 0;
    for (let k = 0; k < poly.length; ++k) {
        const p = poly[k], q = poly[(k + 1) % poly.length];
        area += p[0] * q[1] - q[0] * p[1];
    }
    return Math.abs(0.5 * area);
}

describe('IntrAlignedBox2OrientedBox2', () => {
    const ti = new IntrAlignedBox2OrientedBox2TI();

    it('matches the aligned-aligned answer when the oriented box is axis aligned', () => {
        const b0 = alignedBox(0, 0, 2, 2);
        // Center (3,1), extents (1,1): touches box0 along x = 2.
        expect(ti.test(b0, orientedBox(3, 1, 0, 1, 1)).intersect).toBe(true);
        // Center (3.5,1): separated along box0.axis[0].
        const r = ti.test(b0, orientedBox(3.5, 1, 0, 1, 1));
        expect(r.intersect).toBe(false);
        expect(r.separating).toBe(0);
    });

    it('reports separation along box0.axis[1]', () => {
        const r = ti.test(alignedBox(0, 0, 2, 2), orientedBox(1, 5, 0, 1, 1));
        expect(r.intersect).toBe(false);
        expect(r.separating).toBe(1);
    });

    it('reports the first separating axis found, in upstream order', () => {
        // A 45-degree square far along the diagonal is separated by several
        // axes; upstream tests box0.axis[0] first, so 0 is reported.
        const b0 = alignedBox(-1, -1, 1, 1);
        const r = ti.test(b0, orientedBox(3, 3, Math.PI / 4, 1, 1));
        expect(r.intersect).toBe(false);
        expect(r.separating).toBe(0);
    });

    it('detects the diamond/square touching configuration exactly', () => {
        // Unit square [-1,1]^2 and a 45-degree square of extent 1 centered on
        // the x axis. Its leftmost vertex is at cx - sqrt(2). Contact occurs
        // when cx - sqrt(2) = 1.
        const b0 = alignedBox(-1, -1, 1, 1);
        const touch = 1 + Math.SQRT2;
        expect(ti.test(b0, orientedBox(touch, 0, Math.PI / 4, 1, 1)).intersect)
            .toBe(true);
        expect(ti.test(b0, orientedBox(touch + 1e-9, 0, Math.PI / 4, 1, 1)).intersect)
            .toBe(false);
        expect(ti.test(b0, orientedBox(touch - 1e-9, 0, Math.PI / 4, 1, 1)).intersect)
            .toBe(true);
    });

    it('needs an oriented-box axis to separate a diagonal configuration', () => {
        // A long thin 45-degree box placed off the corner of the square: the
        // projections onto both aligned axes overlap, so only a box1 axis can
        // report separation.
        const b0 = alignedBox(-1, -1, 1, 1);
        const box1 = orientedBox(0, 0, Math.PI / 4, 0.1, 10);
        expect(ti.test(b0, box1).intersect).toBe(true);
        const shifted = orientedBox(2.5, 2.5, Math.PI / 4, 0.1, 10);
        const r = ti.test(b0, shifted);
        expect(r.intersect).toBe(false);
        expect(r.separating).toBe(2);
    });

    it('handles degenerate (zero-extent) oriented boxes', () => {
        const b0 = alignedBox(0, 0, 2, 2);
        // A point at the box corner.
        expect(ti.test(b0, orientedBox(2, 2, 0.7, 0, 0)).intersect).toBe(true);
        // A point outside.
        expect(ti.test(b0, orientedBox(2.5, 2, 0.7, 0, 0)).intersect).toBe(false);
        // A segment (one zero extent) crossing the box.
        expect(ti.test(b0, orientedBox(1, 1, 0.3, 5, 0)).intersect).toBe(true);
    });

    it('agrees with a polygon-clipping overlap oracle on random configurations', () => {
        const rand = makeRandom(31337);
        let numIntersect = 0, numSeparate = 0, numSkipped = 0;
        for (let trial = 0; trial < 2000; ++trial) {
            const x0 = 4 * rand() - 2, y0 = 4 * rand() - 2;
            const b0 = alignedBox(x0, y0, x0 + 0.2 + 2 * rand(), y0 + 0.2 + 2 * rand());
            const b1 = orientedBox(4 * rand() - 2, 4 * rand() - 2,
                2 * Math.PI * rand(), 0.1 + rand(), 0.1 + rand());

            const area = overlapArea(b0, b1);
            const intersect = ti.test(b0, b1).intersect;
            if (area > 0 && area < 1e-9) {
                // Nearly tangential: the oracle cannot decide reliably.
                ++numSkipped;
                continue;
            }
            expect(intersect).toBe(area > 0);
            if (intersect) {
                ++numIntersect;
            } else {
                ++numSeparate;
            }
        }
        expect(numIntersect).toBeGreaterThan(100);
        expect(numSeparate).toBeGreaterThan(100);
        expect(numSkipped).toBeLessThan(20);
    });
});

describe('IntrAlignedBox2OrientedBox2 verification', () => {
    const query = new IntrAlignedBox2OrientedBox2TI();
    const obbQuery = new IntrOrientedBox2OrientedBox2TI();

    const boxArb = fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4),
        positive(4), positive(4))
        .map(([x, y, w, h]) => alignedBox(x, y, x + w, y + h));

    const obbArb = fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4),
        wellScaled(-Math.PI, Math.PI), positive(3), positive(3))
        .map(([cx, cy, a, e0, e1]) => orientedBox(cx, cy, a, e0, e1));

    // The aligned box in centered form, as an oriented box with the standard
    // axes. Every product the two queries form is then bit-identical, so the
    // cross-check below is an exact equality.
    function asOriented(box: AlignedBox): OrientedBox {
        const { center, extent } = box.getCenteredForm();
        return OrientedBox.fromCenterAxisExtent(center,
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])], extent);
    }

    it('agrees with the oriented-oriented query on the same boxes', () => {
        check(fc.tuple(boxArb, obbArb), ([box0, box1]) => {
            const a = query.test(box0, box1);
            const b = obbQuery.test(asOriented(box0), box1);
            expect(a.intersect).toBe(b.intersect);
            if (!a.intersect) {
                expect(a.separating).toBe(b.separating);
            }
        });
    });

    it('the reported separating axis really separates the two boxes', () => {
        check(fc.tuple(boxArb, obbArb), ([box0, box1]) => {
            const r = query.test(box0, box1);
            if (r.intersect) {
                return;
            }
            const axis = r.separating < 2
                ? Vector.fromArray(r.separating === 0 ? [1, 0] : [0, 1])
                : box1.axis[r.separating - 2];
            const project = (vs: Vector[]): [number, number] => {
                let lo = Infinity, hi = -Infinity;
                for (const v of vs) {
                    const d = v.values[0] * axis.values[0]
                        + v.values[1] * axis.values[1];
                    lo = Math.min(lo, d);
                    hi = Math.max(hi, d);
                }
                return [lo, hi];
            };
            const p0 = project(box0.getVertices());
            const p1 = project(box1.getVertices());
            // The query uses a strict '>' on the radius sum, so a reported
            // separation must be a strict gap. The projections are formed
            // from the same inputs by a different route, so a relative
            // tolerance of 1e-12 absorbs the reassociation only.
            const tol = 1e-12 * (1 + Math.abs(p0[0]) + Math.abs(p0[1])
                + Math.abs(p1[0]) + Math.abs(p1[1]));
            expect(p0[1] < p1[0] + tol || p1[1] < p0[0] + tol).toBe(true);
        });
    });

    it('a point common to both boxes forces intersect = true', () => {
        const rnd = seededRandom(0x5ab21c9);
        check(fc.tuple(boxArb, obbArb), ([box0, box1]) => {
            const r = query.test(box0, box1);
            if (r.intersect) {
                return;
            }
            // Sample the oriented box and test containment in the aligned box.
            const c = box1.center.values;
            const u0 = box1.axis[0].values, u1 = box1.axis[1].values;
            const e = box1.extent.values;
            for (let k = 0; k < 200; ++k) {
                const s0 = (2 * rnd() - 1) * e[0];
                const s1 = (2 * rnd() - 1) * e[1];
                const px = c[0] + s0 * u0[0] + s1 * u1[0];
                const py = c[1] + s0 * u0[1] + s1 * u1[1];
                const inside = px >= box0.min.values[0] && px <= box0.max.values[0]
                    && py >= box0.min.values[1] && py <= box0.max.values[1];
                expect(inside).toBe(false);
            }
        }, 60);
    }, 30000);

    it('an axis-aligned oriented box reduces to the aligned-aligned query', () => {
        const b0 = alignedBox(0, 0, 2, 2);
        // Zero rotation: the oriented box has the standard axes.
        expect(query.test(b0, orientedBox(3, 1, 0, 1, 1)).intersect).toBe(true);
        expect(query.test(b0, orientedBox(3.5, 1, 0, 1, 1)).intersect).toBe(false);
        // Touching along a face is an intersection (closed solids).
        const touch = query.test(b0, orientedBox(3, 1, 0, 1, 1));
        expect(touch.intersect).toBe(true);
    });
});
