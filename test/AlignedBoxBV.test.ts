import { describe, it, expect } from 'vitest';
import { AlignedBoxBV } from '../src/AlignedBoxBV.js';
import { AlignedBox } from '../src/AlignedBox.js';
import { Vector, add, length, mul, normalize, sub } from '../src/Vector.js';
import {
    alignedBox, check, expectVectorClose, fc, finite, unitVector,
    wellScaledVector
} from './helpers/arbitraries.js';

function V(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function unit(x: number, y: number, z: number): Vector {
    const v = V(x, y, z);
    normalize(v);
    return v;
}

function makeBV(min: [number, number, number],
    max: [number, number, number]): AlignedBoxBV {
    return AlignedBoxBV.fromBox(AlignedBox.fromMinMax(V(...min), V(...max)));
}

describe('AlignedBoxBV', () => {
    it('default-constructs the [-1,1]^3 box of AlignedBox3', () => {
        // Upstream is 'AlignedBoxBV() : box{}'. AlignedBox has a
        // user-provided default constructor, so value-initialization runs it
        // (there is no zero-initialization first) and the default box is
        // min = (-1,-1,-1), max = (1,1,1). The port previously produced the
        // degenerate box at the origin.
        const bv = new AlignedBoxBV();
        expect(bv.box.min.values).toEqual([-1, -1, -1]);
        expect(bv.box.max.values).toEqual([1, 1, 1]);
        // The default AlignedBoxBV agrees with the default AlignedBox3 and
        // with the OrientedBoxBV default, which encloses [-1,1]^3 too.
        const box = new AlignedBox(3);
        expect(bv.box.min.values).toEqual(box.min.values);
        expect(bv.box.max.values).toEqual(box.max.values);
    });

    it('copies the box in fromBox', () => {
        const box = AlignedBox.fromMinMax(V(-1, -2, -3), V(4, 5, 6));
        const bv = AlignedBoxBV.fromBox(box);
        box.min.values[0] = 100;
        expect(bv.box.min.values).toEqual([-1, -2, -3]);
        expect(bv.box.max.values).toEqual([4, 5, 6]);
    });

    it('splits along the axis of largest extent, breaking ties toward x', () => {
        // Extents (3,1,2): the largest is along x.
        let s = makeBV([-3, -1, -2], [3, 1, 2]).getSplittingAxis();
        expect(s.origin.values).toEqual([0, 0, 0]);
        expect(s.direction.values).toEqual([1, 0, 0]);

        // Extents (1,3,2): the largest is along y. The center is not zero.
        s = makeBV([0, 0, 0], [2, 6, 4]).getSplittingAxis();
        expect(s.origin.values).toEqual([1, 3, 2]);
        expect(s.direction.values).toEqual([0, 1, 0]);

        // Extents (1,2,3): the largest is along z.
        s = makeBV([-1, -2, -3], [1, 2, 3]).getSplittingAxis();
        expect(s.direction.values).toEqual([0, 0, 1]);

        // A cube: the strict comparisons keep the first maximum, x.
        s = makeBV([-1, -1, -1], [1, 1, 1]).getSplittingAxis();
        expect(s.direction.values).toEqual([1, 0, 0]);

        // Ties between y and z: y is chosen because z is not strictly larger.
        s = makeBV([-1, -2, -2], [1, 2, 2]).getSplittingAxis();
        expect(s.direction.values).toEqual([0, 1, 0]);
    });

    it('tests line intersection', () => {
        const bv = makeBV([-1, -1, -1], [1, 1, 1]);

        // A line through the box.
        expect(AlignedBoxBV.intersectLine(V(0, 0, -5), unit(0, 0, 1), bv))
            .toBe(true);
        // The same line's opposite direction: a line is unbounded, so it
        // still intersects.
        expect(AlignedBoxBV.intersectLine(V(0, 0, -5), unit(0, 0, -1), bv))
            .toBe(true);
        // A line missing the box.
        expect(AlignedBoxBV.intersectLine(V(5, 5, 0), unit(0, 0, 1), bv))
            .toBe(false);
        // A line grazing the corner (1,1,1).
        expect(AlignedBoxBV.intersectLine(V(1, 1, 1), unit(1, -1, 0), bv))
            .toBe(true);
    });

    it('tests ray intersection', () => {
        const bv = makeBV([-1, -1, -1], [1, 1, 1]);

        // The ray points at the box.
        expect(AlignedBoxBV.intersectRay(V(0, 0, -5), unit(0, 0, 1), bv))
            .toBe(true);
        // The ray points away from the box.
        expect(AlignedBoxBV.intersectRay(V(0, 0, -5), unit(0, 0, -1), bv))
            .toBe(false);
        // The ray origin is inside the box.
        expect(AlignedBoxBV.intersectRay(V(0, 0, 0), unit(1, 2, 3), bv))
            .toBe(true);
        // The ray misses the box.
        expect(AlignedBoxBV.intersectRay(V(5, 5, -5), unit(0, 0, 1), bv))
            .toBe(false);
    });

    it('tests segment intersection', () => {
        const bv = makeBV([-1, -1, -1], [1, 1, 1]);

        // The segment crosses the box.
        expect(AlignedBoxBV.intersectSegment(V(0, 0, -5), V(0, 0, 5), bv))
            .toBe(true);
        // The segment stops short of the box.
        expect(AlignedBoxBV.intersectSegment(V(0, 0, -5), V(0, 0, -2), bv))
            .toBe(false);
        // The segment is entirely inside the box.
        expect(AlignedBoxBV.intersectSegment(V(-0.5, 0, 0), V(0.5, 0, 0), bv))
            .toBe(true);
        // The segment touches the face x = 1 at one endpoint.
        expect(AlignedBoxBV.intersectSegment(V(1, 0, 0), V(3, 0, 0), bv))
            .toBe(true);
        // A segment that misses.
        expect(AlignedBoxBV.intersectSegment(V(2, 2, 2), V(3, 3, 3), bv))
            .toBe(false);
    });

    it('agrees with brute-force sampling on random segments', () => {
        let seed = 20240613;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const R = (): number => 6 * rand() - 3;

        const bv = makeBV([-1, -0.5, -2], [1.5, 2, 0.5]);
        const box = bv.box;
        const inBox = (p: Vector): boolean => {
            for (let k = 0; k < 3; ++k) {
                if (p.values[k] < box.min.values[k] ||
                    p.values[k] > box.max.values[k]) {
                    return false;
                }
            }
            return true;
        };

        let numHits = 0;
        let numMisses = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const p0 = V(R(), R(), R());
            const p1 = V(R(), R(), R());
            const intersect = AlignedBoxBV.intersectSegment(p0, p1, bv);

            // Sample the segment densely. Any sampled point inside the box
            // proves an intersection; the converse is not guaranteed, so the
            // sampling can only confirm hits.
            let sampledHit = false;
            const delta = sub(p1, p0);
            for (let i = 0; i <= 400; ++i) {
                if (inBox(add(p0, mul(delta, i / 400)))) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(intersect).toBe(true);
                ++numHits;
            }
            else if (!intersect) {
                ++numMisses;
            }

            // The ray from p0 through p1 must intersect whenever the segment
            // does, and the line must intersect whenever the ray does.
            const direction = sub(p1, p0);
            normalize(direction);
            const rayHit = AlignedBoxBV.intersectRay(p0, direction, bv);
            const lineHit = AlignedBoxBV.intersectLine(p0, direction, bv);
            if (intersect) {
                expect(rayHit).toBe(true);
            }
            if (rayHit) {
                expect(lineHit).toBe(true);
            }
        }
        expect(numHits).toBeGreaterThan(5);
        expect(numMisses).toBeGreaterThan(5);
    });
});

// ---------------------------------------------------------------------------
// V43 verification: the splitting axis and the three linear-component
// predicates against an independent slab clipper.
// ---------------------------------------------------------------------------

// Clip the parameter interval [t0, t1] of P + t*D against the slabs of the
// box [lo, hi]. This is an independent implementation of the standard
// slab test, not a restatement of the ported IntrLine3AlignedBox3 code.
function slabHit(P: Vector, D: Vector, t0: number, t1: number,
    lo: Vector, hi: Vector): boolean {
    let tmin = t0;
    let tmax = t1;
    for (let k = 0; k < 3; ++k) {
        const d = D.get(k);
        if (d !== 0) {
            const a = (lo.get(k) - P.get(k)) / d;
            const b = (hi.get(k) - P.get(k)) / d;
            tmin = Math.max(tmin, Math.min(a, b));
            tmax = Math.min(tmax, Math.max(a, b));
            if (tmin > tmax) {
                return false;
            }
        } else if (P.get(k) < lo.get(k) || P.get(k) > hi.get(k)) {
            return false;
        }
    }
    return true;
}

// The slab answer is only trustworthy away from tangency, so ask it on a
// slightly shrunk and a slightly expanded box; when the two agree, the
// configuration is robust and the query must return that answer.
function robustSlab(box: AlignedBox, P: Vector, D: Vector, t0: number,
    t1: number): boolean | undefined {
    const eps = 1e-9;
    const shrunkLo = new Vector(3), shrunkHi = new Vector(3);
    const grownLo = new Vector(3), grownHi = new Vector(3);
    for (let k = 0; k < 3; ++k) {
        shrunkLo.set(k, box.min.get(k) + eps);
        shrunkHi.set(k, box.max.get(k) - eps);
        grownLo.set(k, box.min.get(k) - eps);
        grownHi.set(k, box.max.get(k) + eps);
    }
    const a = slabHit(P, D, t0, t1, shrunkLo, shrunkHi);
    const b = slabHit(P, D, t0, t1, grownLo, grownHi);
    return a === b ? a : undefined;
}

describe('AlignedBoxBV verification', () => {
    // Any box, including degenerate ones, for the arithmetic properties.
    const bvArb = alignedBox(3, -5, 5).map(box => AlignedBoxBV.fromBox(box));

    // A box with extents bounded away from zero, for the intersection
    // properties. The SAT-based segment/box test treats measure-zero contact
    // with a flat box as separation (see docs/API.md), and alignedBox() emits
    // subnormal coordinates, so with a degenerate box "does the linear
    // component meet the box" is ill-posed rather than wrong.
    const fatBvArb = fc.tuple(wellScaledVector(3, -5, 5),
        fc.array(finite(0.5, 4), { minLength: 3, maxLength: 3 }))
        .map(([center, e]) => {
            const extent = Vector.fromArray(e);
            return AlignedBoxBV.fromBox(AlignedBox.fromMinMax(
                sub(center, extent), add(center, extent)));
        });

    it('splits through the center along the axis of largest extent', () => {
        check(bvArb, bv => {
            const { origin, direction } = bv.getSplittingAxis();
            const { center, extent } = bv.box.getCenteredForm();
            expectVectorClose(origin, center, 0, 0);
            // Independent argmax with ties broken toward the smaller index.
            let maxIndex = 0;
            for (let k = 1; k < 3; ++k) {
                if (extent.get(k) > extent.get(maxIndex)) {
                    maxIndex = k;
                }
            }
            expectVectorClose(direction, Vector.unit(3, maxIndex), 0, 0);
            expect(direction.values.filter(x => x !== 0).length).toBe(1);
        });
    });

    it('returns fresh vectors from getSplittingAxis', () => {
        check(bvArb, bv => {
            const s = bv.getSplittingAxis();
            s.origin.set(0, s.origin.get(0) + 1000);
            s.direction.set(0, s.direction.get(0) + 1000);
            const t = bv.getSplittingAxis();
            expect(t.origin.get(0)).not.toBe(s.origin.get(0));
            expect(t.direction.get(0)).not.toBe(s.direction.get(0));
        });
    });

    it('agrees with the slab clipper for lines, rays and segments', () => {
        check(fc.tuple(fatBvArb, wellScaledVector(3, -8, 8), unitVector(3),
            finite(0.5, 12)), ([bv, P, D, len]) => {
            const Q = add(P, mul(D, len));
            const line = robustSlab(bv.box, P, D, -Infinity, Infinity);
            if (line !== undefined) {
                expect(AlignedBoxBV.intersectLine(P, D, bv)).toBe(line);
            }
            const ray = robustSlab(bv.box, P, D, 0, Infinity);
            if (ray !== undefined) {
                expect(AlignedBoxBV.intersectRay(P, D, bv)).toBe(ray);
            }
            // The segment predicate takes endpoints, so its direction is
            // Q - P and its parameter interval is [0, 1].
            const seg = robustSlab(bv.box, P, sub(Q, P), 0, 1);
            if (seg !== undefined) {
                expect(AlignedBoxBV.intersectSegment(P, Q, bv)).toBe(seg);
            }
        });
    });

    it('is monotone in segment, ray and line', () => {
        check(fc.tuple(fatBvArb, wellScaledVector(3, -8, 8), unitVector(3),
            finite(0.5, 12)), ([bv, P, D, len]) => {
            const Q = add(P, mul(D, len));
            const seg = AlignedBoxBV.intersectSegment(P, Q, bv);
            const ray = AlignedBoxBV.intersectRay(P, D, bv);
            const line = AlignedBoxBV.intersectLine(P, D, bv);
            // segment subset of ray subset of line.
            if (seg) { expect(ray).toBe(true); }
            if (ray) { expect(line).toBe(true); }
        });
    });

    it('hits when an endpoint is inside the box', () => {
        check(fc.tuple(fatBvArb, wellScaledVector(3, -0.9, 0.9),
            wellScaledVector(3, -8, 8)), ([bv, t, Q]) => {
            // A point strictly inside the box.
            const { center, extent } = bv.box.getCenteredForm();
            const P = add(center, Vector.fromArray([
                t.get(0) * extent.get(0), t.get(1) * extent.get(1),
                t.get(2) * extent.get(2)]));
            expect(AlignedBoxBV.intersectSegment(P, Q, bv)).toBe(true);
            const D = sub(Q, P);
            if (length(D) > 1e-6) {
                normalize(D);
                expect(AlignedBoxBV.intersectRay(P, D, bv)).toBe(true);
                expect(AlignedBoxBV.intersectLine(P, D, bv)).toBe(true);
            }
        });
    });
});
