import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox.js';
import { Vector, add, mul } from '../src/Vector.js';
import {
    IntrAlignedBox2Circle2TI,
    IntrAlignedBox2Circle2FI,
    IntrAlignedBox2Circle2FIResultType as Type
} from '../src/IntrAlignedBox2Circle2.js';
import { dot, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, positive, seededRandom,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function box(minx: number, miny: number, maxx: number, maxy: number): AlignedBox {
    return AlignedBox.fromMinMax(v2(minx, miny), v2(maxx, maxy));
}

function circle(cx: number, cy: number, r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v2(cx, cy), r);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The distance from the circle center to the box at time t, where the box
// moves with velocity bv and the circle center with velocity cv.
function sqrDistanceAtTime(b: AlignedBox, bv: Vector, c: Hypersphere,
    cv: Vector, t: number): number {
    const movedBox = AlignedBox.fromMinMax(add(b.min, mul(bv, t)),
        add(b.max, mul(bv, t)));
    const movedCenter = add(c.center, mul(cv, t));
    return new DistPointAlignedBox().compute(movedCenter, movedBox).sqrDistance;
}

const unitBox = box(-1, -1, 1, 1);

describe('IntrAlignedBox2Circle2TI', () => {
    const query = new IntrAlignedBox2Circle2TI();

    it('detects a circle inside the box', () => {
        expect(query.test(unitBox, circle(0, 0, 0.25)).intersect).toBe(true);
    });

    it('detects a box inside the circle', () => {
        expect(query.test(unitBox, circle(0, 0, 10)).intersect).toBe(true);
    });

    it('detects face overlap and separation', () => {
        expect(query.test(unitBox, circle(1.5, 0, 1)).intersect).toBe(true);
        expect(query.test(unitBox, circle(2.5, 0, 1)).intersect).toBe(false);
    });

    it('treats tangency as intersection', () => {
        expect(query.test(unitBox, circle(2, 0, 1)).intersect).toBe(true);
        // Corner tangency: distance from (1,1) to the center is exactly 1.
        const s = Math.SQRT1_2;
        expect(query.test(unitBox, circle(1 + s, 1 + s, 1)).intersect).toBe(true);
    });

    it('separates a circle near a box corner', () => {
        // The corner (1,1) is at distance sqrt(2) from (2,2), which exceeds 1.
        expect(query.test(unitBox, circle(2, 2, 1)).intersect).toBe(false);
        expect(query.test(unitBox, circle(2, 2, 1.5)).intersect).toBe(true);
    });
});

describe('IntrAlignedBox2Circle2FI', () => {
    const query = new IntrAlignedBox2Circle2FI();

    it('reports initial overlap for a circle inside the box', () => {
        const result = query.find(unitBox, v2(0, 0), circle(0.25, -0.5, 0.1),
            v2(0, 0));
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([0.25, -0.5]);
    });

    it('reports initial overlap for a circle crossing a face', () => {
        const result = query.find(unitBox, v2(0, 0), circle(1.5, 0.25, 1),
            v2(0, 0));
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([1, 0.25]);
    });

    it('reports initial overlap for a circle covering a vertex', () => {
        const result = query.find(unitBox, v2(0, 0), circle(1.5, 1.5, 1),
            v2(0, 0));
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactPoint.values).toEqual([1, 1]);
    });

    it('reports no contact for a separated, non-moving circle', () => {
        const result = query.find(unitBox, v2(0, 0), circle(5, 0, 1), v2(0, 0));
        expect(result.intersectionType).toBe(Type.noContact);
        expect(result.contactTime).toBe(0);
    });

    it('reports contact at t = 0 for tangency along a face', () => {
        const result = query.find(unitBox, v2(0, 0), circle(2, 0.5, 1),
            v2(-1, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([1, 0.5]);
    });

    it('reports contact at t = 0 for tangency at a vertex', () => {
        // The 3-4-5 triangle makes the tangency exact in binary floating
        // point: the corner (1,1) is at squared distance 3^2 + 4^2 = 25 from
        // the center (4,5), which equals radius^2 for radius 5. Upstream
        // distinguishes overlap (sqrDistance < sqrRadius) from contact
        // (equality), so this is the +1 case with contactTime 0.
        const result = query.find(unitBox, v2(0, 0), circle(4, 5, 5),
            v2(-1, -1));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([1, 1]);
    });

    it('reports overlap when a vertex is strictly inside the circle', () => {
        const result = query.find(unitBox, v2(0, 0), circle(4, 5, 5.5),
            v2(-1, -1));
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactPoint.values).toEqual([1, 1]);
    });

    it('computes the contact for approach along the +x face', () => {
        const result = query.find(unitBox, v2(0, 0), circle(5, 0.25, 1),
            v2(-1, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(0.25, 12);
    });

    it('computes the contact for approach along the -y face', () => {
        const result = query.find(unitBox, v2(0, 0), circle(-0.5, -4, 0.5),
            v2(0, 2));
        expect(result.intersectionType).toBe(Type.contact);
        // The circle boundary reaches y = -1 after moving 2.5 units, so
        // t = 2.5/2 = 1.25.
        expect(result.contactTime).toBeCloseTo(1.25, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(-0.5, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(-1, 12);
    });

    it('computes the contact for approach toward a vertex', () => {
        // The center starts 3 units (along the diagonal) from the corner
        // (1,1) and moves toward it at unit speed, so contact occurs at
        // t = 3 - radius = 2 at the corner.
        const s = Math.SQRT1_2;
        const result = query.find(unitBox, v2(0, 0),
            circle(1 + 3 * s, 1 + 3 * s, 1), v2(-s, -s));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(2, 10);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(1, 12);
    });

    it('computes the contact for a vertex approach that grazes a face', () => {
        // The circle is diagonally away from the corner but moves in -x
        // only, so it eventually strikes the +x face (or the corner region).
        const result = query.find(unitBox, v2(0, 0), circle(6, 0.5, 1),
            v2(-2, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(2, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(0.5, 12);
    });

    it('reports no contact for a receding circle', () => {
        const result = query.find(unitBox, v2(0, 0), circle(5, 0, 1), v2(1, 0));
        expect(result.intersectionType).toBe(Type.noContact);
    });

    it('reports no contact for a circle passing by the box', () => {
        // The circle travels along x = 5, never nearing the box.
        const result = query.find(unitBox, v2(0, 0), circle(5, -10, 1),
            v2(0, 1));
        expect(result.intersectionType).toBe(Type.noContact);
    });

    it('uses the relative velocity of the two objects', () => {
        // The box chases the circle: with both moving at the same velocity
        // the relative velocity is zero and there is no contact.
        const stationary = query.find(unitBox, v2(1, 1), circle(5, 0, 1),
            v2(1, 1));
        expect(stationary.intersectionType).toBe(Type.noContact);

        // The box moves toward a stationary circle: same answer as a circle
        // moving toward a stationary box with the opposite velocity.
        const boxMoves = query.find(unitBox, v2(1, 0), circle(5, 0.25, 1),
            v2(0, 0));
        const circleMoves = query.find(unitBox, v2(0, 0), circle(5, 0.25, 1),
            v2(-1, 0));
        expect(boxMoves.intersectionType).toBe(circleMoves.intersectionType);
        expect(boxMoves.contactTime).toBeCloseTo(circleMoves.contactTime, 12);
    });

    it('handles boxes that are not centered at the origin', () => {
        const b = box(3, 5, 5, 7);  // center (4,6), extent (1,1)
        const result = query.find(b, v2(0, 0), circle(9, 6.25, 1), v2(-1, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(5, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(6.25, 12);
    });

    it('is symmetric under reflection of the configuration', () => {
        const base = query.find(unitBox, v2(0, 0), circle(5, 2, 0.75),
            v2(-1, -0.5));
        const reflected = query.find(unitBox, v2(0, 0), circle(-5, -2, 0.75),
            v2(1, 0.5));
        expect(reflected.intersectionType).toBe(base.intersectionType);
        expect(reflected.contactTime).toBeCloseTo(base.contactTime, 12);
        expect(reflected.contactPoint.values[0])
            .toBeCloseTo(-base.contactPoint.values[0], 12);
        expect(reflected.contactPoint.values[1])
            .toBeCloseTo(-base.contactPoint.values[1], 12);
    });

    it('agrees with a sampled simulation (randomized)', () => {
        const rand = makeRandom(13572468);
        let numContact = 0, numOverlap = 0, numNone = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const ex = 0.25 + rand(), ey = 0.25 + rand();
            const bcx = 2 * rand() - 1, bcy = 2 * rand() - 1;
            const b = box(bcx - ex, bcy - ey, bcx + ex, bcy + ey);
            const radius = 0.1 + rand();
            const c = circle(8 * rand() - 4, 8 * rand() - 4, radius);
            const bv = v2(2 * rand() - 1, 2 * rand() - 1);
            const cv = v2(2 * rand() - 1, 2 * rand() - 1);
            const result = query.find(b, bv, c, cv);
            const rsqr = radius * radius;

            if (result.intersectionType === Type.initiallyOverlapping) {
                ++numOverlap;
                expect(sqrDistanceAtTime(b, bv, c, cv, 0))
                    .toBeLessThanOrEqual(rsqr + 1e-12);
                expect(result.contactTime).toBe(0);
            }
            else if (result.intersectionType === Type.contact) {
                ++numContact;
                const t = result.contactTime;
                expect(t).toBeGreaterThanOrEqual(0);
                // At the contact time the distance from the moved center to
                // the moved box equals the radius.
                expect(Math.sqrt(sqrDistanceAtTime(b, bv, c, cv, t)))
                    .toBeCloseTo(radius, 8);
                // No earlier contact: sample the interval [0,t).
                for (let k = 0; k < 20; ++k) {
                    const s = t * k / 20;
                    expect(Math.sqrt(sqrDistanceAtTime(b, bv, c, cv, s)))
                        .toBeGreaterThan(radius - 1e-9);
                }
            }
            else {
                ++numNone;
                // No contact at any sampled future time, out to a horizon
                // beyond which the objects only separate.
                for (let k = 0; k <= 200; ++k) {
                    const s = 20 * k / 200;
                    expect(Math.sqrt(sqrDistanceAtTime(b, bv, c, cv, s)))
                        .toBeGreaterThan(radius - 1e-9);
                }
            }
        }
        expect(numContact).toBeGreaterThan(10);
        expect(numOverlap).toBeGreaterThan(10);
        expect(numNone).toBeGreaterThan(10);
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrAlignedBox2Circle2.h.
// ---------------------------------------------------------------------------

describe('IntrAlignedBox2Circle2 verification', () => {
    const ti = new IntrAlignedBox2Circle2TI();
    const fi = new IntrAlignedBox2Circle2FI();

    const arbBox = fc.tuple(wellScaledVector(2, -4, 4),
        fc.array(positive(3, 0.2), { minLength: 2, maxLength: 2 }))
        .map(([c, e]) => AlignedBox.fromMinMax(
            Vector.fromArray([c.get(0) - e[0], c.get(1) - e[1]]),
            Vector.fromArray([c.get(0) + e[0], c.get(1) + e[1]])));
    const arbCircle = fc.tuple(wellScaledVector(2, -6, 6), positive(2, 0.2))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

    // The signed gap between the moving circle and the moving box: negative
    // when they overlap, zero on contact. Computed from the closed-form
    // point-box distance, independent of the query's Voronoi-region analysis.
    function gap(b: AlignedBox, bv: Vector, c: Hypersphere, cv: Vector,
        t: number): number {
        let sum = 0;
        for (let i = 0; i < 2; ++i) {
            const lo = b.min.get(i) + t * bv.get(i);
            const hi = b.max.get(i) + t * bv.get(i);
            const p = c.center.get(i) + t * cv.get(i);
            const d = Math.max(lo - p, 0, p - hi);
            sum += d * d;
        }
        return Math.sqrt(sum) - c.radius;
    }

    it('the TI query matches the closed-form point-box distance', () => {
        check(fc.tuple(arbBox, arbCircle), ([b, c]) => {
            const expected = gap(b, Vector.zero(2), c, Vector.zero(2), 0) <= 0;
            const actual = ti.test(b, c).intersect;
            // Only assert away from exact tangency, where the two formulas can
            // differ by an ulp (the query compares squared quantities).
            if (Math.abs(gap(b, Vector.zero(2), c, Vector.zero(2), 0))
                > 1e-9 * (1 + c.radius)) {
                expect(actual).toBe(expected);
            }
        });
    });

    // Configurations aimed at the box so that a useful fraction of the draws
    // produce contacts.
    const arbMoving = fc.tuple(arbBox, arbCircle, unitVector(2),
        positive(2, 0.2))
        .map(([b, c, jitter, speed]) => {
            const center = mul(0.5, add(b.min, b.max));
            const aim = add(sub(center, c.center), mul(2, jitter));
            const len = Math.sqrt(dot(aim, aim));
            const cv = len > 1e-6 ? mul(speed / len, aim) : mul(speed, jitter);
            return { b, c, cv };
        });

    it('an initial overlap is reported with contact time zero', () => {
        check(arbMoving, ({ b, c, cv }) => {
            const g = gap(b, Vector.zero(2), c, cv, 0);
            if (g >= 0) {
                return;
            }
            const res = fi.find(b, Vector.zero(2), c, cv);
            expect(res.intersectionType).toBe(Type.initiallyOverlapping);
            expect(res.contactTime).toBe(0);
            // The reported point is inside the box and inside the circle.
            for (let i = 0; i < 2; ++i) {
                expect(res.contactPoint.get(i))
                    .toBeGreaterThanOrEqual(b.min.get(i) - 1e-9);
                expect(res.contactPoint.get(i))
                    .toBeLessThanOrEqual(b.max.get(i) + 1e-9);
            }
        });
    });

    it('a reported contact time is the first time the objects touch', () => {
        check(arbMoving, ({ b, c, cv }) => {
            const res = fi.find(b, Vector.zero(2), c, cv);
            if (res.intersectionType !== Type.contact
                || res.contactTime === 0) {
                return;
            }
            const t = res.contactTime;
            expect(Number.isFinite(t)).toBe(true);
            expect(t).toBeGreaterThan(0);
            expectClose(gap(b, Vector.zero(2), c, cv, t), 0, 1e-7, 1e-7);
            for (let k = 1; k < 64; ++k) {
                expect(gap(b, Vector.zero(2), c, cv, (t * k) / 64))
                    .toBeGreaterThan(-1e-7);
            }
            // The contact point is on the box boundary and at distance
            // radius from the circle center at the contact time.
            const p = res.contactPoint;
            const moved = sub(p, add(c.center, mul(t, cv)));
            expectClose(Math.sqrt(dot(moved, moved)), c.radius, 1e-6, 1e-6);
            for (let i = 0; i < 2; ++i) {
                expect(p.get(i)).toBeGreaterThanOrEqual(b.min.get(i) - 1e-7);
                expect(p.get(i)).toBeLessThanOrEqual(b.max.get(i) + 1e-7);
            }
        });
    });

    it('a reported no-contact is confirmed by sampling the motion', () => {
        check(arbMoving, ({ b, c, cv }) => {
            const res = fi.find(b, Vector.zero(2), c, cv);
            if (res.intersectionType !== Type.noContact) {
                return;
            }
            for (let k = 0; k <= 400; ++k) {
                expect(gap(b, Vector.zero(2), c, cv, (20 * k) / 400))
                    .toBeGreaterThan(-1e-7);
            }
            expect(res.contactTime).toBe(0);
            expect(res.contactPoint.equals(Vector.zero(2))).toBe(true);
        }, 80);
    });

    it('the result depends only on the relative velocity', () => {
        check(fc.tuple(arbMoving, wellScaledVector(2, -2, 2)),
            ([{ b, c, cv }, drift]) => {
                const a = fi.find(b, Vector.zero(2), c, cv);
                const d = fi.find(b, drift, c, add(cv, drift));
                expect(d.intersectionType).toBe(a.intersectionType);
                expectClose(d.contactTime, a.contactTime, 1e-9, 1e-9);
            });
    });

    it('is equivariant under translation and under reflection in the axes',
        () => {
            check(fc.tuple(arbMoving, wellScaledVector(2, -5, 5)),
                ([{ b, c, cv }, shift]) => {
                    const a = fi.find(b, Vector.zero(2), c, cv);

                    const bT = AlignedBox.fromMinMax(add(b.min, shift),
                        add(b.max, shift));
                    const cT = Hypersphere.fromCenterRadius(
                        add(c.center, shift), c.radius);
                    const t = fi.find(bT, Vector.zero(2), cT, cv);
                    expect(t.intersectionType).toBe(a.intersectionType);
                    expectClose(t.contactTime, a.contactTime, 1e-9, 1e-9);
                    if (a.intersectionType !== Type.noContact) {
                        expectVectorClose(t.contactPoint,
                            add(a.contactPoint, shift), 1e-8, 1e-8);
                    }

                    // Reflect x -> -x. The query mirrors the circle center
                    // into the first quadrant, so this exercises the sign
                    // bookkeeping.
                    const flip = (v: Vector): Vector =>
                        Vector.fromArray([-v.get(0), v.get(1)]);
                    const bF = AlignedBox.fromMinMax(
                        Vector.fromArray([-b.max.get(0), b.min.get(1)]),
                        Vector.fromArray([-b.min.get(0), b.max.get(1)]));
                    const cF = Hypersphere.fromCenterRadius(flip(c.center),
                        c.radius);
                    const f = fi.find(bF, Vector.zero(2), cF, flip(cv));
                    expect(f.intersectionType).toBe(a.intersectionType);
                    expectClose(f.contactTime, a.contactTime, 1e-9, 1e-9);
                    if (a.intersectionType !== Type.noContact) {
                        expectVectorClose(f.contactPoint, flip(a.contactPoint),
                            1e-8, 1e-8);
                    }
                });
        });

    it('reports no NaN in any field', () => {
        check(arbMoving, ({ b, c, cv }) => {
            const res = fi.find(b, Vector.zero(2), c, cv);
            expect(Number.isNaN(res.contactTime)).toBe(false);
            expect(Number.isNaN(res.contactPoint.get(0))).toBe(false);
            expect(Number.isNaN(res.contactPoint.get(1))).toBe(false);
        });
    });

    it('rejects inputs of the wrong dimension', () => {
        const b3 = AlignedBox.fromMinMax(Vector.zero(3), Vector.zero(3));
        const c2 = circle(0, 0, 1);
        expect(() => ti.test(b3, c2)).toThrow('mismatched sizes');
        expect(() => fi.find(b3, Vector.zero(3), c2, Vector.zero(2)))
            .toThrow('mismatched sizes');
    });
});

// ---------------------------------------------------------------------------
// Minkowski-sum first-contact oracle (sibling check of the
// IntrAlignedBox3Sphere3 dynamic-query fix, #490).
//
// The moving circle touches the box exactly when the ray C + t*V -- C the
// circle center relative to the box center, V the relative velocity -- meets
// the Minkowski sum of the box with the disk of radius r. In 2D that sum is
// bounded by the four box edges pushed out by r and the four circles of
// radius r at the box vertices. Each of those eight pieces is a subset of the
// sum, so a ray starting outside meets each piece no earlier than it enters
// the sum, and it enters the sum on one of them: the smallest nonnegative
// parameter over all pieces is the first contact time. The oracle uses none
// of the query's Voronoi-region case analysis.
//
// The 3D sibling accepted the first rounded-edge probe that reported a hit,
// which is only an upper bound for the entry time; the 2D case analysis
// instead *selects* one piece, because from a point outside a convex planar
// region the boundary pieces are met in a single angular order. Each of
// upstream's tests is the sign of Cross(V, P - C) for P one of the junctions
// where an arc meets a flat edge (P = K +/- radius * e_i), so the decision
// tree is a search on that order and lands on the piece the ray actually
// enters. These properties pin that selection against the oracle.
// ---------------------------------------------------------------------------

describe('IntrAlignedBox2Circle2 Minkowski-sum oracle', () => {
    const fi = new IntrAlignedBox2Circle2FI();

    type OracleResult = { type: Type, t: number, point: number[] };

    // First contact of the ray C + t*V with the Minkowski sum of the box
    // [-K,K] and the disk of radius 'radius'.
    function minkowskiFirstContact(K: readonly number[],
        C: readonly number[], radius: number,
        V: readonly number[]): OracleResult {
        let sqrLen = 0;
        for (let i = 0; i < 2; ++i) {
            const d = Math.max(Math.abs(C[i]) - K[i], 0);
            sqrLen += d * d;
        }
        if (Math.sqrt(sqrLen) <= radius) {
            return {
                type: Type.initiallyOverlapping, t: 0, point: [C[0], C[1]]
            };
        }

        let best = Number.POSITIVE_INFINITY;
        let bestPoint = [0, 0];
        const consider = (t: number, q: number[]): void => {
            if (t >= 0 && t < best) {
                best = t;
                bestPoint = q;
            }
        };

        // The four edges, offset by 'radius' along their outward normals.
        for (let i = 0; i < 2; ++i) {
            const j = 1 - i;
            for (const sgn of [-1, 1]) {
                if (V[i] === 0) {
                    continue;
                }
                const t = (sgn * (K[i] + radius) - C[i]) / V[i];
                if (t < 0) {
                    continue;
                }
                const pj = C[j] + t * V[j];
                if (Math.abs(pj) <= K[j]) {
                    const q = [0, 0];
                    q[i] = sgn * K[i];
                    q[j] = pj;
                    consider(t, q);
                }
            }
        }

        // The four circles of radius 'radius' at the box vertices.
        const a2 = V[0] * V[0] + V[1] * V[1];
        if (a2 > 0) {
            for (const s0 of [-1, 1]) {
                for (const s1 of [-1, 1]) {
                    const q = [s0 * K[0], s1 * K[1]];
                    const e0 = C[0] - q[0];
                    const e1 = C[1] - q[1];
                    const a1 = V[0] * e0 + V[1] * e1;
                    const a0 = e0 * e0 + e1 * e1 - radius * radius;
                    const discr = a1 * a1 - a2 * a0;
                    if (discr < 0) {
                        continue;
                    }
                    consider((-a1 - Math.sqrt(discr)) / a2, q);
                }
            }
        }

        if (!Number.isFinite(best)) {
            return { type: Type.noContact, t: 0, point: [0, 0] };
        }
        return { type: Type.contact, t: best, point: bestPoint };
    }

    // The signed gap between the moving circle and the box [-K,K].
    function gapOf(K: readonly number[], C: readonly number[], radius: number,
        V: readonly number[], t: number): number {
        let sqrLen = 0;
        for (let i = 0; i < 2; ++i) {
            const d = Math.max(Math.abs(C[i] + t * V[i]) - K[i], 0);
            sqrLen += d * d;
        }
        return Math.sqrt(sqrLen) - radius;
    }

    // The gap is convex in t (the distance from a point to a convex set along
    // a line is convex), so a ternary search finds its minimum.
    function minimumGap(K: readonly number[], C: readonly number[],
        radius: number, V: readonly number[], tHi: number): number {
        let lo = 0;
        let hi = tHi;
        for (let i = 0; i < 100; ++i) {
            const m0 = lo + (hi - lo) / 3;
            const m1 = hi - (hi - lo) / 3;
            if (gapOf(K, C, radius, V, m0) <= gapOf(K, C, radius, V, m1)) {
                hi = m1;
            }
            else {
                lo = m0;
            }
        }
        return gapOf(K, C, radius, V, 0.5 * (lo + hi));
    }

    // Compare one configuration -- box [-K,K] translated to 'origin', circle
    // center origin+C, relative velocity V -- against the oracle. Returns
    // false when the configuration is skipped: initially overlapping, at
    // rest, or tangential. A tangential configuration is one whose minimum
    // gap is within rounding of zero; the query and the oracle round the
    // discriminants of the same quadratics differently there and can disagree
    // on whether the objects touch at all.
    function checkAgainstOracle(K: number[], C: number[], radius: number,
        V: number[], origin: number[] = [0, 0]): boolean {
        const speed = Math.hypot(V[0], V[1]);
        if (speed === 0 || gapOf(K, C, radius, V, 0) <= 0) {
            return false;
        }
        const size = Math.hypot(K[0], K[1]) + radius + Math.hypot(C[0], C[1]);
        const oracle = minkowskiFirstContact(K, C, radius, V);
        const tHi = oracle.type === Type.contact
            ? 2 * oracle.t + size / speed
            : 100 * size / speed;
        if (Math.abs(minimumGap(K, C, radius, V, tHi)) < 1e-9 * size) {
            return false;
        }

        const b = box(origin[0] - K[0], origin[1] - K[1],
            origin[0] + K[0], origin[1] + K[1]);
        const c = circle(origin[0] + C[0], origin[1] + C[1], radius);
        const res = fi.find(b, v2(0, 0), c, v2(V[0], V[1]));
        expect(res.intersectionType).toBe(oracle.type);
        if (oracle.type === Type.contact) {
            // The query and the oracle solve the same quadratics, so they
            // agree to the rounding of their coefficients; a time tolerance
            // is a length tolerance divided by the speed.
            expectClose(res.contactTime, oracle.t, 1e-9 * size / speed, 1e-9);
            // The contact point is on the box and at distance 'radius' from
            // the circle center at the contact time. (It is not compared with
            // the oracle's point directly: when the ray enters within a
            // whisker of an arc/edge junction, two pieces give contact times
            // that differ only to second order while their points differ to
            // first order, so such a comparison would be flaky where the
            // times are not.)
            const px = origin[0] + C[0] + res.contactTime * V[0];
            const py = origin[1] + C[1] + res.contactTime * V[1];
            expectClose(Math.hypot(res.contactPoint.get(0) - px,
                res.contactPoint.get(1) - py), radius, 1e-7 * size, 1e-7);
            for (let i = 0; i < 2; ++i) {
                expect(res.contactPoint.get(i))
                    .toBeGreaterThanOrEqual(b.min.get(i) - 1e-9 * size);
                expect(res.contactPoint.get(i))
                    .toBeLessThanOrEqual(b.max.get(i) + 1e-9 * size);
            }
        }
        return true;
    }

    // x advanced by 'ulps' representable doubles away from zero.
    function nudge(x: number, ulps: number): number {
        if (ulps === 0 || x === 0 || !Number.isFinite(x)) {
            return x;
        }
        const view = new DataView(new ArrayBuffer(8));
        view.setFloat64(0, x);
        const bits = view.getBigUint64(0);
        view.setBigUint64(0, x > 0 ? bits + BigInt(ulps) : bits - BigInt(ulps));
        return view.getFloat64(0);
    }

    it('agrees with the oracle on fast-check configurations', () => {
        const arbConfig = fc.tuple(
            fc.array(positive(3, 0.05), { minLength: 2, maxLength: 2 }),
            positive(3, 0.05), unitVector(2), positive(6, 1.05),
            unitVector(2), positive(3, 0.1), wellScaledVector(2, -5, 5));
        let tested = 0;
        check(arbConfig,
            ([K, radius, dir, far, jitter, speed, origin]) => {
                // Start outside the rounded box and aim the motion back at
                // it, with a jitter that makes the corner and the face cases
                // both common.
                const out = (Math.hypot(K[0], K[1]) + radius) * far;
                const C = [dir.get(0) * out, dir.get(1) * out];
                const aim = [-C[0] + 2 * jitter.get(0) * out / far,
                    -C[1] + 2 * jitter.get(1) * out / far];
                const len = Math.hypot(aim[0], aim[1]);
                if (len < 1e-6) {
                    return;
                }
                const V = [aim[0] * speed / len, aim[1] * speed / len];
                if (checkAgainstOracle(K, C, radius, V,
                    [origin.get(0), origin.get(1)])) {
                    ++tested;
                }
            }, 400);
        expect(tested).toBeGreaterThan(250);
    }, 30000);

    it('agrees with the oracle on corner-aimed motions across scales', () => {
        // Motions aimed at a box corner, with a jitter of the order of the
        // corner region, so the ray often just clips or just misses the
        // rounded corner. That is the family the 3D defects came from. The
        // generator spans three decades of scale.
        const rnd = seededRandom(0x2f6e2b1);
        let tested = 0;
        let atVertex = 0;
        let onEdge = 0;
        for (let iter = 0; iter < 6000; ++iter) {
            const scale = Math.exp(6 * rnd() - 3);
            const K = [scale * (0.05 + 2 * rnd()), scale * (0.05 + 2 * rnd())];
            const radius = scale * (0.02 + 1.5 * rnd());
            const angle = 2 * Math.PI * rnd();
            const dist = scale * (2 + 10 * rnd());
            const C = [dist * Math.cos(angle), dist * Math.sin(angle)];
            const target = [K[0] * (rnd() < 0.5 ? -1 : 1),
                K[1] * (rnd() < 0.5 ? -1 : 1)];
            const aim = [
                target[0] - C[0] + (radius + K[0]) * 0.25 * (2 * rnd() - 1),
                target[1] - C[1] + (radius + K[1]) * 0.25 * (2 * rnd() - 1)];
            const len = Math.hypot(aim[0], aim[1]);
            if (len < 1e-12) {
                continue;
            }
            const speed = scale * (0.1 + 3 * rnd());
            const V = [aim[0] * speed / len, aim[1] * speed / len];
            if (!checkAgainstOracle(K, C, radius, V)) {
                continue;
            }
            ++tested;
            const res = fi.find(box(-K[0], -K[1], K[0], K[1]), v2(0, 0),
                circle(C[0], C[1], radius), v2(V[0], V[1]));
            if (res.intersectionType === Type.contact) {
                const corner = Math.abs(Math.abs(res.contactPoint.get(0))
                    - K[0]) < 1e-12 * scale
                    && Math.abs(Math.abs(res.contactPoint.get(1)) - K[1])
                    < 1e-12 * scale;
                if (corner) {
                    ++atVertex;
                }
                else {
                    ++onEdge;
                }
            }
        }
        expect(tested).toBeGreaterThan(3000);
        // Both the rounded-vertex and the flat-edge branches are reached.
        expect(atVertex).toBeGreaterThan(200);
        expect(onEdge).toBeGreaterThan(200);
    }, 30000);

    it('agrees with the oracle when aimed at an arc/edge junction', () => {
        // The ray is aimed exactly at one of the eight junction points where
        // a rounded corner meets a flat edge, then its direction is perturbed
        // by a few ulps. Those points are where upstream's sign tests change
        // branch, so this is where a case analysis that selects the wrong
        // piece -- or has no case for a piece -- shows up.
        const rnd = seededRandom(0x13579bd);
        let tested = 0;
        for (let iter = 0; iter < 6000; ++iter) {
            const scale = Math.exp(6 * rnd() - 3);
            const K = [scale * (0.05 + 2 * rnd()), scale * (0.05 + 2 * rnd())];
            const radius = scale * (0.02 + 1.5 * rnd());
            const sx = rnd() < 0.5 ? -1 : 1;
            const sy = rnd() < 0.5 ? -1 : 1;
            const junction = rnd() < 0.5
                ? [sx * (K[0] + radius), sy * K[1]]
                : [sx * K[0], sy * (K[1] + radius)];
            const angle = 2 * Math.PI * rnd();
            const dist = scale * (2 + 10 * rnd());
            const C = [dist * Math.cos(angle), dist * Math.sin(angle)];
            const aim = [junction[0] - C[0], junction[1] - C[1]];
            const len = Math.hypot(aim[0], aim[1]);
            if (len < 1e-12) {
                continue;
            }
            const speed = scale * (0.1 + 3 * rnd());
            const ulps = Math.floor(rnd() * 9) - 4;
            const V = [nudge(aim[0] * speed / len, ulps),
                nudge(aim[1] * speed / len, -ulps)];
            if (checkAgainstOracle(K, C, radius, V)) {
                ++tested;
            }
        }
        expect(tested).toBeGreaterThan(3000);
    }, 30000);

    it('agrees with the oracle from inside the corner square', () => {
        // The circle center lies in [K, K + radius]^2 but outside the quarter
        // disk: upstream's VertexSeparated, which probes the rounded vertex
        // only. (The 3D sibling had to probe the three rounded edges meeting
        // the vertex as well; in 2D a ray from the corner square can reach a
        // flat edge only through the arc, because the segment to any point of
        // the edge crosses the axis through the vertex strictly inside the
        // vertex circle.)
        const rnd = seededRandom(0x7fedcba);
        let tested = 0;
        let contacts = 0;
        for (let iter = 0; iter < 6000; ++iter) {
            const scale = Math.exp(6 * rnd() - 3);
            const K = [scale * (0.05 + 2 * rnd()), scale * (0.05 + 2 * rnd())];
            const radius = scale * (0.02 + 1.5 * rnd());
            let dx = 0;
            let dy = 0;
            for (let k = 0; k < 64; ++k) {
                dx = radius * rnd();
                dy = radius * rnd();
                if (dx * dx + dy * dy > radius * radius) {
                    break;
                }
            }
            if (dx * dx + dy * dy <= radius * radius) {
                continue;
            }
            const C = [(rnd() < 0.5 ? -1 : 1) * (K[0] + dx),
                (rnd() < 0.5 ? -1 : 1) * (K[1] + dy)];
            const speed = scale * (0.1 + 3 * rnd());
            const angle = 2 * Math.PI * rnd();
            const V = [speed * Math.cos(angle), speed * Math.sin(angle)];
            if (!checkAgainstOracle(K, C, radius, V)) {
                continue;
            }
            ++tested;
            if (fi.find(box(-K[0], -K[1], K[0], K[1]), v2(0, 0),
                circle(C[0], C[1], radius), v2(V[0], V[1])).intersectionType
                === Type.contact) {
                ++contacts;
            }
        }
        expect(tested).toBeGreaterThan(3000);
        expect(contacts).toBeGreaterThan(300);
    }, 30000);

    it('agrees with the oracle for extreme extent and radius ratios', () => {
        // Boxes whose extents differ by up to four decades, and radii from
        // far below to far above the extents: a thin box is where a rounded
        // corner can span a whole side of the sum.
        const rnd = seededRandom(0x0abcdef);
        let tested = 0;
        for (let iter = 0; iter < 6000; ++iter) {
            const scale = Math.exp(6 * rnd() - 3);
            const K = [scale * Math.exp(8 * rnd() - 4) * 0.1,
                scale * Math.exp(8 * rnd() - 4) * 0.1];
            const radius = scale * Math.exp(8 * rnd() - 4) * 0.1;
            const angle = 2 * Math.PI * rnd();
            const dist = (Math.hypot(K[0], K[1]) + radius) * (1.05 + 8 * rnd());
            const C = [dist * Math.cos(angle), dist * Math.sin(angle)];
            const target = [K[0] * (rnd() < 0.5 ? -1 : 1),
                K[1] * (rnd() < 0.5 ? -1 : 1)];
            const aim = [
                target[0] - C[0] + (radius + K[0]) * 0.5 * (2 * rnd() - 1),
                target[1] - C[1] + (radius + K[1]) * 0.5 * (2 * rnd() - 1)];
            const len = Math.hypot(aim[0], aim[1]);
            if (len < 1e-12) {
                continue;
            }
            const speed = scale * (0.1 + 3 * rnd());
            const V = [aim[0] * speed / len, aim[1] * speed / len];
            if (checkAgainstOracle(K, C, radius, V)) {
                ++tested;
            }
        }
        expect(tested).toBeGreaterThan(3000);
    }, 30000);

    it('agrees with the oracle on integer configurations', () => {
        // Small integers are exact in binary64, so upstream's sign tests hit
        // their >= 0 and <= 0 boundaries exactly rather than by a rounding
        // accident; this pins the tie-breaking of the case analysis.
        const rnd = seededRandom(0x2468ace);
        let tested = 0;
        for (let iter = 0; iter < 6000; ++iter) {
            const m = 6;
            const K = [1 + Math.floor(rnd() * m), 1 + Math.floor(rnd() * m)];
            const radius = 1 + Math.floor(rnd() * m);
            const C = [
                (rnd() < 0.5 ? -1 : 1)
                * (Math.floor(rnd() * 2 * m) - m + K[0] + radius),
                (rnd() < 0.5 ? -1 : 1)
                * (Math.floor(rnd() * 2 * m) - m + K[1] + radius)];
            const V = [Math.floor(rnd() * 9) - 4, Math.floor(rnd() * 9) - 4];
            if (checkAgainstOracle(K, C, radius, V)) {
                ++tested;
            }
        }
        expect(tested).toBeGreaterThan(1000);
    }, 30000);

    it('agrees with the oracle on degenerate boxes, radii and velocities',
        () => {
            // Zero radius (a moving point), zero extents (a segment or a
            // point box), axis-aligned velocities, and circle centers placed
            // exactly on the Voronoi-region boundaries of the case analysis.
            const rnd = seededRandom(0x555aaa1);
            let tested = 0;
            for (let iter = 0; iter < 6000; ++iter) {
                const scale = Math.exp(6 * rnd() - 3);
                let K = [scale * (0.05 + 2 * rnd()),
                    scale * (0.05 + 2 * rnd())];
                let radius = scale * (0.02 + 1.5 * rnd());
                const pick = Math.floor(rnd() * 4);
                if (pick === 0) {
                    radius = 0;
                }
                if (pick === 1) {
                    K = [0, 0];
                }
                if (pick === 2) {
                    K = [K[0], 0];
                }
                const angle = 2 * Math.PI * rnd();
                const dist = scale * (2 + 10 * rnd());
                let C = [dist * Math.cos(angle), dist * Math.sin(angle)];
                const snap = Math.floor(rnd() * 4);
                if (snap === 0) {
                    C = [C[0], K[1] + radius];
                }
                if (snap === 1) {
                    C = [K[0] + radius, C[1]];
                }
                if (snap === 2) {
                    C = [C[0], K[1]];
                }
                const target = [K[0] * (rnd() < 0.5 ? -1 : 1),
                    K[1] * (rnd() < 0.5 ? -1 : 1)];
                let aim = [
                    target[0] - C[0] + (radius + K[0]) * 0.5 * (2 * rnd() - 1),
                    target[1] - C[1] + (radius + K[1]) * 0.5 * (2 * rnd() - 1)];
                const axis = Math.floor(rnd() * 4);
                if (axis === 0) {
                    aim = [aim[0], 0];
                }
                if (axis === 1) {
                    aim = [0, aim[1]];
                }
                const len = Math.hypot(aim[0], aim[1]);
                if (len < 1e-12) {
                    continue;
                }
                const speed = scale * (0.1 + 3 * rnd());
                const V = [aim[0] * speed / len, aim[1] * speed / len];
                if (checkAgainstOracle(K, C, radius, V)) {
                    ++tested;
                }
            }
            expect(tested).toBeGreaterThan(2000);
        }, 30000);

    // --- worked configurations at the junctions ---------------------------
    // Box [-1,1]^2 with radius 1. The rounded box has junctions at (2,1) and
    // (1,2) in the first quadrant, and all of the following rays move along
    // (-1,-1) directions, so the arithmetic is exact in binary64.

    it('takes the rounded vertex when aimed exactly at a junction', () => {
        // Aimed exactly at the junction (1,2) of the top edge and the corner
        // arc: the sign test r*V[0] - Dot(V,Perp(delta)) is exactly zero.
        const a = fi.find(unitBox, v2(0, 0), circle(3, 4, 1), v2(-1, -1));
        expect(a.intersectionType).toBe(Type.contact);
        expect(a.contactTime).toBe(2);
        expect(a.contactPoint.values).toEqual([1, 1]);
        expect(checkAgainstOracle([1, 1], [3, 4], 1, [-1, -1])).toBe(true);

        // Aimed exactly at the junction (2,1) of the right edge and the same
        // arc: -r*V[1] - Dot(V,Perp(delta)) is exactly zero.
        const b = fi.find(unitBox, v2(0, 0), circle(4, 3, 1), v2(-1, -1));
        expect(b.intersectionType).toBe(Type.contact);
        expect(b.contactTime).toBe(2);
        expect(b.contactPoint.values).toEqual([1, 1]);
        expect(checkAgainstOracle([1, 1], [4, 3], 1, [-1, -1])).toBe(true);
    });

    it('takes the flat edge past the junction and the far arc past its end',
        () => {
            // Passes inside the junction (1,2) and lands on the top edge at
            // (0,2), so the contact point is (0,1) on the box.
            const a = fi.find(unitBox, v2(0, 0), circle(3, 4, 1),
                v2(-1.5, -1));
            expect(a.intersectionType).toBe(Type.contact);
            expect(a.contactTime).toBe(2);
            expect(a.contactPoint.values).toEqual([0, 1]);
            expect(checkAgainstOracle([1, 1], [3, 4], 1, [-1.5, -1]))
                .toBe(true);

            // Passes beyond the far end (-1,2) of the top edge, so the first
            // contact is on the arc at the vertex (-1,1).
            const b = fi.find(unitBox, v2(0, 0), circle(3, 4, 1),
                v2(-2.05, -1));
            expect(b.intersectionType).toBe(Type.contact);
            expect(b.contactPoint.values).toEqual([-1, 1]);
            expect(checkAgainstOracle([1, 1], [3, 4], 1, [-2.05, -1]))
                .toBe(true);

            // A little further and the ray misses the rounded box entirely.
            const c = fi.find(unitBox, v2(0, 0), circle(3, 4, 1),
                v2(-2.2, -1));
            expect(c.intersectionType).toBe(Type.noContact);
            expect(checkAgainstOracle([1, 1], [3, 4], 1, [-2.2, -1]))
                .toBe(true);
        });
});
