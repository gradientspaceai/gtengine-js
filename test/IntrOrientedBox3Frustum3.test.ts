import { describe, it, expect } from 'vitest';
import { Frustum3 } from '../src/Frustum3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import { IntrOrientedBox3Frustum3TI } from '../src/IntrOrientedBox3Frustum3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// An axis-aligned oriented box.
function alignedBox(center: number[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(vec(center),
        [vec([1, 0, 0]), vec([0, 1, 0]), vec([0, 0, 1])], vec(extent));
}

// An oriented box whose first axis is the given direction.
function tiltedBox(center: number[], axis0: number[], extent: number[]):
    OrientedBox {
    const a0 = vec(axis0);
    normalize(a0);
    const basis = [a0, new Vector(3), new Vector(3)];
    computeOrthogonalComplement3(1, basis);
    return OrientedBox.fromCenterAxisExtent(vec(center), basis, vec(extent));
}

// Independent point-in-frustum test.
function insideFrustum(f: Frustum3, x: Vector): boolean {
    const diff = sub(x, f.origin);
    const d = dot(diff, f.dVector);
    if (d < f.dMin || d > f.dMax) {
        return false;
    }
    const scale = d / f.dMin;
    return Math.abs(dot(diff, f.rVector)) <= f.rBound * scale
        && Math.abs(dot(diff, f.uVector)) <= f.uBound * scale;
}

describe('IntrOrientedBox3Frustum3', () => {
    const ti = new IntrOrientedBox3Frustum3TI();

    // The default frustum: origin (0,0,0), view direction (0,0,1),
    // dMin = 1, dMax = 2, uBound = rBound = 1.
    const frustum = new Frustum3();

    it('accepts a small box inside the frustum', () => {
        const b = alignedBox([0, 0, 1.5], [0.1, 0.1, 0.1]);
        expect(ti.test(b, frustum).intersect).toBe(true);
    });

    it('rejects a box behind the eye point', () => {
        const b = alignedBox([0, 0, -5], [1, 1, 1]);
        expect(ti.test(b, frustum).intersect).toBe(false);
    });

    it('rejects a box beyond the far plane', () => {
        const b = alignedBox([0, 0, 10], [1, 1, 1]);
        expect(ti.test(b, frustum).intersect).toBe(false);
    });

    it('rejects a box entirely to the right of the frustum', () => {
        const b = alignedBox([10, 0, 1.5], [1, 1, 1]);
        expect(ti.test(b, frustum).intersect).toBe(false);
    });

    it('rejects a box entirely above the frustum', () => {
        const b = alignedBox([0, 10, 1.5], [1, 1, 1]);
        expect(ti.test(b, frustum).intersect).toBe(false);
    });

    it('accepts a box straddling the near plane', () => {
        const b = alignedBox([0, 0, 1], [0.5, 0.5, 0.5]);
        expect(ti.test(b, frustum).intersect).toBe(true);
    });

    it('accepts a large box that encloses the frustum', () => {
        const b = alignedBox([0, 0, 0], [20, 20, 20]);
        expect(ti.test(b, frustum).intersect).toBe(true);
    });

    it('accepts a tilted box poking into the frustum side', () => {
        const b = tiltedBox([1.6, 0, 1.5], [1, 0, 1], [0.6, 0.2, 0.2]);
        expect(ti.test(b, frustum).intersect).toBe(true);
    });

    it('handles a degenerate box with zero extents as a point', () => {
        const insidePoint = alignedBox([0, 0, 1.5], [0, 0, 0]);
        expect(ti.test(insidePoint, frustum).intersect).toBe(true);
        const outsidePoint = alignedBox([0, 0, 0.5], [0, 0, 0]);
        expect(ti.test(outsidePoint, frustum).intersect).toBe(false);
    });

    it('reports an intersection whenever a sampled point lies in both', () => {
        const rand = makeRandom(1414213);
        const f = Frustum3.fromParameters(
            vec([0.5, -1, 0.25]), vec([0, 0, 1]), vec([0, 1, 0]),
            vec([1, 0, 0]), 1, 4, 0.8, 1.2);
        let hits = 0;
        for (let trial = 0; trial < 250; ++trial) {
            const b = tiltedBox(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 1],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [0.2 + rand(), 0.2 + rand(), 0.2 + rand()]);
            const result = ti.test(b, f).intersect;

            // Sample the box on a grid; any sample inside the frustum forces
            // the query to report an intersection.
            let found = false;
            const n = 12;
            for (let i = 0; i <= n && !found; ++i) {
                const a0 = -b.extent.values[0] + (2 * b.extent.values[0] * i) / n;
                for (let j = 0; j <= n && !found; ++j) {
                    const a1 = -b.extent.values[1]
                        + (2 * b.extent.values[1] * j) / n;
                    for (let k = 0; k <= n; ++k) {
                        const a2 = -b.extent.values[2]
                            + (2 * b.extent.values[2] * k) / n;
                        const p = add(b.center, add(mul(a0, b.axis[0]),
                            add(mul(a1, b.axis[1]), mul(a2, b.axis[2]))));
                        if (insideFrustum(f, p)) {
                            found = true;
                            break;
                        }
                    }
                }
            }

            if (found) {
                ++hits;
                expect(result).toBe(true);
            }
        }
        expect(hits).toBeGreaterThan(20);
    });

    it('reports an intersection when a frustum vertex is inside the box', () => {
        const rand = makeRandom(99991);
        const vertices = frustum.computeVertices();
        let hits = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const b = tiltedBox(
                [4 * rand() - 2, 4 * rand() - 2, 4 * rand()],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [0.1 + rand(), 0.1 + rand(), 0.1 + rand()]);
            let found = false;
            for (const v of vertices) {
                const diff = sub(v, b.center);
                if (Math.abs(dot(diff, b.axis[0])) <= b.extent.values[0]
                    && Math.abs(dot(diff, b.axis[1])) <= b.extent.values[1]
                    && Math.abs(dot(diff, b.axis[2])) <= b.extent.values[2]) {
                    found = true;
                    break;
                }
            }
            if (found) {
                ++hits;
                expect(ti.test(b, frustum).intersect).toBe(true);
            }
        }
        expect(hits).toBeGreaterThan(5);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, wellScaledVector, seededRandom
} from './helpers/arbitraries.js';
import { length as vlength } from '../src/Vector.js';

const frustumArb = fc.tuple(wellScaledVector(3, -3, 3), rotationFrame(3),
    fc.double({ min: 0.3, max: 2, noNaN: true }),
    fc.double({ min: 1.5, max: 6, noNaN: true }),
    fc.double({ min: 0.2, max: 2, noNaN: true }),
    fc.double({ min: 0.2, max: 2, noNaN: true }))
    .map(([e, R, dMin, extra, u, r]) => Frustum3.fromParameters(e, R[0], R[1],
        R[2], dMin, dMin + extra, u, r));

const boxArbV31 = fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
    fc.double({ min: 0.1, max: 2, noNaN: true }),
    fc.double({ min: 0.1, max: 2, noNaN: true }),
    fc.double({ min: 0.1, max: 2, noNaN: true }))
    .map(([c, R, e0, e1, e2]) => OrientedBox.fromCenterAxisExtent(c, R,
        Vector.fromArray([e0, e1, e2])));

// Membership in the solid frustum, from the definition in Frustum3.h: a point
// X is inside when n <= Dot(D,X-E) <= f and the (r,u) bounds scale linearly
// with the depth.
function inFrustum(f: Frustum3, p: Vector, tol: number): boolean {
    const q = sub(p, f.origin);
    const d = dot(f.dVector, q);
    if (d < f.dMin - tol || d > f.dMax + tol) {
        return false;
    }
    const scale = d / f.dMin;
    return Math.abs(dot(f.rVector, q)) <= f.rBound * scale + tol
        && Math.abs(dot(f.uVector, q)) <= f.uBound * scale + tol;
}

describe('IntrOrientedBox3Frustum3 verification', () => {
    const tiq = new IntrOrientedBox3Frustum3TI();

    it('reports an intersection whenever a sampled common point exists', () => {
        const rnd = seededRandom(0x11c93f6a);
        for (let trial = 0; trial < 150; ++trial) {
            const a = rnd() * 2 * Math.PI;
            const R = [
                Vector.fromArray([Math.cos(a), Math.sin(a), 0]),
                Vector.fromArray([-Math.sin(a), Math.cos(a), 0]),
                Vector.fromArray([0, 0, 1])];
            const f = Frustum3.fromParameters(
                Vector.fromArray([rnd() * 4 - 2, rnd() * 4 - 2,
                    rnd() * 4 - 2]),
                R[2], R[1], R[0], 0.5 + rnd(), 2 + rnd() * 3,
                0.3 + rnd(), 0.3 + rnd());
            const b = OrientedBox.fromCenterAxisExtent(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3,
                    rnd() * 6 - 3]),
                [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                    Vector.fromArray([0, 0, 1])],
                Vector.fromArray([0.2 + rnd(), 0.2 + rnd(), 0.2 + rnd()]));
            const got = tiq.test(b, f).intersect;
            // Sample the box for a point inside the frustum.
            let common = false;
            for (let i = 0; i <= 10 && !common; ++i) {
                for (let j = 0; j <= 10 && !common; ++j) {
                    for (let k = 0; k <= 10; ++k) {
                        const p = add(b.center, add(
                            mul((2 * i / 10 - 1) * b.extent.values[0],
                                b.axis[0]),
                            add(mul((2 * j / 10 - 1) * b.extent.values[1],
                                b.axis[1]),
                            mul((2 * k / 10 - 1) * b.extent.values[2],
                                b.axis[2]))));
                        if (inFrustum(f, p, 0)) {
                            common = true;
                            break;
                        }
                    }
                }
            }
            if (!common) {
                // Sample the frustum for a point inside the box.
                for (let i = 0; i <= 10 && !common; ++i) {
                    const d = f.dMin + ((f.dMax - f.dMin) * i) / 10;
                    const scale = d / f.dMin;
                    for (let j = 0; j <= 10 && !common; ++j) {
                        for (let k = 0; k <= 10; ++k) {
                            const p = add(f.origin, add(mul(d, f.dVector),
                                add(mul((2 * j / 10 - 1) * f.rBound * scale,
                                    f.rVector),
                                mul((2 * k / 10 - 1) * f.uBound * scale,
                                    f.uVector))));
                            const q = sub(p, b.center);
                            if (Math.abs(dot(q, b.axis[0]))
                                <= b.extent.values[0]
                                && Math.abs(dot(q, b.axis[1]))
                                <= b.extent.values[1]
                                && Math.abs(dot(q, b.axis[2]))
                                <= b.extent.values[2]) {
                                common = true;
                                break;
                            }
                        }
                    }
                }
            }
            if (common) {
                expect(got).toBe(true);
            }
        }
    }, 30000);

    it('a box centred on a sampled frustum point always intersects', () => {
        check(fc.tuple(frustumArb, rotationFrame(3),
            fc.double({ min: 0, max: 1, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }),
            fc.double({ min: 0.05, max: 1, noNaN: true })),
            ([f, R, td, tr, tu, e]) => {
                const d = f.dMin + (f.dMax - f.dMin) * td;
                const scale = d / f.dMin;
                const p = add(f.origin, add(mul(d, f.dVector),
                    add(mul(tr * f.rBound * scale, f.rVector),
                        mul(tu * f.uBound * scale, f.uVector))));
                expect(inFrustum(f, p, 1e-9)).toBe(true);
                const b = OrientedBox.fromCenterAxisExtent(p, R,
                    Vector.fromArray([e, e, e]));
                expect(tiq.test(b, f).intersect).toBe(true);
            });
    });

    it('a box far outside never intersects', () => {
        check(fc.tuple(frustumArb, boxArbV31, rotationFrame(3)),
            ([f, b, R]) => {
                const reach = f.dMax + f.rBound * f.dMax / f.dMin
                    + f.uBound * f.dMax / f.dMin + vlength(b.extent) + 1;
                // Move the box behind the frustum origin, opposite dVector.
                const far = OrientedBox.fromCenterAxisExtent(
                    sub(f.origin, mul(2 * reach, f.dVector)), R, b.extent);
                expect(tiq.test(far, f).intersect).toBe(false);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(frustumArb, boxArbV31, rotationFrame(3),
            wellScaledVector(3, -3, 3)), ([f, b, R, tr]) => {
            const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
            const xf = (p: Vector): Vector => add(tr, rot(p));
            const f2 = Frustum3.fromParameters(xf(f.origin), rot(f.dVector),
                rot(f.uVector), rot(f.rVector), f.dMin, f.dMax, f.uBound,
                f.rBound);
            const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                b.axis.map(x => rot(x)), b.extent);
            const a = tiq.test(b, f).intersect;
            const c = tiq.test(b2, f2).intersect;
            if (a === c) {
                return;
            }
            // Only a configuration that grazes the frustum boundary may flip.
            // Verify by shrinking and growing the box: if the two answers
            // bracket the disagreement, the box really is on the boundary.
            const shrink = OrientedBox.fromCenterAxisExtent(b.center, b.axis,
                mul(0.98, b.extent));
            const grow = OrientedBox.fromCenterAxisExtent(b.center, b.axis,
                mul(1.02, b.extent));
            expect(tiq.test(shrink, f).intersect)
                .not.toBe(tiq.test(grow, f).intersect);
        });
    });

    it('the default frustum and a unit box at the near plane intersect', () => {
        const f = new Frustum3();
        const b = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([0, 0, 1.5]),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([0.25, 0.25, 0.25]));
        expect(tiq.test(b, f).intersect).toBe(true);
        // A box behind the eye point misses.
        const behind = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([0, 0, -3]), b.axis, b.extent);
        expect(tiq.test(behind, f).intersect).toBe(false);
    });
});
