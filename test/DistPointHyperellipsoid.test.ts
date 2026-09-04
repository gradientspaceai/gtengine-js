import { describe, expect, it } from 'vitest';
import { DistPointHyperellipsoid } from '../src/DistPointHyperellipsoid.js';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ellipsoid(center: number[], axis: number[][],
    extent: number[]): Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(v(...center),
        axis.map(a => v(...a)), v(...extent));
}

function axisAligned(extent: number[]): Hyperellipsoid {
    const n = extent.length;
    const axes: number[][] = [];
    for (let i = 0; i < n; ++i) {
        const a = new Array<number>(n).fill(0);
        a[i] = 1;
        axes.push(a);
    }
    const center = new Array<number>(n).fill(0);
    return ellipsoid(center, axes, extent);
}

// The value of the hyperellipsoid equation at a point in the coordinate
// system of the hyperellipsoid; it is 1 exactly on the surface.
function equationValue(p: Vector, e: Hyperellipsoid): number {
    const delta = sub(p, e.center);
    let sum = 0;
    for (let i = 0; i < e.dimension; ++i) {
        const y = dot(delta, e.axis[i]) / e.extent.values[i];
        sum += y * y;
    }
    return sum;
}

// Minimum distance from a 2D point to an axis-aligned ellipse, computed by a
// dense angular sampling followed by local refinement.
function bruteForce2D(p: Vector, e0: number, e1: number): number {
    const at = (t: number): number =>
        length(sub(p, v(e0 * Math.cos(t), e1 * Math.sin(t))));
    let best = Number.MAX_VALUE;
    let bt = 0;
    const n = 4096;
    for (let i = 0; i < n; ++i) {
        const t = 2 * Math.PI * i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    let h = 2 * Math.PI / n;
    for (let pass = 0; pass < 80; ++pass) {
        for (const s of [1, -1]) {
            const d = at(bt + s * h);
            if (d < best) {
                best = d;
                bt = bt + s * h;
            }
        }
        h *= 0.7;
    }
    return best;
}

describe('DistPointHyperellipsoid', () => {
    const query = new DistPointHyperellipsoid();

    it('reduces to the circle distance when the extents are equal', () => {
        const circle = axisAligned([2, 2]);
        const p = v(5, 0);
        const result = query.compute(p, circle);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
    });

    it('reduces to the sphere distance in 3D', () => {
        const sphere = axisAligned([3, 3, 3]);
        const p = v(0, 0, 10);
        const result = query.compute(p, sphere);
        expect(result.distance).toBeCloseTo(7, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(3, 10);
    });

    it('returns zero distance for a point on the ellipse', () => {
        const e = axisAligned([2, 1]);
        const p = v(2, 0);
        const result = query.compute(p, e);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
    });

    it('returns the minor semi-axis for the ellipse center', () => {
        // The center is equidistant from the two ends of the minor axis; the
        // algorithm reports one of them.
        const e = axisAligned([3, 1]);
        const result = query.compute(v(0, 0), e);
        expect(result.distance).toBeCloseTo(1, 10);
        expect(Math.abs(result.closest[1].values[1])).toBeCloseTo(1, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 10);
    });

    it('returns the smallest semi-axis for the ellipsoid center', () => {
        const e = axisAligned([4, 3, 2]);
        const result = query.compute(v(0, 0, 0), e);
        expect(result.distance).toBeCloseTo(2, 10);
    });

    it('handles a point inside but off center', () => {
        const e = axisAligned([2, 1]);
        const result = query.compute(v(0.5, 0), e);
        // The closest point is on the minor axis at (0,+-1)? No: the closest
        // ellipse point to an interior point near the center is found by the
        // bisection; verify by the equation and by a brute-force minimum.
        expect(equationValue(result.closest[1], e)).toBeCloseTo(1, 9);
        expect(result.distance).toBeCloseTo(bruteForce2D(v(0.5, 0), 2, 1), 8);
    });

    it('agrees with computeAxisAligned for an axis-aligned hyperellipsoid',
        () => {
            const e = axisAligned([3, 1.5, 0.75]);
            const p = v(1, -2, 4);
            const r0 = query.compute(p, e);
            const r1 = query.computeAxisAligned(p, e.extent);
            expect(r1.distance).toBeCloseTo(r0.distance, 12);
            for (let i = 0; i < 3; ++i) {
                expect(r1.closest[1].values[i])
                    .toBeCloseTo(r0.closest[1].values[i], 12);
            }
        });

    it('handles a rotated and translated ellipse', () => {
        const c = Math.SQRT1_2;
        const e = ellipsoid([1, 2], [[c, c], [-c, c]], [3, 1]);
        // A point on the surface: center + 3 * axis[0].
        const onSurface = add(e.center, mul(3, e.axis[0]));
        const r0 = query.compute(onSurface, e);
        expect(r0.distance).toBeCloseTo(0, 10);

        // A point far along axis[1].
        const outside = add(e.center, mul(5, e.axis[1]));
        const r1 = query.compute(outside, e);
        expect(r1.distance).toBeCloseTo(4, 10);
    });

    it('agrees with a dense sampling for random 2D queries', () => {
        let seed = 112233445;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const e0 = 0.3 + 3 * rand();
            const e1 = 0.3 + 3 * rand();
            const e = axisAligned([e0, e1]);
            const p = v(8 * rand() - 4, 8 * rand() - 4);
            const result = query.compute(p, e);
            expect(result.distance).toBeCloseTo(bruteForce2D(p, e0, e1), 6);
            expect(equationValue(result.closest[1], e)).toBeCloseTo(1, 8);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 10);
        }
    });

    it('produces surface points whose normal is parallel to the offset',
        () => {
            let seed = 998877665;
            const rand = (): number => {
                seed = (1103515245 * seed + 12345) % 2147483648;
                return seed / 2147483648;
            };
            for (let trial = 0; trial < 200; ++trial) {
                const e = axisAligned([0.4 + 2 * rand(), 0.4 + 2 * rand(),
                    0.4 + 2 * rand()]);
                const p = v(8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4);
                const result = query.compute(p, e);
                const x = result.closest[1];
                expect(equationValue(x, e)).toBeCloseTo(1, 7);

                // The gradient of the hyperellipsoid equation at x is
                // 2*(x[i]/e[i]^2); it must be parallel to p - x.
                const grad = v(
                    x.values[0] / (e.extent.values[0] * e.extent.values[0]),
                    x.values[1] / (e.extent.values[1] * e.extent.values[1]),
                    x.values[2] / (e.extent.values[2] * e.extent.values[2]));
                const diff = sub(p, x);
                const lg = length(grad);
                const ld = length(diff);
                if (lg > 1e-8 && ld > 1e-6) {
                    const cosine = dot(grad, diff) / (lg * ld);
                    expect(Math.abs(cosine)).toBeCloseTo(1, 5);
                }
                expect(ld).toBeCloseTo(result.distance, 10);
            }
        });

    it('throws for a dimension mismatch', () => {
        const e = axisAligned([1, 2]);
        expect(() => query.compute(v(1, 2, 3), e)).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistPointHyperellipsoid.h.
// ---------------------------------------------------------------------------

// Extents are bounded away from zero: the bisector divides by e[i] and by the
// smallest extent, and SqrDistanceSpecial divides by e[i]^2 - e[N-1]^2.
const extentArb = (n: number): fc.Arbitrary<number[]> =>
    fc.array(positive(4, 0.2), { minLength: n, maxLength: n });

// rotationFrame only builds frames of dimension 2 and 3; higher dimensions
// use the coordinate frame, which still exercises the permutation and
// reflection bookkeeping of SqrDistance.
const frameArb = (n: number): fc.Arbitrary<Vector[]> => (n <= 3
    ? rotationFrame(n)
    : fc.constant(Array.from({ length: n }, (_, i) => Vector.unit(n, i))));

const hyperellipsoidArb = (n: number): fc.Arbitrary<Hyperellipsoid> =>
    fc.tuple(wellScaledVector(n, -5, 5), frameArb(n), extentArb(n))
        .map(([c, axes, e]) => Hyperellipsoid.fromCenterAxisExtent(c, axes,
            Vector.fromArray(e)));

// Upstream's Bisector brackets the root by smin = z[last] - 1 where z[last]
// is the query-point coordinate along the smallest-extent axis divided by
// that extent. When z[last] is tiny but nonzero, smin rounds to exactly -1,
// the root is closer to -1 than the double-precision spacing there, and the
// final division by (s + pSqr[i]) cancels catastrophically -- the closest
// point leaves the surface and the distance can even be Infinity. See the
// upstream-bug note in src/DistPointHyperellipsoid.ts. Properties that
// compare against an independent computation skip that band; a coordinate of
// exactly zero takes a different branch and is correct, so it is kept.
function isWellConditioned(h: Hyperellipsoid, p: Vector): boolean {
    let jmin = 0;
    for (let i = 1; i < h.dimension; ++i) {
        if (h.extent.values[i] < h.extent.values[jmin]) {
            jmin = i;
        }
    }
    const z = Math.abs(dot(sub(p, h.center), h.axis[jmin]))
        / h.extent.values[jmin];
    return z === 0 || z > 1e-5;
}

// Coordinates of X in the hyperellipsoid frame.
function localCoords(h: Hyperellipsoid, x: Vector): number[] {
    const delta = sub(x, h.center);
    const y: number[] = [];
    for (let i = 0; i < h.dimension; ++i) {
        y.push(dot(delta, h.axis[i]));
    }
    return y;
}

// First-order optimality: at the closest surface point X the residual P - X
// must be parallel to the surface normal, whose frame coordinates are
// y[i]/e[i]^2. Returns the largest 2x2 "cross product" of the two vectors,
// scaled by their magnitudes so the test is relative.
function optimalityResidual(h: Hyperellipsoid, p: Vector,
    x: Vector): number {
    const n = h.dimension;
    const y = localCoords(h, x);
    const z = localCoords(h, p);
    const grad: number[] = [];
    const res: number[] = [];
    for (let i = 0; i < n; ++i) {
        grad.push(y[i] / (h.extent.values[i] * h.extent.values[i]));
        res.push(z[i] - y[i]);
    }
    let gn = 0;
    let rn = 0;
    for (let i = 0; i < n; ++i) {
        gn += grad[i] * grad[i];
        rn += res[i] * res[i];
    }
    gn = Math.sqrt(gn);
    rn = Math.sqrt(rn);
    if (gn === 0 || rn === 0) {
        return 0;
    }
    let worst = 0;
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            const c = Math.abs(res[i] * grad[j] - res[j] * grad[i])
                / (gn * rn);
            if (c > worst) {
                worst = c;
            }
        }
    }
    return worst;
}

// Brute-force minimum distance from a point to a 3D ellipsoid surface, over
// the (theta,phi) parameterization, with a local pattern-search refinement.
function bruteForce3D(p: Vector, h: Hyperellipsoid): number {
    const e = h.extent.values;
    const at = (t: number, u: number): number => {
        const st = Math.sin(t);
        const x = add(h.center, add(mul(e[0] * st * Math.cos(u), h.axis[0]),
            add(mul(e[1] * st * Math.sin(u), h.axis[1]),
                mul(e[2] * Math.cos(t), h.axis[2]))));
        return length(sub(p, x));
    };
    const nt = 90;
    const nu = 180;
    let best = Number.MAX_VALUE;
    let bt = 0;
    let bu = 0;
    for (let i = 0; i <= nt; ++i) {
        const t = Math.PI * i / nt;
        for (let j = 0; j < nu; ++j) {
            const u = 2 * Math.PI * j / nu;
            const d = at(t, u);
            if (d < best) {
                best = d;
                bt = t;
                bu = u;
            }
        }
    }
    let ht = Math.PI / nt;
    let hu = 2 * Math.PI / nu;
    for (let pass = 0; pass < 120; ++pass) {
        for (const [dt, du] of [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1],
            [1, -1], [-1, 1], [-1, -1]]) {
            const d = at(bt + dt * ht, bu + du * hu);
            if (d < best) {
                best = d;
                bt += dt * ht;
                bu += du * hu;
            }
        }
        ht *= 0.75;
        hu *= 0.75;
    }
    return best;
}

describe('DistPointHyperellipsoid verification', () => {
    const query = new DistPointHyperellipsoid();

    for (const n of [2, 3, 4]) {
        it(`result is self consistent in ${n}D and closest[1] is on the surface`,
            () => {
                check(fc.tuple(wellScaledVector(n, -8, 8),
                    hyperellipsoidArb(n)), ([p, h]) => {
                    if (!isWellConditioned(h, p)) {
                        return;
                    }
                    const r = query.compute(p, h);
                    expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                        1e-12);
                    expectVectorClose(r.closest[0], p, 0, 0);
                    expectClose(r.distance,
                        length(sub(r.closest[0], r.closest[1])), 1e-8, 1e-8);
                    // The bisection stops at the double-precision
                    // resolution of the root s. Near the center of the
                    // hyperellipsoid the closest point is ill-conditioned in
                    // s, so the surface equation is satisfied to about 1e-8
                    // rather than to machine precision.
                    expectClose(equationValue(r.closest[1], h), 1, 1e-7,
                        1e-7);
                }, 120);
            }, 30000);
    }

    it('satisfies the first-order optimality condition (2D and 3D)', () => {
        check(fc.tuple(fc.constantFrom(2, 3),
            fc.array(wellScaled(-8, 8), { minLength: 4, maxLength: 4 })),
        ([n, raw]) => {
            const h = Hyperellipsoid.fromCenterAxisExtent(new Vector(n),
                Array.from({ length: n }, (_, i) => Vector.unit(n, i)),
                Vector.fromArray(Array.from({ length: n },
                    (_, i) => 0.5 + Math.abs(raw[i]) % 3)));
            const p = Vector.fromArray(raw.slice(0, n));
            if (!isWellConditioned(h, p)) {
                return;
            }
            const r = query.compute(p, h);
            expect(optimalityResidual(h, p, r.closest[1]))
                .toBeLessThanOrEqual(1e-6);
        }, 200);
    }, 30000);

    it('matches a brute-force minimization over a 2D ellipse', () => {
        check(fc.tuple(wellScaledVector(2, -8, 8), positive(4, 0.2),
            positive(4, 0.2)), ([p, e0, e1]) => {
            const h = axisAligned([e0, e1]);
            if (!isWellConditioned(h, p)) {
                return;
            }
            expectClose(query.compute(p, h).distance, bruteForce2D(p, e0, e1),
                1e-6, 1e-6);
        }, 120);
    }, 30000);

    it('matches a brute-force minimization over a 3D ellipsoid', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), hyperellipsoidArb(3)),
            ([p, h]) => {
                if (!isWellConditioned(h, p)) {
                    return;
                }
                expectClose(query.compute(p, h).distance, bruteForce3D(p, h),
                    1e-5, 1e-5);
            }, 25);
    }, 30000);

    it('is not larger than the distance to any sampled surface point', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), hyperellipsoidArb(3),
            wellScaledVector(3, -1, 1)), ([p, h, u]) => {
            if (!isWellConditioned(h, p)) {
                return;
            }
            const d = query.compute(p, h).distance;
            // Project u onto the surface through the frame coordinates.
            let sum = 0;
            for (let i = 0; i < 3; ++i) {
                sum += u.values[i] * u.values[i];
            }
            if (sum < 1e-6) {
                return;
            }
            const s = 1 / Math.sqrt(sum);
            let q = h.center.clone();
            for (let i = 0; i < 3; ++i) {
                q = add(q, mul(s * u.values[i] * h.extent.values[i],
                    h.axis[i]));
            }
            expect(d).toBeLessThanOrEqual(length(sub(p, q)) + 1e-8);
        }, 120);
    }, 30000);

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), hyperellipsoidArb(3),
            rotationFrame(3), wellScaledVector(3, -4, 4)),
        ([p, h, R, tr]) => {
            const rot = (x: Vector): Vector => {
                let y = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    y = add(y, mul(x.values[i], R[i]));
                }
                return y;
            };
            const moved = Hyperellipsoid.fromCenterAxisExtent(
                add(rot(h.center), tr), h.axis.map(a => rot(a)), h.extent);
            if (!isWellConditioned(h, p)) {
                return;
            }
            // Near the center of a near-spheroid the minimizer is barely
            // determined (a whole circle of surface points is almost
            // equidistant), so the frame-dependent bisection can land on
            // different points; skip that regime.
            if (length(sub(p, h.center))
                < 0.05 * Math.min(...h.extent.values)) {
                return;
            }
            const r0 = query.compute(p, h);
            const r1 = query.compute(add(rot(p), tr), moved);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            // The closest point comes out of a bisection whose iterates
            // depend on the frame, and near a spheroid's symmetry axis the
            // minimizer is barely determined, so the point comparison needs a
            // much looser tolerance than the distance comparison.
            expectVectorClose(add(rot(r0.closest[1]), tr), r1.closest[1],
                1e-4, 1e-4);
        }, 120);
    }, 30000);

    it('computeAxisAligned agrees with compute on the canonical frame', () => {
        check(fc.tuple(fc.constantFrom(2, 3, 4),
            fc.array(wellScaled(-8, 8), { minLength: 4, maxLength: 4 }),
            fc.array(positive(4, 0.2), { minLength: 4, maxLength: 4 })),
        ([n, raw, ext]) => {
            const extent = Vector.fromArray(ext.slice(0, n));
            const h = Hyperellipsoid.fromCenterAxisExtent(new Vector(n),
                Array.from({ length: n }, (_, i) => Vector.unit(n, i)),
                extent);
            const p = Vector.fromArray(raw.slice(0, n));
            const r0 = query.compute(p, h);
            const r1 = query.computeAxisAligned(p, extent);
            // Both paths run the identical code, so they must agree bit for
            // bit even in the ill-conditioned band; guard only against the
            // Infinity case, where the comparison is meaningless.
            if (!Number.isFinite(r0.distance)) {
                expect(Number.isFinite(r1.distance)).toBe(false);
                return;
            }
            expectClose(r0.distance, r1.distance, 1e-12, 1e-12);
            expectVectorClose(r0.closest[1], r1.closest[1], 1e-12, 1e-12);
        }, 200);
    });

    it('places the center at the smallest semi-axis', () => {
        // From the center the closest surface point lies along the shortest
        // axis, at distance min_i e[i].
        check(hyperellipsoidArb(3), h => {
            const r = query.compute(h.center, h);
            const emin = Math.min(...h.extent.values);
            expectClose(r.distance, emin, 1e-9, 1e-9);
        }, 120);
    });

    it('returns zero for points already on the surface', () => {
        check(fc.tuple(hyperellipsoidArb(3), wellScaledVector(3, -1, 1)),
            ([h, u]) => {
                let sum = 0;
                for (let i = 0; i < 3; ++i) {
                    sum += u.values[i] * u.values[i];
                }
                if (sum < 1e-4) {
                    return;
                }
                const s = 1 / Math.sqrt(sum);
                let q = h.center.clone();
                for (let i = 0; i < 3; ++i) {
                    q = add(q, mul(s * u.values[i] * h.extent.values[i],
                        h.axis[i]));
                }
                expect(query.compute(q, h).distance).toBeLessThanOrEqual(1e-8);
            }, 120);
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), hyperellipsoidArb(3)),
            ([p, h]) => {
                const p0 = p.clone();
                const c = h.center.clone();
                const e = h.extent.clone();
                const a = h.axis.map(x => x.clone());
                const r = query.compute(p, h);
                expect(p.values).toEqual(p0.values);
                expect(h.center.values).toEqual(c.values);
                expect(h.extent.values).toEqual(e.values);
                h.axis.forEach((x, i) =>
                    expect(x.values).toEqual(a[i].values));
                r.closest[0].values[0] = 555;
                expect(p.values).toEqual(p0.values);
            }, 120);
    });
    it('is correct at the exact center but degrades in a narrow band', () => {
        // Documented upstream defect: a query point whose smallest-extent
        // coordinate is exactly zero is handled by the special branch and is
        // exact, but a tiny nonzero coordinate loses the bisection bracket.
        // This test pins the correct behavior on the exact axis hyperplane
        // and records the failure mode so a future upstream fix is noticed.
        const h = axisAligned([4, 3, 2]);
        const onPlane = query.compute(v(1.5, 0.5, 0), h);
        expect(equationValue(onPlane.closest[1], h)).toBeCloseTo(1, 12);
        const nearPlane = query.compute(v(1.5, 0.5, 1e-13), h);
        // Upstream (and hence the port) leaves the surface here.
        expect(Math.abs(equationValue(nearPlane.closest[1], h) - 1))
            .toBeGreaterThan(1e-6);
        // And at 1e-16 the divisor underflows to zero.
        expect(query.compute(v(1.5, 0.5, 1e-16), h).distance)
            .toBe(Number.POSITIVE_INFINITY);
    });
});
