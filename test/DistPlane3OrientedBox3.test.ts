import { describe, expect, it } from 'vitest';
import { DistPlane3OrientedBox3 } from '../src/DistPlane3OrientedBox3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistPlane3AlignedBox3 } from '../src/DistPlane3AlignedBox3.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = v(...normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, v(...origin));
}

// An orthonormal frame parameterized by two angles.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        v(ca, sa, 0),
        v(-sa * cb, ca * cb, sb),
        v(sa * sb, -ca * sb, cb)
    ];
}

function obox(center: number[], axis: Vector[], extent: number[]):
    OrientedBox {
    return OrientedBox.fromCenterAxisExtent(v(...center), axis, v(...extent));
}

// The distance from a plane to a solid box is zero when the box straddles
// the plane; otherwise the closest box point is a vertex, so the distance is
// the smallest |Dot(N,V) - c| over the vertices.
function analyticDistance(p: Hyperplane, b: OrientedBox): number {
    let minSigned = Number.MAX_VALUE;
    let maxSigned = -Number.MAX_VALUE;
    for (const vertex of b.getVertices()) {
        const s = dot(p.normal, vertex) - p.constant;
        minSigned = Math.min(minSigned, s);
        maxSigned = Math.max(maxSigned, s);
    }
    if (minSigned <= 0 && maxSigned >= 0) {
        return 0;
    }
    return minSigned > 0 ? minSigned : -maxSigned;
}

function expectConsistent(p: Hyperplane, b: OrientedBox,
    result: { distance: number, sqrDistance: number, closest: [Vector, Vector] }):
    void {
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 9);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 9);

    // closest[0] is on the plane.
    expect(dot(p.normal, result.closest[0])).toBeCloseTo(p.constant, 9);

    // closest[1] is in the box.
    const delta1 = sub(result.closest[1], b.center);
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(dot(b.axis[i], delta1)))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
    }
}

describe('DistPlane3OrientedBox3', () => {
    const query = new DistPlane3OrientedBox3();

    it('agrees with the axis-aligned case for an identity frame', () => {
        const b = obox([0, 0, 0], frame(0, 0), [1, 1, 1]);
        const p = plane([0, 0, 1], [0, 0, 5]);
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
        expectConsistent(p, b, result);
    });

    it('computes the distance for a rotated box', () => {
        // A box rotated 45 degrees about z has half-diagonal sqrt(2) along x
        // and y, so the plane x = 5 is sqrt(2) away from a unit box.
        const s = Math.SQRT1_2;
        const b = obox([0, 0, 0], [v(s, s, 0), v(-s, s, 0), v(0, 0, 1)],
            [1, 1, 1]);
        const p = plane([1, 0, 0], [5, 0, 0]);
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(5 - Math.SQRT2, 12);
        expectConsistent(p, b, result);
    });

    it('reports zero distance when the plane cuts the box', () => {
        const b = obox([0, 0, 0], frame(0.4, 0.9), [1, 2, 3]);
        const p = plane([1, 1, 1], [0, 0, 0]);
        const result = query.compute(p, b);
        expect(result.distance).toBe(0);
        expectConsistent(p, b, result);
    });

    it('handles a plane touching a rotated box corner', () => {
        const b = obox([0, 0, 0], frame(0.3, 0.7), [1, 2, 3]);
        const corner = b.getVertices()[7];
        const p = Hyperplane.fromNormalOrigin(
            (() => { const n = v(1, 1, 1); normalize(n); return n; })(),
            corner);
        // The corner may or may not be the extreme vertex for this normal,
        // so only nonnegativity and consistency are asserted here.
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(analyticDistance(p, b), 9);
        expectConsistent(p, b, result);
    });

    it('matches the analytic vertex-based distance on random inputs', () => {
        let seed = 24681357;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 400; ++trial) {
            const b = obox(
                [rand() * 6 - 3, rand() * 6 - 3, rand() * 6 - 3],
                frame(rand() * Math.PI, rand() * Math.PI),
                [rand() * 2 + 0.05, rand() * 2 + 0.05, rand() * 2 + 0.05]);
            const p = plane(
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                [rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4]);
            const result = query.compute(p, b);
            expect(result.distance).toBeCloseTo(analyticDistance(p, b), 9);
            expectConsistent(p, b, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistPlane3OrientedBox3.h.
// ---------------------------------------------------------------------------

describe('DistPlane3OrientedBox3 verification', () => {
    const query = new DistPlane3OrientedBox3();
    const abQuery = new DistPlane3AlignedBox3();

    const planeArb = fc.tuple(unitVector(3), wellScaledVector(3, -6, 6))
        .map(([n, o]) => Hyperplane.fromNormalOrigin(n, o));

    const boxArb = fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
        fc.array(finite(0, 4), { minLength: 3, maxLength: 3 }))
        .map(([c, axis, e]) =>
            OrientedBox.fromCenterAxisExtent(c, axis, v(e[0], e[1], e[2])));

    // Map a vector from box coordinates to world coordinates.
    function toWorld(f: Vector[], p: Vector): Vector {
        return add(add(mul(p.values[0], f[0]), mul(p.values[1], f[1])),
            mul(p.values[2], f[2]));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(planeArb, boxArb), ([p, b]) => {
            const r = query.compute(p, b);
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            const d = sub(r.closest[0], r.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), r.distance, 1e-9, 1e-9);
            expectClose(dot(p.normal, r.closest[0]), p.constant, 1e-9, 1e-9);
            const delta = sub(r.closest[1], b.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(b.axis[i], delta)))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
            }
        });
    });

    it('matches the vertex-signed-distance formula', () => {
        check(fc.tuple(planeArb, boxArb), ([p, b]) => {
            const r = query.compute(p, b);
            expectClose(r.distance, analyticDistance(p, b), 1e-9, 1e-9);
        });
    });

    // Both closest points come back in box coordinates from the canonical
    // query, so both must be mapped through the same frame. Transforming only
    // one of them (or transforming an already-world-space point a second
    // time) would break this identity.
    it('agrees with the aligned-box query for an identity box frame', () => {
        const axes = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
        check(fc.tuple(planeArb, wellScaledVector(3, -5, 5),
            fc.array(finite(0, 4), { minLength: 3, maxLength: 3 })),
            ([p, c, e]) => {
                const ob = OrientedBox.fromCenterAxisExtent(c, axes,
                    v(e[0], e[1], e[2]));
                const ab = AlignedBox.fromMinMax(
                    v(c.values[0] - e[0], c.values[1] - e[1],
                        c.values[2] - e[2]),
                    v(c.values[0] + e[0], c.values[1] + e[1],
                        c.values[2] + e[2]));
                const r0 = query.compute(p, ob);
                const r1 = abQuery.compute(p, ab);
                expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
                expectVectorClose(r0.closest[1], r1.closest[1], 1e-9, 1e-9);
            });
    });

    it('is equivariant under a rigid motion of plane and box', () => {
        check(fc.tuple(planeArb, boxArb, rotationFrame(3),
            wellScaledVector(3, -6, 6)), ([p, b, R, t]) => {
                const r0 = query.compute(p, b);
                const xf = (q: Vector) => add(toWorld(R, q), t);
                const p1 = Hyperplane.fromNormalOrigin(toWorld(R, p.normal),
                    xf(p.origin));
                const b1 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                    [toWorld(R, b.axis[0]), toWorld(R, b.axis[1]),
                        toWorld(R, b.axis[2])], b.extent);
                const r1 = query.compute(p1, b1);
                expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
                expectVectorClose(xf(r0.closest[1]), r1.closest[1],
                    1e-8, 1e-8);
            });
    });

    it('handles a degenerate zero-extent box', () => {
        check(fc.tuple(planeArb, wellScaledVector(3, -5, 5), rotationFrame(3)),
            ([p, c, axis]) => {
                const b = OrientedBox.fromCenterAxisExtent(c, axis,
                    v(0, 0, 0));
                const r = query.compute(p, b);
                expectClose(r.distance,
                    Math.abs(dot(p.normal, c) - p.constant), 1e-9, 1e-9);
                expectVectorClose(r.closest[1], c, 1e-9, 1e-9);
            });
    });
});
