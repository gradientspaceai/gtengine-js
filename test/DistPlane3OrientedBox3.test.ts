import { describe, expect, it } from 'vitest';
import { DistPlane3OrientedBox3 } from '../src/DistPlane3OrientedBox3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, normalize, sub } from '../src/Vector.js';

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
