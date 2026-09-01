import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere';
import { Triangle } from '../src/Triangle';
import { Vector, add, dot, mul, sub } from '../src/Vector';
import { DistPointTriangle } from '../src/DistPointTriangle';
import {
    IntrSphere3Triangle3FI,
    IntrSphere3Triangle3FIResultType
} from '../src/IntrSphere3Triangle3';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function sphere(c: number[], r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(c), r);
}

function triangle(a: number[], b: number[], c: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(a), Vector.fromArray(b),
        Vector.fromArray(c));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The squared distance from the moving sphere center to the (moving) triangle
// at time t.
function sqrDistAtTime(sph: Hypersphere, sphVel: Vector, tri: Triangle,
    triVel: Vector, t: number): number {
    const c = add(sph.center, mul(t, sphVel));
    const moved = Triangle.fromVertices(
        add(tri.v[0], mul(t, triVel)),
        add(tri.v[1], mul(t, triVel)),
        add(tri.v[2], mul(t, triVel)));
    return new DistPointTriangle().compute(c, moved).sqrDistance;
}

describe('IntrSphere3Triangle3', () => {
    const fi = new IntrSphere3Triangle3FI();
    const T = IntrSphere3Triangle3FIResultType;
    // A triangle in the z = 0 plane.
    const tri = triangle([0, 0, 0], [4, 0, 0], [0, 4, 0]);
    const zero = vec(0, 0, 0);

    it('reports an initial overlap', () => {
        const s = sphere([1, 1, 0.5], 1);
        const result = fi.find(s, vec(0, 0, -1), tri, zero);
        expect(result.intersectionType).toBe(T.initiallyOverlapping);
        expect(result.contactTime).toBe(0);
        // The contact point is the triangle point closest to the center.
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[2]).toBeCloseTo(0, 12);
    });

    it('reports initial tangential contact as type +1 at time 0', () => {
        // The distance from the center to the triangle is exactly the radius.
        const s = sphere([1, 1, 1], 1);
        const result = fi.find(s, vec(0, 0, -1), tri, zero);
        expect(result.intersectionType).toBe(T.contact);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values[2]).toBeCloseTo(0, 12);
    });

    it('reports no contact for a separated, non-moving sphere', () => {
        const s = sphere([1, 1, 5], 1);
        const result = fi.find(s, zero, tri, zero);
        expect(result.intersectionType).toBe(T.noContact);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([0, 0, 0]);
    });

    it('finds the face contact for a sphere dropping onto the interior', () => {
        // The center starts 5 above the face and moves down at unit speed;
        // contact when the center is at height 1.
        const s = sphere([1, 1, 5], 1);
        const result = fi.find(s, vec(0, 0, -1), tri, zero);
        expect(result.intersectionType).toBe(T.contact);
        expect(result.contactTime).toBeCloseTo(4, 12);
        // Upstream reports the sphere CENTER at the time of contact, not the
        // point where the surfaces touch (see the PR notes).
        expect(result.contactPoint.values[2]).toBeCloseTo(1, 12);
        // The distance to the triangle at the reported time is the radius.
        expect(Math.sqrt(sqrDistAtTime(s, vec(0, 0, -1), tri, zero,
            result.contactTime))).toBeCloseTo(1, 9);
    });

    it('finds the same contact when the triangle moves instead', () => {
        const s = sphere([1, 1, 5], 1);
        const result = fi.find(s, zero, tri, vec(0, 0, 1));
        expect(result.intersectionType).toBe(T.contact);
        expect(result.contactTime).toBeCloseTo(4, 12);
        expect(Math.sqrt(sqrDistAtTime(s, zero, tri, vec(0, 0, 1),
            result.contactTime))).toBeCloseTo(1, 9);
    });

    it('finds the contact from the negative side of the plane', () => {
        const s = sphere([1, 1, -5], 1);
        const result = fi.find(s, vec(0, 0, 1), tri, zero);
        expect(result.intersectionType).toBe(T.contact);
        expect(result.contactTime).toBeCloseTo(4, 12);
    });

    it('reports no contact when the sphere moves away from the plane', () => {
        const s = sphere([1, 1, 5], 1);
        const result = fi.find(s, vec(0, 0, 1), tri, zero);
        expect(result.intersectionType).toBe(T.noContact);
    });

    it('finds a vertex (sphere-wedge) contact', () => {
        // The sphere approaches the vertex (0,0,0) along the -x axis from
        // outside the triangle; contact when the center is at distance 1.
        const s = sphere([-6, 0, 0], 1);
        const v = vec(1, 0, 0);
        const result = fi.find(s, v, tri, zero);
        expect(result.intersectionType).toBe(T.contact);
        expect(result.contactTime).toBeCloseTo(5, 12);
        expect(Math.sqrt(sqrDistAtTime(s, v, tri, zero, result.contactTime)))
            .toBeCloseTo(1, 9);
    });

    it('finds an edge (half-cylinder) contact', () => {
        // Approach the edge from (0,0,0) to (4,0,0) from -y, at x = 2.
        const s = sphere([2, -6, 0], 1);
        const v = vec(0, 1, 0);
        const result = fi.find(s, v, tri, zero);
        expect(result.intersectionType).toBe(T.contact);
        expect(result.contactTime).toBeCloseTo(5, 12);
        expect(Math.sqrt(sqrDistAtTime(s, v, tri, zero, result.contactTime)))
            .toBeCloseTo(1, 9);
    });

    it('reports no contact when the sphere passes wide of the triangle', () => {
        const s = sphere([-6, -6, 0], 0.5);
        const result = fi.find(s, vec(1, 0, 0), tri, zero);
        expect(result.intersectionType).toBe(T.noContact);
    });

    it('reports no contact when the relative velocity is zero', () => {
        const s = sphere([1, 1, 5], 1);
        const v = vec(0, 0, -1);
        const result = fi.find(s, v, tri, v);
        expect(result.intersectionType).toBe(T.noContact);
    });

    it('agrees with a dense time sampling on random configurations', () => {
        const rnd = makeRandom(9001);
        let earlyContact = 0;
        let missedContact = 0;
        let badDistance = 0;
        let contacts = 0;
        const tMax = 8;
        const steps = 2000;

        for (let trial = 0; trial < 150; ++trial) {
            const t3 = Triangle.fromVertices(
                vec(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1),
                vec(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1),
                vec(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1));
            const e0 = sub(t3.v[1], t3.v[0]);
            const e1 = sub(t3.v[2], t3.v[0]);
            const n = vec(
                e0.values[1] * e1.values[2] - e0.values[2] * e1.values[1],
                e0.values[2] * e1.values[0] - e0.values[0] * e1.values[2],
                e0.values[0] * e1.values[1] - e0.values[1] * e1.values[0]);
            if (Math.sqrt(dot(n, n)) < 0.5) {
                continue;  // near-degenerate triangle
            }

            const s = Hypersphere.fromCenterRadius(
                vec(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4),
                0.2 + 0.8 * rnd());
            const centroid = mul(1 / 3, add(t3.v[0], add(t3.v[1], t3.v[2])));
            // Aim roughly at the triangle, then jitter.
            const v = add(mul(0.25, sub(centroid, s.center)),
                vec(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5));
            const triVel = vec(0.2 * rnd() - 0.1, 0.2 * rnd() - 0.1,
                0.2 * rnd() - 0.1);
            if (dot(sub(v, triVel), sub(v, triVel)) < 1e-6) {
                continue;
            }

            const result = fi.find(s, v, t3, triVel);
            const rsqr = s.radius * s.radius;

            // Brute force: the earliest sampled time at which the distance
            // from the moving center to the moving triangle is <= radius.
            let sampledTime = -1;
            for (let k = 0; k <= steps; ++k) {
                const t = (tMax * k) / steps;
                if (sqrDistAtTime(s, v, t3, triVel, t) <= rsqr) {
                    sampledTime = t;
                    break;
                }
            }

            if (result.intersectionType === T.noContact) {
                // A missed contact is only acceptable if the sampled contact
                // is beyond the sampled window; a sampled hit well inside the
                // window means the query missed it.
                if (sampledTime >= 0 && sampledTime < tMax - 1e-9) {
                    ++missedContact;
                }
                continue;
            }

            ++contacts;
            const tbar = result.contactTime;
            if (tbar < -1e-9) {
                ++earlyContact;
                continue;
            }
            // At the reported time, the sphere touches the triangle (or, for
            // an initial overlap, is closer than the radius).
            const d = Math.sqrt(sqrDistAtTime(s, v, t3, triVel, tbar));
            if (result.intersectionType === T.initiallyOverlapping) {
                if (d > s.radius + 1e-9) {
                    ++badDistance;
                }
            }
            else if (Math.abs(d - s.radius) > 1e-6) {
                ++badDistance;
            }
            // The reported point is the sphere center at contact time.
            const expected = add(s.center, mul(tbar, v));
            const diff = sub(result.contactPoint, expected);
            if (result.intersectionType === T.contact && tbar > 0 &&
                Math.sqrt(dot(diff, diff)) > 1e-9) {
                ++badDistance;
            }
            // The contact must be no later than the first sampled contact
            // (within one sampling step).
            if (sampledTime >= 0 && tbar > sampledTime + tMax / steps + 1e-6) {
                ++earlyContact;
            }
        }

        expect(contacts).toBeGreaterThan(15);
        expect([earlyContact, missedContact, badDistance]).toEqual([0, 0, 0]);
    });
});
