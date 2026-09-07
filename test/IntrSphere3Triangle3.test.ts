import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import {
    IntrSphere3Triangle3FI,
    IntrSphere3Triangle3FIResultType
} from '../src/IntrSphere3Triangle3.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrSphere3Triangle3.h.
// ---------------------------------------------------------------------------

describe('IntrSphere3Triangle3 verification', () => {
    const fi = new IntrSphere3Triangle3FI();
    const T = IntrSphere3Triangle3FIResultType;
    const dpt = new DistPointTriangle();

    const arbTriangle = fc.tuple(wellScaledVector(3, -3, 3),
        wellScaledVector(3, -3, 3), wellScaledVector(3, -3, 3))
        .filter(([a, b, c]) => {
            const e0 = sub(b, a), e1 = sub(c, a);
            const n = cross(e0, e1);
            return dot(n, n) > 0.25;
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));
    const arbSphere = fc.tuple(wellScaledVector(3, -6, 6), positive(2, 0.25))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

    // The distance from the moving sphere center to the moving triangle at
    // the given time.
    function gap(sph: Hypersphere, sv: Vector, tri: Triangle, tv: Vector,
        time: number): number {
        const c = add(sph.center, mul(time, sv));
        const moved = Triangle.fromVertices(
            add(tri.v[0], mul(time, tv)),
            add(tri.v[1], mul(time, tv)),
            add(tri.v[2], mul(time, tv)));
        return dpt.compute(c, moved).distance - sph.radius;
    }

    // A configuration whose sphere is aimed at the triangle centroid, so that
    // a useful fraction of the draws produce contacts.
    const arbMoving = fc.tuple(arbSphere, arbTriangle, unitVector(3),
        positive(2, 0.2))
        .map(([sph, tri, jitter, speed]) => {
            const centroid = mul(1 / 3,
                add(tri.v[0], add(tri.v[1], tri.v[2])));
            const aim = add(sub(centroid, sph.center), mul(1.5, jitter));
            const len = Math.sqrt(dot(aim, aim));
            const sv = len > 1e-6 ? mul(speed / len, aim) : mul(speed, jitter);
            return { sph, tri, sv };
        });

    it('reports an initial overlap when the sphere already meets the triangle',
        () => {
            check(fc.tuple(arbSphere, arbTriangle, wellScaledVector(3, -2, 2)),
                ([sph, tri, sv]) => {
                    const ptRes = dpt.compute(sph.center, tri);
                    const rsqr = sph.radius * sph.radius;
                    if (ptRes.sqrDistance > rsqr) {
                        return;
                    }
                    const res = fi.find(sph, sv, tri, Vector.zero(3));
                    expect(res.contactTime).toBe(0);
                    expect(res.intersectionType).toBe(
                        ptRes.sqrDistance < rsqr ? T.initiallyOverlapping
                            : T.contact);
                    // The contact point is the triangle point closest to the
                    // sphere center.
                    expectVectorClose(res.contactPoint, ptRes.closest[1],
                        1e-12, 1e-12);
                });
        });

    it('reports no contact for separated objects with zero relative velocity',
        () => {
            check(fc.tuple(arbSphere, arbTriangle, wellScaledVector(3, -2, 2)),
                ([sph, tri, v]) => {
                    if (dpt.compute(sph.center, tri).sqrDistance
                        <= sph.radius * sph.radius) {
                        return;
                    }
                    // Equal velocities: the relative velocity is exactly zero.
                    const res = fi.find(sph, v, tri, v.clone());
                    expect(res.intersectionType).toBe(T.noContact);
                    expect(res.contactTime).toBe(0);
                    expect(res.contactPoint.equals(Vector.zero(3))).toBe(true);
                });
        });

    it('a reported contact time is a tangency and no earlier time touches',
        () => {
            check(arbMoving, ({ sph, tri, sv }) => {
                const res = fi.find(sph, sv, tri, Vector.zero(3));
                if (res.intersectionType !== T.contact
                    || res.contactTime === 0) {
                    return;
                }
                const t = res.contactTime;
                expect(Number.isFinite(t)).toBe(true);
                expect(t).toBeGreaterThan(0);
                // At the contact time the sphere touches the triangle.
                expectClose(gap(sph, sv, tri, Vector.zero(3), t), 0,
                    1e-6, 1e-6);
                // No strictly earlier time has the sphere overlapping the
                // triangle.
                for (let k = 1; k < 64; ++k) {
                    const earlier = (t * k) / 64;
                    expect(gap(sph, sv, tri, Vector.zero(3), earlier))
                        .toBeGreaterThan(-1e-7);
                }
                // The reported contact point is the position of the sphere
                // center at the contact time (an upstream quirk: for the
                // initial-overlap case the contact point is a triangle point
                // instead).
                expectVectorClose(res.contactPoint,
                    add(sph.center, mul(t, sv)), 1e-12, 1e-12);
            });
        });

    it('a reported no-contact is confirmed by sampling the motion', () => {
        check(arbMoving, ({ sph, tri, sv }) => {
            const res = fi.find(sph, sv, tri, Vector.zero(3));
            if (res.intersectionType !== T.noContact) {
                return;
            }
            for (let k = 0; k <= 400; ++k) {
                const t = (20 * k) / 400;
                expect(gap(sph, sv, tri, Vector.zero(3), t))
                    .toBeGreaterThan(-1e-7);
            }
        }, 60);
    });

    it('the contact time depends only on the relative velocity', () => {
        check(fc.tuple(arbMoving, wellScaledVector(3, -2, 2)),
            ([{ sph, tri, sv }, drift]) => {
                const a = fi.find(sph, sv, tri, Vector.zero(3));
                const b = fi.find(sph, add(sv, drift), tri, drift);
                expect(b.intersectionType).toBe(a.intersectionType);
                expectClose(b.contactTime, a.contactTime, 1e-9, 1e-9);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(arbMoving, rotationFrame(3), wellScaledVector(3, -3, 3)),
            ([{ sph, tri, sv }, frame, shift]) => {
                const rot = (v: Vector): Vector => Vector.fromArray([
                    dot(frame[0], v), dot(frame[1], v), dot(frame[2], v)]);
                const map = (v: Vector): Vector => add(rot(v), shift);
                const a = fi.find(sph, sv, tri, Vector.zero(3));
                const b = fi.find(
                    Hypersphere.fromCenterRadius(map(sph.center), sph.radius),
                    rot(sv),
                    Triangle.fromVertices(map(tri.v[0]), map(tri.v[1]),
                        map(tri.v[2])),
                    Vector.zero(3));
                expect(b.intersectionType).toBe(a.intersectionType);
                if (a.intersectionType !== T.noContact) {
                    expectClose(b.contactTime, a.contactTime, 1e-7, 1e-7);
                    expectVectorClose(b.contactPoint, map(a.contactPoint),
                        1e-6, 1e-6);
                }
            });
    });

    it('reports no NaN when a contact is found', () => {
        check(arbMoving, ({ sph, tri, sv }) => {
            const res = fi.find(sph, sv, tri, Vector.zero(3));
            expect(Number.isNaN(res.contactTime)).toBe(false);
            for (let i = 0; i < 3; ++i) {
                expect(Number.isNaN(res.contactPoint.get(i))).toBe(false);
            }
        });
    });

    it('uses the upstream zero-divisor semantics for a degenerate triangle'
        + ' edge', () => {
        // The "triangle" degenerates to the segment from (0,0,0) to (1,0,0),
        // so E[2] = v0 - v2 is the zero vector and sqrLenE[2] is zero.
        // Upstream computes E[2] * nu[2] / sqrLenE[2] with the Vector
        // operator/, which yields the ZERO vector for a zero divisor, so the
        // first half-cylinder test degenerates into the sphere test at the
        // vertex v0 and fires with tbar = 5 - sqrt(3)/2. Dividing the scalars
        // first would give NaN, all comparisons would be false, and the query
        // would fall through to the (0,1)-cylinder and report 4 instead.
        const tri = triangle([0, 0, 0], [1, 0, 0], [0, 0, 0]);
        const sph = sphere([0.5, 5, 0], 1);
        const res = fi.find(sph, vec(0, -1, 0), tri, vec(0, 0, 0));
        expect(res.intersectionType).toBe(T.contact);
        expectClose(res.contactTime, 5 - Math.sqrt(0.75));
    });
});
