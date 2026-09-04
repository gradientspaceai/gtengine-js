import { describe, it, expect } from 'vitest';
import {
    getContainerCapsule3,
    inContainerCapsule3,
    inContainerSphereCapsule3,
    inContainerCapsuleCapsule3,
    mergeContainersCapsule3
} from '../src/ContCapsule3.js';
import { Capsule, type Capsule3 } from '../src/Capsule.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Segment } from '../src/Segment.js';
import { DistPointSegment } from '../src/DistPointSegment.js';
import { Vector, add, length, mul, sub } from '../src/Vector.js';
import { Line } from '../src/Line.js';
import { DistPointLine } from '../src/DistPointLine.js';
import {
    check, expectClose, fc, rotationFrame, seededRandom, unitVector,
    wellScaledVector
} from './helpers/arbitraries.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function capsule(p0: Vector, p1: Vector, radius: number): Capsule3 {
    return Capsule.fromSegmentRadius(Segment.fromEndpoints(p0, p1), radius);
}

function distanceToSegment(point: Vector, seg: Segment): number {
    return new DistPointSegment().compute(point, seg).distance;
}

function makeRandom(seed: number): () => number {
    let s = seed;
    return () => {
        s = (s * 1103515245 + 12345) % 2147483648;
        return s / 2147483648;
    };
}

describe('getContainerCapsule3', () => {
    it('contains every input point (random clouds)', () => {
        const rand = makeRandom(2024);
        for (let trial = 0; trial < 25; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 40; ++i) {
                points.push(v(
                    12 * (rand() - 0.5),
                    2 * (rand() - 0.5),
                    2 * (rand() - 0.5)));
            }
            const c = getContainerCapsule3(points);
            expect(c.dimension).toBe(3);
            for (const p of points) {
                expect(distanceToSegment(p, c.segment))
                    .toBeLessThanOrEqual(c.radius + 1e-9);
            }
        }
    });

    it('recovers a capsule-like axis for points on a cylinder wall', () => {
        const points: Vector[] = [];
        for (const x of [-3, 3]) {
            for (let k = 0; k < 16; ++k) {
                const a = (2 * Math.PI * k) / 16;
                points.push(v(x, Math.cos(a), Math.sin(a)));
            }
        }
        const c = getContainerCapsule3(points);
        expect(Math.abs(c.segment.p[0].values[0])).toBeCloseTo(3, 8);
        expect(Math.abs(c.segment.p[1].values[0])).toBeCloseTo(3, 8);
        expect(c.radius).toBeCloseTo(1, 10);
        // The caps are pulled in as far as possible; the extreme rings lie on
        // the hemispherical caps, so the segment endpoints reach the rings.
        for (const p of points) {
            expect(inContainerCapsule3(p, c)).toBe(true);
        }
    });

    it('degenerates to a sphere for points on a sphere', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            for (let j = 1; j < 8; ++j) {
                const theta = (2 * Math.PI * i) / 8;
                const phi = (Math.PI * j) / 8;
                points.push(v(
                    Math.sin(phi) * Math.cos(theta),
                    Math.sin(phi) * Math.sin(theta),
                    Math.cos(phi)));
            }
        }
        const c = getContainerCapsule3(points);
        expect(c.radius).toBeCloseTo(1, 8);
        const cf = c.segment.getCenteredForm();
        expect(cf.extent).toBeLessThan(1e-8);
    });

    it('handles coincident points (degenerate capsule)', () => {
        const p = v(1, -2, 4);
        const c = getContainerCapsule3([p, p.clone(), p.clone()]);
        expect(c.radius).toBeCloseTo(0, 12);
        expect(c.segment.p[0].values).toEqual([1, -2, 4]);
        expect(c.segment.p[1].values).toEqual([1, -2, 4]);
    });

    it('throws for an empty point set and for non-3D points', () => {
        expect(() => getContainerCapsule3([])).toThrow();
        expect(() => getContainerCapsule3([Vector.fromArray([1, 2])]))
            .toThrow();
    });
});

describe('inContainerCapsule3', () => {
    const c = capsule(v(-2, 0, 0), v(2, 0, 0), 1);

    it('accepts interior and boundary points', () => {
        expect(inContainerCapsule3(v(0, 0, 0), c)).toBe(true);
        expect(inContainerCapsule3(v(0, 1, 0), c)).toBe(true);
        expect(inContainerCapsule3(v(3, 0, 0), c)).toBe(true);
        expect(inContainerCapsule3(v(-3, 0, 0), c)).toBe(true);
    });

    it('rejects exterior points', () => {
        expect(inContainerCapsule3(v(0, 1.0001, 0), c)).toBe(false);
        expect(inContainerCapsule3(v(3.0001, 0, 0), c)).toBe(false);
        expect(inContainerCapsule3(v(2, 1, 1), c)).toBe(false);
    });
});

describe('inContainerSphereCapsule3', () => {
    const c = capsule(v(-2, 0, 0), v(2, 0, 0), 2);

    it('accepts a contained sphere', () => {
        const s = Hypersphere.fromCenterRadius(v(0, 0, 0), 1);
        expect(inContainerSphereCapsule3(s, c)).toBe(true);
    });

    it('accepts a tangent sphere', () => {
        const s = Hypersphere.fromCenterRadius(v(0, 1, 0), 1);
        expect(inContainerSphereCapsule3(s, c)).toBe(true);
    });

    it('rejects a sphere poking out', () => {
        const s = Hypersphere.fromCenterRadius(v(0, 1.5, 0), 1);
        expect(inContainerSphereCapsule3(s, c)).toBe(false);
    });

    it('rejects a sphere larger than the capsule radius', () => {
        const s = Hypersphere.fromCenterRadius(v(0, 0, 0), 3);
        expect(inContainerSphereCapsule3(s, c)).toBe(false);
    });
});

describe('inContainerCapsuleCapsule3', () => {
    const outer = capsule(v(-4, 0, 0), v(4, 0, 0), 2);

    it('accepts a capsule inside another', () => {
        const inner = capsule(v(-1, 0, 0), v(1, 0, 0), 1);
        expect(inContainerCapsuleCapsule3(inner, outer)).toBe(true);
    });

    it('accepts a capsule equal to the container', () => {
        expect(inContainerCapsuleCapsule3(outer, outer)).toBe(true);
    });

    it('rejects a capsule that pokes out at one end', () => {
        const inner = capsule(v(-1, 0, 0), v(6, 0, 0), 1);
        expect(inContainerCapsuleCapsule3(inner, outer)).toBe(false);
    });

    it('rejects a fatter capsule', () => {
        const inner = capsule(v(-1, 0, 0), v(1, 0, 0), 3);
        expect(inContainerCapsuleCapsule3(inner, outer)).toBe(false);
    });
});

describe('mergeContainersCapsule3', () => {
    it('returns the container when one capsule contains the other', () => {
        const big = capsule(v(-4, 0, 0), v(4, 0, 0), 2);
        const small = capsule(v(-1, 0, 0), v(1, 0, 0), 1);
        const m0 = mergeContainersCapsule3(small, big);
        expect(m0.radius).toBeCloseTo(2, 12);
        expect(m0.equals(big)).toBe(true);

        const m1 = mergeContainersCapsule3(big, small);
        expect(m1.equals(big)).toBe(true);
    });

    it('returns a copy, not an alias, of the containing capsule', () => {
        const big = capsule(v(-4, 0, 0), v(4, 0, 0), 2);
        const small = capsule(v(-1, 0, 0), v(1, 0, 0), 1);
        const merge = mergeContainersCapsule3(small, big);
        merge.radius = 99;
        expect(big.radius).toBe(2);
    });

    it('contains both collinear input capsules', () => {
        const c0 = capsule(v(-4, 0, 0), v(0, 0, 0), 1);
        const c1 = capsule(v(0, 0, 0), v(4, 0, 0), 1);
        const merge = mergeContainersCapsule3(c0, c1);
        expect(inContainerCapsuleCapsule3(c0, merge)).toBe(true);
        expect(inContainerCapsuleCapsule3(c1, merge)).toBe(true);
        expect(merge.radius).toBeCloseTo(1, 10);
    });

    it('contains both inputs for random capsule pairs', () => {
        const rand = makeRandom(555);
        const randomCapsule = (): Capsule3 => {
            const p0 = v(6 * (rand() - 0.5), 6 * (rand() - 0.5), 6 * (rand() - 0.5));
            const p1 = v(6 * (rand() - 0.5), 6 * (rand() - 0.5), 6 * (rand() - 0.5));
            return capsule(p0, p1, 0.25 + rand());
        };

        for (let trial = 0; trial < 50; ++trial) {
            const c0 = randomCapsule();
            const c1 = randomCapsule();
            const merge = mergeContainersCapsule3(c0, c1);
            for (const c of [c0, c1]) {
                for (const end of c.segment.p) {
                    expect(distanceToSegment(end, merge.segment) + c.radius)
                        .toBeLessThanOrEqual(merge.radius + 1e-8);
                }
            }
        }
    });

    it('is symmetric in radius for a symmetric configuration', () => {
        const c0 = capsule(v(-3, 1, 0), v(3, 1, 0), 1);
        const c1 = capsule(v(-3, -1, 0), v(3, -1, 0), 1);
        const m01 = mergeContainersCapsule3(c0, c1);
        const m10 = mergeContainersCapsule3(c1, c0);
        expect(m01.radius).toBeCloseTo(m10.radius, 12);
        expect(m01.radius).toBeCloseTo(2, 12);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContCapsule3.h semantics.
// ---------------------------------------------------------------------------

describe('ContCapsule3 verification', () => {
    // Elongated clouds: t * axis + small transverse noise. The principal
    // direction is then well separated from the other two, so the
    // least-squares line (an eigenvector problem) is well conditioned and
    // equivariance properties hold to a tight tolerance.
    // The parameters t are spread evenly over [-5, 5] rather than drawn at
    // random, so the variance along the axis (about 8) always dominates the
    // transverse variance (at most 0.09) and the principal direction is never
    // ambiguous. Random t values can all collapse to zero, which makes the
    // eigenvector arbitrary and every equivariance property meaningless.
    const elongatedCloud = fc.tuple(unitVector(3), wellScaledVector(3, -4, 4),
        fc.array(wellScaledVector(3, -0.3, 0.3),
            { minLength: 4, maxLength: 12 }))
        .map(([axis, origin, noises]) =>
            noises.map((noise, i) => {
                const t = -5 + (10 * i) / (noises.length - 1);
                return add(add(origin, mul(t, axis)), noise);
            }));

    const rigid = (frame: Vector[], t: Vector) => (p: Vector): Vector =>
        add(add(add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])),
            mul(p.get(2), frame[2])), t);

    // The design claim: every input point is inside the fitted capsule.
    it('the fitted capsule contains every input point', () => {
        check(elongatedCloud, (points: Vector[]) => {
            const capsule = getContainerCapsule3(points);
            for (const p of points) {
                expect(distanceToSegment(p, capsule.segment))
                    .toBeLessThanOrEqual(capsule.radius + 1e-9);
            }
        });
    });

    // The radius is the largest distance from the points to the *infinite*
    // fitted axis, not to the capsule segment; cross-check with DistPointLine
    // on the line through the capsule segment.
    it('the radius is the largest distance to the fitted axis line', () => {
        check(elongatedCloud, (points: Vector[]) => {
            const capsule = getContainerCapsule3(points);
            const cf = capsule.segment.getCenteredForm();
            if (cf.extent === 0) {
                // Degenerate cloud (all points coincident): the segment is a
                // single point and there is no axis direction to test.
                return;
            }
            const axisLine = Line.fromOriginDirection(cf.center, cf.direction);
            let maxDist = 0;
            for (const p of points) {
                maxDist = Math.max(maxDist,
                    new DistPointLine().compute(p, axisLine).distance);
            }
            expectClose(capsule.radius, maxDist, 1e-9, 1e-9);
        });
    });

    // Rigid motions: the capsule of the transformed cloud is the transform of
    // the capsule. The fit is an eigen decomposition, so the axis direction
    // may come back negated and the two segment endpoints swapped.
    it('is equivariant under rigid motions', () => {
        check(fc.tuple(elongatedCloud, rotationFrame(3), wellScaledVector(3)),
            ([points, frame, t]: [Vector[], Vector[], Vector]) => {
                const xform = rigid(frame, t);
                const c0 = getContainerCapsule3(points);
                const c1 = getContainerCapsule3(points.map(xform));
                expectClose(c1.radius, c0.radius, 1e-8, 1e-8);
                const want = [xform(c0.segment.p[0]), xform(c0.segment.p[1])];
                const got = [c1.segment.p[0], c1.segment.p[1]];
                const dist = (a: Vector, b: Vector): number => length(sub(a, b));
                const same = Math.max(dist(want[0], got[0]), dist(want[1], got[1]));
                const swapped = Math.max(dist(want[0], got[1]), dist(want[1], got[0]));
                expect(Math.min(same, swapped)).toBeLessThanOrEqual(1e-7);
            });
    });

    // inContainerSphereCapsule3 must agree with sampling the sphere surface
    // against inContainerCapsule3 (away from tangency, where rounding rules).
    it('sphere containment agrees with sampling the sphere', () => {
        const rand = seededRandom(0x5eedca);
        for (let trial = 0; trial < 300; ++trial) {
            const p0 = Vector.fromArray(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3]);
            const p1 = Vector.fromArray(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3]);
            const cap = capsule(p0, p1, 0.5 + 2 * rand());
            const sphere = Hypersphere.fromCenterRadius(
                Vector.fromArray([6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3]),
                2 * rand());
            const inside = inContainerSphereCapsule3(sphere, cap);
            // Slack of the analytic test; skip near-tangent configurations.
            const slack = (cap.radius - sphere.radius)
                - distanceToSegment(sphere.center, cap.segment);
            if (Math.abs(slack) < 1e-6) {
                continue;
            }
            expect(inside).toBe(slack > 0);
            if (inside) {
                for (let k = 0; k < 40; ++k) {
                    const z = 2 * rand() - 1;
                    const a = 2 * Math.PI * rand();
                    const s = Math.sqrt(Math.max(0, 1 - z * z));
                    const q = add(sphere.center, mul(sphere.radius,
                        Vector.fromArray(
                            [s * Math.cos(a), s * Math.sin(a), z])));
                    expect(inContainerCapsule3(q, cap)).toBe(true);
                }
            }
        }
    }, 30000);

    // A capsule is inside another exactly when the spheres at its two
    // endpoints are; sample the contained capsule to confirm.
    it('capsule containment agrees with sampling the contained capsule', () => {
        const rand = seededRandom(0x5eedcb);
        for (let trial = 0; trial < 300; ++trial) {
            const mk = (scale: number): Capsule3 => capsule(
                Vector.fromArray([scale * (2 * rand() - 1), scale * (2 * rand() - 1),
                    scale * (2 * rand() - 1)]),
                Vector.fromArray([scale * (2 * rand() - 1), scale * (2 * rand() - 1),
                    scale * (2 * rand() - 1)]),
                0.2 + scale * rand());
            const inner = mk(1);
            const outer = mk(3);
            if (!inContainerCapsuleCapsule3(inner, outer)) {
                continue;
            }
            // Sample the inner capsule's surface: p(s) + r * unit.
            for (let k = 0; k < 40; ++k) {
                const s = rand();
                const axisPoint = add(mul(1 - s, inner.segment.p[0]),
                    mul(s, inner.segment.p[1]));
                const z = 2 * rand() - 1;
                const a = 2 * Math.PI * rand();
                const t = Math.sqrt(Math.max(0, 1 - z * z));
                const q = add(axisPoint, mul(inner.radius,
                    Vector.fromArray([t * Math.cos(a), t * Math.sin(a), z])));
                expect(distanceToSegment(q, outer.segment))
                    .toBeLessThanOrEqual(outer.radius + 1e-9);
            }
        }
    }, 30000);

    // MergeContainers claims a capsule containing both inputs (though not
    // necessarily of least volume). Sample both input capsules and confirm.
    it('merge contains sampled surface points of both inputs', () => {
        const rand = seededRandom(0x5eedcc);
        for (let trial = 0; trial < 200; ++trial) {
            const mk = (): Capsule3 => capsule(
                Vector.fromArray([6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3]),
                Vector.fromArray([6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3]),
                0.2 + 1.5 * rand());
            const c0 = mk(), c1 = mk();
            const merge = mergeContainersCapsule3(c0, c1);
            for (const input of [c0, c1]) {
                for (let k = 0; k < 30; ++k) {
                    const s = rand();
                    const axisPoint = add(mul(1 - s, input.segment.p[0]),
                        mul(s, input.segment.p[1]));
                    const z = 2 * rand() - 1;
                    const a = 2 * Math.PI * rand();
                    const t = Math.sqrt(Math.max(0, 1 - z * z));
                    const q = add(axisPoint, mul(input.radius,
                        Vector.fromArray([t * Math.cos(a), t * Math.sin(a), z])));
                    expect(distanceToSegment(q, merge.segment))
                        .toBeLessThanOrEqual(merge.radius + 1e-9);
                }
            }
        }
    }, 30000);

    // Degenerate clouds: all points equal gives a zero-radius, zero-extent
    // capsule centered on the point.
    it('collapses to a point for coincident inputs', () => {
        check(wellScaledVector(3), (p: Vector) => {
            const cap = getContainerCapsule3([p, p.clone(), p.clone()]);
            expect(cap.radius).toBeLessThanOrEqual(1e-12);
            expectClose(length(sub(cap.segment.p[0], p)), 0, 1e-12, 1e-12);
            expectClose(length(sub(cap.segment.p[1], p)), 0, 1e-12, 1e-12);
        });
    });
});
