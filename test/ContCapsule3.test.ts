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
import { Vector } from '../src/Vector.js';

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
