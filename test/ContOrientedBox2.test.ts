import { describe, it, expect } from 'vitest';
import {
    getContainerOrientedBox2,
    inContainerOrientedBox2,
    mergeContainersOrientedBox2
} from '../src/ContOrientedBox2.js';
import { OrientedBox, type OrientedBox2 } from '../src/OrientedBox.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { ApprQuery } from '../src/ApprQuery.js';

function v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function box2(cx: number, cy: number, angle: number, e0: number,
    e1: number): OrientedBox2 {
    const c = Math.cos(angle), s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(
        v(cx, cy), [v(c, s), v(-s, c)], v(e0, e1));
}

// Independent containment check with a tolerance.
function containedWithin(point: Vector, box: OrientedBox2,
    tol: number): boolean {
    const diff = sub(point, box.center);
    for (let i = 0; i < 2; ++i) {
        if (Math.abs(dot(diff, box.axis[i])) > box.extent.values[i] + tol) {
            return false;
        }
    }
    return true;
}

function makeRandom(seed: number): () => number {
    let s = seed;
    return () => {
        s = (s * 1103515245 + 12345) % 2147483648;
        return s / 2147483648;
    };
}

describe('getContainerOrientedBox2', () => {
    it('recovers the exact extents of an axis-aligned rectangle', () => {
        // A 6 x 2 rectangle centered at (1,-1); include interior points so
        // the covariance is well conditioned.
        const points: Vector[] = [];
        for (let i = 0; i <= 6; ++i) {
            for (let j = 0; j <= 2; ++j) {
                points.push(v(1 - 3 + i, -1 - 1 + j));
            }
        }

        const box = getContainerOrientedBox2(points);
        expect(box).not.toBeNull();
        const b = box as OrientedBox2;
        expect(b.center.values[0]).toBeCloseTo(1, 10);
        expect(b.center.values[1]).toBeCloseTo(-1, 10);

        // ApprGaussian2 returns eigenvalues in increasing order, so axis[0]
        // is the low-variance direction (y here) and axis[1] the x direction.
        expect(Math.abs(b.axis[0].values[1])).toBeCloseTo(1, 10);
        expect(Math.abs(b.axis[1].values[0])).toBeCloseTo(1, 10);
        expect(b.extent.values[0]).toBeCloseTo(1, 10);
        expect(b.extent.values[1]).toBeCloseTo(3, 10);

        for (const p of points) {
            expect(containedWithin(p, b, 1e-10)).toBe(true);
        }
    });

    it('recovers a rotated rectangle', () => {
        const angle = 0.7;
        const c = Math.cos(angle), s = Math.sin(angle);
        const points: Vector[] = [];
        for (let i = -4; i <= 4; ++i) {
            for (let j = -1; j <= 1; ++j) {
                const x = i * 0.5, y = j * 0.5;
                points.push(v(c * x - s * y + 3, s * x + c * y - 2));
            }
        }

        const box = getContainerOrientedBox2(points) as OrientedBox2;
        expect(box.center.values[0]).toBeCloseTo(3, 10);
        expect(box.center.values[1]).toBeCloseTo(-2, 10);
        // Extents in increasing eigenvalue order: 0.5 across, 2 along.
        expect(box.extent.values[0]).toBeCloseTo(0.5, 10);
        expect(box.extent.values[1]).toBeCloseTo(2, 10);
        // The long axis is the rotated x direction.
        expect(Math.abs(dot(box.axis[1], v(c, s)))).toBeCloseTo(1, 10);
    });

    it('contains every input point (random clouds)', () => {
        const rand = makeRandom(4242);
        for (let trial = 0; trial < 25; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 30; ++i) {
                points.push(v(8 * (rand() - 0.5), 3 * (rand() - 0.5)));
            }
            const box = getContainerOrientedBox2(points) as OrientedBox2;
            expect(box).not.toBeNull();
            for (const p of points) {
                expect(containedWithin(p, box, 1e-10)).toBe(true);
            }
        }
    });

    it('handles collinear points (one zero extent)', () => {
        const points = [v(0, 0), v(1, 1), v(2, 2), v(3, 3)];
        const box = getContainerOrientedBox2(points) as OrientedBox2;
        expect(box.center.values[0]).toBeCloseTo(1.5, 10);
        expect(box.center.values[1]).toBeCloseTo(1.5, 10);
        expect(box.extent.values[0]).toBeCloseTo(0, 10);
        expect(box.extent.values[1]).toBeCloseTo(1.5 * Math.SQRT2, 10);
    });

    it('handles coincident points (degenerate box)', () => {
        const points = [v(2, 5), v(2, 5), v(2, 5)];
        const box = getContainerOrientedBox2(points) as OrientedBox2;
        expect(box.center.values[0]).toBeCloseTo(2, 12);
        expect(box.center.values[1]).toBeCloseTo(5, 12);
        expect(box.extent.values[0]).toBeCloseTo(0, 12);
        expect(box.extent.values[1]).toBeCloseTo(0, 12);
    });

    it('throws for an empty point set', () => {
        expect(() => getContainerOrientedBox2([])).toThrow();
    });

    it('returns null when the Gaussian fit reports failure', () => {
        // Index validation is off by default, so force a failure by turning
        // it on and passing a point count below the fitter minimum.
        ApprQuery.validateIndices = true;
        try {
            expect(getContainerOrientedBox2([v(1, 2)])).toBeNull();
        } finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('throws for non-2D points', () => {
        expect(() => getContainerOrientedBox2([Vector.fromArray([1, 2, 3])]))
            .toThrow();
    });

    it('does not alias the fitter parameters', () => {
        const points = [v(0, 0), v(2, 0), v(2, 1), v(0, 1)];
        const box = getContainerOrientedBox2(points) as OrientedBox2;
        const again = getContainerOrientedBox2(points) as OrientedBox2;
        expect(box).not.toBe(again);
        expect(box.center.values).toEqual(again.center.values);
    });
});

describe('inContainerOrientedBox2', () => {
    const box = box2(1, 2, Math.PI / 4, 2, 1);

    it('accepts the center', () => {
        expect(inContainerOrientedBox2(box.center, box)).toBe(true);
    });

    it('accepts the vertices of an axis-aligned box (boundary included)', () => {
        // The upstream test is |Dot(diff, axis[i])| > extent[i], which is
        // exact only when the projections are exact; use an axis-aligned box
        // so that the vertex projections are exactly the extents.
        const aligned = box2(1, 2, 0, 2, 1);
        for (const vertex of aligned.getVertices()) {
            expect(inContainerOrientedBox2(vertex, aligned)).toBe(true);
        }
    });

    it('accepts the vertices of a rotated box up to rounding', () => {
        for (const vertex of box.getVertices()) {
            expect(containedWithin(vertex, box, 1e-12)).toBe(true);
        }
    });

    it('rejects points just outside each face', () => {
        const s = Math.SQRT1_2;
        // Along axis[0] = (s,s) at distance 2.001 from the center.
        expect(inContainerOrientedBox2(
            v(1 + 2.001 * s, 2 + 2.001 * s), box)).toBe(false);
        // Along axis[1] = (-s,s) at distance 1.001 from the center.
        expect(inContainerOrientedBox2(
            v(1 - 1.001 * s, 2 + 1.001 * s), box)).toBe(false);
    });

    it('throws for mismatched dimensions', () => {
        expect(() => inContainerOrientedBox2(Vector.fromArray([1, 2, 3]), box))
            .toThrow();
    });
});

describe('mergeContainersOrientedBox2', () => {
    it('reproduces a box merged with itself', () => {
        const box = box2(1, -2, 0.3, 2, 1);
        const merge = mergeContainersOrientedBox2(box, box);
        expect(merge.center.values[0]).toBeCloseTo(1, 10);
        expect(merge.center.values[1]).toBeCloseTo(-2, 10);
        const extents = [merge.extent.values[0], merge.extent.values[1]]
            .sort((a, b) => a - b);
        expect(extents[0]).toBeCloseTo(1, 10);
        expect(extents[1]).toBeCloseTo(2, 10);
    });

    it('merges two axis-aligned boxes into their bounding box', () => {
        const b0 = box2(0, 0, 0, 1, 1);
        const b1 = box2(4, 0, 0, 1, 1);
        const merge = mergeContainersOrientedBox2(b0, b1);
        expect(merge.center.values[0]).toBeCloseTo(2, 10);
        expect(merge.center.values[1]).toBeCloseTo(0, 10);
        expect(merge.extent.values[0]).toBeCloseTo(3, 10);
        expect(merge.extent.values[1]).toBeCloseTo(1, 10);
    });

    it('contains all vertices of both inputs (random pairs)', () => {
        const rand = makeRandom(777);
        for (let trial = 0; trial < 50; ++trial) {
            const b0 = box2(6 * (rand() - 0.5), 6 * (rand() - 0.5),
                2 * Math.PI * rand(), 0.2 + 2 * rand(), 0.2 + 2 * rand());
            const b1 = box2(6 * (rand() - 0.5), 6 * (rand() - 0.5),
                2 * Math.PI * rand(), 0.2 + 2 * rand(), 0.2 + 2 * rand());
            const merge = mergeContainersOrientedBox2(b0, b1);
            for (const b of [b0, b1]) {
                for (const vertex of b.getVertices()) {
                    expect(containedWithin(vertex, merge, 1e-9)).toBe(true);
                }
            }
        }
    });

    it('handles anti-aligned axes by negating the second box axis', () => {
        const b0 = box2(0, 0, 0, 2, 1);
        // Rotated by pi, so its axis[0] = (-1,0) points opposite to b0's.
        const b1 = box2(0, 0, Math.PI, 2, 1);
        const merge = mergeContainersOrientedBox2(b0, b1);
        expect(Math.abs(merge.axis[0].values[0])).toBeCloseTo(1, 10);
        expect(merge.extent.values[0]).toBeCloseTo(2, 10);
        expect(merge.extent.values[1]).toBeCloseTo(1, 10);
    });

    it('produces orthonormal axes', () => {
        const b0 = box2(1, 2, 0.4, 1, 3);
        const b1 = box2(-2, 1, 1.9, 2, 0.5);
        const merge = mergeContainersOrientedBox2(b0, b1);
        expect(dot(merge.axis[0], merge.axis[0])).toBeCloseTo(1, 12);
        expect(dot(merge.axis[1], merge.axis[1])).toBeCloseTo(1, 12);
        expect(dot(merge.axis[0], merge.axis[1])).toBeCloseTo(0, 12);
    });

    it('throws for non-2D boxes', () => {
        const b3 = new OrientedBox(3);
        const b2 = box2(0, 0, 0, 1, 1);
        expect(() => mergeContainersOrientedBox2(b3, b2)).toThrow();
    });
});
