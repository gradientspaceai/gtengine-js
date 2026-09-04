import { describe, it, expect } from 'vitest';
import {
    getContainerOrientedBox2,
    inContainerOrientedBox2,
    mergeContainersOrientedBox2
} from '../src/ContOrientedBox2.js';
import { OrientedBox, type OrientedBox2 } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, fc, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';
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

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContOrientedBox2.h semantics.
// ---------------------------------------------------------------------------

describe('ContOrientedBox2 verification', () => {
    const box2Arb = fc.tuple(wellScaledVector(2, -6, 6), rotationFrame(2),
        fc.double({ min: 0.1, max: 4, noNaN: true }),
        fc.double({ min: 0.1, max: 4, noNaN: true }))
        .map(([c, frame, e0, e1]) => OrientedBox.fromCenterAxisExtent(
            c, frame, Vector.fromArray([e0, e1])));

    // Anisotropic lattice cloud: the two principal variances (about 4 and 0.7)
    // are well separated, so the fitted frame is unambiguous.
    const baseGrid: Vector[] = [];
    for (let i = -3; i <= 3; ++i) {
        for (let j = -1; j <= 1; ++j) {
            baseGrid.push(v(i, j));
        }
    }

    // GetContainer adjusts the fitted box so that it is tight along its own
    // axes: every input point projects inside [-e_j, e_j] and both bounds are
    // attained.
    it('is tight along its own axes and contains every point', () => {
        check(fc.array(wellScaledVector(2, -6, 6), { minLength: 1, maxLength: 14 }),
            (points: Vector[]) => {
                const box = getContainerOrientedBox2(points);
                if (box === null) {
                    return;
                }
                for (const p of points) {
                    expect(containedWithin(p, box, 1e-9)).toBe(true);
                }
                for (let j = 0; j < 2; ++j) {
                    let lo = Infinity, hi = -Infinity;
                    for (const p of points) {
                        const s = dot(sub(p, box.center), box.axis[j]);
                        lo = Math.min(lo, s);
                        hi = Math.max(hi, s);
                    }
                    expectClose(hi, box.extent.get(j), 1e-9, 1e-9);
                    expectClose(lo, -box.extent.get(j), 1e-9, 1e-9);
                }
            });
    });

    it('keeps the box axes orthonormal', () => {
        check(fc.array(wellScaledVector(2, -6, 6), { minLength: 1, maxLength: 14 }),
            (points: Vector[]) => {
                const box = getContainerOrientedBox2(points);
                if (box === null) {
                    return;
                }
                expectClose(dot(box.axis[0], box.axis[0]), 1, 1e-12, 1e-12);
                expectClose(dot(box.axis[1], box.axis[1]), 1, 1e-12, 1e-12);
                expectClose(dot(box.axis[0], box.axis[1]), 0, 1e-12, 1e-12);
            });
    });

    // Rigid motions: the extents are invariant and the center follows.
    it('is equivariant under rigid motions', () => {
        const reference = getContainerOrientedBox2(baseGrid)!;
        check(fc.tuple(rotationFrame(2), wellScaledVector(2)),
            ([frame, t]: [Vector[], Vector]) => {
                const xform = (p: Vector): Vector =>
                    add(add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])), t);
                const moved = getContainerOrientedBox2(baseGrid.map(xform))!;
                const want = [reference.extent.get(0), reference.extent.get(1)]
                    .sort((a, b) => a - b);
                const got = [moved.extent.get(0), moved.extent.get(1)]
                    .sort((a, b) => a - b);
                expectClose(got[0], want[0], 1e-9, 1e-9);
                expectClose(got[1], want[1], 1e-9, 1e-9);
                expect(length(sub(moved.center, xform(reference.center))))
                    .toBeLessThanOrEqual(1e-8);
            });
    });

    // inContainer against an independent projection test.
    it('inContainer agrees with a brute-force projection test', () => {
        check(fc.tuple(box2Arb, wellScaledVector(2, -12, 12)),
            ([box, p]: [OrientedBox2, Vector]) => {
                const diff = sub(p, box.center);
                const s0 = dot(diff, box.axis[0]);
                const s1 = dot(diff, box.axis[1]);
                const slack = Math.min(box.extent.get(0) - Math.abs(s0),
                    box.extent.get(1) - Math.abs(s1));
                if (Math.abs(slack) < 1e-9) {
                    return;
                }
                expect(inContainerOrientedBox2(p, box)).toBe(slack > 0);
            });
    });

    // MergeContainers claims a box containing both inputs (not necessarily of
    // least area). Confirm on random pairs, and check tightness along the
    // merged axes.
    it('merge contains both inputs and is tight along its own axes', () => {
        check(fc.tuple(box2Arb, box2Arb),
            ([b0, b1]: [OrientedBox2, OrientedBox2]) => {
                const merge = mergeContainersOrientedBox2(b0, b1);
                const vertices = [...b0.getVertices(), ...b1.getVertices()];
                for (const vertex of vertices) {
                    expect(containedWithin(vertex, merge, 1e-9)).toBe(true);
                }
                expectClose(dot(merge.axis[0], merge.axis[0]), 1, 1e-12, 1e-12);
                expectClose(dot(merge.axis[1], merge.axis[1]), 1, 1e-12, 1e-12);
                expectClose(dot(merge.axis[0], merge.axis[1]), 0, 1e-12, 1e-12);
                for (let j = 0; j < 2; ++j) {
                    let lo = Infinity, hi = -Infinity;
                    for (const vertex of vertices) {
                        const s = dot(sub(vertex, merge.center), merge.axis[j]);
                        lo = Math.min(lo, s);
                        hi = Math.max(hi, s);
                    }
                    expectClose(hi, merge.extent.get(j), 1e-9, 1e-9);
                    expectClose(lo, -merge.extent.get(j), 1e-9, 1e-9);
                }
            });
    });

    // Merging a box with itself reproduces it (up to rounding): the averaged
    // axis is the box axis and the projections recover the extents.
    it('merging a box with itself reproduces it', () => {
        check(box2Arb, (b: OrientedBox2) => {
            const merge = mergeContainersOrientedBox2(b, b);
            expect(length(sub(merge.center, b.center))).toBeLessThanOrEqual(1e-9);
            expectClose(merge.extent.get(0), b.extent.get(0), 1e-9, 1e-9);
            expectClose(merge.extent.get(1), b.extent.get(1), 1e-9, 1e-9);
        });
    });
});
