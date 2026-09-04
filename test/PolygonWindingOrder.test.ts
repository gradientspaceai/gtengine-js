import { describe, it, expect } from 'vitest';
import { PolygonWindingOrder } from '../src/PolygonWindingOrder.js';
import { Vector } from '../src/Vector.js';
import { check, fc, latticeVector } from './helpers/arbitraries.js';

function poly(points: [number, number][]): Vector[] {
    return points.map(p => Vector.fromArray(p));
}

function reversed(points: [number, number][]): [number, number][] {
    return points.slice().reverse();
}

// An independent winding test: the sign of twice the signed area computed by
// the shoelace formula. It is positive for counterclockwise polygons.
function shoelaceIsCCW(points: [number, number][]): boolean {
    let sum = 0;
    const n = points.length;
    for (let i = n - 1, j = 0; j < n; i = j++) {
        sum += points[i][0] * points[j][1] - points[j][0] * points[i][1];
    }
    return sum > 0;
}

describe('PolygonWindingOrder', () => {
    const query = new PolygonWindingOrder();

    it('classifies a counterclockwise triangle', () => {
        const p: [number, number][] = [[0, 0], [1, 0], [0, 1]];
        expect(query.isCounterClockwise(poly(p))).toBe(true);
        expect(query.isCounterClockwise(poly(reversed(p)))).toBe(false);
    });

    it('classifies a counterclockwise square', () => {
        const p: [number, number][] = [[0, 0], [2, 0], [2, 2], [0, 2]];
        expect(query.isCounterClockwise(poly(p))).toBe(true);
        expect(query.isCounterClockwise(poly(reversed(p)))).toBe(false);
    });

    it('is invariant to a cyclic rotation of the vertex list', () => {
        const p: [number, number][] = [[0, 0], [3, 0], [3, 1], [1, 1], [1, 2], [0, 2]];
        for (let k = 0; k < p.length; ++k) {
            const rotated = p.slice(k).concat(p.slice(0, k));
            expect(query.isCounterClockwise(poly(rotated))).toBe(true);
            expect(query.isCounterClockwise(poly(reversed(rotated)))).toBe(false);
        }
    });

    it('classifies a concave polygon whose lower-left vertex is convex', () => {
        // An L-shaped polygon in counterclockwise order. The lower-left
        // vertex (0,0) is convex, as it must be for any simple polygon.
        const p: [number, number][] = [[0, 0], [4, 0], [4, 1], [1, 1], [1, 3], [0, 3]];
        expect(shoelaceIsCCW(p)).toBe(true);
        expect(query.isCounterClockwise(poly(p))).toBe(true);
        expect(query.isCounterClockwise(poly(reversed(p)))).toBe(false);
    });

    it('classifies a star-shaped concave polygon', () => {
        const p: [number, number][] = [
            [0, 0], [2, 1], [4, 0], [3, 2], [4, 4], [2, 3], [0, 4], [1, 2]
        ];
        expect(shoelaceIsCCW(p)).toBe(true);
        expect(query.isCounterClockwise(poly(p))).toBe(true);
        expect(query.isCounterClockwise(poly(reversed(p)))).toBe(false);
    });

    it('handles a polygon with several vertices at the minimum x', () => {
        // Three vertices share x = 0; the lexicographic minimum is (0,0).
        const p: [number, number][] = [[0, 0], [2, 0], [2, 3], [0, 3], [0, 2], [0, 1]];
        expect(shoelaceIsCCW(p)).toBe(true);
        expect(query.isCounterClockwise(poly(p))).toBe(true);
        expect(query.isCounterClockwise(poly(reversed(p)))).toBe(false);
    });

    it('handles negative coordinates', () => {
        const p: [number, number][] = [[-3, -3], [-1, -3], [-1, -1], [-3, -1]];
        expect(query.isCounterClockwise(poly(p))).toBe(true);
        expect(query.isCounterClockwise(poly(reversed(p)))).toBe(false);
    });

    it('agrees with the shoelace formula on many convex polygons', () => {
        let seed = 13579;
        const next = () => {
            seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF;
            return seed / 0x7FFFFFFF;
        };
        for (let trial = 0; trial < 200; ++trial) {
            // Points on a circle, sorted by angle, form a convex polygon in
            // counterclockwise order.
            const n = 3 + Math.floor(next() * 6);
            const angles: number[] = [];
            for (let i = 0; i < n; ++i) {
                angles.push(next() * 2 * Math.PI);
            }
            angles.sort((a, b) => a - b);
            const cx = next() * 10 - 5;
            const cy = next() * 10 - 5;
            const r = 1 + next() * 4;
            const p: [number, number][] = angles.map(a =>
                [cx + r * Math.cos(a), cy + r * Math.sin(a)] as [number, number]);
            // Skip near-degenerate configurations where two angles nearly
            // coincide, which can make the shoelace sign unreliable.
            let degenerate = false;
            for (let i = 1; i < n; ++i) {
                if (angles[i] - angles[i - 1] < 1e-3) {
                    degenerate = true;
                }
            }
            if (degenerate) {
                continue;
            }
            expect(query.isCounterClockwise(poly(p))).toBe(shoelaceIsCCW(p));
            expect(query.isCounterClockwise(poly(reversed(p)))).toBe(!shoelaceIsCCW(p));
        }
    });

    it('returns false for a degenerate collinear "polygon"', () => {
        // The turn at the lower-left vertex has zero dot-perp, so the query
        // reports "not counterclockwise", matching the upstream 'dotPerp > 0'.
        const p: [number, number][] = [[0, 0], [1, 0], [2, 0]];
        expect(query.isCounterClockwise(poly(p))).toBe(false);
    });

    it('handles a self-touching polygon that pinches at a non-minimal vertex', () => {
        // A "bowtie-free" figure that touches itself at (2,2) but whose
        // lower-left vertex (0,0) is an ordinary convex corner. The query
        // only inspects that corner, so the classification follows the
        // local turn.
        const p: [number, number][] = [
            [0, 0], [4, 0], [2, 2], [4, 4], [0, 4], [2, 2]
        ];
        expect(query.isCounterClockwise(poly(p))).toBe(true);
        expect(query.isCounterClockwise(poly(reversed(p)))).toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// A polygon that is star-shaped about the origin, hence simple and wound
// counterclockwise: the vertices sit at distinct angles on a 10-degree grid,
// sorted increasingly, with integer radii. The four axis directions are always
// present, so no angular gap reaches 180 degrees and the origin is strictly
// inside the kernel. The 10-degree separation and radii of at least 1 keep the
// signed area comfortably away from zero, so the sign computed by the port at
// the lower-left vertex and the sign of the shoelace sum cannot disagree
// through round-off.
const starPolygon: fc.Arbitrary<[number, number][]> =
    fc.uniqueArray(fc.integer({ min: 0, max: 35 }), { maxLength: 8 })
        .map(extra => Array.from(new Set([0, 9, 18, 27, ...extra]))
            .sort((a, b) => a - b))
        .chain(steps => fc.tuple(
            fc.constant(steps),
            fc.array(fc.integer({ min: 1, max: 5 }),
                { minLength: steps.length, maxLength: steps.length })))
        .map(([steps, radii]) => steps.map((s, i): [number, number] => {
            const a = (s * 10) * Math.PI / 180;
            return [radii[i] * Math.cos(a), radii[i] * Math.sin(a)];
        }));

describe('PolygonWindingOrder verification', () => {
    const query = new PolygonWindingOrder();

    it('agrees with the shoelace sign on star-shaped polygons', () => {
        check(starPolygon, p => {
            expect(query.isCounterClockwise(poly(p))).toBe(shoelaceIsCCW(p));
            const r = reversed(p);
            expect(query.isCounterClockwise(poly(r))).toBe(shoelaceIsCCW(r));
        });
    });

    it('is invariant to every cyclic rotation of the vertex list', () => {
        check(starPolygon, p => {
            const expected = query.isCounterClockwise(poly(p));
            for (let k = 1; k < p.length; ++k) {
                const rotated = p.slice(k).concat(p.slice(0, k));
                expect(query.isCounterClockwise(poly(rotated))).toBe(expected);
            }
        });
    });

    it('is invariant to translation, which moves the lower-left vertex with it', () => {
        check(fc.tuple(starPolygon, latticeVector(2, -20, 20)), ([p, t]) => {
            const moved = p.map(([x, y]): [number, number] =>
                [x + t.values[0], y + t.values[1]]);
            expect(query.isCounterClockwise(poly(moved)))
                .toBe(query.isCounterClockwise(poly(p)));
        });
    });

    it('is invariant to a positive uniform scale and to a point reflection', () => {
        // A point reflection is a rotation by pi, so it preserves orientation
        // even though it changes which vertex is lexicographically smallest.
        check(fc.tuple(starPolygon, fc.integer({ min: 1, max: 8 })),
            ([p, s]) => {
                const expected = query.isCounterClockwise(poly(p));
                const scaled = p.map(([x, y]): [number, number] => [s * x, s * y]);
                expect(query.isCounterClockwise(poly(scaled))).toBe(expected);
                const reflected = p.map(([x, y]): [number, number] => [-x, -y]);
                expect(query.isCounterClockwise(poly(reflected))).toBe(expected);
            });
    });

    it('reverses when the polygon is mirrored in a coordinate axis', () => {
        // Mirroring reverses the orientation of a simple polygon, so the
        // answer must flip (the polygon is never degenerate here).
        check(starPolygon, p => {
            const expected = query.isCounterClockwise(poly(p));
            const mirrored = p.map(([x, y]): [number, number] => [-x, y]);
            expect(query.isCounterClockwise(poly(mirrored))).toBe(!expected);
        });
    });

    it('picks the lexicographically smallest vertex, as Vector.lessThan orders', () => {
        // Upstream selects the "lower left" vertex with std::array's
        // lexicographic operator<, which compares x first and y only on a tie.
        // The port must therefore prefer the smaller x, not the smaller y.
        check(starPolygon, p => {
            let best = 0;
            for (let i = 1; i < p.length; ++i) {
                if (p[i][0] < p[best][0] ||
                    (p[i][0] === p[best][0] && p[i][1] < p[best][1])) {
                    best = i;
                }
            }
            const n = p.length;
            const v = p[best];
            const next = p[(best + 1) % n];
            const prev = p[(best + n - 1) % n];
            const dotPerp = (next[0] - v[0]) * (prev[1] - v[1])
                - (next[1] - v[1]) * (prev[0] - v[0]);
            expect(query.isCounterClockwise(poly(p))).toBe(dotPerp > 0);
        });
    });

    it('returns false for degenerate inputs with no turn at the minimum', () => {
        // A single vertex, a repeated vertex and a collinear chain all give
        // dotPerp == 0, and upstream returns 'dotPerp > 0'.
        check(fc.tuple(latticeVector(2, -5, 5), fc.integer({ min: 1, max: 6 })),
            ([v, n]) => {
                const p: [number, number][] = [];
                for (let i = 0; i < n; ++i) {
                    p.push([v.values[0], v.values[1]]);
                }
                expect(query.isCounterClockwise(poly(p))).toBe(false);
            });
        check(fc.tuple(latticeVector(2, -5, 5), latticeVector(2, -5, 5),
            fc.integer({ min: 2, max: 6 })), ([a, d, n]) => {
            if (d.values[0] === 0 && d.values[1] === 0) { return; }
            const p: [number, number][] = [];
            for (let i = 0; i < n; ++i) {
                p.push([a.values[0] + i * d.values[0],
                    a.values[1] + i * d.values[1]]);
            }
            expect(query.isCounterClockwise(poly(p))).toBe(false);
        });
    });
});
