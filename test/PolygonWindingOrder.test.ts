import { describe, it, expect } from 'vitest';
import { PolygonWindingOrder } from '../src/PolygonWindingOrder';
import { Vector } from '../src/Vector';

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
