import { describe, expect, it } from 'vitest';
import { PolylineOffset } from '../src/PolylineOffset.js';
import { Vector, dot, length as vectorLength, normalize, sub } from '../src/Vector.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}


function expectPolyline(actual: readonly Vector[], expected: number[][],
    digits = 12): void {
    expect(actual.length).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(actual[i].get(0)).toBeCloseTo(expected[i][0], digits);
        expect(actual[i].get(1)).toBeCloseTo(expected[i][1], digits);
    }
}

describe('PolylineOffset construction', () => {
    it('computes the segment directions and right normals of an open polyline', () => {
        const vertices = [v2(0, 0), v2(2, 0), v2(2, 2)];
        const query = new PolylineOffset(vertices, true);
        // One direction/normal per segment; an open polyline has n-1 segments.
        expect(query.getDirections().length).toBe(2);
        expect(query.getNormals().length).toBe(2);
        expectPolyline(query.getDirections(), [[1, 0], [0, 1]]);
        // Perp(x,y) = (y,-x), so the normal points to the right of the segment.
        expectPolyline(query.getNormals(), [[0, -1], [1, 0]]);
    });

    it('computes a closing segment for a closed polyline', () => {
        const vertices = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
        const query = new PolylineOffset(vertices, false);
        expect(query.getDirections().length).toBe(4);
        expectPolyline(query.getDirections(),
            [[1, 0], [0, 1], [-1, 0], [0, -1]]);
        expectPolyline(query.getNormals(),
            [[0, -1], [1, 0], [0, 1], [-1, 0]]);
    });

    it('rejects too few vertices', () => {
        expect(() => new PolylineOffset([v2(0, 0)], true))
            .toThrowError('Invalid number of polyline vertices.');
        expect(() => new PolylineOffset([], true))
            .toThrowError('Invalid number of polyline vertices.');
        expect(() => new PolylineOffset([v2(0, 0), v2(1, 0)], false))
            .toThrowError('Invalid number of polyline vertices.');
        expect(() => new PolylineOffset([], false))
            .toThrowError('Invalid number of polyline vertices.');
        // The minimum counts are accepted.
        expect(() => new PolylineOffset([v2(0, 0), v2(1, 0)], true)).not.toThrow();
        expect(() => new PolylineOffset([v2(0, 0), v2(1, 0), v2(0, 1)], false))
            .not.toThrow();
    });
});

describe('PolylineOffset execute validation', () => {
    const query = new PolylineOffset([v2(0, 0), v2(1, 0), v2(1, 1)], true);

    it('requires a positive offset distance', () => {
        expect(() => query.execute(0, true, false))
            .toThrowError('The offset distance must be positive.');
        expect(() => query.execute(-1, true, false))
            .toThrowError('The offset distance must be positive.');
    });

    it('requires at least one side to be requested', () => {
        expect(() => query.execute(1, false, false))
            .toThrowError('Expecting a directive to compute an offset polyline.');
    });

    it('leaves the polyline that was not requested empty', () => {
        const right = query.execute(1, true, false);
        expect(right.rightPolyline.length).toBe(3);
        expect(right.leftPolyline).toEqual([]);

        const left = query.execute(1, false, true);
        expect(left.rightPolyline).toEqual([]);
        expect(left.leftPolyline.length).toBe(3);
    });
});

describe('PolylineOffset open polylines', () => {
    it('offsets an axis-aligned L-shape by hand-computable amounts', () => {
        // V = (0,0), (2,0), (2,2). Directions (1,0) and (0,1) with right
        // normals (0,-1) and (1,0). At the corner the bisector is
        // (0,-1) + (1,0) = (1,-1) and 1 + Dot(N0,N1) = 1.
        const query = new PolylineOffset([v2(0, 0), v2(2, 0), v2(2, 2)], true);
        const result = query.execute(1, true, true);
        expectPolyline(result.rightPolyline, [[0, -1], [3, -1], [3, 2]]);
        expectPolyline(result.leftPolyline, [[0, 1], [1, 1], [1, 2]]);
    });

    it('offsets a collinear polyline to a parallel polyline', () => {
        // Parallel consecutive segments are allowed; the bisector at the
        // interior vertex is 2*N and 1 + Dot(N,N) = 2, so the offset is d*N.
        const query = new PolylineOffset(
            [v2(0, 0), v2(1, 0), v2(2, 0), v2(5, 0)], true);
        const result = query.execute(0.5, true, true);
        expectPolyline(result.rightPolyline,
            [[0, -0.5], [1, -0.5], [2, -0.5], [5, -0.5]]);
        expectPolyline(result.leftPolyline,
            [[0, 0.5], [1, 0.5], [2, 0.5], [5, 0.5]]);
    });

    it('offsets a single segment to two parallel segments', () => {
        const query = new PolylineOffset([v2(0, 0), v2(0, 3)], true);
        const result = query.execute(2, true, true);
        // Direction (0,1), right normal Perp(0,1) = (1,0).
        expectPolyline(result.rightPolyline, [[2, 0], [2, 3]]);
        expectPolyline(result.leftPolyline, [[-2, 0], [-2, 3]]);
    });

    it('offsets a right-angle corner of a 3-4-5 triangle leg', () => {
        // A 90-degree turn: (0,0) -> (3,0) -> (3,4) is the same corner shape
        // as the L, scaled; check the interior offset vertex explicitly.
        const query = new PolylineOffset([v2(0, 0), v2(3, 0), v2(3, 4)], true);
        const result = query.execute(0.25, true, true);
        expectPolyline(result.rightPolyline,
            [[0, -0.25], [3.25, -0.25], [3.25, 4]]);
        expectPolyline(result.leftPolyline,
            [[0, 0.25], [2.75, 0.25], [2.75, 4]]);
    });

    it('the left offset is the reflection of the right offset', () => {
        const vertices = [v2(0, 0), v2(2, 1), v2(3, -1), v2(6, 2)];
        const query = new PolylineOffset(vertices, true);
        const result = query.execute(0.4, true, true);
        for (let i = 0; i < vertices.length; ++i) {
            const r = sub(result.rightPolyline[i], vertices[i]);
            const l = sub(result.leftPolyline[i], vertices[i]);
            expect(r.get(0)).toBeCloseTo(-l.get(0), 12);
            expect(r.get(1)).toBeCloseTo(-l.get(1), 12);
        }
    });
});

describe('PolylineOffset closed polylines', () => {
    it('offsets a counterclockwise unit square outward and inward', () => {
        // For a counterclockwise polygon the right normals point outward.
        const square = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
        const query = new PolylineOffset(square, false);

        const outward = query.execute(1, true, false);
        expectPolyline(outward.rightPolyline,
            [[-1, -1], [2, -1], [2, 2], [-1, 2]]);

        const inward = query.execute(0.25, false, true);
        expectPolyline(inward.leftPolyline,
            [[0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]]);
    });

    it('offsets a clockwise square in the opposite senses', () => {
        // Reversing the orientation swaps which side is the outside.
        const square = [v2(0, 0), v2(0, 1), v2(1, 1), v2(1, 0)];
        const query = new PolylineOffset(square, false);
        const result = query.execute(1, true, true);
        // The right normals now point inward, so the "right" offset shrinks
        // past the center and folds over.
        expectPolyline(result.leftPolyline,
            [[-1, -1], [-1, 2], [2, 2], [2, -1]]);
    });

    it('offsets a regular n-gon to a regular n-gon of larger circumradius', () => {
        // Offsetting a regular n-gon outward by d scales the circumradius from
        // R to R + d/cos(pi/n).
        for (const n of [3, 5, 6, 12]) {
            const R = 2;
            const d = 0.35;
            const vertices: Vector[] = [];
            for (let k = 0; k < n; ++k) {
                const a = (2 * Math.PI * k) / n;
                vertices.push(v2(R * Math.cos(a), R * Math.sin(a)));
            }
            const query = new PolylineOffset(vertices, false);
            const result = query.execute(d, true, false);
            const expectedR = R + d / Math.cos(Math.PI / n);
            expect(result.rightPolyline.length).toBe(n);
            for (let k = 0; k < n; ++k) {
                const p = result.rightPolyline[k];
                const a = (2 * Math.PI * k) / n;
                expect(vectorLength(p)).toBeCloseTo(expectedR, 11);
                // Same angular position as the source vertex.
                expect(p.get(0)).toBeCloseTo(expectedR * Math.cos(a), 11);
                expect(p.get(1)).toBeCloseTo(expectedR * Math.sin(a), 11);
            }
        }
    });

    it('a triangle is the minimum closed polyline', () => {
        const tri = [v2(0, 0), v2(4, 0), v2(0, 3)];
        const query = new PolylineOffset(tri, false);
        const result = query.execute(0.5, true, true);
        expect(result.rightPolyline.length).toBe(3);
        expect(result.leftPolyline.length).toBe(3);
    });
});

describe('PolylineOffset randomized cross-check', () => {
    // Each offset vertex must lie at signed distance d (measured along the
    // right normal) from the lines of both segments that meet at it. This is
    // the defining property of the offset and is independent of the bisector
    // formula the code uses.
    function checkDistances(vertices: readonly Vector[], isOpen: boolean,
        d: number): void {
        const query = new PolylineOffset(vertices, isOpen);
        const normals = query.getNormals();
        const result = query.execute(d, true, true);
        const n = vertices.length;

        for (let i = 0; i < n; ++i) {
            // Segment indices adjacent to vertex i.
            const adjacent: number[] = [];
            if (i > 0) {
                adjacent.push(i - 1);
            }
            else if (!isOpen) {
                adjacent.push(normals.length - 1);
            }
            if (i < normals.length && (isOpen ? i < n - 1 : true)) {
                adjacent.push(i);
            }

            for (const s of adjacent) {
                const base = vertices[s];
                const N = normals[s];
                expect(dot(sub(result.rightPolyline[i], base), N))
                    .toBeCloseTo(d, 10);
                expect(dot(sub(result.leftPolyline[i], base), N))
                    .toBeCloseTo(-d, 10);
            }
        }
    }

    it('holds for random open polylines', () => {
        let seed = 20260830;
        const rand = (): number => {
            seed = (1664525 * seed + 1013904223) >>> 0;
            return seed / 4294967296;
        };

        for (let trial = 0; trial < 25; ++trial) {
            const n = 3 + Math.floor(6 * rand());
            const vertices: Vector[] = [];
            let x = 0;
            let y = 0;
            for (let i = 0; i < n; ++i) {
                vertices.push(v2(x, y));
                // Advance with a bounded turn so consecutive directions are
                // never antiparallel (the singular case upstream disallows).
                x += 1 + rand();
                y += 2 * rand() - 1;
            }
            checkDistances(vertices, true, 0.1 + rand());
        }
    });

    it('holds for random convex closed polylines', () => {
        let seed = 13579246;
        const rand = (): number => {
            seed = (1664525 * seed + 1013904223) >>> 0;
            return seed / 4294967296;
        };

        for (let trial = 0; trial < 25; ++trial) {
            const n = 3 + Math.floor(8 * rand());
            // A convex polygon: sorted angles on an ellipse, counterclockwise.
            const angles: number[] = [];
            for (let i = 0; i < n; ++i) {
                angles.push(2 * Math.PI * rand());
            }
            angles.sort((a, b) => a - b);
            const a = 1 + rand();
            const b = 1 + rand();
            const vertices: Vector[] = angles.map(
                t => v2(a * Math.cos(t), b * Math.sin(t)));

            // Skip nearly coincident vertices and near-180-degree turns; both
            // make the bisector formula ill-conditioned (upstream documents
            // the exact 180-degree turn as a singularity).
            let ok = true;
            const dirs: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                const j = (i + 1) % n;
                const e = sub(vertices[j], vertices[i]);
                if (vectorLength(e) < 1e-2) {
                    ok = false;
                    break;
                }
                normalize(e);
                dirs.push(e);
            }
            if (ok) {
                for (let i = 0; i < n; ++i) {
                    if (1 + dot(dirs[i], dirs[(i + 1) % n]) < 0.1) {
                        ok = false;
                    }
                }
            }
            if (!ok) {
                continue;
            }

            checkDistances(vertices, false, 0.05 + 0.2 * rand());
        }
    });
});
