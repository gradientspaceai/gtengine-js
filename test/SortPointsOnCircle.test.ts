import { describe, it, expect } from 'vitest';
import { SortPointsOnCircle } from '../src/SortPointsOnCircle.js';

type Point2 = [number, number];

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Points at the 8 compass directions (multiples of 45 degrees) on the unit
// circle around 'center', listed in index order E, NE, N, NW, W, SW, S, SE.
function compassPoints(center: Point2): Point2[] {
    const s = Math.SQRT1_2;
    const offsets: Point2[] = [
        [1, 0],    // 0: angle 0
        [s, s],    // 1: angle pi/4
        [0, 1],    // 2: angle pi/2
        [-s, s],   // 3: angle 3*pi/4
        [-1, 0],   // 4: angle pi
        [-s, -s],  // 5: angle -3*pi/4
        [0, -1],   // 6: angle -pi/2
        [s, -s]    // 7: angle -pi/4
    ];
    return offsets.map(([x, y]) => [center[0] + x, center[1] + y] as Point2);
}

describe('SortPointsOnCircle', () => {
    // Both algorithms sort by signed angle in ascending order, with the
    // angles measured from the reference ray direction D: negative angles
    // in (-pi, 0) first, then angle 0, then positive angles up to pi.
    const sorters = [
        ['byAngle', SortPointsOnCircle.byAngle.bind(SortPointsOnCircle)],
        ['byGeometry', SortPointsOnCircle.byGeometry.bind(SortPointsOnCircle)]
    ] as const;

    for (const [name, sort] of sorters) {
        describe(name, () => {
            it('sorts all quadrants CCW with D = (1,0) about the origin', () => {
                const P = compassPoints([0, 0]);
                // Ascending signed angle: -3pi/4, -pi/2, -pi/4, 0, pi/4,
                // pi/2, 3pi/4, pi.
                expect(sort(P, [0, 0], [1, 0], true))
                    .toEqual([5, 6, 7, 0, 1, 2, 3, 4]);
            });

            it('subtracts the center before sorting', () => {
                const C: Point2 = [10, -3];
                const P = compassPoints(C);
                expect(sort(P, C, [1, 0], true))
                    .toEqual([5, 6, 7, 0, 1, 2, 3, 4]);
            });

            // Signed-zero quirk (faithful to upstream): for the two tests
            // below, the point on the ray opposite to D gets transformed
            // coordinates (-r, -0), and atan2(-0, -r) = -pi, so byAngle
            // places it first instead of last (angle +pi). byGeometry
            // compares -0 == 0 and places it last. The two orderings are
            // the same cyclic order; only the starting point differs.
            it('clockwise sorting reverses the angular sweep', () => {
                const P = compassPoints([0, 0]);
                // With sortCCW = false the clockwise angles are positive, so
                // ascending order is CCW angles 3pi/4, pi/2, pi/4 (clockwise
                // negatives first), then 0, then -pi/4, -pi/2, -3pi/4, and
                // the boundary point W at +/-pi.
                const expected = name === 'byAngle'
                    ? [4, 3, 2, 1, 0, 7, 6, 5]
                    : [3, 2, 1, 0, 7, 6, 5, 4];
                expect(sort(P, [0, 0], [1, 0], false)).toEqual(expected);
            });

            it('rotating the reference ray D rotates the starting point', () => {
                const P = compassPoints([0, 0]);
                // D = (0,1): world angle theta maps to signed angle
                // theta - pi/2, so N is angle 0 and S is the boundary at
                // +/-pi.
                const expected = name === 'byAngle'
                    ? [6, 7, 0, 1, 2, 3, 4, 5]
                    : [7, 0, 1, 2, 3, 4, 5, 6];
                expect(sort(P, [0, 0], [0, 1], true)).toEqual(expected);
            });

            it('D need not be unit length', () => {
                const P = compassPoints([0, 0]);
                expect(sort(P, [0, 0], [7, 0], true))
                    .toEqual([5, 6, 7, 0, 1, 2, 3, 4]);
            });

            it('the point on the reference ray has angle 0; the opposite ray sorts last', () => {
                const P: Point2[] = [
                    [-2, 0],  // 0: on the opposite ray, angle pi
                    [3, 0],   // 1: on the reference ray, angle 0
                    [0, -1],  // 2: angle -pi/2
                    [0, 2]    // 3: angle pi/2
                ];
                expect(sort(P, [0, 0], [1, 0], true)).toEqual([2, 1, 3, 0]);
            });

            it('ties on the same ray are broken by distance to the center', () => {
                const P: Point2[] = [
                    [2, 2],   // 0: angle pi/4, radius 2*sqrt(2)
                    [4, 0],   // 1: angle 0, farthest on the reference ray
                    [1, 1],   // 2: angle pi/4, radius sqrt(2)
                    [1, 0],   // 3: angle 0, nearest on the reference ray
                    [3, 3],   // 4: angle pi/4, radius 3*sqrt(2)
                    [2, 0]    // 5: angle 0, middle of the reference ray
                ];
                expect(sort(P, [0, 0], [1, 0], true))
                    .toEqual([3, 5, 1, 2, 0, 4]);
            });

            it('ties on the opposite (angle pi) ray are broken by distance', () => {
                const P: Point2[] = [
                    [-3, 0],  // 0: angle pi, radius 3
                    [1, 0],   // 1: angle 0
                    [-1, 0],  // 2: angle pi, radius 1
                    [-2, 0]   // 3: angle pi, radius 2
                ];
                expect(sort(P, [0, 0], [1, 0], true)).toEqual([1, 2, 3, 0]);
            });

            it('handles a single point and an empty list', () => {
                expect(sort([], [0, 0], [1, 0], true)).toEqual([]);
                expect(sort([[5, 5]], [1, 1], [1, 0], true)).toEqual([0]);
            });
        });
    }

    it('byAngle and byGeometry agree on random point sets', () => {
        const rand = makeRandom(0x50c1);
        for (let trial = 0; trial < 50; ++trial) {
            const C: Point2 = [4 * rand() - 2, 4 * rand() - 2];
            const angle = 2 * Math.PI * rand();
            const D: Point2 = [Math.cos(angle), Math.sin(angle)];
            const sortCCW = rand() < 0.5;
            const P: Point2[] = [];
            const count = 3 + Math.floor(rand() * 30);
            for (let i = 0; i < count; ++i) {
                P.push([C[0] + 10 * (rand() - 0.5), C[1] + 10 * (rand() - 0.5)]);
            }
            const indicesA = SortPointsOnCircle.byAngle(P, C, D, sortCCW);
            const indicesG = SortPointsOnCircle.byGeometry(P, C, D, sortCCW);
            expect(indicesG).toEqual(indicesA);
        }
    });

    it('random CCW orderings are circular shifts of the world-space atan2 order', () => {
        const rand = makeRandom(0xc1c1e);
        for (let trial = 0; trial < 20; ++trial) {
            const C: Point2 = [rand(), rand()];
            const angle = 2 * Math.PI * rand();
            const D: Point2 = [Math.cos(angle), Math.sin(angle)];
            const count = 5 + Math.floor(rand() * 10);
            const P: Point2[] = [];
            for (let i = 0; i < count; ++i) {
                const theta = 2 * Math.PI * rand();
                const r = 0.5 + rand();
                P.push([C[0] + r * Math.cos(theta), C[1] + r * Math.sin(theta)]);
            }
            const indices = SortPointsOnCircle.byAngle(P, C, D, true);
            // Independent computation: sorting by world-space angle produces
            // the same cyclic order, up to rotation of the starting index.
            const world = P.map((p, i) => ({
                i,
                a: Math.atan2(p[1] - C[1], p[0] - C[0])
            })).sort((u, v) => u.a - v.a).map((o) => o.i);
            const start = world.indexOf(indices[0]);
            const rotated = world.slice(start).concat(world.slice(0, start));
            expect(indices).toEqual(rotated);
        }
    });
});
