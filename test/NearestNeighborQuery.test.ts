import { describe, it, expect } from 'vitest';
import {
    NearestNeighborQuery,
    PositionSite,
    PositionDirectionSite
} from '../src/NearestNeighborQuery.js';
import type {
    NearestNeighborNode,
    NearestNeighborSortedPoint
} from '../src/NearestNeighborQuery.js';
import { Vector } from '../src/Vector.js';
import { check, fc, scaled, wellScaledVector } from './helpers/arbitraries.js';

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function makeSites(points: number[][]): PositionSite[] {
    return points.map((p) => new PositionSite(Vector.fromArray(p)));
}

function randomPoints(rand: () => number, count: number, dimension: number,
    scale: number): number[][] {
    const points: number[][] = [];
    for (let i = 0; i < count; ++i) {
        const p: number[] = [];
        for (let d = 0; d < dimension; ++d) {
            p.push(scale * (2 * rand() - 1));
        }
        points.push(p);
    }
    return points;
}

function sqrDistance(p: number[], q: number[]): number {
    let sum = 0;
    for (let d = 0; d < p.length; ++d) {
        const diff = p[d] - q[d];
        sum += diff * diff;
    }
    return sum;
}

// Independent brute-force reference: all site indices within 'radius' of
// 'point', sorted ascending by (squared distance, index).
function bruteForceNeighbors(points: number[][], point: number[],
    radius: number): { sqrLength: number, index: number }[] {
    const sqrRadius = radius * radius;
    const result: { sqrLength: number, index: number }[] = [];
    for (let i = 0; i < points.length; ++i) {
        const sqrLength = sqrDistance(points[i], point);
        if (sqrLength <= sqrRadius) {
            result.push({ sqrLength, index: i });
        }
    }
    result.sort((a, b) => (a.sqrLength - b.sqrLength) || (a.index - b.index));
    return result;
}

// Verify the kd-tree structure produced by the upstream construction order:
// node ranges follow the halfNumSites recursion, and the nth_element
// partition guarantees hold for every internal node.
function verifyTree(query: NearestNeighborQuery, numSites: number): void {
    const nodes = query.getNodes();
    const sorted = query.getSortedPoints();

    const checkNode = (nodeIndex: number, offset: number, count: number): void => {
        const node: NearestNeighborNode = nodes[nodeIndex];
        expect(node.numSites).toBe(count);
        if (node.siteOffset === -1) {
            // Internal node: left range [offset, offset+half), right range
            // [offset+half, offset+count), all left values <= split and all
            // right values >= split along the node axis.
            const half = Math.trunc(count / 2);
            for (let i = offset; i < offset + half; ++i) {
                expect(sorted[i].position.get(node.axis)).toBeLessThanOrEqual(node.split);
            }
            for (let i = offset + half; i < offset + count; ++i) {
                expect(sorted[i].position.get(node.axis)).toBeGreaterThanOrEqual(node.split);
            }
            // Note: sorted[offset + half] itself may have been rearranged by
            // the recursive builds of the children, so only the partition
            // inequalities are checked here (they are preserved because the
            // child rearrangements stay within [offset, offset+half) and
            // [offset+half, offset+count)).
            checkNode(node.left, offset, half);
            checkNode(node.right, offset + half, count - half);
        } else {
            // Leaf node.
            expect(node.siteOffset).toBe(offset);
            expect(node.left).toBe(-1);
            expect(node.right).toBe(-1);
            expect(node.axis).toBe(-1);
            expect(node.split).toBe(Number.MAX_VALUE);
        }
    };

    checkNode(0, 0, numSites);

    // The sorted points are a permutation of the original site indices.
    const indices = sorted.map((sp: NearestNeighborSortedPoint) => sp.index).sort((a, b) => a - b);
    expect(indices).toEqual(Array.from({ length: numSites }, (_, i) => i));
}

describe('NearestNeighborQuery', () => {
    describe('construction validation', () => {
        it('rejects invalid max levels', () => {
            const sites = makeSites([[0, 0]]);
            expect(() => new NearestNeighborQuery(sites, 1, 0)).toThrow('Invalid max level.');
            expect(() => new NearestNeighborQuery(sites, 1, 33)).toThrow('Invalid max level.');
        });

        it('rejects an empty site list', () => {
            expect(() => new NearestNeighborQuery([], 1, 8)).toThrow('Empty point list.');
        });

        it('stores the construction parameters', () => {
            const sites = makeSites([[0, 0], [1, 1], [2, 2]]);
            const query = new NearestNeighborQuery(sites, 2, 8);
            expect(query.getMaxLeafSize()).toBe(2);
            expect(query.getMaxLevel()).toBe(8);
        });
    });

    describe('tree structure', () => {
        it('a small site set fits in the root leaf', () => {
            const sites = makeSites([[0, 0], [1, 1]]);
            const query = new NearestNeighborQuery(sites, 4, 8);
            expect(query.getNumNodes()).toBe(1);
            expect(query.getDepth()).toBe(0);
            expect(query.getLargestNodeSize()).toBe(2);
            verifyTree(query, 2);
        });

        it('splits along axes in order 0,1,...,N-1 with median offsets', () => {
            const rand = makeRandom(0x5eed);
            const points = randomPoints(rand, 100, 2, 10);
            const query = new NearestNeighborQuery(makeSites(points), 1, 16);
            verifyTree(query, 100);
            // Root splits on axis 0, its children on axis 1.
            const nodes = query.getNodes();
            expect(nodes[0].axis).toBe(0);
            expect(nodes[nodes[0].left].axis).toBe(1);
            expect(nodes[nodes[0].right].axis).toBe(1);
            // A binary tree grown two nodes at a time has an odd node count.
            expect(query.getNumNodes() % 2).toBe(1);
            expect(query.getDepth()).toBeGreaterThan(0);
        });

        it('maxLevel bounds the depth and forces larger leaves', () => {
            const rand = makeRandom(0xcafe);
            const points = randomPoints(rand, 200, 2, 10);
            const query = new NearestNeighborQuery(makeSites(points), 1, 3);
            // Splits occur only while level <= maxLevel, so nodes exist at
            // levels 0..maxLevel+1.
            expect(query.getDepth()).toBeLessThanOrEqual(4);
            expect(query.getLargestNodeSize()).toBeGreaterThan(1);
            verifyTree(query, 200);
        });

        it('builds a valid tree in 3D', () => {
            const rand = makeRandom(0x3d3d);
            const points = randomPoints(rand, 60, 3, 5);
            const query = new NearestNeighborQuery(makeSites(points), 2, 10);
            verifyTree(query, 60);
        });
    });

    describe('findNeighbors', () => {
        it('rejects an invalid maximum number of neighbors', () => {
            const sites = makeSites([[0, 0]]);
            const query = new NearestNeighborQuery(sites, 1, 8);
            expect(() => query.findNeighbors(Vector.fromArray([0, 0]), 1, 0))
                .toThrow('Invalid maximum number of neighbors.');
        });

        it('finds exact neighbors on an integer grid', () => {
            // 5x5 grid, index = x + 5*y.
            const points: number[][] = [];
            for (let y = 0; y < 5; ++y) {
                for (let x = 0; x < 5; ++x) {
                    points.push([x, y]);
                }
            }
            const query = new NearestNeighborQuery(makeSites(points), 1, 8);
            const neighbors = query.findNeighbors(Vector.fromArray([2, 2]), 1, 25);
            // Within radius 1 of (2,2): itself at distance 0 and the four
            // axis neighbors at distance 1. Heap pop order is nonincreasing
            // (squared distance, index), so ties come out by decreasing
            // index.
            expect(neighbors).toEqual([17, 13, 11, 7, 12]);
        });

        it('returns an empty array when nothing is in range', () => {
            const sites = makeSites([[10, 10], [-10, -10]]);
            const query = new NearestNeighborQuery(sites, 1, 8);
            expect(query.findNeighbors(Vector.fromArray([0, 0]), 1, 4)).toEqual([]);
        });

        it('handles duplicate points and radius zero', () => {
            const points = [[1, 1], [1, 1], [3, 4], [1, 1], [1, 1], [1, 1], [7, -2]];
            const query = new NearestNeighborQuery(makeSites(points), 1, 8);
            const at = Vector.fromArray([1, 1]);
            // All five duplicates are at distance zero.
            const all = query.findNeighbors(at, 0, 10);
            expect(all.slice().sort((a, b) => a - b)).toEqual([0, 1, 3, 4, 5]);
            // The neighbor count is capped by maxNeighbors.
            const capped = query.findNeighbors(at, 0, 3);
            expect(capped.length).toBe(3);
            const unique = new Set(capped);
            expect(unique.size).toBe(3);
            for (const index of capped) {
                expect([0, 1, 3, 4, 5]).toContain(index);
            }
        });

        it('matches brute force exactly when capacity is unbounded (2D randomized)', () => {
            const rand = makeRandom(0x2d2d);
            for (let trial = 0; trial < 30; ++trial) {
                const numPoints = 1 + Math.floor(rand() * 60);
                const points = randomPoints(rand, numPoints, 2, 10);
                const query = new NearestNeighborQuery(makeSites(points),
                    1 + Math.floor(rand() * 4), 1 + Math.floor(rand() * 16));
                verifyTree(query, numPoints);
                for (let q = 0; q < 10; ++q) {
                    const target = randomPoints(rand, 1, 2, 12)[0];
                    const radius = 6 * rand();
                    const expected = bruteForceNeighbors(points, target, radius);
                    const neighbors = query.findNeighbors(
                        Vector.fromArray(target), radius, numPoints);
                    // With capacity for every candidate, the heap retains
                    // all of them and pops in descending (distance, index)
                    // order, which is the reverse of the brute-force order.
                    expect(neighbors).toEqual(
                        expected.map((entry) => entry.index).reverse());
                }
            }
        });

        it('returns the k smallest distances when capacity is limited (2D randomized)', () => {
            const rand2 = makeRandom(0x7a7a);
            for (let trial = 0; trial < 30; ++trial) {
                const numPoints = 5 + Math.floor(rand2() * 60);
                const points = randomPoints(rand2, numPoints, 2, 10);
                const query = new NearestNeighborQuery(makeSites(points), 2, 12);
                for (let q = 0; q < 10; ++q) {
                    const target = randomPoints(rand2, 1, 2, 12)[0];
                    const radius = 8 * rand2();
                    const maxNeighbors = 1 + Math.floor(rand2() * 8);
                    const expected = bruteForceNeighbors(points, target, radius);
                    const neighbors = query.findNeighbors(
                        Vector.fromArray(target), radius, maxNeighbors);

                    // Count: min(maxNeighbors, number in range).
                    expect(neighbors.length).toBe(Math.min(maxNeighbors, expected.length));

                    // Distinct indices, all within the radius, and their
                    // distances are exactly the k smallest distances.
                    expect(new Set(neighbors).size).toBe(neighbors.length);
                    const returnedDistances = neighbors
                        .map((index) => sqrDistance(points[index], target))
                        .sort((a, b) => a - b);
                    const expectedDistances = expected
                        .slice(0, neighbors.length)
                        .map((entry) => entry.sqrLength);
                    expect(returnedDistances).toEqual(expectedDistances);

                    // Output order: nonincreasing squared distance.
                    for (let i = 0; i + 1 < neighbors.length; ++i) {
                        expect(sqrDistance(points[neighbors[i]], target))
                            .toBeGreaterThanOrEqual(
                                sqrDistance(points[neighbors[i + 1]], target));
                    }
                }
            }
        });

        it('matches brute force in 3D with a deep tree', () => {
            const rand = makeRandom(0x3dbf);
            const numPoints = 128;
            const points = randomPoints(rand, numPoints, 3, 5);
            const query = new NearestNeighborQuery(makeSites(points), 1, 20);
            verifyTree(query, numPoints);
            for (let q = 0; q < 20; ++q) {
                const target = randomPoints(rand, 1, 3, 6)[0];
                const radius = 4 * rand();
                const expected = bruteForceNeighbors(points, target, radius);
                const neighbors = query.findNeighbors(
                    Vector.fromArray(target), radius, numPoints);
                expect(neighbors).toEqual(
                    expected.map((entry) => entry.index).reverse());
            }
        });

        it('finds the single site itself', () => {
            const query = new NearestNeighborQuery(makeSites([[2, 3]]), 1, 8);
            expect(query.findNeighbors(Vector.fromArray([2, 3]), 0.5, 4)).toEqual([0]);
        });
    });

    describe('sites', () => {
        it('PositionSite copies its input and returns a copy', () => {
            const p = Vector.fromArray([1, 2]);
            const site = new PositionSite(p);
            p.set(0, 99);
            expect(site.getPosition().values).toEqual([1, 2]);
            const returned = site.getPosition();
            returned.set(1, -1);
            expect(site.position.values).toEqual([1, 2]);
        });

        it('PositionDirectionSite carries a direction and reports position', () => {
            const site = new PositionDirectionSite(
                Vector.fromArray([1, 2]), Vector.fromArray([0, 1]));
            expect(site.getPosition().values).toEqual([1, 2]);
            expect(site.direction.values).toEqual([0, 1]);
        });
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md).
//
// Although the upstream TODO calls the query "approximate", the pruning tests
// are conservative: an internal node's left subtree holds only keys <= split
// and its right subtree only keys >= split, so any site within 'radius' of the
// query point survives both descent tests. The query is therefore exact for
// the sites inside the radius, and the properties below compare it against
// brute force. The distances are computed as Dot(diff, diff) in both places
// (same summation order), so the comparisons are exact.
// ---------------------------------------------------------------------------

// A random configuration: dimension, sites, leaf size, level cap, query point.
const nnConfiguration = fc.record({
    dimension: fc.integer({ min: 1, max: 3 }),
    numSites: fc.integer({ min: 1, max: 40 }),
    maxLeafSize: fc.integer({ min: 1, max: 5 }),
    maxLevel: fc.integer({ min: 1, max: 8 }),
    coords: fc.array(scaled(-10, 10, 4096), { minLength: 123, maxLength: 123 }),
    query: fc.array(scaled(-12, 12, 4096), { minLength: 3, maxLength: 3 }),
    radius: scaled(0, 14, 64),
    maxNeighbors: fc.integer({ min: 1, max: 8 })
}).map(({ dimension, numSites, maxLeafSize, maxLevel, coords, query, radius,
    maxNeighbors }) => {
    const points: number[][] = [];
    for (let i = 0; i < numSites; ++i) {
        const p: number[] = [];
        for (let d = 0; d < dimension; ++d) { p.push(coords[(3 * i + d) % coords.length]); }
        points.push(p);
    }
    return {
        points, maxLeafSize, maxLevel, radius, maxNeighbors,
        query: query.slice(0, dimension)
    };
});

describe('NearestNeighborQuery verification', () => {
    it('findNeighbors returns the k nearest sites inside the radius', () => {
        check(nnConfiguration, ({ points, maxLeafSize, maxLevel, query, radius,
            maxNeighbors }) => {
            const q = new NearestNeighborQuery(makeSites(points), maxLeafSize, maxLevel);
            const neighbors = q.findNeighbors(Vector.fromArray(query), radius,
                maxNeighbors);
            const reference = bruteForceNeighbors(points, query, radius);
            const expectedCount = Math.min(maxNeighbors, reference.length);
            expect(neighbors.length).toBe(expectedCount);

            // The reported neighbors are inside the radius and their distances
            // are exactly the smallest 'expectedCount' brute-force distances.
            const reported = neighbors.map(i => ({
                index: i, sqrLength: sqrDistance(points[i], query)
            }));
            expect([...reported].map(r => r.sqrLength).sort((a, b) => a - b))
                .toEqual(reference.slice(0, expectedCount).map(r => r.sqrLength));
            // The output order is the heap-pop order: nonincreasing distance,
            // ties broken by decreasing site index.
            for (let i = 1; i < reported.length; ++i) {
                const previous = reported[i - 1], current = reported[i];
                expect(previous.sqrLength > current.sqrLength
                    || (previous.sqrLength === current.sqrLength
                        && previous.index > current.index)).toBe(true);
            }
            // When the cutoff is not a tie, the retained set is determined and
            // must be exactly the brute-force prefix.
            const tie = expectedCount < reference.length
                && reference[expectedCount - 1].sqrLength
                    === reference[expectedCount].sqrLength;
            if (!tie) {
                expect([...neighbors].sort((a, b) => a - b))
                    .toEqual(reference.slice(0, expectedCount)
                        .map(r => r.index).sort((a, b) => a - b));
            }
        });
    });

    it('a radius covering the whole set returns the k globally nearest sites', () => {
        check(nnConfiguration, ({ points, maxLeafSize, maxLevel, query, maxNeighbors }) => {
            const q = new NearestNeighborQuery(makeSites(points), maxLeafSize, maxLevel);
            const neighbors = q.findNeighbors(Vector.fromArray(query), 1000, maxNeighbors);
            expect(neighbors.length).toBe(Math.min(maxNeighbors, points.length));
            const worstReported = Math.max(...neighbors.map(i =>
                sqrDistance(points[i], query)));
            // Nothing outside the reported set is closer than its worst
            // member. (Which member of a tie is retained depends on the
            // traversal order, so the set itself is compared by distance.)
            const returned = new Set(neighbors);
            for (let i = 0; i < points.length; ++i) {
                if (!returned.has(i)) {
                    expect(sqrDistance(points[i], query))
                        .toBeGreaterThanOrEqual(worstReported);
                }
            }
        });
    });

    it('radius zero returns only sites coincident with the query point', () => {
        check(nnConfiguration, ({ points, maxLeafSize, maxLevel, maxNeighbors }) => {
            const q = new NearestNeighborQuery(makeSites(points), maxLeafSize, maxLevel);
            // Query at an existing site: it must be found at distance 0.
            const target = points[points.length - 1];
            const neighbors = q.findNeighbors(Vector.fromArray(target), 0, maxNeighbors);
            expect(neighbors.length).toBeGreaterThanOrEqual(1);
            for (const i of neighbors) {
                expect(sqrDistance(points[i], target)).toBe(0);
            }
        });
    });

    it('the tree satisfies the upstream construction invariants', () => {
        check(nnConfiguration, ({ points, maxLeafSize, maxLevel }) => {
            const q = new NearestNeighborQuery(makeSites(points), maxLeafSize, maxLevel);
            verifyTree(q, points.length);
            // Build splits while level <= maxLevel, so the deepest level is
            // maxLevel + 1.
            expect(q.getDepth()).toBeLessThanOrEqual(maxLevel + 1);
            const nodes = q.getNodes();
            let internal = 0, leaves = 0, largest = 0, covered = 0;
            for (const node of nodes) {
                if (node.siteOffset === -1) {
                    ++internal;
                } else {
                    ++leaves;
                    largest = Math.max(largest, node.numSites);
                    covered += node.numSites;
                }
            }
            expect(nodes.length).toBe(2 * internal + 1);
            expect(leaves).toBe(internal + 1);
            expect(covered).toBe(points.length);
            expect(q.getLargestNodeSize()).toBe(largest);
            expect(q.getNumNodes()).toBe(nodes.length);
        });
    });

    it('sites are copied, so later edits do not disturb the tree', () => {
        check(nnConfiguration, ({ points, maxLeafSize, maxLevel }) => {
            const sites = makeSites(points);
            const q = new NearestNeighborQuery(sites, maxLeafSize, maxLevel);
            const before = q.getSortedPoints().map(sp => [...sp.position.values]);
            // PositionSite copies its input and getPosition returns a copy, so
            // mutating a site cannot reach the tree.
            for (const site of sites) { site.position.set(0, 1e6); }
            expect(q.getSortedPoints().map(sp => [...sp.position.values]))
                .toEqual(before);
        });
    });

    it('rejects invalid construction and query arguments', () => {
        check(fc.integer({ min: -4, max: 40 }), (maxLevel) => {
            const sites = makeSites([[0, 0], [1, 1], [2, 2]]);
            if (maxLevel > 0 && maxLevel <= 32) {
                expect(() => new NearestNeighborQuery(sites, 1, maxLevel)).not.toThrow();
            } else {
                expect(() => new NearestNeighborQuery(sites, 1, maxLevel))
                    .toThrow('Invalid max level.');
            }
        });
        // An empty site list reaches Build, whose "Empty point list" assert
        // fires (upstream behaves the same way).
        expect(() => new NearestNeighborQuery([], 1, 8)).toThrow('Empty point list.');
        const q = new NearestNeighborQuery(makeSites([[0, 0], [3, 4]]), 1, 8);
        expect(() => q.findNeighbors(Vector.fromArray([0, 0]), 1, 0))
            .toThrow('Invalid maximum number of neighbors.');
    });

    it('PositionSite and PositionDirectionSite copy their inputs', () => {
        check(fc.tuple(wellScaledVector(3), wellScaledVector(3)), ([p, d]) => {
            const site = new PositionSite(p);
            const pd = new PositionDirectionSite(p, d);
            expect(site.position).not.toBe(p);
            expect(site.getPosition().values).toEqual(p.values);
            expect(site.getPosition()).not.toBe(site.position);
            expect(pd.position.values).toEqual(p.values);
            expect(pd.direction.values).toEqual(d.values);
            expect(pd.position).not.toBe(p);
            expect(pd.direction).not.toBe(d);
        });
    });
});
