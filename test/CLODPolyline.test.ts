import { describe, it, expect } from 'vitest';
import { CLODPolyline } from '../src/CLODPolyline';
import { Vector } from '../src/Vector';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function vec(...values: number[]): Vector {
    const v = new Vector(values.length);
    for (let i = 0; i < values.length; ++i) {
        v.values[i] = values[i];
    }
    return v;
}

function key(v: Vector): string {
    return v.values.map(x => x.toFixed(12)).join(',');
}

// An independent brute-force computation of the vertex weights: the distance
// from V[z] to the segment <V[m],V[p]> divided by the length of that segment.
function bruteForceWeight(m: number, z: number, p: number,
    vertices: readonly Vector[]): number {
    const a = vertices[m].values;
    const b = vertices[p].values;
    const c = vertices[z].values;
    const n = a.length;
    let lenSqr = 0;
    for (let k = 0; k < n; ++k) {
        lenSqr += (b[k] - a[k]) * (b[k] - a[k]);
    }
    const len = Math.sqrt(lenSqr);
    if (len === 0) {
        return Number.MAX_VALUE;
    }
    let t = 0;
    for (let k = 0; k < n; ++k) {
        t += (c[k] - a[k]) * (b[k] - a[k]);
    }
    t = Math.max(0, Math.min(1, t / lenSqr));
    let distSqr = 0;
    for (let k = 0; k < n; ++k) {
        const d = c[k] - (a[k] + t * (b[k] - a[k]));
        distSqr += d * d;
    }
    return Math.sqrt(distSqr) / len;
}

function bruteForceWeights(vertices: readonly Vector[],
    closed: boolean): number[] {
    const n = vertices.length;
    const w = new Array<number>(n).fill(0);
    if (closed) {
        w[0] = bruteForceWeight(n - 1, 0, 1, vertices);
        w[n - 1] = bruteForceWeight(n - 2, n - 1, 0, vertices);
    }
    else {
        w[0] = Number.MAX_VALUE;
        w[n - 1] = Number.MAX_VALUE;
    }
    for (let z = 1; z < n - 1; ++z) {
        w[z] = bruteForceWeight(z - 1, z, z + 1, vertices);
    }
    return w;
}

// Check that the active edges form a single closed cycle (closed polyline) or
// a single open path (open polyline) over the vertices {0..L-1}.
function checkTopology(polyline: CLODPolyline): void {
    const L = polyline.getLevelOfDetail();
    const numEdges = polyline.getNumEdges();
    const edges = polyline.getEdges();
    const closed = polyline.getClosed();
    expect(numEdges).toBe(closed ? L : L - 1);

    const degree = new Array<number>(L).fill(0);
    const adjacency: number[][] = [];
    for (let i = 0; i < L; ++i) {
        adjacency.push([]);
    }
    for (let e = 0; e < numEdges; ++e) {
        const a = edges[2 * e];
        const b = edges[2 * e + 1];
        expect(a).toBeGreaterThanOrEqual(0);
        expect(a).toBeLessThan(L);
        expect(b).toBeGreaterThanOrEqual(0);
        expect(b).toBeLessThan(L);
        expect(a).not.toBe(b);
        ++degree[a];
        ++degree[b];
        adjacency[a].push(b);
        adjacency[b].push(a);
    }

    const degreeCount = new Map<number, number>();
    for (const d of degree) {
        degreeCount.set(d, (degreeCount.get(d) ?? 0) + 1);
    }
    if (closed) {
        expect(degreeCount.get(2)).toBe(L);
    }
    else {
        expect(degreeCount.get(1) ?? 0).toBe(2);
        expect(degreeCount.get(2) ?? 0).toBe(L - 2);
    }

    // The graph is connected (a single polyline, not several pieces).
    const visited = new Array<boolean>(L).fill(false);
    const stack = [0];
    visited[0] = true;
    let count = 1;
    while (stack.length > 0) {
        const v = stack.pop() as number;
        for (const w of adjacency[v]) {
            if (!visited[w]) {
                visited[w] = true;
                ++count;
                stack.push(w);
            }
        }
    }
    expect(count).toBe(L);
}

describe('CLODPolyline', () => {
    const openVertices = [
        vec(0, 0), vec(1, 0.5), vec(2, 0.1), vec(3, 1.7), vec(4, 0.3),
        vec(5, 2.2), vec(6, 0.05), vec(7, 1.1)
    ];

    const closedVertices = [
        vec(0, 0), vec(2, 0.2), vec(4, 0), vec(4.5, 2), vec(4, 4),
        vec(2, 3.7), vec(0, 4), vec(-0.5, 2)
    ];

    it('rejects too-few vertices', () => {
        expect(() => new CLODPolyline([vec(0, 0)], false))
            .toThrow(/Invalid inputs/);
        expect(() => new CLODPolyline([vec(0, 0), vec(1, 1)], true))
            .toThrow(/Invalid inputs/);
    });

    it('handles the minimal open polyline (two vertices)', () => {
        const p = new CLODPolyline([vec(0, 0), vec(3, 4)], false);
        expect(p.getNumVertices()).toBe(2);
        expect(p.getNumEdges()).toBe(1);
        expect(p.getMinLevelOfDetail()).toBe(2);
        expect(p.getMaxLevelOfDetail()).toBe(2);
        expect(p.getEdges().slice(0, 2)).toEqual([0, 1]);
        expect(p.getClosed()).toBe(false);
    });

    it('handles the minimal closed polyline (a triangle)', () => {
        const p = new CLODPolyline([vec(0, 0), vec(1, 0), vec(0, 1)], true);
        expect(p.getNumVertices()).toBe(3);
        expect(p.getNumEdges()).toBe(3);
        expect(p.getMinLevelOfDetail()).toBe(3);
        expect(p.getMaxLevelOfDetail()).toBe(3);
        expect(p.getEdges().slice(0, 6)).toEqual([0, 1, 1, 2, 2, 0]);
        expect(p.getClosed()).toBe(true);
    });

    it('reports the LOD range for an open polyline', () => {
        const p = new CLODPolyline(openVertices, false);
        expect(p.getMinLevelOfDetail()).toBe(2);
        expect(p.getMaxLevelOfDetail()).toBe(openVertices.length);
        expect(p.getLevelOfDetail()).toBe(openVertices.length);
        expect(p.getNumEdges()).toBe(openVertices.length - 1);
    });

    it('reports the LOD range for a closed polyline', () => {
        const p = new CLODPolyline(closedVertices, true);
        expect(p.getMinLevelOfDetail()).toBe(3);
        expect(p.getMaxLevelOfDetail()).toBe(closedVertices.length);
        expect(p.getLevelOfDetail()).toBe(closedVertices.length);
        expect(p.getNumEdges()).toBe(closedVertices.length);
    });

    it('stores a permutation of the input vertices', () => {
        for (const [vertices, closed] of
            [[openVertices, false], [closedVertices, true]] as const) {
            const p = new CLODPolyline(vertices, closed);
            const got = p.getVertices().map(key).sort();
            const want = vertices.map(key).sort();
            expect(got).toEqual(want);
        }
    });

    it('does not alias the caller vertices', () => {
        const input = openVertices.map(v => v.clone());
        const p = new CLODPolyline(input, false);
        const before = p.getVertices().map(key);
        input[0].values[0] = 1000;
        expect(p.getVertices().map(key)).toEqual(before);
    });

    it('orders the vertices by decreasing collapse weight', () => {
        // The upstream algorithm computes all weights once and pops them from
        // a min-heap; the stored vertex order is therefore the input sorted
        // by decreasing weight, the first entry being the most important.
        for (const [vertices, closed] of
            [[openVertices, false], [closedVertices, true]] as const) {
            const p = new CLODPolyline(vertices, closed);
            const weights = bruteForceWeights(vertices, closed);
            const stored = p.getVertices();
            const orderedWeights = stored.map(v => {
                const index = vertices.findIndex(q => key(q) === key(v));
                return weights[index];
            });
            for (let i = 1; i < orderedWeights.length; ++i) {
                expect(orderedWeights[i]).toBeLessThanOrEqual(
                    orderedWeights[i - 1]);
            }
        }
    });

    it('keeps the endpoints of an open polyline at the highest LOD priority',
        () => {
            const p = new CLODPolyline(openVertices, false);
            const stored = p.getVertices().map(key);
            const first = key(openVertices[0]);
            const last = key(openVertices[openVertices.length - 1]);
            expect(stored.slice(0, 2).sort()).toEqual([first, last].sort());
        });

    it('produces a valid polyline topology at every open LOD level', () => {
        const p = new CLODPolyline(openVertices, false);
        for (let L = p.getMaxLevelOfDetail(); L >= p.getMinLevelOfDetail(); --L) {
            p.setLevelOfDetail(L);
            expect(p.getLevelOfDetail()).toBe(L);
            checkTopology(p);
        }
    });

    it('produces a valid polyline topology at every closed LOD level', () => {
        const p = new CLODPolyline(closedVertices, true);
        for (let L = p.getMaxLevelOfDetail(); L >= p.getMinLevelOfDetail(); --L) {
            p.setLevelOfDetail(L);
            expect(p.getLevelOfDetail()).toBe(L);
            checkTopology(p);
        }
    });

    it('restores the full edge array when the LOD is raised again', () => {
        for (const [vertices, closed] of
            [[openVertices, false], [closedVertices, true]] as const) {
            const p = new CLODPolyline(vertices, closed);
            const full = p.getEdges().slice();
            const fullNumEdges = p.getNumEdges();
            p.setLevelOfDetail(p.getMinLevelOfDetail());
            p.setLevelOfDetail(p.getMaxLevelOfDetail());
            expect(p.getNumEdges()).toBe(fullNumEdges);
            expect(p.getEdges().slice(0, 2 * fullNumEdges))
                .toEqual(full.slice(0, 2 * fullNumEdges));
            checkTopology(p);
        }
    });

    it('reduces an open polyline to the segment between its endpoints', () => {
        // This pins down the fix for the upstream ComputeEdges defect: the
        // single edge that survives at the minimum level of detail joins the
        // two original endpoints (permuted indices 0 and 1).
        const p = new CLODPolyline(openVertices, false);
        p.setLevelOfDetail(2);
        expect(p.getNumEdges()).toBe(1);
        expect(p.getEdges().slice(0, 2).sort()).toEqual([0, 1]);
        const kept = p.getVertices().slice(0, 2).map(key).sort();
        const ends = [key(openVertices[0]),
            key(openVertices[openVertices.length - 1])].sort();
        expect(kept).toEqual(ends);
    });

    it('ignores out-of-range level-of-detail requests', () => {
        const p = new CLODPolyline(openVertices, false);
        p.setLevelOfDetail(5);
        expect(p.getLevelOfDetail()).toBe(5);
        p.setLevelOfDetail(1);
        expect(p.getLevelOfDetail()).toBe(5);
        p.setLevelOfDetail(1000);
        expect(p.getLevelOfDetail()).toBe(5);
    });

    it('reaches the level of detail in one step or many', () => {
        const oneStep = new CLODPolyline(closedVertices, true);
        const stepwise = new CLODPolyline(closedVertices, true);
        oneStep.setLevelOfDetail(4);
        for (let L = closedVertices.length - 1; L >= 4; --L) {
            stepwise.setLevelOfDetail(L);
        }
        expect(oneStep.getNumEdges()).toBe(stepwise.getNumEdges());
        expect(oneStep.getEdges().slice(0, 2 * oneStep.getNumEdges()))
            .toEqual(stepwise.getEdges().slice(0, 2 * stepwise.getNumEdges()));
    });

    it('handles nearly collinear vertices (near-zero weights)', () => {
        const vertices = [vec(0, 0), vec(1, 1e-12), vec(2, 0), vec(3, 1e-12),
            vec(4, 0), vec(5, 3)];
        const p = new CLODPolyline(vertices, false);
        for (let L = 6; L >= 2; --L) {
            p.setLevelOfDetail(L);
            checkTopology(p);
        }
        // The steeply-offset vertex (5,3) is an endpoint and is retained; the
        // remaining survivor at the lowest level is the other endpoint.
        p.setLevelOfDetail(2);
        const kept = p.getVertices().slice(0, 2).map(key).sort();
        expect(kept).toEqual([key(vec(0, 0)), key(vec(5, 3))].sort());
    });

    it('handles duplicate vertices (zero-length spans)', () => {
        const vertices = [vec(0, 0), vec(1, 1), vec(1, 1), vec(2, 0),
            vec(3, 2)];
        const p = new CLODPolyline(vertices, false);
        for (let L = 5; L >= 2; --L) {
            p.setLevelOfDetail(L);
            checkTopology(p);
        }
    });

    it('works in 3D and for random polylines', () => {
        const rand = makeRandom(0xc10d);
        for (let trial = 0; trial < 8; ++trial) {
            const n = 5 + trial;
            const vertices: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                vertices.push(vec(10 * rand(), 10 * rand(), 10 * rand()));
            }
            for (const closed of [false, true]) {
                const p = new CLODPolyline(vertices, closed);
                for (let L = p.getMaxLevelOfDetail();
                    L >= p.getMinLevelOfDetail(); --L) {
                    p.setLevelOfDetail(L);
                    checkTopology(p);
                }
                // Vertices retained at level L are the first L stored ones,
                // a nested sequence as L decreases.
                p.setLevelOfDetail(p.getMaxLevelOfDetail());
                expect(p.getVertices().map(key).sort())
                    .toEqual(vertices.map(key).sort());
            }
        }
    });
});
