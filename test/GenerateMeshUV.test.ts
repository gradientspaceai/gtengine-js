import { describe, it, expect } from 'vitest';
import { GenerateMeshUV } from '../src/GenerateMeshUV';
import { Vector } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

interface Mesh {
    vertices: Vector[];
    indices: number[];
}

// An n-by-n grid of vertices on [0,1]^2 in the plane z = 0, triangulated with
// counterclockwise triangles. The vertex at row r and column c has index
// r*n + c and position (c*h, r*h, 0) with h = 1/(n-1). This mesh has
// rectangle topology, so its boundary is a single closed polyline.
function makeGrid(n: number, warp = 0): Mesh {
    const h = 1 / (n - 1);
    const vertices: Vector[] = [];
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            // 'warp' lifts the interior out of the plane to make a curved
            // surface patch; the boundary stays planar.
            const interior = (r > 0 && r < n - 1 && c > 0 && c < n - 1);
            const z = interior ? warp * Math.sin(Math.PI * c * h) * Math.sin(Math.PI * r * h) : 0;
            vertices.push(v3(c * h, r * h, z));
        }
    }
    const indices: number[] = [];
    for (let r = 0; r + 1 < n; ++r) {
        for (let c = 0; c + 1 < n; ++c) {
            const v00 = r * n + c;
            const v10 = r * n + c + 1;
            const v01 = (r + 1) * n + c;
            const v11 = (r + 1) * n + c + 1;
            indices.push(v00, v10, v11);
            indices.push(v00, v11, v01);
        }
    }
    return { vertices, indices };
}

function makeTCoords(count: number): Vector[] {
    const tcoords: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        tcoords.push(new Vector(2));
    }
    return tcoords;
}

function run(mesh: Mesh, numIterations: number, useSquareTopology: boolean,
    progress: ((iteration: number) => void) | null = null): Vector[] {
    const tcoords = makeTCoords(mesh.vertices.length);
    const generator = new GenerateMeshUV(0, progress);
    generator.generate(numIterations, useSquareTopology, mesh.vertices.length,
        mesh.vertices, mesh.indices.length, mesh.indices, tcoords);
    return tcoords;
}

// A probe that exposes the protected vertex-graph members so that tests can
// verify the mean value system is actually satisfied at convergence.
class GenerateMeshUVProbe extends GenerateMeshUV {
    // The maximum over the interior vertices of the distance between the
    // vertex's texture coordinate and the mean value average of its
    // neighbors' texture coordinates.
    maxResidual(): number {
        let maxLen = 0;
        for (let i = this.mNumBoundaryEdges; i < this.mNumVertices; ++i) {
            const v0 = this.mOrderedVertices[i];
            const range0 = this.mVertexGraph[v0].range0;
            const range1 = this.mVertexGraph[v0].range1;
            let t0 = 0, t1 = 0, weightSum = 0;
            for (let j = 0; j < range1; ++j) {
                const data = this.mVertexGraphData[range0 + j];
                weightSum += data.second;
                t0 += data.second * this.mTCoords[data.first].values[0];
                t1 += data.second * this.mTCoords[data.first].values[1];
            }
            const d0 = t0 / weightSum - this.mTCoords[v0].values[0];
            const d1 = t1 / weightSum - this.mTCoords[v0].values[1];
            maxLen = Math.max(maxLen, Math.hypot(d0, d1));
        }
        return maxLen;
    }

    numBoundaryEdges(): number {
        return this.mNumBoundaryEdges;
    }

    orderedVertices(): number[] {
        return this.mOrderedVertices;
    }

    distanceOf(v: number): number {
        return this.mVertexGraph[v].distance;
    }

    // Every mean value weight of an interior edge must be positive.
    minInteriorWeight(): number {
        let minWeight = Number.MAX_VALUE;
        for (const edge of this.mInteriorEdges) {
            for (const [a, b] of [[edge.V[0], edge.V[1]], [edge.V[1], edge.V[0]]]) {
                const range0 = this.mVertexGraph[a].range0;
                const range1 = this.mVertexGraph[a].range1;
                for (let j = 0; j < range1; ++j) {
                    const data = this.mVertexGraphData[range0 + j];
                    if (data.first === b) {
                        minWeight = Math.min(minWeight, data.second);
                    }
                }
            }
        }
        return minWeight;
    }
}

describe('GenerateMeshUV', () => {
    it('maps a planar unit-square grid boundary exactly onto [0,1]^2', () => {
        // n - 1 is a power of two, so the arc-length normalization of the
        // boundary is exact in binary floating point and the four forced
        // corners land exactly on the grid corners.
        const n = 9;
        const mesh = makeGrid(n);
        const tcoords = run(mesh, 2, true);
        const h = 1 / (n - 1);

        // The boundary walk starts at vertex 0, which is the corner (0,0), and
        // the arc-length parametrization of a uniform grid boundary lands the
        // four forced corners exactly on the grid corners. Therefore every
        // boundary texture coordinate equals the vertex's (x,y).
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                if (r === 0 || r === n - 1 || c === 0 || c === n - 1) {
                    const v = r * n + c;
                    expect(tcoords[v].values[0]).toBeCloseTo(c * h, 12);
                    expect(tcoords[v].values[1]).toBeCloseTo(r * h, 12);
                }
            }
        }
    });

    it('reproduces the identity map on a planar grid (linear precision of ' +
        'mean value coordinates)', () => {
        const n = 5;
        const mesh = makeGrid(n);
        const tcoords = run(mesh, 2000, true);
        const h = 1 / (n - 1);

        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                const v = r * n + c;
                expect(tcoords[v].values[0]).toBeCloseTo(c * h, 8);
                expect(tcoords[v].values[1]).toBeCloseTo(r * h, 8);
            }
        }
    });

    it('keeps the square boundary on the square when rounding shifts a forced ' +
        'corner (upstream quirk)', () => {
        // For n = 6 the boundary arc length totals 4.000000000000001, so
        // std::lower_bound / lowerBound places the first forced corner one
        // boundary vertex later than the exact arithmetic would. Upstream has
        // no tolerance here; the port preserves the behavior. The result is
        // still a monotone parametrization of the square boundary, but one
        // vertex is displaced by a grid step.
        const n = 6;
        const tcoords = run(makeGrid(n), 2, true);
        const onSquare = (t: Vector): boolean => {
            const [u, v] = t.values;
            return (Math.abs(u) < 1e-12 || Math.abs(u - 1) < 1e-12 ||
                Math.abs(v) < 1e-12 || Math.abs(v - 1) < 1e-12) &&
                u >= -1e-12 && u <= 1 + 1e-12 && v >= -1e-12 && v <= 1 + 1e-12;
        };
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                if (r === 0 || r === n - 1 || c === 0 || c === n - 1) {
                    expect(onSquare(tcoords[r * n + c])).toBe(true);
                }
            }
        }
        // The displaced vertex: (1, 1/5) geometrically, but uv (1, 0).
        expect(tcoords[1 * n + 5].values[0]).toBeCloseTo(1, 12);
        expect(tcoords[1 * n + 5].values[1]).toBeCloseTo(0, 12);
    });

    it('places every disk-topology boundary coordinate on the uv-circle and ' +
        'every interior coordinate strictly inside', () => {
        const mesh = makeGrid(6, 0.35);
        const n = 6;
        const tcoords = run(mesh, 2000, false);

        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                const v = r * n + c;
                const dx = tcoords[v].values[0] - 0.5;
                const dy = tcoords[v].values[1] - 0.5;
                const radius = Math.hypot(dx, dy);
                if (r === 0 || r === n - 1 || c === 0 || c === n - 1) {
                    expect(radius).toBeCloseTo(0.5, 12);
                } else {
                    // The interior vertices are convex combinations of their
                    // neighbors, hence lie in the convex hull of the boundary
                    // polygon, which is inside the uv-disk.
                    expect(radius).toBeLessThan(0.5);
                    expect(tcoords[v].values[0]).toBeGreaterThan(0);
                    expect(tcoords[v].values[0]).toBeLessThan(1);
                    expect(tcoords[v].values[1]).toBeGreaterThan(0);
                    expect(tcoords[v].values[1]).toBeLessThan(1);
                }
            }
        }
    });

    it('converges to a solution of the mean value system with positive weights',
        () => {
            const mesh = makeGrid(7, 0.4);
            const tcoords = makeTCoords(mesh.vertices.length);
            const probe = new GenerateMeshUVProbe(0, null);
            probe.generate(3000, false, mesh.vertices.length, mesh.vertices,
                mesh.indices.length, mesh.indices, tcoords);

            expect(probe.maxResidual()).toBeLessThan(1e-9);
            expect(probe.minInteriorWeight()).toBeGreaterThan(0);
        });

    it('orders the vertices by topological distance from the boundary', () => {
        const n = 7;
        const mesh = makeGrid(n);
        const tcoords = makeTCoords(mesh.vertices.length);
        const probe = new GenerateMeshUVProbe(0, null);
        probe.generate(2, true, mesh.vertices.length, mesh.vertices,
            mesh.indices.length, mesh.indices, tcoords);

        // A closed boundary polyline has as many vertices as edges.
        expect(probe.numBoundaryEdges()).toBe(4 * (n - 1));

        const ordered = probe.orderedVertices();
        expect(ordered.length).toBe(n * n);
        expect(new Set(ordered).size).toBe(n * n);

        // The distances along the ordering are nondecreasing, the first
        // numBoundaryEdges entries are the distance-0 (boundary) vertices, and
        // a vertex's distance is its Chebyshev-like ring index in the grid.
        let previous = -1;
        for (const v of ordered) {
            const d = probe.distanceOf(v);
            expect(d).toBeGreaterThanOrEqual(previous);
            previous = d;
        }
        for (let i = 0; i < probe.numBoundaryEdges(); ++i) {
            expect(probe.distanceOf(ordered[i])).toBe(0);
        }
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                const ring = Math.min(r, c, n - 1 - r, n - 1 - c);
                expect(probe.distanceOf(r * n + c)).toBe(ring);
            }
        }
    });

    it('rounds an odd iteration count up and reports progress per iteration',
        () => {
            const mesh = makeGrid(4);
            const seen: number[] = [];
            run(mesh, 7, true, i => { seen.push(i); });
            // 7 is rounded up to 8.
            expect(seen).toEqual([1, 2, 3, 4, 5, 6, 7, 8]);

            const none: number[] = [];
            run(mesh, 0, true, i => { none.push(i); });
            expect(none).toEqual([]);
        });

    it('is deterministic and independent of the generator instance', () => {
        const mesh = makeGrid(5, 0.2);
        const a = run(mesh, 200, false);
        const b = run(mesh, 200, false);
        const generator = new GenerateMeshUV();
        const c = makeTCoords(mesh.vertices.length);
        generator.generate(200, false, mesh.vertices.length, mesh.vertices,
            mesh.indices.length, mesh.indices, c);
        // Reuse the same generator; the state must be fully reinitialized.
        const d = makeTCoords(mesh.vertices.length);
        generator.generate(200, false, mesh.vertices.length, mesh.vertices,
            mesh.indices.length, mesh.indices, d);

        for (let i = 0; i < a.length; ++i) {
            expect(b[i].values[0]).toBe(a[i].values[0]);
            expect(b[i].values[1]).toBe(a[i].values[1]);
            expect(c[i].values[0]).toBe(a[i].values[0]);
            expect(c[i].values[1]).toBe(a[i].values[1]);
            expect(d[i].values[0]).toBe(a[i].values[0]);
            expect(d[i].values[1]).toBe(a[i].values[1]);
        }
    });

    it('handles the degenerate single-triangle mesh in disk topology', () => {
        const mesh: Mesh = {
            vertices: [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)],
            indices: [0, 1, 2]
        };
        const tcoords = run(mesh, 4, false);

        // Three boundary vertices, no interior vertices; each texture
        // coordinate is on the uv-circle and they are pairwise distinct.
        for (let i = 0; i < 3; ++i) {
            expect(Math.hypot(tcoords[i].values[0] - 0.5,
                tcoords[i].values[1] - 0.5)).toBeCloseTo(0.5, 12);
        }
        expect(tcoords[0].values[0]).toBeCloseTo(1, 12);
        expect(tcoords[0].values[1]).toBeCloseTo(0.5, 12);
        const distinct = new Set(tcoords.map(t => `${t.values[0].toFixed(9)},${t.values[1].toFixed(9)}`));
        expect(distinct.size).toBe(3);
    });

    it('preserves the boundary while a warped interior stays in the unit ' +
        'square for square topology', () => {
        const n = 5;
        const mesh = makeGrid(n, 0.6);
        const tcoords = run(mesh, 1000, true);
        const h = 1 / (n - 1);

        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                const v = r * n + c;
                if (r === 0 || r === n - 1 || c === 0 || c === n - 1) {
                    expect(tcoords[v].values[0]).toBeCloseTo(c * h, 12);
                    expect(tcoords[v].values[1]).toBeCloseTo(r * h, 12);
                } else {
                    expect(tcoords[v].values[0]).toBeGreaterThan(0);
                    expect(tcoords[v].values[0]).toBeLessThan(1);
                    expect(tcoords[v].values[1]).toBeGreaterThan(0);
                    expect(tcoords[v].values[1]).toBeLessThan(1);
                }
            }
        }

        // The warped patch is symmetric about the diagonal x = y, so the
        // texture coordinates must be symmetric too.
        for (let r = 1; r < n - 1; ++r) {
            for (let c = 1; c < n - 1; ++c) {
                expect(tcoords[r * n + c].values[0]).toBeCloseTo(
                    tcoords[c * n + r].values[1], 9);
            }
        }
    });
});
