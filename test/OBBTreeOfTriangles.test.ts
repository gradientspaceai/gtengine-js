import { describe, it, expect } from 'vitest';
import { IntrLine3Triangle3FI } from '../src/IntrLine3Triangle3.js';
import { IntrRay3Triangle3FI } from '../src/IntrRay3Triangle3.js';
import { IntrSegment3Triangle3FI } from '../src/IntrSegment3Triangle3.js';
import { Line } from '../src/Line.js';
import { OBBNode, OBBTree } from '../src/OBBTree.js';
import {
    OBBTreeOfTriangles, OBBTreeOfTrianglesIntersection
} from '../src/OBBTreeOfTriangles.js';
import { Ray } from '../src/Ray.js';
import { Segment } from '../src/Segment.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, div, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

type Tri = [number, number, number];

interface Mesh {
    vertices: Vector[];
    triangles: Tri[];
}

// A triangulated height field over [0,1]^2 with (n+1)^2 vertices and 2*n*n
// triangles. The height function makes the mesh genuinely three-dimensional,
// so the covariance matrices of the centroids are nondegenerate.
function makeHeightField(n: number): Mesh {
    const vertices: Vector[] = [];
    for (let i = 0; i <= n; ++i) {
        const x = i / n;
        for (let j = 0; j <= n; ++j) {
            const y = j / n;
            const z = 0.3 * Math.sin(3.0 * x + 1.0) * Math.cos(2.0 * y - 0.5)
                + 0.1 * x * y;
            vertices.push(v3(x, y, z));
        }
    }

    const triangles: Tri[] = [];
    const idx = (i: number, j: number): number => i * (n + 1) + j;
    for (let i = 0; i < n; ++i) {
        for (let j = 0; j < n; ++j) {
            triangles.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
            triangles.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
        }
    }

    return { vertices: vertices, triangles: triangles };
}

// A closed octahedron centered at the origin.
function makeOctahedron(radius: number): Mesh {
    const r = radius;
    const vertices = [
        v3(+r, 0, 0), v3(-r, 0, 0),
        v3(0, +r, 0), v3(0, -r, 0),
        v3(0, 0, +r), v3(0, 0, -r)
    ];
    const triangles: Tri[] = [
        [0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
        [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]
    ];
    return { vertices: vertices, triangles: triangles };
}

function meshTriangle(mesh: Mesh, t: number): Triangle {
    const tri = mesh.triangles[t];
    return Triangle.fromVertices(mesh.vertices[tri[0]], mesh.vertices[tri[1]],
        mesh.vertices[tri[2]]);
}

interface Hit {
    triangleIndex: number;
    parameter: number;
    point: Vector;
}

// An independent brute-force computation of the intersections of a linear
// component with every triangle of the mesh, sorted the way execute() sorts.
function bruteForce(mesh: Mesh, queryType: number, P: Vector, Q: Vector): Hit[] {
    const hits: Hit[] = [];
    for (let t = 0; t < mesh.triangles.length; ++t) {
        const triangle = meshTriangle(mesh, t);
        if (queryType === OBBTreeOfTriangles.LINE_QUERY) {
            const result = new IntrLine3Triangle3FI().find(
                Line.fromOriginDirection(P, Q), triangle);
            if (result.intersect) {
                hits.push({
                    triangleIndex: t, parameter: result.parameter,
                    point: result.point
                });
            }
        } else if (queryType === OBBTreeOfTriangles.RAY_QUERY) {
            const result = new IntrRay3Triangle3FI().find(
                Ray.fromOriginDirection(P, Q), triangle);
            if (result.intersect) {
                hits.push({
                    triangleIndex: t, parameter: result.parameter,
                    point: result.point
                });
            }
        } else {
            const result = new IntrSegment3Triangle3FI().find(
                Segment.fromEndpoints(P, Q), triangle);
            if (result.intersect) {
                // Convert the centered-form parameter s to t in [0,1].
                hits.push({
                    triangleIndex: t,
                    parameter: result.parameter / length(sub(Q, P)) + 0.5,
                    point: result.point
                });
            }
        }
    }

    hits.sort((a, b) => {
        if (a.parameter !== b.parameter) {
            return a.parameter < b.parameter ? -1 : +1;
        }
        return a.triangleIndex - b.triangleIndex;
    });
    return hits;
}

function expectHitsEqual(actual: OBBTreeOfTrianglesIntersection[],
    expected: Hit[], tol: number = 1e-12): void {
    expect(actual.length).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(actual[i].triangleIndex).toBe(expected[i].triangleIndex);
        expect(actual[i].parameter).toBeCloseTo(expected[i].parameter, 12);
        for (let k = 0; k < 3; ++k) {
            expect(Math.abs(actual[i].point.get(k) - expected[i].point.get(k)))
                .toBeLessThan(tol);
        }
    }
}

// Verify that every reachable node's box contains all vertices of all
// triangles the node represents.
function checkContainment(tree: OBBTreeOfTriangles, mesh: Mesh,
    tol: number = 1e-9): number {
    const nodes = tree.getNodes();
    const partition = tree.getPartition();
    let visited = 0;

    const visit = (nodeIndex: number): void => {
        ++visited;
        const node = nodes[nodeIndex];
        expect(node.minIndex).not.toBe(OBBNode.invalid);
        expect(node.maxIndex).not.toBe(OBBNode.invalid);
        for (let i = node.minIndex; i <= node.maxIndex; ++i) {
            const tri = mesh.triangles[partition[i]];
            for (let k = 0; k < 3; ++k) {
                const diff = sub(mesh.vertices[tri[k]], node.box.center);
                for (let j = 0; j < 3; ++j) {
                    const d = dot(diff, node.box.axis[j]);
                    expect(Math.abs(d)).toBeLessThanOrEqual(
                        node.box.extent.get(j) + tol);
                }
            }
        }
        if (node.leftChild !== OBBNode.invalid) {
            visit(node.leftChild);
            visit(node.rightChild);
        }
    };

    visit(0);
    return visited;
}

// A small deterministic pseudorandom generator (xorshift32 scaled to [0,1)).
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return (): number => {
        s ^= s << 13; s >>>= 0;
        s ^= s >>> 17;
        s ^= s << 5; s >>>= 0;
        return s / 4294967296;
    };
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

describe('OBBTreeOfTriangles construction', () => {
    it('rejects invalid input', () => {
        const tree = new OBBTreeOfTriangles();
        expect(() => tree.createFromTriangles([v3(0, 0, 0), v3(1, 0, 0)],
            [[0, 1, 1]])).toThrow('Invalid input.');
        expect(() => tree.createFromTriangles(
            [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)], [])).toThrow('Invalid input.');
    });

    it('copies the input and computes centroids', () => {
        const mesh = makeOctahedron(2);
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles);

        expect(tree.getVertices().length).toBe(mesh.vertices.length);
        expect(tree.getTriangles().length).toBe(mesh.triangles.length);
        for (let t = 0; t < mesh.triangles.length; ++t) {
            expect([...tree.getTriangles()[t]]).toEqual(mesh.triangles[t]);
        }

        // The stored vertices are copies; mutating the input must not change
        // the tree.
        const before = tree.getVertices()[0].get(0);
        mesh.vertices[0].set(0, 100);
        expect(tree.getVertices()[0].get(0)).toBe(before);
        mesh.vertices[0].set(0, 2);

        const centroids = tree.getCentroids();
        expect(centroids.length).toBe(mesh.triangles.length);
        for (let t = 0; t < mesh.triangles.length; ++t) {
            const tri = mesh.triangles[t];
            const expected = div(add(add(mesh.vertices[tri[0]],
                mesh.vertices[tri[1]]), mesh.vertices[tri[2]]), 3);
            for (let k = 0; k < 3; ++k) {
                expect(centroids[t].get(k)).toBeCloseTo(expected.get(k), 14);
            }
        }
    });

    it('builds a balanced tree whose partition is a permutation', () => {
        const mesh = makeHeightField(4); // 32 triangles
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles);

        // 32 triangles: the height is log2(32) = 5.
        expect(tree.getHeight()).toBe(5);
        expect(tree.getNodes().length).toBe(2 ** 6 - 1);

        const partition = tree.getPartition();
        expect(partition.length).toBe(mesh.triangles.length);
        expect([...partition].sort((a, b) => a - b)).toEqual(
            mesh.triangles.map((_, i) => i));

        // Each leaf of a full-height tree holds a single triangle and every
        // node is visited by the containment walk.
        const visited = checkContainment(tree, mesh);
        expect(visited).toBe(2 * mesh.triangles.length - 1);
    });

    it('honors an explicit height and clamps the leaf ranges', () => {
        const mesh = makeHeightField(4);
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles, 2);

        expect(tree.getHeight()).toBe(2);
        checkContainment(tree, mesh);

        // The four leaves partition the 32 triangles into contiguous ranges
        // of 8 each.
        const nodes = tree.getNodes();
        const leaves: OBBNode[] = [];
        const gather = (index: number): void => {
            const node = nodes[index];
            if (node.leftChild !== OBBNode.invalid) {
                gather(node.leftChild);
                gather(node.rightChild);
            } else {
                leaves.push(node);
            }
        };
        gather(0);
        expect(leaves.length).toBe(4);
        let total = 0;
        for (const leaf of leaves) {
            total += leaf.maxIndex - leaf.minIndex + 1;
        }
        expect(total).toBe(mesh.triangles.length);

        // A height larger than the number of levels the centroids can fill
        // still produces a valid tree; the extra levels are simply unused.
        // (The base class clamps a height larger than 31 to 31, but such a
        // tree preallocates 2^32-1 nodes, so it is not exercised here.)
        const deep = new OBBTreeOfTriangles();
        deep.createFromTriangles(mesh.vertices, mesh.triangles, 8);
        expect(deep.getHeight()).toBe(8);
        checkContainment(deep, mesh);
    });

    it('creates a single-node tree for one triangle', () => {
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(vertices, [[0, 1, 2]]);

        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        const node = tree.getNodes()[0];
        expect(node.minIndex).toBe(0);
        expect(node.maxIndex).toBe(0);
        expect(node.leftChild).toBe(OBBNode.invalid);
        expect(node.rightChild).toBe(OBBNode.invalid);

        // The root is a leaf, so its box is the degenerate triangle box.
        expect(node.box.extent.get(2)).toBe(0);
    });

    it('exposes the query-type constants and the base fullHeight sentinel', () => {
        expect(OBBTreeOfTriangles.LINE_QUERY).toBe(0);
        expect(OBBTreeOfTriangles.RAY_QUERY).toBe(1);
        expect(OBBTreeOfTriangles.SEGMENT_QUERY).toBe(2);
        expect(OBBTree.fullHeight).toBe(Number.MAX_SAFE_INTEGER);
        expect(new OBBTreeOfTriangles()).toBeInstanceOf(OBBTree);
    });
});

// ---------------------------------------------------------------------------
// Leaf boxes
// ---------------------------------------------------------------------------

describe('OBBTreeOfTriangles leaf boxes', () => {
    it('builds a degenerate orthonormal frame around each triangle', () => {
        const mesh = makeOctahedron(1.5);
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles);

        const nodes = tree.getNodes();
        const partition = tree.getPartition();
        let numLeaves = 0;
        for (const node of nodes) {
            if (node.minIndex === OBBNode.invalid ||
                node.leftChild !== OBBNode.invalid) {
                continue;
            }
            ++numLeaves;
            const t = partition[node.minIndex];
            const tri = mesh.triangles[t];
            const V = [mesh.vertices[tri[0]], mesh.vertices[tri[1]],
                mesh.vertices[tri[2]]];

            // The box center is the triangle centroid.
            const centroid = div(add(add(V[0], V[1]), V[2]), 3);
            for (let k = 0; k < 3; ++k) {
                expect(node.box.center.get(k)).toBeCloseTo(centroid.get(k), 14);
            }

            // The axes are orthonormal.
            for (let i = 0; i < 3; ++i) {
                expect(length(node.box.axis[i])).toBeCloseTo(1, 14);
                for (let j = i + 1; j < 3; ++j) {
                    expect(Math.abs(dot(node.box.axis[i], node.box.axis[j])))
                        .toBeLessThan(1e-13);
                }
            }

            // axis[0] is the normalized edge V1-V0 and axis[2] is a unit
            // triangle normal.
            const edge10 = sub(V[1], V[0]);
            normalize(edge10);
            for (let k = 0; k < 3; ++k) {
                expect(node.box.axis[0].get(k)).toBeCloseTo(edge10.get(k), 14);
            }
            const normal = cross(sub(V[1], V[0]), sub(V[2], V[0]));
            normalize(normal);
            for (let k = 0; k < 3; ++k) {
                expect(node.box.axis[2].get(k)).toBeCloseTo(normal.get(k), 13);
            }

            // The third extent is zero and the box contains the triangle.
            expect(node.box.extent.get(2)).toBe(0);
            for (const vertex of V) {
                const diff = sub(vertex, node.box.center);
                expect(Math.abs(dot(diff, node.box.axis[0])))
                    .toBeLessThanOrEqual(node.box.extent.get(0) + 1e-12);
                expect(Math.abs(dot(diff, node.box.axis[1])))
                    .toBeLessThanOrEqual(node.box.extent.get(1) + 1e-12);
                expect(Math.abs(dot(diff, node.box.axis[2])))
                    .toBeLessThan(1e-12);
            }

            // At least one extent is attained (the extents are tight).
            let tight = false;
            for (const vertex of V) {
                const diff = sub(vertex, node.box.center);
                if (Math.abs(Math.abs(dot(diff, node.box.axis[0]))
                    - node.box.extent.get(0)) < 1e-12) {
                    tight = true;
                }
            }
            expect(tight).toBe(true);
        }
        expect(numLeaves).toBe(mesh.triangles.length);
    });
});

// ---------------------------------------------------------------------------
// Intersection queries
// ---------------------------------------------------------------------------

describe('OBBTreeOfTriangles execute', () => {
    it('matches brute force for a ray through an octahedron', () => {
        const mesh = makeOctahedron(1);
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles);

        // A ray along +x from outside enters one face and exits another.
        const P = v3(-3, 0.1, 0.05);
        const Q = v3(1, 0, 0);
        const hits = tree.execute(OBBTreeOfTriangles.RAY_QUERY, P, Q);
        expect(hits.length).toBe(2);
        expectHitsEqual(hits, bruteForce(mesh,
            OBBTreeOfTriangles.RAY_QUERY, P, Q));

        // The parameters are increasing and the reported points lie on the
        // ray.
        expect(hits[0].parameter).toBeLessThan(hits[1].parameter);
        for (const hit of hits) {
            const expected = add(P, mul(hit.parameter, Q));
            for (let k = 0; k < 3; ++k) {
                expect(hit.point.get(k)).toBeCloseTo(expected.get(k), 12);
            }
        }
    });

    it('returns no hits when the linear component misses the mesh', () => {
        const mesh = makeOctahedron(1);
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles);

        // A ray pointing away from the mesh.
        expect(tree.execute(OBBTreeOfTriangles.RAY_QUERY, v3(-3, 0, 0),
            v3(-1, 0, 0)).length).toBe(0);

        // A line that passes far from the mesh.
        expect(tree.execute(OBBTreeOfTriangles.LINE_QUERY, v3(0, 0, 5),
            v3(1, 0, 0)).length).toBe(0);

        // A segment that stops short of the mesh.
        expect(tree.execute(OBBTreeOfTriangles.SEGMENT_QUERY, v3(-3, 0, 0),
            v3(-2, 0, 0)).length).toBe(0);
    });

    it('reports segment parameters in [0,1] and points on the segment', () => {
        const mesh = makeOctahedron(1);
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles);

        const P = v3(-2, 0.2, -0.1);
        const Q = v3(2, 0.2, -0.1);
        const hits = tree.execute(OBBTreeOfTriangles.SEGMENT_QUERY, P, Q);
        expect(hits.length).toBe(2);
        expectHitsEqual(hits, bruteForce(mesh,
            OBBTreeOfTriangles.SEGMENT_QUERY, P, Q));

        for (const hit of hits) {
            expect(hit.parameter).toBeGreaterThanOrEqual(0);
            expect(hit.parameter).toBeLessThanOrEqual(1);
            const expected = add(mul(1 - hit.parameter, P), mul(hit.parameter, Q));
            for (let k = 0; k < 3; ++k) {
                expect(hit.point.get(k)).toBeCloseTo(expected.get(k), 12);
            }
        }

        // The two endpoints of the mesh crossing are symmetric about the
        // segment midpoint, so the parameters straddle 0.5.
        expect(hits[0].parameter).toBeLessThan(0.5);
        expect(hits[1].parameter).toBeGreaterThan(0.5);
    });

    it('reports both triangles for a hit on a shared edge (upstream #167 fix)', () => {
        // Two triangles in the plane z = 0 sharing the edge from (1,0,0) to
        // (0,1,0). Upstream collects hits in a std::set ordered by parameter
        // alone, so one of these two coincident hits would be dropped.
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0)];
        const triangles: Tri[] = [[0, 1, 2], [1, 3, 2]];
        const mesh: Mesh = { vertices: vertices, triangles: triangles };
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(vertices, triangles);

        const P = v3(0.5, 0.5, 1);
        const Q = v3(0, 0, -1);
        const hits = tree.execute(OBBTreeOfTriangles.RAY_QUERY, P, Q);

        const expected = bruteForce(mesh, OBBTreeOfTriangles.RAY_QUERY, P, Q);
        expect(expected.length).toBe(2);
        expect(hits.length).toBe(2);
        expect(hits[0].parameter).toBe(hits[1].parameter);
        expect(hits[0].parameter).toBeCloseTo(1, 14);
        // The tie is broken by the triangle index, so the order is
        // deterministic.
        expect(hits.map((h) => h.triangleIndex)).toEqual([0, 1]);
        expectHitsEqual(hits, expected);
    });

    it('orders coincident hits by triangle index', () => {
        const a = new OBBTreeOfTrianglesIntersection(7, v3(0, 0, 0), 2);
        const b = new OBBTreeOfTrianglesIntersection(3, v3(0, 0, 0), 2);
        const c = new OBBTreeOfTrianglesIntersection(9, v3(0, 0, 0), 1);
        expect(b.lessThan(a)).toBe(true);
        expect(a.lessThan(b)).toBe(false);
        expect(c.lessThan(a)).toBe(true);
        expect(a.lessThan(c)).toBe(false);
        expect(a.lessThan(a)).toBe(false);

        const d = new OBBTreeOfTrianglesIntersection();
        expect(d.triangleIndex).toBe(OBBTreeOfTrianglesIntersection.invalid);
        expect(d.parameter).toBe(0);
        expect(d.point.size).toBe(3);
    });

    it('matches brute force for many random linear components', () => {
        const mesh = makeHeightField(5); // 50 triangles
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles);
        checkContainment(tree, mesh);

        const rand = makeRandom(0x51ee7c0d);
        const numTrials = 120;
        let numHitting = 0;

        for (let trial = 0; trial < numTrials; ++trial) {
            // An origin above the height field and a target inside it, so
            // that many of the queries hit the mesh.
            const P = v3(0.5 + 1.6 * (rand() - 0.5), 0.5 + 1.6 * (rand() - 0.5),
                1.0 + rand());
            const target = v3(rand(), rand(), 0.3 * (rand() - 0.5));
            const D = sub(target, P);
            normalize(D);

            for (const queryType of [OBBTreeOfTriangles.LINE_QUERY,
                OBBTreeOfTriangles.RAY_QUERY]) {
                const hits = tree.execute(queryType, P, D);
                const expected = bruteForce(mesh, queryType, P, D);
                numHitting += expected.length > 0 ? 1 : 0;
                expectHitsEqual(hits, expected, 1e-10);
            }

            const Q = add(P, mul(4, D));
            const segHits = tree.execute(OBBTreeOfTriangles.SEGMENT_QUERY, P, Q);
            expectHitsEqual(segHits,
                bruteForce(mesh, OBBTreeOfTriangles.SEGMENT_QUERY, P, Q), 1e-10);
        }

        // The trials are meaningful only if a good fraction of them hit.
        expect(numHitting).toBeGreaterThan(numTrials);
    });

    it('matches brute force on truncated trees of every height', () => {
        const mesh = makeHeightField(4); // 32 triangles
        const rand = makeRandom(0x2b1c9f47);

        // Precompute a set of rays that hit the mesh.
        const rays: { P: Vector; D: Vector }[] = [];
        while (rays.length < 20) {
            const P = v3(rand(), rand(), 1.5);
            const target = v3(rand(), rand(), -0.2);
            const D = sub(target, P);
            normalize(D);
            rays.push({ P: P, D: D });
        }

        for (let height = 0; height <= 6; ++height) {
            const tree = new OBBTreeOfTriangles();
            tree.createFromTriangles(mesh.vertices, mesh.triangles, height);
            expect(tree.getHeight()).toBe(height);
            checkContainment(tree, mesh);

            for (const ray of rays) {
                expectHitsEqual(
                    tree.execute(OBBTreeOfTriangles.RAY_QUERY, ray.P, ray.D),
                    bruteForce(mesh, OBBTreeOfTriangles.RAY_QUERY, ray.P, ray.D),
                    1e-10);
            }
        }
    });

    it('gives the same hits for a line as for the two opposing rays', () => {
        const mesh = makeOctahedron(1);
        const tree = new OBBTreeOfTriangles();
        tree.createFromTriangles(mesh.vertices, mesh.triangles);

        const P = v3(0.05, 0.07, 0.03);
        const D = v3(1, 0, 0);
        const lineHits = tree.execute(OBBTreeOfTriangles.LINE_QUERY, P, D);
        const forward = tree.execute(OBBTreeOfTriangles.RAY_QUERY, P, D);
        const backward = tree.execute(OBBTreeOfTriangles.RAY_QUERY, P,
            mul(-1, D));
        expect(lineHits.length).toBe(2);
        expect(forward.length).toBe(1);
        expect(backward.length).toBe(1);

        const lineIndices = lineHits.map((h) => h.triangleIndex).sort(
            (a, b) => a - b);
        const rayIndices = [forward[0].triangleIndex,
            backward[0].triangleIndex].sort((a, b) => a - b);
        expect(lineIndices).toEqual(rayIndices);

        // The line parameter of the backward hit is the negative of the ray
        // parameter.
        const negative = lineHits.find(
            (h) => h.triangleIndex === backward[0].triangleIndex);
        expect(negative).toBeDefined();
        expect(negative!.parameter).toBeCloseTo(-backward[0].parameter, 12);
    });
});
