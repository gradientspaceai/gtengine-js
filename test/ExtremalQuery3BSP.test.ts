import { describe, it, expect } from 'vitest';
import { ExtremalQuery3BSP } from '../src/ExtremalQuery3BSP';
import { ExtremalQuery3PRJ } from '../src/ExtremalQuery3PRJ';
import { Polyhedron3 } from '../src/Polyhedron3';
import { Vector, dot, normalize } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A deterministic pseudorandom generator, so failures reproduce.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

// ---------------------------------------------------------------------------
// Reference implementation: the maximum of Dot(direction, V) over all the
// vertices used by the polyhedron. The BSP query and the brute-force query
// need not agree on the *index* when several vertices are equally extreme
// (which happens for directions lying on an arc of the Gauss map), so the
// checks below compare the extreme dot value, which is unique.
// ---------------------------------------------------------------------------

function extremeValues(polytope: Polyhedron3, direction: Vector):
    { maxValue: number, minValue: number } {
    const vertices = polytope.getVertices();
    let maxValue = -Number.MAX_VALUE;
    let minValue = Number.MAX_VALUE;
    for (const i of polytope.getUniqueIndices()) {
        const d = dot(direction, vertices[i]);
        maxValue = Math.max(maxValue, d);
        minValue = Math.min(minValue, d);
    }
    return { maxValue: maxValue, minValue: minValue };
}

// ---------------------------------------------------------------------------
// Simplicial convex polyhedra. The BSP algorithm builds the Gauss map from the
// face normals of the triangles, so it requires strictly convex polytopes
// whose triangles are not coplanar with their neighbors (a triangulated cube
// face would produce two identical normals whose cross product is zero).
// ---------------------------------------------------------------------------

// The regular tetrahedron.
const tetraVertices: Vector[] = [
    v3(1, 1, 1), v3(1, -1, -1), v3(-1, 1, -1), v3(-1, -1, 1)
];
const tetraIndices: number[] = [
    0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2
];

// The regular octahedron with vertices at +/- the coordinate axes.
const octaVertices: Vector[] = [
    v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
    v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)
];
const octaIndices: number[] = [
    0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
    2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5
];

// The regular icosahedron.
function makeIcosahedron(): { vertices: Vector[], indices: number[] } {
    const t = (1 + Math.sqrt(5)) / 2;
    const raw: [number, number, number][] = [
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]
    ];
    const vertices = raw.map(p => {
        const v = v3(p[0], p[1], p[2]);
        normalize(v);
        return v;
    });
    const indices: number[] = [
        0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11,
        1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
        3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9,
        4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1
    ];
    return { vertices: vertices, indices: indices };
}

// One level of triangle subdivision of a mesh inscribed in a sphere, with the
// new vertices projected back to the sphere. The result is strictly convex and
// simplicial.
function subdivideOnSphere(vertices: Vector[], indices: number[]):
    { vertices: Vector[], indices: number[] } {
    const outVertices = vertices.map(v => v.clone());
    const midpoints = new Map<string, number>();
    const midpoint = (i0: number, i1: number): number => {
        const key = (i0 < i1 ? `${i0},${i1}` : `${i1},${i0}`);
        const found = midpoints.get(key);
        if (found !== undefined) {
            return found;
        }
        const a = outVertices[i0];
        const b = outVertices[i1];
        const m = v3(0.5 * (a.get(0) + b.get(0)), 0.5 * (a.get(1) + b.get(1)),
            0.5 * (a.get(2) + b.get(2)));
        normalize(m);
        const index = outVertices.length;
        outVertices.push(m);
        midpoints.set(key, index);
        return index;
    };

    const outIndices: number[] = [];
    for (let t = 0; t < indices.length; t += 3) {
        const i0 = indices[t + 0];
        const i1 = indices[t + 1];
        const i2 = indices[t + 2];
        const m01 = midpoint(i0, i1);
        const m12 = midpoint(i1, i2);
        const m20 = midpoint(i2, i0);
        outIndices.push(i0, m01, m20, i1, m12, m01, i2, m20, m12, m01, m12, m20);
    }
    return { vertices: outVertices, indices: outIndices };
}

// Verify that every vertex is on the inner side of every face plane.
function isConvex(vertices: Vector[], indices: number[]): boolean {
    for (let t = 0; t < indices.length; t += 3) {
        const a = vertices[indices[t + 0]];
        const b = vertices[indices[t + 1]];
        const c = vertices[indices[t + 2]];
        const e1 = [b.get(0) - a.get(0), b.get(1) - a.get(1), b.get(2) - a.get(2)];
        const e2 = [c.get(0) - a.get(0), c.get(1) - a.get(1), c.get(2) - a.get(2)];
        const n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0]
        ];
        for (const v of vertices) {
            const d = n[0] * (v.get(0) - a.get(0)) + n[1] * (v.get(1) - a.get(1))
                + n[2] * (v.get(2) - a.get(2));
            if (d > 1e-12) {
                return false;
            }
        }
    }
    return true;
}

function makePolyhedron(vertices: Vector[], indices: number[]): Polyhedron3 {
    const polytope = new Polyhedron3(vertices, indices.length, indices, true);
    expect(polytope.isValid()).toBe(true);
    return polytope;
}

// Check that the query returns genuinely extreme vertices for the direction.
function checkDirection(query: ExtremalQuery3BSP, polytope: Polyhedron3,
    direction: Vector, tolerance: number): void {
    const result = query.getExtremeVertices(direction);
    const vertices = polytope.getVertices();
    const expected = extremeValues(polytope, direction);

    expect(result.positiveDirection).toBeGreaterThanOrEqual(0);
    expect(result.negativeDirection).toBeGreaterThanOrEqual(0);
    expect(dot(direction, vertices[result.positiveDirection]))
        .toBeCloseTo(expected.maxValue, tolerance);
    expect(dot(direction, vertices[result.negativeDirection]))
        .toBeCloseTo(expected.minValue, tolerance);
}

describe('ExtremalQuery3BSP', () => {
    it('builds a nonempty BSP tree for a tetrahedron', () => {
        const polytope = makePolyhedron(tetraVertices, tetraIndices);
        const query = new ExtremalQuery3BSP(polytope);

        // Four triangles give six mesh edges, hence six spherical-edge arcs.
        // Every vertex of a tetrahedron has exactly three adjacent triangles,
        // so no bisector arcs are created (the recursion requires a
        // separation greater than one and different from numTriangles - 1).
        expect(query.getNumNodes()).toBeGreaterThan(0);
        expect(query.getTreeDepth()).toBeGreaterThanOrEqual(2);
        expect(query.getPolytope()).toBe(polytope);
        expect(query.getFaceNormals().length).toBe(4);
    });

    it('finds the extreme vertices of a tetrahedron along the vertex directions', () => {
        const polytope = makePolyhedron(tetraVertices, tetraIndices);
        const query = new ExtremalQuery3BSP(polytope);

        // Along the direction of a vertex, that vertex is the unique maximum
        // and the opposite face is the minimum, so the indices are exact.
        for (let i = 0; i < tetraVertices.length; ++i) {
            const direction = tetraVertices[i].clone();
            normalize(direction);
            const result = query.getExtremeVertices(direction);
            expect(result.positiveDirection).toBe(i);
            expect(dot(direction, tetraVertices[result.negativeDirection]))
                .toBeCloseTo(-1 / Math.sqrt(3), 10);
        }
    });

    it('finds the extreme vertices of an octahedron along the axis directions', () => {
        const polytope = makePolyhedron(octaVertices, octaIndices);
        const query = new ExtremalQuery3BSP(polytope);

        const axes: [Vector, number, number][] = [
            [v3(1, 0, 0), 0, 1],
            [v3(-1, 0, 0), 1, 0],
            [v3(0, 1, 0), 2, 3],
            [v3(0, -1, 0), 3, 2],
            [v3(0, 0, 1), 4, 5],
            [v3(0, 0, -1), 5, 4]
        ];
        for (const [direction, pos, neg] of axes) {
            const result = query.getExtremeVertices(direction);
            expect(result.positiveDirection).toBe(pos);
            expect(result.negativeDirection).toBe(neg);
        }
    });

    it('agrees with the brute-force query on random directions (tetrahedron, octahedron)', () => {
        // These two polytopes are small enough that the upstream BSP
        // construction is exact; see the KNOWN UPSTREAM DEFECT note in
        // src/ExtremalQuery3BSP.ts and the icosahedron test below.
        const cases: [Vector[], number[]][] = [
            [tetraVertices, tetraIndices],
            [octaVertices, octaIndices]
        ];
        const rnd = makeRandom(20260901);
        for (const [vertices, indices] of cases) {
            const polytope = makePolyhedron(vertices, indices);
            const bsp = new ExtremalQuery3BSP(polytope);
            const prj = new ExtremalQuery3PRJ(polytope);
            for (let k = 0; k < 500; ++k) {
                const direction = v3(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
                if (dot(direction, direction) < 1e-6) {
                    continue;
                }
                normalize(direction);
                checkDirection(bsp, polytope, direction, 10);

                // The brute-force projection query returns the same extreme
                // values (it need not return the same indices, because ties
                // are broken differently).
                const reference = prj.getExtremeVertices(direction);
                const result = bsp.getExtremeVertices(direction);
                expect(dot(direction, vertices[result.positiveDirection]))
                    .toBeCloseTo(dot(direction, vertices[reference.positiveDirection]), 10);
                expect(dot(direction, vertices[result.negativeDirection]))
                    .toBeCloseTo(dot(direction, vertices[reference.negativeDirection]), 10);
            }
        }
    });

    it('records the upstream BSP defect on larger polytopes', () => {
        // The upstream construction is not exact for polytopes with more
        // complex Gauss maps; see the KNOWN UPSTREAM DEFECT note in
        // src/ExtremalQuery3BSP.ts. This test pins the observed behavior so a
        // change in the construction is noticed: the answers are always
        // vertices of the polytope, they are correct for the large majority of
        // directions, but they are not always correct.
        const ico = makeIcosahedron();
        const sphere = subdivideOnSphere(ico.vertices, ico.indices);
        expect(isConvex(sphere.vertices, sphere.indices)).toBe(true);

        const rnd = makeRandom(31337);
        for (const mesh of [ico, sphere]) {
            const polytope = makePolyhedron(mesh.vertices, mesh.indices);
            const query = new ExtremalQuery3BSP(polytope);
            expect(query.getNumNodes()).toBeGreaterThan(0);

            let wrong = 0;
            let total = 0;
            for (let k = 0; k < 400; ++k) {
                const direction = v3(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
                if (dot(direction, direction) < 1e-6) {
                    continue;
                }
                normalize(direction);
                const result = query.getExtremeVertices(direction);

                // The answers are always valid vertex indices.
                expect(result.positiveDirection).toBeGreaterThanOrEqual(0);
                expect(result.positiveDirection).toBeLessThan(mesh.vertices.length);
                expect(result.negativeDirection).toBeGreaterThanOrEqual(0);
                expect(result.negativeDirection).toBeLessThan(mesh.vertices.length);

                const expected = extremeValues(polytope, direction);
                ++total;
                const okPos = Math.abs(dot(direction, mesh.vertices[result.positiveDirection])
                    - expected.maxValue) <= 1e-9;
                const okNeg = Math.abs(dot(direction, mesh.vertices[result.negativeDirection])
                    - expected.minValue) <= 1e-9;
                if (!okPos || !okNeg) {
                    ++wrong;
                }
            }
            expect(total).toBeGreaterThan(300);
            expect(wrong / total).toBeLessThan(0.05);
        }
    });

    it('is deterministic: two queries over the same polytope agree exactly', () => {
        const ico = makeIcosahedron();
        const polytope = makePolyhedron(ico.vertices, ico.indices);
        const q0 = new ExtremalQuery3BSP(polytope);
        const q1 = new ExtremalQuery3BSP(polytope);
        expect(q1.getNumNodes()).toBe(q0.getNumNodes());
        expect(q1.getTreeDepth()).toBe(q0.getTreeDepth());

        const rnd = makeRandom(4242);
        for (let k = 0; k < 100; ++k) {
            const direction = v3(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            normalize(direction);
            expect(q1.getExtremeVertices(direction)).toEqual(
                q0.getExtremeVertices(direction));
        }
    });

    it('returns antipodal answers for a centrally symmetric polytope', () => {
        // The octahedron is centrally symmetric, so reversing the direction
        // swaps the two extreme vertices.
        const polytope = makePolyhedron(octaVertices, octaIndices);
        const query = new ExtremalQuery3BSP(polytope);

        const rnd = makeRandom(99);
        for (let k = 0; k < 100; ++k) {
            const direction = v3(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(direction, direction) < 1e-6) {
                continue;
            }
            normalize(direction);
            const forward = query.getExtremeVertices(direction);
            const backward = query.getExtremeVertices(
                v3(-direction.get(0), -direction.get(1), -direction.get(2)));
            expect(backward.negativeDirection).toBe(forward.positiveDirection);
            expect(backward.positiveDirection).toBe(forward.negativeDirection);
        }
    });

    it('rejects a polytope whose triangles do not form a manifold mesh', () => {
        // Repeating a triangle makes the mesh insertion fail; upstream would
        // silently associate the null triangle with a face normal.
        const indices = tetraIndices.concat([0, 2, 1]);
        const polytope = new Polyhedron3(tetraVertices, indices.length, indices, true);
        expect(polytope.isValid()).toBe(true);
        expect(() => new ExtremalQuery3BSP(polytope)).toThrow();
    });
});
