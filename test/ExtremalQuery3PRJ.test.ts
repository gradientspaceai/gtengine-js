import { describe, it, expect } from 'vitest';
import { ExtremalQuery3PRJ } from '../src/ExtremalQuery3PRJ.js';
import { Polyhedron3 } from '../src/Polyhedron3.js';
import { Vector, dot } from '../src/Vector.js';
import { check, expectClose, fc, rotationFrame, scaled, unitVector } from './helpers/arbitraries.js';

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
// Reference implementation: the argmax/argmin of the dot product over the
// unique vertices, with the same tie rule as the port (strict comparisons over
// ascending unique indices, so the smallest index wins a tie). Note that this
// uses the raw dot product, not the centroid-relative one, which is the point
// of the check: subtracting the centroid must not change the answer.
// ---------------------------------------------------------------------------

function bruteForce(polytope: Polyhedron3, direction: Vector):
    { positiveDirection: number, negativeDirection: number } {
    const vertices = polytope.getVertices();
    let positiveDirection = -1;
    let negativeDirection = -1;
    let maxValue = -Number.MAX_VALUE;
    let minValue = Number.MAX_VALUE;
    for (const i of polytope.getUniqueIndices()) {
        const d = dot(direction, vertices[i]);
        if (d < minValue) {
            minValue = d;
            negativeDirection = i;
        }
        if (d > maxValue) {
            maxValue = d;
            positiveDirection = i;
        }
    }
    return {
        positiveDirection: positiveDirection,
        negativeDirection: negativeDirection
    };
}

// ---------------------------------------------------------------------------
// Convex polyhedra
// ---------------------------------------------------------------------------

// The unit cube [0,1]^3, faces counterclockwise from outside.
const cubeVertices: Vector[] = [
    v3(0, 0, 0), v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0),
    v3(0, 0, 1), v3(1, 0, 1), v3(1, 1, 1), v3(0, 1, 1)
];

const cubeIndices: number[] = [
    0, 3, 2, 0, 2, 1,
    4, 5, 6, 4, 6, 7,
    0, 1, 5, 0, 5, 4,
    3, 7, 6, 3, 6, 2,
    0, 4, 7, 0, 7, 3,
    1, 2, 6, 1, 6, 5
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

const tetraVertices: Vector[] = [
    v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)
];

const tetraIndices: number[] = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];

function makePolyhedron(vertices: Vector[], indices: number[]): Polyhedron3 {
    return new Polyhedron3(vertices, indices.length, indices, true);
}

function makeCube(): Polyhedron3 {
    return makePolyhedron(cubeVertices, cubeIndices);
}

function makeOcta(): Polyhedron3 {
    return makePolyhedron(octaVertices, octaIndices);
}

function makeTetra(): Polyhedron3 {
    return makePolyhedron(tetraVertices, tetraIndices);
}

// An icosahedron, optionally subdivided: every vertex lies on the unit sphere,
// so the point set is in convex position and the triangulation is the convex
// hull. Subdividing once gives 42 vertices and 80 faces.
function makeGeodesic(subdivisions: number): { vertices: Vector[], indices: number[] } {
    const t = (1 + Math.sqrt(5)) / 2;
    const raw: number[][] = [
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]
    ];
    const vertices: Vector[] = raw.map(p => {
        const len = Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        return v3(p[0] / len, p[1] / len, p[2] / len);
    });
    let faces: number[][] = [
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]
    ];

    const midpointCache = new Map<string, number>();
    const midpoint = (a: number, b: number): number => {
        const key = a < b ? a + '_' + b : b + '_' + a;
        const found = midpointCache.get(key);
        if (found !== undefined) {
            return found;
        }
        const p = vertices[a];
        const q = vertices[b];
        const m = [
            0.5 * (p.get(0) + q.get(0)),
            0.5 * (p.get(1) + q.get(1)),
            0.5 * (p.get(2) + q.get(2))
        ];
        const len = Math.sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
        const index = vertices.length;
        vertices.push(v3(m[0] / len, m[1] / len, m[2] / len));
        midpointCache.set(key, index);
        return index;
    };

    for (let s = 0; s < subdivisions; ++s) {
        const next: number[][] = [];
        for (const f of faces) {
            const a = midpoint(f[0], f[1]);
            const b = midpoint(f[1], f[2]);
            const c = midpoint(f[2], f[0]);
            next.push([f[0], a, c], [f[1], b, a], [f[2], c, b], [a, b, c]);
        }
        faces = next;
    }

    const indices: number[] = [];
    for (const f of faces) {
        indices.push(f[0], f[1], f[2]);
    }
    return { vertices: vertices, indices: indices };
}

// A random rotation, scaling and translation applied to a vertex array.
// Affine maps with an invertible linear part preserve convexity, so the result
// is still a convex polyhedron with the same face triangulation.
function affine(vertices: Vector[], rand: () => number): Vector[] {
    const ca = Math.cos(2 * Math.PI * rand());
    const sa = Math.sin(2 * Math.PI * rand());
    const cb = Math.cos(2 * Math.PI * rand());
    const sb = Math.sin(2 * Math.PI * rand());
    const sx = 0.5 + 2 * rand();
    const sy = 0.5 + 2 * rand();
    const sz = 0.5 + 2 * rand();
    const tx = 10 * rand() - 5;
    const ty = 10 * rand() - 5;
    const tz = 10 * rand() - 5;
    return vertices.map(p => {
        const x = sx * p.get(0);
        const y = sy * p.get(1);
        const z = sz * p.get(2);
        // Rotate about z then about x.
        const x1 = ca * x - sa * y;
        const y1 = sa * x + ca * y;
        const y2 = cb * y1 - sb * z;
        const z2 = sb * y1 + cb * z;
        return v3(x1 + tx, y2 + ty, z2 + tz);
    });
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

describe('ExtremalQuery3PRJ construction', () => {
    it('inherits the base-class members', () => {
        const cube = makeCube();
        const query = new ExtremalQuery3PRJ(cube);
        expect(query.getPolytope()).toBe(cube);
        expect(query.getFaceNormals().length).toBe(12);
    });

    it('throws for an invalid polyhedron', () => {
        // Fewer than 12 indices: Polyhedron3 construction fails and the
        // vertex pool is null, so the base-class constructor throws.
        const bad = new Polyhedron3(cubeVertices, 6, [0, 1, 2, 0, 2, 3], true);
        expect(bad.isValid()).toBe(false);
        expect(() => new ExtremalQuery3PRJ(bad)).toThrow();
    });
});

describe('ExtremalQuery3PRJ known values', () => {
    it('finds the extreme cube vertices along the coordinate axes', () => {
        const query = new ExtremalQuery3PRJ(makeCube());

        // Four cube vertices share the extreme x; strict comparisons give the
        // first one visited, and the unique indices ascend, so the smallest
        // index wins. x = 1 at indices 1, 2, 5, 6; x = 0 at 0, 3, 4, 7.
        let r = query.getExtremeVertices(v3(1, 0, 0));
        expect(r.positiveDirection).toBe(1);
        expect(r.negativeDirection).toBe(0);

        // y = 1 at indices 2, 3, 6, 7; y = 0 at 0, 1, 4, 5.
        r = query.getExtremeVertices(v3(0, 1, 0));
        expect(r.positiveDirection).toBe(2);
        expect(r.negativeDirection).toBe(0);

        // z = 1 at indices 4, 5, 6, 7; z = 0 at 0, 1, 2, 3.
        r = query.getExtremeVertices(v3(0, 0, 1));
        expect(r.positiveDirection).toBe(4);
        expect(r.negativeDirection).toBe(0);

        // Negating the direction swaps the two results.
        r = query.getExtremeVertices(v3(-1, 0, 0));
        expect(r.positiveDirection).toBe(0);
        expect(r.negativeDirection).toBe(1);
    });

    it('finds the unique extreme cube vertices along a body diagonal', () => {
        const query = new ExtremalQuery3PRJ(makeCube());
        const r = query.getExtremeVertices(v3(1, 1, 1));
        // (1,1,1) is index 6 and (0,0,0) is index 0; both are unique extremes.
        expect(r.positiveDirection).toBe(6);
        expect(r.negativeDirection).toBe(0);
    });

    it('finds the extreme octahedron vertices', () => {
        const query = new ExtremalQuery3PRJ(makeOcta());

        let r = query.getExtremeVertices(v3(1, 0, 0));
        expect(r.positiveDirection).toBe(0);
        expect(r.negativeDirection).toBe(1);

        r = query.getExtremeVertices(v3(0, 1, 0));
        expect(r.positiveDirection).toBe(2);
        expect(r.negativeDirection).toBe(3);

        r = query.getExtremeVertices(v3(0, 0, 1));
        expect(r.positiveDirection).toBe(4);
        expect(r.negativeDirection).toBe(5);

        // A direction near +x but tilted toward +y still picks (1,0,0), since
        // its projection 1 beats (0,1,0)'s projection 0.1.
        r = query.getExtremeVertices(v3(1, 0.1, 0));
        expect(r.positiveDirection).toBe(0);
        expect(r.negativeDirection).toBe(1);
    });

    it('finds the extreme tetrahedron vertices', () => {
        const query = new ExtremalQuery3PRJ(makeTetra());
        const r = query.getExtremeVertices(v3(1, 1, 1));
        // (1,0,0), (0,1,0), (0,0,1) all project to 1: the smallest index wins.
        expect(r.positiveDirection).toBe(1);
        expect(r.negativeDirection).toBe(0);
    });

    it('is scale invariant in the direction', () => {
        const query = new ExtremalQuery3PRJ(makeCube());
        const base = query.getExtremeVertices(v3(0.3, -0.7, 0.2));
        for (const s of [1e-8, 0.5, 1, 3, 1e8]) {
            const r = query.getExtremeVertices(v3(0.3 * s, -0.7 * s, 0.2 * s));
            expect(r).toEqual(base);
        }
    });
});

describe('ExtremalQuery3PRJ degenerate inputs', () => {
    it('returns the first unique index for the zero direction', () => {
        // Every projection is 0. The first vertex visited sets both extremes,
        // and no strict comparison ever succeeds afterwards.
        const query = new ExtremalQuery3PRJ(makeCube());
        const r = query.getExtremeVertices(v3(0, 0, 0));
        const first = makeCube().getUniqueIndices()[0];
        expect(r.positiveDirection).toBe(first);
        expect(r.negativeDirection).toBe(first);
    });

    it('ignores vertex-pool entries the polyhedron does not reference', () => {
        // Polyhedron3 supports a pool shared by several polyhedra. The extra
        // points are far away in +x and -x, so a naive loop over the whole
        // pool would return them.
        const pool = cubeVertices.concat([v3(1000, 0, 0), v3(-1000, 0, 0)]);
        const polytope = makePolyhedron(pool, cubeIndices);
        expect(polytope.getUniqueIndices().length).toBe(8);
        const query = new ExtremalQuery3PRJ(polytope);
        const r = query.getExtremeVertices(v3(1, 0, 0));
        expect(r.positiveDirection).toBe(1);
        expect(r.negativeDirection).toBe(0);
    });

    it('handles a polyhedron far from the origin', () => {
        // The projections are taken relative to the centroid, so a large
        // translation must not change which vertices are extreme.
        const shift = 1e6;
        const shifted = cubeVertices.map(p =>
            v3(p.get(0) + shift, p.get(1) + shift, p.get(2) + shift));
        const query = new ExtremalQuery3PRJ(makePolyhedron(shifted, cubeIndices));
        const r = query.getExtremeVertices(v3(1, 0, 0));
        expect(r.positiveDirection).toBe(1);
        expect(r.negativeDirection).toBe(0);
    });

    it('returns valid indices for every direction it is given', () => {
        const query = new ExtremalQuery3PRJ(makeCube());
        const rand = makeRandom(1234);
        for (let trial = 0; trial < 200; ++trial) {
            const d = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            const r = query.getExtremeVertices(d);
            expect(r.positiveDirection).toBeGreaterThanOrEqual(0);
            expect(r.positiveDirection).toBeLessThan(cubeVertices.length);
            expect(r.negativeDirection).toBeGreaterThanOrEqual(0);
            expect(r.negativeDirection).toBeLessThan(cubeVertices.length);
        }
    });
});

describe('ExtremalQuery3PRJ against brute force', () => {
    it('matches argmax/argmin of the raw dot product on the platonic solids', () => {
        const solids = [makeCube(), makeOcta(), makeTetra()];
        const rand = makeRandom(20250902);
        for (const polytope of solids) {
            const query = new ExtremalQuery3PRJ(polytope);
            for (let trial = 0; trial < 300; ++trial) {
                const d = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
                expect(query.getExtremeVertices(d)).toEqual(bruteForce(polytope, d));
            }
        }
    });

    it('matches brute force on random convex polyhedra', () => {
        const rand = makeRandom(777);
        const base = makeGeodesic(1);
        expect(base.vertices.length).toBe(42);
        expect(base.indices.length).toBe(3 * 80);

        for (let mesh = 0; mesh < 20; ++mesh) {
            const vertices = affine(base.vertices, rand);
            const polytope = makePolyhedron(vertices, base.indices);
            expect(polytope.getUniqueIndices().length).toBe(42);
            const query = new ExtremalQuery3PRJ(polytope);
            for (let trial = 0; trial < 40; ++trial) {
                const d = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
                expect(query.getExtremeVertices(d)).toEqual(bruteForce(polytope, d));
            }
        }
    });

    it('reports a supporting vertex: no vertex projects further', () => {
        const rand = makeRandom(4242);
        const geo = makeGeodesic(0);
        const polytope = makePolyhedron(affine(geo.vertices, rand), geo.indices);
        const vertices = polytope.getVertices();
        const query = new ExtremalQuery3PRJ(polytope);
        for (let trial = 0; trial < 300; ++trial) {
            const d = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            const r = query.getExtremeVertices(d);
            const maxDot = dot(d, vertices[r.positiveDirection]);
            const minDot = dot(d, vertices[r.negativeDirection]);
            for (const i of polytope.getUniqueIndices()) {
                const value = dot(d, vertices[i]);
                expect(value).toBeLessThanOrEqual(maxDot);
                expect(value).toBeGreaterThanOrEqual(minDot);
            }
        }
    });

    it('is antisymmetric: negating the direction swaps the extremes', () => {
        const rand = makeRandom(31415);
        const geo = makeGeodesic(1);
        const polytope = makePolyhedron(affine(geo.vertices, rand), geo.indices);
        const query = new ExtremalQuery3PRJ(polytope);
        for (let trial = 0; trial < 200; ++trial) {
            const d = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            const r = query.getExtremeVertices(d);
            const rn = query.getExtremeVertices(v3(-d.get(0), -d.get(1), -d.get(2)));
            // With distinct projections (the generic case for random
            // directions on a strictly convex point set) the swap is exact.
            expect(rn.positiveDirection).toBe(r.negativeDirection);
            expect(rn.negativeDirection).toBe(r.positiveDirection);
        }
    });

    it('picks the vertex the face normals support', () => {
        // For a convex polyhedron, the outward normal of a face is a direction
        // in which every vertex of that face is extreme. The query must return
        // one of them.
        const cube = makeCube();
        const query = new ExtremalQuery3PRJ(cube);
        const normals = query.getFaceNormals();
        const indices = cube.getIndices();
        for (let t = 0; t < normals.length; ++t) {
            const r = query.getExtremeVertices(normals[t]);
            const face = [indices[3 * t], indices[3 * t + 1], indices[3 * t + 2]];
            const vertices = cube.getVertices();
            const supported = dot(normals[t], vertices[face[0]]);
            expect(dot(normals[t], vertices[r.positiveDirection]))
                .toBeCloseTo(supported, 12);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). ExtremalQuery3PRJ is the
// brute-force query, so the reference is the argmax/argmin of the *raw* dot
// product (the port projects relative to the vertex average, a shift that must
// not change the answer). Near-ties are excluded from the index comparison:
// when two vertices project to within round-off of each other the two
// computations may legitimately choose different vertices, so those draws only
// assert equality of the extreme values.
// ---------------------------------------------------------------------------

const prjSolid = fc.record({
    which: fc.integer({ min: 0, max: 3 }),
    frame: rotationFrame(3),
    scale: fc.tuple(scaled(0.5, 3, 16), scaled(0.5, 3, 16), scaled(0.5, 3, 16)),
    translate: fc.tuple(scaled(-8, 8, 32), scaled(-8, 8, 32), scaled(-8, 8, 32))
}).map(({ which, frame, scale, translate }) => {
    const geodesic = makeGeodesic(1);
    const base = [
        { vertices: cubeVertices, indices: cubeIndices },
        { vertices: octaVertices, indices: octaIndices },
        { vertices: tetraVertices, indices: tetraIndices },
        { vertices: geodesic.vertices, indices: geodesic.indices }
    ][which];
    const map = (p: Vector): Vector => {
        const s = [scale[0] * p.get(0), scale[1] * p.get(1), scale[2] * p.get(2)];
        return v3(
            frame[0].get(0) * s[0] + frame[1].get(0) * s[1] + frame[2].get(0) * s[2] + translate[0],
            frame[0].get(1) * s[0] + frame[1].get(1) * s[1] + frame[2].get(1) * s[2] + translate[1],
            frame[0].get(2) * s[0] + frame[1].get(2) * s[1] + frame[2].get(2) * s[2] + translate[2]);
    };
    const indices = base.indices.slice();
    return makePolyhedron(base.vertices.map(map), indices);
});

// The gap between the best and the second-best projection, relative to the
// spread of the projections; small values mean the argmax is ill-conditioned.
function projectionGap(polytope: Polyhedron3, direction: Vector): number {
    const vertices = polytope.getVertices();
    const values = polytope.getUniqueIndices()
        .map(i => dot(direction, vertices[i])).sort((a, b) => a - b);
    const spread = Math.max(1, values[values.length - 1] - values[0]);
    const top = values[values.length - 1] - values[values.length - 2];
    const bottom = values[1] - values[0];
    return Math.min(top, bottom) / spread;
}

describe('ExtremalQuery3PRJ verification', () => {
    it('matches the raw-dot-product argmax/argmin', () => {
        check(fc.tuple(prjSolid, unitVector(3)), ([polytope, direction]) => {
            const query = new ExtremalQuery3PRJ(polytope);
            const result = query.getExtremeVertices(direction);
            const vertices = polytope.getVertices();
            const reference = bruteForce(polytope, direction);
            const hi = dot(direction, vertices[result.positiveDirection]);
            const lo = dot(direction, vertices[result.negativeDirection]);
            const refHi = dot(direction, vertices[reference.positiveDirection]);
            const refLo = dot(direction, vertices[reference.negativeDirection]);
            expectClose(hi, refHi, 1e-9, 1e-12);
            expectClose(lo, refLo, 1e-9, 1e-12);
            if (projectionGap(polytope, direction) > 1e-8) {
                expect(result).toEqual(reference);
            }
        });
    });

    it('reports supporting vertices for every direction', () => {
        check(fc.tuple(prjSolid, unitVector(3)), ([polytope, direction]) => {
            const query = new ExtremalQuery3PRJ(polytope);
            const { positiveDirection, negativeDirection } =
                query.getExtremeVertices(direction);
            const unique = polytope.getUniqueIndices();
            expect(unique).toContain(positiveDirection);
            expect(unique).toContain(negativeDirection);
            const vertices = polytope.getVertices();
            const hi = dot(direction, vertices[positiveDirection]);
            const lo = dot(direction, vertices[negativeDirection]);
            // The scale of the projections bounds the round-off of the
            // centroid-relative evaluation.
            const scale = Math.max(1, ...unique.map(i =>
                Math.abs(dot(direction, vertices[i]))));
            for (const i of unique) {
                const d = dot(direction, vertices[i]);
                expect(d).toBeLessThanOrEqual(hi + 1e-12 * scale);
                expect(d).toBeGreaterThanOrEqual(lo - 1e-12 * scale);
            }
        });
    });

    it('negating the direction swaps the two extremes', () => {
        check(fc.tuple(prjSolid, unitVector(3)), ([polytope, direction]) => {
            const query = new ExtremalQuery3PRJ(polytope);
            const forward = query.getExtremeVertices(direction);
            const backward = query.getExtremeVertices(
                v3(-direction.get(0), -direction.get(1), -direction.get(2)));
            if (projectionGap(polytope, direction) > 1e-8) {
                expect(backward.positiveDirection).toBe(forward.negativeDirection);
                expect(backward.negativeDirection).toBe(forward.positiveDirection);
            }
        });
    });

    it('is invariant under translation of the polytope', () => {
        check(fc.tuple(prjSolid, unitVector(3),
            fc.tuple(scaled(-3, 3, 16), scaled(-3, 3, 16), scaled(-3, 3, 16))),
        ([polytope, direction, t]) => {
            if (projectionGap(polytope, direction) <= 1e-7) { return; }
            const moved = makePolyhedron(
                polytope.getVertices().map(p => v3(p.get(0) + t[0],
                    p.get(1) + t[1], p.get(2) + t[2])),
                polytope.getIndices().slice());
            expect(new ExtremalQuery3PRJ(moved).getExtremeVertices(direction))
                .toEqual(new ExtremalQuery3PRJ(polytope).getExtremeVertices(direction));
        });
    });

    it('is invariant under positive scaling of the direction', () => {
        check(fc.tuple(prjSolid, unitVector(3), scaled(0.25, 8, 16)),
            ([polytope, direction, s]) => {
                const query = new ExtremalQuery3PRJ(polytope);
                if (projectionGap(polytope, direction) <= 1e-7) { return; }
                expect(query.getExtremeVertices(
                    v3(s * direction.get(0), s * direction.get(1), s * direction.get(2))))
                    .toEqual(query.getExtremeVertices(direction));
            });
    });

    it('resolves exact ties to the smallest unique index', () => {
        // A square prism: the four top vertices tie for the +z direction and
        // the four bottom vertices tie for -z. Strict comparisons over the
        // ascending unique indices keep the first vertex visited.
        const query = new ExtremalQuery3PRJ(makeCube());
        const result = query.getExtremeVertices(v3(0, 0, 1));
        expect(result.positiveDirection).toBe(4);   // (0,0,1), first z = 1
        expect(result.negativeDirection).toBe(0);   // (0,0,0), first z = 0
        // The zero direction makes every projection 0, so both extremes are
        // the first unique index (only the strict < and > can fire, and the
        // initial sentinels are +-MAX_VALUE).
        const zero = query.getExtremeVertices(v3(0, 0, 0));
        expect(zero.positiveDirection).toBe(0);
        expect(zero.negativeDirection).toBe(0);
    });
});
