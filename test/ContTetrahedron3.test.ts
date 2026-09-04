import { describe, it, expect } from 'vitest';
import { inContainerTetrahedron3 } from '../src/ContTetrahedron3.js';
import { Tetrahedron3 } from '../src/Tetrahedron3.js';
import { Vector, add, mul } from '../src/Vector.js';
import { check, fc, latticeVector } from './helpers/arbitraries.js';
import { exactDyadic, orient3 } from './helpers/exact.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Independent containment test via barycentric coordinates computed from
// Cramer's rule on the edge basis.
function referenceInTetra(p: Vector, t: Tetrahedron3): boolean | null {
    const [a, b, c, d] = t.v;
    const m = [
        [b.values[0] - a.values[0], c.values[0] - a.values[0],
            d.values[0] - a.values[0]],
        [b.values[1] - a.values[1], c.values[1] - a.values[1],
            d.values[1] - a.values[1]],
        [b.values[2] - a.values[2], c.values[2] - a.values[2],
            d.values[2] - a.values[2]]
    ];
    const det = (m0: number[][]) =>
        m0[0][0] * (m0[1][1] * m0[2][2] - m0[1][2] * m0[2][1])
        - m0[0][1] * (m0[1][0] * m0[2][2] - m0[1][2] * m0[2][0])
        + m0[0][2] * (m0[1][0] * m0[2][1] - m0[1][1] * m0[2][0]);
    const D = det(m);
    if (D === 0) {
        return null;
    }
    const r = [p.values[0] - a.values[0], p.values[1] - a.values[1],
        p.values[2] - a.values[2]];
    const bary: number[] = [];
    for (let col = 0; col < 3; ++col) {
        const mc = m.map(row => row.slice());
        for (let row = 0; row < 3; ++row) {
            mc[row][col] = r[row];
        }
        bary.push(det(mc) / D);
    }
    const b0 = 1 - bary[0] - bary[1] - bary[2];
    const all = [b0, bary[0], bary[1], bary[2]];
    // Near-boundary points are ambiguous under round-off.
    if (all.some(x => Math.abs(x) < 1e-9)) {
        return null;
    }
    return all.every(x => x >= 0);
}

// The canonical tetrahedron used by Tetrahedron3's default constructor:
// (0,0,0), (1,0,0), (0,1,0), (0,0,1).
const canonical = new Tetrahedron3();

describe('ContTetrahedron3', () => {
    it('contains the vertices, the centroid and edge midpoints', () => {
        for (const vertex of canonical.v) {
            expect(inContainerTetrahedron3(vertex, canonical)).toBe(true);
        }
        let centroid = new Vector(3);
        for (const vertex of canonical.v) {
            centroid = add(centroid, mul(0.25, vertex));
        }
        expect(inContainerTetrahedron3(centroid, canonical)).toBe(true);
        expect(inContainerTetrahedron3(v3(0.5, 0.5, 0), canonical)).toBe(true);
        expect(inContainerTetrahedron3(v3(0, 0.5, 0.5), canonical)).toBe(true);
    });

    it('excludes points outside each face', () => {
        expect(inContainerTetrahedron3(v3(-0.01, 0.1, 0.1), canonical))
            .toBe(false);
        expect(inContainerTetrahedron3(v3(0.1, -0.01, 0.1), canonical))
            .toBe(false);
        expect(inContainerTetrahedron3(v3(0.1, 0.1, -0.01), canonical))
            .toBe(false);
        // Outside the slanted face x + y + z = 1.
        expect(inContainerTetrahedron3(v3(0.4, 0.4, 0.4), canonical))
            .toBe(false);
        expect(inContainerTetrahedron3(v3(2, 2, 2), canonical)).toBe(false);
    });

    it('is unchanged by a rigid transform of the tetrahedron', () => {
        // Rotate about z by 0.7 and translate.
        const ca = Math.cos(0.7), sa = Math.sin(0.7);
        const xform = (p: Vector) => v3(
            ca * p.values[0] - sa * p.values[1] + 3,
            sa * p.values[0] + ca * p.values[1] - 1,
            p.values[2] + 2);
        const tetra = Tetrahedron3.fromArray(canonical.v.map(xform));
        expect(inContainerTetrahedron3(xform(v3(0.2, 0.2, 0.2)), tetra))
            .toBe(true);
        expect(inContainerTetrahedron3(xform(v3(0.6, 0.6, 0.6)), tetra))
            .toBe(false);
    });

    it('reports containment for a degenerate (flat) tetrahedron only on '
        + 'the plane', () => {
        const flat = Tetrahedron3.fromVertices(
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0));
        // All triple products vanish for in-plane points, so <= 0 holds.
        expect(inContainerTetrahedron3(v3(0.3, 0.3, 0), flat)).toBe(true);
        // A point off the plane fails on one of the faces.
        const above = inContainerTetrahedron3(v3(0.3, 0.3, 1), flat);
        const below = inContainerTetrahedron3(v3(0.3, 0.3, -1), flat);
        expect(above && below).toBe(false);
    });

    it('rejects non-3D points', () => {
        expect(() => inContainerTetrahedron3(Vector.fromArray([0, 0]),
            canonical)).toThrow();
    });

    it('agrees with barycentric coordinates on random points and tetrahedra',
        () => {
            let seed = 987654321;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) % 2147483648;
                return seed / 2147483648;
            };
            const rv = () => v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2);

            let checks = 0;
            let insideCount = 0;
            for (let trial = 0; trial < 60; ++trial) {
                const a = rv(), b = rv(), c = rv(), d = rv();
                // inContainerTetrahedron3 assumes the ordering for which the
                // interior triple products are nonpositive; use the reference
                // barycentric test to identify the expected answer and skip
                // the mirrored orientation, which upstream does not support.
                const tetra = Tetrahedron3.fromVertices(a, b, c, d);
                const interior = add(add(mul(0.25, a), mul(0.25, b)),
                    add(mul(0.25, c), mul(0.25, d)));
                if (!inContainerTetrahedron3(interior, tetra)) {
                    continue; // opposite orientation
                }
                for (let k = 0; k < 40; ++k) {
                    const p = rv();
                    const expected = referenceInTetra(p, tetra);
                    if (expected === null) {
                        continue;
                    }
                    expect(inContainerTetrahedron3(p, tetra)).toBe(expected);
                    ++checks;
                    if (expected) {
                        ++insideCount;
                    }
                }
            }
            expect(checks).toBeGreaterThan(500);
            expect(insideCount).toBeGreaterThan(0);
        });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContTetrahedron3.h semantics.
// ---------------------------------------------------------------------------

describe('ContTetrahedron3 verification', () => {
    // Lattice tetrahedra and lattice query points, so the four triple scalar
    // products can be evaluated exactly with bigint arithmetic.
    const latticeTetra = fc.tuple(latticeVector(3, -5, 5), latticeVector(3, -5, 5),
        latticeVector(3, -5, 5), latticeVector(3, -5, 5))
        .map(vs => Tetrahedron3.fromArray(vs));

    // Exact sign of DotCross(b - a, c - a, p - a) via the bigint orient3
    // predicate: det[b-a, c-a, p-a] is exactly that triple scalar product.
    const exactSign = (a: Vector, b: Vector, c: Vector, p: Vector): number => {
        const scaled = exactDyadic([...a.values, ...b.values, ...c.values,
            ...p.values]);
        return orient3(scaled.slice(0, 3), scaled.slice(3, 6),
            scaled.slice(6, 9), scaled.slice(9, 12));
    };

    // Upstream evaluates DotCross(edge20, edge10, diffP0) for face <0,2,1>,
    // which is det[v2-v0, v1-v0, p-v0], and rejects when it is positive; the
    // other three faces are analogous. Cross-check every branch against the
    // exact predicate.
    it('agrees with the exact triple-scalar-product signs', () => {
        check(fc.tuple(latticeTetra, latticeVector(3, -7, 7)),
            ([tetra, p]: [Tetrahedron3, Vector]) => {
                const t = tetra.v;
                const signs = [
                    exactSign(t[0], t[2], t[1], p),
                    exactSign(t[0], t[1], t[3], p),
                    exactSign(t[0], t[3], t[2], p),
                    exactSign(t[1], t[2], t[3], p)
                ];
                const expected = signs.every(s => s <= 0);
                expect(inContainerTetrahedron3(p, tetra)).toBe(expected);
            });
    });

    // The query is short-circuited face by face; confirm that the result does
    // not depend on which face rejects first by testing every cyclic
    // relabelling of the vertices that preserves the orientation. The four
    // triple products are the same set for each relabelling.
    it('is invariant under orientation-preserving relabellings', () => {
        check(fc.tuple(latticeTetra, latticeVector(3, -7, 7)),
            ([tetra, p]: [Tetrahedron3, Vector]) => {
                const t = tetra.v;
                const base = inContainerTetrahedron3(p, tetra);
                // Even permutations of (0,1,2,3) preserve orientation.
                const even = [[1, 2, 0, 3], [2, 0, 1, 3], [0, 3, 1, 2],
                    [3, 1, 0, 2], [1, 0, 3, 2], [2, 1, 3, 0]];
                for (const perm of even) {
                    const relabelled = Tetrahedron3.fromArray(
                        perm.map(i => t[i]));
                    expect(inContainerTetrahedron3(p, relabelled)).toBe(base);
                }
            });
    });

    // Reversing the orientation of the tetrahedron (an odd permutation) flips
    // every triple product, so only points on the boundary of a degenerate
    // tetrahedron can be contained by both orderings.
    it('an orientation reversal only keeps boundary points', () => {
        check(fc.tuple(latticeTetra, latticeVector(3, -7, 7)),
            ([tetra, p]: [Tetrahedron3, Vector]) => {
                const t = tetra.v;
                const flipped = Tetrahedron3.fromArray(
                    [t[1], t[0], t[2], t[3]]);
                if (inContainerTetrahedron3(p, tetra)
                    && inContainerTetrahedron3(p, flipped)) {
                    // Every triple product is zero, so p is on all four face
                    // planes.
                    expect(exactSign(t[0], t[2], t[1], p)).toBe(0);
                    expect(exactSign(t[0], t[1], t[3], p)).toBe(0);
                    expect(exactSign(t[0], t[3], t[2], p)).toBe(0);
                    expect(exactSign(t[1], t[2], t[3], p)).toBe(0);
                }
            });
    });

    // Convexity: if two lattice points are contained, so is their midpoint.
    it('is convex (midpoints of contained points are contained)', () => {
        check(fc.tuple(latticeTetra, latticeVector(3, -7, 7),
            latticeVector(3, -7, 7)),
            ([tetra, p, q]: [Tetrahedron3, Vector, Vector]) => {
                if (!inContainerTetrahedron3(p, tetra)
                    || !inContainerTetrahedron3(q, tetra)) {
                    return;
                }
                const mid = mul(0.5, add(p, q));
                expect(inContainerTetrahedron3(mid, tetra)).toBe(true);
            });
    });

    // The vertices and the centroid are contained, provided the tetrahedron
    // has the vertex ordering the upstream query assumes, which is
    // det[v1-v0, v2-v0, v3-v0] > 0 (the face <0,2,1> test is the negation of
    // that determinant).
    it('contains its own vertices and centroid', () => {
        check(latticeTetra.filter(tetra => {
            const t = tetra.v;
            const scaled = exactDyadic([...t[0].values, ...t[1].values,
                ...t[2].values, ...t[3].values]);
            return orient3(scaled.slice(0, 3), scaled.slice(3, 6),
                scaled.slice(6, 9), scaled.slice(9, 12)) > 0;
        }), (tetra: Tetrahedron3) => {
            for (const vertex of tetra.v) {
                expect(inContainerTetrahedron3(vertex.clone(), tetra))
                    .toBe(true);
            }
            let centroid = new Vector(3);
            for (const vertex of tetra.v) { centroid = add(centroid, vertex); }
            expect(inContainerTetrahedron3(mul(0.25, centroid), tetra))
                .toBe(true);
        });
    });
});
