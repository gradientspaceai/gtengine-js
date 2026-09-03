import { describe, it, expect } from 'vitest';
import { inContainerTetrahedron3 } from '../src/ContTetrahedron3.js';
import { Tetrahedron3 } from '../src/Tetrahedron3.js';
import { Vector, add, mul } from '../src/Vector.js';

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
