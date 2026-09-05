import { describe, it, expect } from 'vitest';
import { IntpLinearNonuniform3 } from '../src/IntpLinearNonuniform3.js';
import type { IntpLinearNonuniform3TetrahedronMesh } from '../src/IntpLinearNonuniform3.js';
import { Vector } from '../src/Vector.js';
import { computeBarycentrics3 } from '../src/Vector3.js';
import { Delaunay3 } from '../src/Delaunay3.js';
import { Delaunay3Mesh } from '../src/Delaunay3Mesh.js';
import { check, expectClose, fc, latticeVector, seededRandom, unitVector }
    from './helpers/arbitraries.js';

// A minimal tetrahedron mesh adapter that satisfies the interface required
// by IntpLinearNonuniform3. The tetrahedra are index quadruples into
// 'vertices'.
class TestMesh implements IntpLinearNonuniform3TetrahedronMesh {
    constructor(
        public vertices: Vector[],
        public tetrahedra: number[][],
        public epsilon: number = 0) {
    }

    getContainingTetrahedron(P: Vector): number {
        for (let t = 0; t < this.tetrahedra.length; ++t) {
            const bary = this.getBarycentrics(t, P);
            // A small tolerance keeps points on the shared faces of the
            // decomposition inside one of the tetrahedra.
            const tol = -1e-12;
            if (bary !== null && bary[0] >= tol && bary[1] >= tol
                && bary[2] >= tol && bary[3] >= tol) {
                return t;
            }
        }
        return -1;
    }

    getTetrahedronIndices(t: number): readonly number[] | null {
        return 0 <= t && t < this.tetrahedra.length ? this.tetrahedra[t] : null;
    }

    getBarycentrics(t: number, P: Vector): readonly number[] | null {
        const tet = this.tetrahedra[t];
        const result = computeBarycentrics3(P, this.vertices[tet[0]],
            this.vertices[tet[1]], this.vertices[tet[2]], this.vertices[tet[3]],
            this.epsilon);
        return result.valid ? result.bary : null;
    }
}

// The unit cube split into five tetrahedra.
function makeCubeMesh(): TestMesh {
    const vertices = [
        Vector.fromArray([0, 0, 0]),  // 0
        Vector.fromArray([1, 0, 0]),  // 1
        Vector.fromArray([1, 1, 0]),  // 2
        Vector.fromArray([0, 1, 0]),  // 3
        Vector.fromArray([0, 0, 1]),  // 4
        Vector.fromArray([1, 0, 1]),  // 5
        Vector.fromArray([1, 1, 1]),  // 6
        Vector.fromArray([0, 1, 1])   // 7
    ];
    const tetrahedra = [
        [0, 1, 2, 5],
        [0, 2, 3, 7],
        [0, 4, 5, 7],
        [2, 5, 6, 7],
        [0, 2, 5, 7]
    ];
    return new TestMesh(vertices, tetrahedra);
}

describe('IntpLinearNonuniform3', () => {
    it('throws when there are no samples', () => {
        expect(() => new IntpLinearNonuniform3(makeCubeMesh(), []))
            .toThrow('Invalid input.');
    });

    it('reproduces the samples at the vertices', () => {
        const mesh = makeCubeMesh();
        const F = [3, -1, 5, 2, 0, 7, -2, 1.5];
        const interp = new IntpLinearNonuniform3(mesh, F);
        for (let i = 0; i < mesh.vertices.length; ++i) {
            const result = interp.evaluate(mesh.vertices[i]);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(F[i], 14);
        }
    });

    it('reproduces an affine function exactly', () => {
        const mesh = makeCubeMesh();
        const f = (x: number, y: number, z: number) => 1 - 2 * x + 0.5 * y + 4 * z;
        const F = mesh.vertices.map(v => f(v.get(0), v.get(1), v.get(2)));
        const interp = new IntpLinearNonuniform3(mesh, F);
        for (let i = 0; i <= 5; ++i) {
            for (let j = 0; j <= 5; ++j) {
                for (let k = 0; k <= 5; ++k) {
                    const x = i / 5, y = j / 5, z = k / 5;
                    const result = interp.evaluate(Vector.fromArray([x, y, z]));
                    expect(result.valid).toBe(true);
                    expect(result.F).toBeCloseTo(f(x, y, z), 12);
                }
            }
        }
    });

    it('is a barycentric combination inside a tetrahedron', () => {
        const vertices = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([0, 1, 0]),
            Vector.fromArray([0, 0, 1])
        ];
        const mesh = new TestMesh(vertices, [[0, 1, 2, 3]]);
        const F = [4, -2, 6, 1];
        const interp = new IntpLinearNonuniform3(mesh, F);
        const P = Vector.fromArray([0.25, 0.25, 0.25]);
        const result = interp.evaluate(P);
        expect(result.valid).toBe(true);
        expect(result.F).toBeCloseTo((F[0] + F[1] + F[2] + F[3]) / 4, 12);
    });

    it('reports points outside the mesh as invalid', () => {
        const mesh = makeCubeMesh();
        const F = [3, -1, 5, 2, 0, 7, -2, 1.5];
        const interp = new IntpLinearNonuniform3(mesh, F);
        for (const P of [[-1, 0.5, 0.5], [2, 0.5, 0.5], [0.5, -1, 0.5],
            [0.5, 0.5, 3], [5, 5, 5]]) {
            expect(interp.evaluate(Vector.fromArray(P)).valid).toBe(false);
        }
    });

    it('reports a degenerate tetrahedron as invalid', () => {
        // The four vertices are coplanar.
        const vertices = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([0, 1, 0]),
            Vector.fromArray([1, 1, 0])
        ];
        const mesh = new TestMesh(vertices, [[0, 1, 2, 3]]);
        const interp = new IntpLinearNonuniform3(mesh, [1, 2, 3, 4]);
        expect(interp.evaluate(Vector.fromArray([0.25, 0.25, 0])).valid).toBe(false);

        const forced = Object.create(mesh) as TestMesh;
        forced.getContainingTetrahedron = () => 0;
        const interp2 = new IntpLinearNonuniform3(forced, [1, 2, 3, 4]);
        expect(interp2.evaluate(Vector.fromArray([0.25, 0.25, 0])).valid).toBe(false);
    });

    it('reports a missing index quadruple as invalid', () => {
        const mesh = makeCubeMesh();
        const broken = Object.create(mesh) as TestMesh;
        broken.getTetrahedronIndices = () => null;
        const interp = new IntpLinearNonuniform3(broken, [3, -1, 5, 2, 0, 7, -2, 1.5]);
        expect(interp.evaluate(Vector.fromArray([0.2, 0.1, 0.05])).valid).toBe(false);
    });

    it('is continuous across the shared face of two tetrahedra', () => {
        const mesh = makeCubeMesh();
        const F = [3, -1, 5, 2, 0, 7, -2, 1.5];
        const interp = new IntpLinearNonuniform3(mesh, F);
        // Tetrahedra 0 and 4 share the face <0, 2, 5>. Sample the face and
        // check that the value is the barycentric blend of F[0], F[2], F[5].
        const v0 = mesh.vertices[0], v2 = mesh.vertices[2], v5 = mesh.vertices[5];
        for (let i = 1; i <= 4; ++i) {
            for (let j = 1; i + j <= 5; ++j) {
                const b0 = i / 6, b1 = j / 6, b2 = 1 - b0 - b1;
                const P = Vector.fromArray([
                    b0 * v0.get(0) + b1 * v2.get(0) + b2 * v5.get(0),
                    b0 * v0.get(1) + b1 * v2.get(1) + b2 * v5.get(1),
                    b0 * v0.get(2) + b1 * v2.get(2) + b2 * v5.get(2)
                ]);
                const result = interp.evaluate(P);
                expect(result.valid).toBe(true);
                expect(result.F).toBeCloseTo(b0 * F[0] + b1 * F[2] + b2 * F[5], 12);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against the upstream
// IntpLinearNonuniform3.h and cross-checks with the real Delaunay3Mesh.
// ---------------------------------------------------------------------------

describe('IntpLinearNonuniform3 verification', () => {
    // Regression for the port defect found in V29: the duck-typed mesh
    // interface required getIndices(t), but Delaunay3Mesh (the mesh the
    // upstream header names as the source of the tetrahedralization) exposes
    // getIndices() for the whole flat index array plus
    // getTetrahedronIndices(t) for one tetrahedron. TypeScript's
    // parameter-arity assignability let the zero-argument accessor satisfy
    // the interface, so every query blended the samples of the FIRST
    // tetrahedron. On the old code the value below was 2.5, not 0.5.
    it('interpolates through the real Delaunay3Mesh', () => {
        const points = [
            Vector.fromArray([0, 0, 0]), Vector.fromArray([2, 0, 0]),
            Vector.fromArray([0, 2, 0]), Vector.fromArray([0, 0, 2]),
            Vector.fromArray([2, 2, 2]), Vector.fromArray([2, 2, 0]),
            Vector.fromArray([0, 2, 2]), Vector.fromArray([2, 0, 2])
        ];
        const delaunay = new Delaunay3();
        expect(delaunay.compute(points)).toBe(true);
        const mesh: IntpLinearNonuniform3TetrahedronMesh =
            new Delaunay3Mesh(delaunay);
        const f = (p: Vector) =>
            1 + 2 * p.values[0] - 3 * p.values[1] + 0.5 * p.values[2];
        const F = points.map(f);
        const interp = new IntpLinearNonuniform3(mesh, F);
        const P = Vector.fromArray([1, 1, 1]);
        const result = interp.evaluate(P);
        expect(result.valid).toBe(true);
        expectClose(result.F, f(P), 1e-12, 1e-12);
    });

    // Linear interpolation of a linear function over a tetrahedralization is
    // exact in exact arithmetic; the lattice generator keeps the vertices and
    // the sample values exactly representable.
    it('reproduces an affine function on a Delaunay tetrahedralization of lattice points',
        () => {
            let numValid = 0;
            check(fc.tuple(
                fc.array(latticeVector(3, -5, 5), { minLength: 8, maxLength: 14 }),
                fc.integer({ min: -4, max: 4 }),
                fc.integer({ min: -4, max: 4 }),
                fc.integer({ min: -4, max: 4 }),
                fc.integer({ min: -4, max: 4 }),
                fc.integer({ min: 0, max: 0xffff })
            ), ([pts, a, b, c, d, seed]) => {
                const delaunay = new Delaunay3();
                if (!delaunay.compute(pts) || delaunay.getDimension() !== 3) {
                    return true;   // degenerate (coplanar or duplicate) input
                }
                const mesh = new Delaunay3Mesh(delaunay);
                const f = (p: Vector) => a + b * p.values[0] + c * p.values[1]
                    + d * p.values[2];
                const F = mesh.getVertices().map(f);
                const interp = new IntpLinearNonuniform3(mesh, F);

                const rand = seededRandom(seed + 3);
                for (let n = 0; n < 15; ++n) {
                    const P = Vector.fromArray([-5 + 10 * rand(),
                        -5 + 10 * rand(), -5 + 10 * rand()]);
                    const result = interp.evaluate(P);
                    if (!result.valid) {
                        continue;
                    }
                    ++numValid;
                    // The barycentrics come from exact rational arithmetic
                    // rounded once, so the error is a few ulps of the largest
                    // term.
                    expectClose(result.F, f(P), 1e-9, 1e-12);
                }
                return true;
            }, 25);
            expect(numValid).toBeGreaterThan(20);
        }, 30000);

    // The interpolation must be exactly the barycentric combination of the
    // samples at the vertices of the containing tetrahedron.
    it('equals the barycentric combination of the containing tetrahedron', () => {
        let numValid = 0;
        check(fc.tuple(
            fc.array(latticeVector(3, -4, 4), { minLength: 8, maxLength: 12 }),
            fc.array(fc.integer({ min: -20, max: 20 }),
                { minLength: 12, maxLength: 12 }),
            fc.integer({ min: 0, max: 0xffff })
        ), ([pts, samples, seed]) => {
            const delaunay = new Delaunay3();
            if (!delaunay.compute(pts) || delaunay.getDimension() !== 3) {
                return true;
            }
            const mesh = new Delaunay3Mesh(delaunay);
            const F = mesh.getVertices().map((_, i) => samples[i % samples.length]);
            const interp = new IntpLinearNonuniform3(mesh, F);

            const rand = seededRandom(seed + 11);
            for (let n = 0; n < 10; ++n) {
                const P = Vector.fromArray([-4 + 8 * rand(), -4 + 8 * rand(),
                    -4 + 8 * rand()]);
                const t = mesh.getContainingTetrahedron(P);
                const result = interp.evaluate(P);
                if (t === -1) {
                    expect(result.valid).toBe(false);
                    continue;
                }
                const bary = mesh.getBarycentrics(t, P);
                const idx = mesh.getTetrahedronIndices(t);
                if (bary === null || idx === null) {
                    expect(result.valid).toBe(false);
                    continue;
                }
                expect(result.valid).toBe(true);
                ++numValid;
                const expected = bary[0] * F[idx[0]] + bary[1] * F[idx[1]]
                    + bary[2] * F[idx[2]] + bary[3] * F[idx[3]];
                expect(result.F).toBe(expected);
            }
            return true;
        }, 25);
        expect(numValid).toBeGreaterThan(20);
    }, 30000);

    // The upstream contract: the result is valid if and only if the point is
    // in the convex hull, so a point far outside is never valid.
    it('reports points far outside the hull as invalid', () => {
        check(fc.tuple(
            fc.array(latticeVector(3, -4, 4), { minLength: 8, maxLength: 12 }),
            unitVector(3)
        ), ([pts, dir]) => {
            const delaunay = new Delaunay3();
            if (!delaunay.compute(pts) || delaunay.getDimension() !== 3) {
                return true;
            }
            const mesh = new Delaunay3Mesh(delaunay);
            const F = mesh.getVertices().map((_, i) => i);
            const interp = new IntpLinearNonuniform3(mesh, F);
            // dir is unit length, so 'far' is 1000 units from the origin and
            // therefore well outside the hull of points in [-4,4]^3.
            const far = Vector.fromArray(dir.values.map(x => 1e3 * x));
            expect(interp.evaluate(far).valid).toBe(false);
            return true;
        }, 100);
    }, 30000);

    // The samples array is aliased, not copied (upstream stores Real const*).
    it('aliases the sample array rather than copying it', () => {
        const mesh = makeCubeMesh();
        const F = [1, 2, 3, 4, 5, 6, 7, 8];
        const interp = new IntpLinearNonuniform3(mesh, F);
        const P = Vector.fromArray([0.2, 0.3, 0.1]);
        const before = interp.evaluate(P);
        expect(before.valid).toBe(true);
        F[0] += 10;
        const after = interp.evaluate(P);
        expect(after.F).not.toBe(before.F);
    });
});
