import { describe, it, expect } from 'vitest';
import { IntpLinearNonuniform2 } from '../src/IntpLinearNonuniform2.js';
import type { IntpLinearNonuniform2TriangleMesh } from '../src/IntpLinearNonuniform2.js';
import { Vector } from '../src/Vector.js';
import { computeBarycentrics2 } from '../src/Vector2.js';
import { Delaunay2 } from '../src/Delaunay2.js';
import { Delaunay2Mesh } from '../src/Delaunay2Mesh.js';
import { check, expectClose, fc, latticeVector, seededRandom, unitVector }
    from './helpers/arbitraries.js';

// A minimal triangle mesh adapter that satisfies the interface required by
// IntpLinearNonuniform2. The triangles are index triples into 'vertices'.
class TestMesh implements IntpLinearNonuniform2TriangleMesh {
    constructor(
        public vertices: Vector[],
        public triangles: number[][],
        public epsilon: number = 0) {
    }

    getContainingTriangle(P: Vector): number {
        for (let t = 0; t < this.triangles.length; ++t) {
            const bary = this.getBarycentrics(t, P);
            if (bary !== null && bary[0] >= 0 && bary[1] >= 0 && bary[2] >= 0) {
                return t;
            }
        }
        return -1;
    }

    getTriangleIndices(t: number): readonly number[] | null {
        return 0 <= t && t < this.triangles.length ? this.triangles[t] : null;
    }

    getBarycentrics(t: number, P: Vector): readonly number[] | null {
        const tri = this.triangles[t];
        const result = computeBarycentrics2(P, this.vertices[tri[0]],
            this.vertices[tri[1]], this.vertices[tri[2]], this.epsilon);
        return result.valid ? result.bary : null;
    }
}

// The unit square split into two triangles.
function makeSquareMesh(): TestMesh {
    const vertices = [
        Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]),
        Vector.fromArray([1, 1]),
        Vector.fromArray([0, 1])
    ];
    return new TestMesh(vertices, [[0, 1, 2], [0, 2, 3]]);
}

describe('IntpLinearNonuniform2', () => {
    it('throws when there are no samples', () => {
        expect(() => new IntpLinearNonuniform2(makeSquareMesh(), []))
            .toThrow('Invalid input.');
    });

    it('reproduces the samples at the vertices', () => {
        const mesh = makeSquareMesh();
        const F = [3, -1, 5, 2];
        const interp = new IntpLinearNonuniform2(mesh, F);
        for (let i = 0; i < mesh.vertices.length; ++i) {
            const result = interp.evaluate(mesh.vertices[i]);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(F[i], 14);
        }
    });

    it('reproduces an affine function exactly', () => {
        const mesh = makeSquareMesh();
        const f = (x: number, y: number) => 2 - 3 * x + 0.5 * y;
        const F = mesh.vertices.map(v => f(v.get(0), v.get(1)));
        const interp = new IntpLinearNonuniform2(mesh, F);
        for (let i = 0; i <= 10; ++i) {
            for (let j = 0; j <= 10; ++j) {
                const x = i / 10, y = j / 10;
                const result = interp.evaluate(Vector.fromArray([x, y]));
                expect(result.valid).toBe(true);
                expect(result.F).toBeCloseTo(f(x, y), 12);
            }
        }
    });

    it('is a barycentric combination inside a triangle', () => {
        const mesh = makeSquareMesh();
        const F = [3, -1, 5, 2];
        const interp = new IntpLinearNonuniform2(mesh, F);
        // The centroid of triangle 0 = <V0, V1, V2>.
        const P = Vector.fromArray([2 / 3, 1 / 3]);
        const result = interp.evaluate(P);
        expect(result.valid).toBe(true);
        expect(result.F).toBeCloseTo((F[0] + F[1] + F[2]) / 3, 12);
    });

    it('reports points outside the mesh as invalid', () => {
        const mesh = makeSquareMesh();
        const F = [3, -1, 5, 2];
        const interp = new IntpLinearNonuniform2(mesh, F);
        for (const P of [[-1, 0.5], [2, 0.5], [0.5, -1], [0.5, 2], [1.5, 1.5]]) {
            const result = interp.evaluate(Vector.fromArray(P));
            expect(result.valid).toBe(false);
        }
    });

    it('reports a degenerate triangle as invalid', () => {
        // The three vertices are collinear, so the barycentric computation
        // fails and the interpolation is not valid.
        const vertices = [
            Vector.fromArray([0, 0]),
            Vector.fromArray([1, 1]),
            Vector.fromArray([2, 2])
        ];
        const mesh = new TestMesh(vertices, [[0, 1, 2]]);
        const interp = new IntpLinearNonuniform2(mesh, [1, 2, 3]);
        // getContainingTriangle returns -1 because the barycentrics fail.
        expect(interp.evaluate(Vector.fromArray([1, 1])).valid).toBe(false);

        // Force the containing-triangle lookup to succeed so that the
        // barycentric failure path inside evaluate is exercised.
        const forced = Object.create(mesh) as TestMesh;
        forced.getContainingTriangle = () => 0;
        const interp2 = new IntpLinearNonuniform2(forced, [1, 2, 3]);
        expect(interp2.evaluate(Vector.fromArray([1, 1])).valid).toBe(false);
    });

    it('reports a missing index triple as invalid', () => {
        const mesh = makeSquareMesh();
        const broken = Object.create(mesh) as TestMesh;
        broken.getTriangleIndices = () => null;
        const interp = new IntpLinearNonuniform2(broken, [3, -1, 5, 2]);
        expect(interp.evaluate(Vector.fromArray([0.5, 0.25])).valid).toBe(false);
    });

    it('is continuous across the shared edge of the two triangles', () => {
        const mesh = makeSquareMesh();
        const F = [3, -1, 5, 2];
        const interp = new IntpLinearNonuniform2(mesh, F);
        // The shared edge is the diagonal from V0 to V2. On the edge the
        // value is the linear blend of F[0] and F[2] regardless of which
        // triangle is used.
        for (let k = 0; k <= 10; ++k) {
            const s = k / 10;
            const result = interp.evaluate(Vector.fromArray([s, s]));
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo((1 - s) * F[0] + s * F[2], 12);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against the upstream
// IntpLinearNonuniform2.h and cross-checks with the real Delaunay2Mesh.
// ---------------------------------------------------------------------------

describe('IntpLinearNonuniform2 verification', () => {
    // Regression for the port defect found in V29: the duck-typed mesh
    // interface required getIndices(t), but Delaunay2Mesh and PlanarMesh (the
    // meshes the upstream header names as the sources of the triangulation)
    // expose getIndices() for the whole flat index array plus
    // getTriangleIndices(t) for one triangle. TypeScript's parameter-arity
    // assignability let the zero-argument accessor satisfy the interface, so
    // every query blended the samples of the FIRST triangle. On the old code
    // the value below was 0.5454..., not 0.
    it('interpolates through the real Delaunay2Mesh', () => {
        const points = [
            Vector.fromArray([0, 0]), Vector.fromArray([2, 0]),
            Vector.fromArray([2, 2]), Vector.fromArray([0, 2]),
            Vector.fromArray([1, 0.9])
        ];
        const delaunay = new Delaunay2();
        expect(delaunay.compute(points)).toBe(true);
        const mesh: IntpLinearNonuniform2TriangleMesh = new Delaunay2Mesh(delaunay);
        const f = (p: Vector) => 1 + 2 * p.values[0] - 3 * p.values[1];
        const F = points.map(f);
        const interp = new IntpLinearNonuniform2(mesh, F);
        const P = Vector.fromArray([1, 1]);
        const result = interp.evaluate(P);
        expect(result.valid).toBe(true);
        expectClose(result.F, f(P), 1e-12, 1e-12);
    });

    // Linear interpolation of a linear function over a triangulation is exact
    // in exact arithmetic. The lattice generator keeps the vertices and the
    // sample values exactly representable, so only the barycentric divide
    // introduces rounding.
    it('reproduces an affine function on a Delaunay triangulation of lattice points',
        () => {
            let numValid = 0;
            check(fc.tuple(
                fc.array(latticeVector(2, -6, 6), { minLength: 6, maxLength: 12 }),
                fc.integer({ min: -4, max: 4 }),
                fc.integer({ min: -4, max: 4 }),
                fc.integer({ min: -4, max: 4 }),
                fc.integer({ min: 0, max: 0xffff })
            ), ([pts, a, b, c, seed]) => {
                const delaunay = new Delaunay2();
                if (!delaunay.compute(pts) || delaunay.getDimension() !== 2) {
                    return true;   // degenerate (collinear or duplicate) input
                }
                const mesh = new Delaunay2Mesh(delaunay);
                const f = (p: Vector) => a + b * p.values[0] + c * p.values[1];
                const F = mesh.getVertices().map(f);
                const interp = new IntpLinearNonuniform2(mesh, F);

                const rand = seededRandom(seed + 1);
                for (let n = 0; n < 20; ++n) {
                    const P = Vector.fromArray(
                        [-6 + 12 * rand(), -6 + 12 * rand()]);
                    const result = interp.evaluate(P);
                    if (!result.valid) {
                        continue;
                    }
                    ++numValid;
                    // The barycentrics come from exact rational arithmetic
                    // rounded once, so the error is a few ulps of the largest
                    // term (|F| <= 4 + 4*6 + 4*6 = 52 here).
                    expectClose(result.F, f(P), 1e-9, 1e-12);
                }
                return true;
            }, 40);
            expect(numValid).toBeGreaterThan(50);
        }, 30000);

    // The interpolation must be exactly the barycentric combination of the
    // samples at the vertices of the containing triangle.
    it('equals the barycentric combination of the containing triangle', () => {
        let numValid = 0;
        check(fc.tuple(
            fc.array(latticeVector(2, -5, 5), { minLength: 6, maxLength: 10 }),
            fc.array(fc.integer({ min: -20, max: 20 }),
                { minLength: 10, maxLength: 10 }),
            fc.integer({ min: 0, max: 0xffff })
        ), ([pts, samples, seed]) => {
            const delaunay = new Delaunay2();
            if (!delaunay.compute(pts) || delaunay.getDimension() !== 2) {
                return true;
            }
            const mesh = new Delaunay2Mesh(delaunay);
            const F = mesh.getVertices().map((_, i) => samples[i % samples.length]);
            const interp = new IntpLinearNonuniform2(mesh, F);

            const rand = seededRandom(seed + 7);
            for (let n = 0; n < 15; ++n) {
                const P = Vector.fromArray([-5 + 10 * rand(), -5 + 10 * rand()]);
                const t = mesh.getContainingTriangle(P);
                const result = interp.evaluate(P);
                if (t === -1) {
                    expect(result.valid).toBe(false);
                    continue;
                }
                const bary = mesh.getBarycentrics(t, P);
                const idx = mesh.getTriangleIndices(t);
                if (bary === null || idx === null) {
                    expect(result.valid).toBe(false);
                    continue;
                }
                expect(result.valid).toBe(true);
                ++numValid;
                const expected = bary[0] * F[idx[0]] + bary[1] * F[idx[1]]
                    + bary[2] * F[idx[2]];
                expect(result.F).toBe(expected);
            }
            return true;
        }, 40);
        expect(numValid).toBeGreaterThan(50);
    }, 30000);

    // The upstream contract: the result is valid if and only if the point is
    // in the convex hull, so a point far outside is never valid.
    it('reports points far outside the hull as invalid', () => {
        check(fc.tuple(
            fc.array(latticeVector(2, -4, 4), { minLength: 6, maxLength: 10 }),
            unitVector(2)
        ), ([pts, dir]) => {
            const delaunay = new Delaunay2();
            if (!delaunay.compute(pts) || delaunay.getDimension() !== 2) {
                return true;
            }
            const mesh = new Delaunay2Mesh(delaunay);
            const F = mesh.getVertices().map((_, i) => i);
            const interp = new IntpLinearNonuniform2(mesh, F);
            // dir is unit length, so 'far' is 1000 units from the origin and
            // therefore well outside the hull of points in [-4,4]^2.
            const far = Vector.fromArray(
                [1e3 * dir.values[0], 1e3 * dir.values[1]]);
            expect(interp.evaluate(far).valid).toBe(false);
            return true;
        });
    });

    // The samples array is aliased, not copied (upstream stores Real const*).
    it('aliases the sample array rather than copying it', () => {
        const mesh = makeSquareMesh();
        const F = [3, -1, 5, 2];
        const interp = new IntpLinearNonuniform2(mesh, F);
        const before = interp.evaluate(Vector.fromArray([0.5, 0.25]));
        F[0] += 10;
        const after = interp.evaluate(Vector.fromArray([0.5, 0.25]));
        expect(after.F).not.toBe(before.F);
    });
});
