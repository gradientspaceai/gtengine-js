import { describe, it, expect } from 'vitest';
import { IntpLinearNonuniform3 } from '../src/IntpLinearNonuniform3';
import type { IntpLinearNonuniform3TetrahedronMesh } from '../src/IntpLinearNonuniform3';
import { Vector } from '../src/Vector';
import { computeBarycentrics3 } from '../src/Vector3';

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

    getIndices(t: number): readonly number[] | null {
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
        broken.getIndices = () => null;
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
