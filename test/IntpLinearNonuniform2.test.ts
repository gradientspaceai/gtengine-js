import { describe, it, expect } from 'vitest';
import { IntpLinearNonuniform2 } from '../src/IntpLinearNonuniform2';
import type { IntpLinearNonuniform2TriangleMesh } from '../src/IntpLinearNonuniform2';
import { Vector } from '../src/Vector';
import { computeBarycentrics2 } from '../src/Vector2';

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

    getIndices(t: number): readonly number[] | null {
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
        broken.getIndices = () => null;
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
