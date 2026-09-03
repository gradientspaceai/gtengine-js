import { describe, it, expect } from 'vitest';
import { Tetrahedron3 } from '../src/Tetrahedron3.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross, dotCross, unitCross } from '../src/Vector3.js';

function V(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The signed volume of the tetrahedron, computed independently of the class.
function signedVolume(t: Tetrahedron3): number {
    return dotCross(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0]),
        sub(t.v[3], t.v[0])) / 6;
}

// The outward unit normal of face 'face', computed independently from the
// face index table and a point-inside test.
function outwardFaceNormal(t: Tetrahedron3, face: number): Vector {
    const indices = Tetrahedron3.getFaceIndices(face);
    const normal = unitCross(sub(t.v[indices[1]], t.v[indices[0]]),
        sub(t.v[indices[2]], t.v[indices[0]]));
    const centroid = t.computeCentroid();
    const inward = dot(normal, sub(centroid, t.v[indices[0]]));
    return (inward <= 0 ? normal : mul(normal, -1));
}

// The unnormalized (area-weighted) face normal from the winding of the face.
function rawFaceNormal(t: Tetrahedron3, face: number): Vector {
    const indices = Tetrahedron3.getFaceIndices(face);
    return cross(sub(t.v[indices[1]], t.v[indices[0]]),
        sub(t.v[indices[2]], t.v[indices[0]]));
}

const s3 = 1 / Math.sqrt(3);

describe('Tetrahedron3', () => {
    it('default-constructs the canonical tetrahedron', () => {
        const t = new Tetrahedron3();
        expect(t.v[0].values).toEqual([0, 0, 0]);
        expect(t.v[1].values).toEqual([1, 0, 0]);
        expect(t.v[2].values).toEqual([0, 1, 0]);
        expect(t.v[3].values).toEqual([0, 0, 1]);
        expect(signedVolume(t)).toBeCloseTo(1 / 6, 15);
    });

    it('copies the vertices in the factories and in clone', () => {
        const p = [V(1, 2, 3), V(4, 5, 6), V(7, 8, 10), V(-1, 0, 2)];
        const t = Tetrahedron3.fromVertices(p[0], p[1], p[2], p[3]);
        p[0].values[0] = 100;
        expect(t.v[0].values).toEqual([1, 2, 3]);

        const u = Tetrahedron3.fromArray(t.v);
        const w = t.clone();
        t.v[2].values[1] = -50;
        expect(u.v[2].values).toEqual([7, 8, 10]);
        expect(w.v[2].values).toEqual([7, 8, 10]);

        expect(() => Tetrahedron3.fromArray([V(0, 0, 0)])).toThrow();
        expect(() => Tetrahedron3.fromArray(
            [V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), Vector.fromArray([0, 0])]))
            .toThrow();
    });

    it('reports the documented face, edge and augmented index tables', () => {
        expect(Tetrahedron3.getFaceIndices(0)).toEqual([0, 2, 1]);
        expect(Tetrahedron3.getFaceIndices(1)).toEqual([0, 1, 3]);
        expect(Tetrahedron3.getFaceIndices(2)).toEqual([0, 3, 2]);
        expect(Tetrahedron3.getFaceIndices(3)).toEqual([1, 2, 3]);
        expect(Tetrahedron3.getAllFaceIndices()).toEqual(
            [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3]);
        expect(Tetrahedron3.getAllEdgeIndices()).toEqual(
            [0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3]);
        for (let e = 0; e < 6; ++e) {
            const pair = Tetrahedron3.getEdgeIndices(e);
            expect(Tetrahedron3.getAllEdgeIndices().slice(2 * e, 2 * e + 2))
                .toEqual([pair[0], pair[1]]);
        }
        expect(() => Tetrahedron3.getFaceIndices(4)).toThrow();
        expect(() => Tetrahedron3.getEdgeIndices(6)).toThrow();
        expect(() => Tetrahedron3.getEdgeAugmented(-1)).toThrow();
        expect(() => Tetrahedron3.getVertexAugmented(4)).toThrow();
    });

    it('has consistent augmented index tables', () => {
        // getEdgeAugmented(e) = {v0,v1,v2,v3} where the edge is {v0,v1} and
        // the triangles sharing the edge are {v0,v2,v1} and {v0,v1,v3}. Both
        // triangles must be faces of the tetrahedron (as index sets).
        const faceSets = [0, 1, 2, 3].map(f =>
            Array.from(Tetrahedron3.getFaceIndices(f)).sort().join(','));
        for (let e = 0; e < 6; ++e) {
            const a = Tetrahedron3.getEdgeAugmented(e);
            expect(new Set(a).size).toBe(4);
            const edge = Tetrahedron3.getEdgeIndices(e);
            expect([a[0], a[1]].sort()).toEqual([edge[0], edge[1]].sort());
            for (const tri of [[a[0], a[2], a[1]], [a[0], a[1], a[3]]]) {
                expect(faceSets).toContain(tri.slice().sort().join(','));
            }
        }

        // getVertexAugmented(v) = {v0,v1,v2,v3} where the vertex is v0 and the
        // triangles sharing it are {v0,v1,v2}, {v0,v2,v3} and {v0,v3,v1}.
        for (let v = 0; v < 4; ++v) {
            const a = Tetrahedron3.getVertexAugmented(v);
            expect(new Set(a).size).toBe(4);
            expect(a[0]).toBe(v);
            for (const tri of [[a[0], a[1], a[2]], [a[0], a[2], a[3]],
                [a[0], a[3], a[1]]]) {
                expect(faceSets).toContain(tri.slice().sort().join(','));
            }
        }
    });

    it('computes the face normals of the canonical tetrahedron', () => {
        const t = new Tetrahedron3();
        expect(t.computeFaceNormal(0).values).toEqual([0, 0, -1]);
        expect(t.computeFaceNormal(1).values).toEqual([0, -1, 0]);
        expect(t.computeFaceNormal(2).values).toEqual([-1, 0, 0]);
        const n3 = t.computeFaceNormal(3);
        expect(n3.values[0]).toBeCloseTo(s3, 15);
        expect(n3.values[1]).toBeCloseTo(s3, 15);
        expect(n3.values[2]).toBeCloseTo(s3, 15);
    });

    it('computes edge and vertex normals as normalized face-normal sums', () => {
        const t = Tetrahedron3.fromVertices(
            V(0.5, -1, 2), V(3, 0.25, 1), V(1, 1, -3), V(-1, 2, 0.5));
        // The vertex ordering is positively oriented, so the winding-based
        // normals of computeFaceNormal/computeEdgeNormal/computeVertexNormal
        // point outward.
        expect(signedVolume(t)).toBeGreaterThan(0);

        // The vertex normal is the negative of the opposite face's normal.
        // The face opposite vertex k is the face whose index set omits k.
        for (let v = 0; v < 4; ++v) {
            let opposite = -1;
            for (let f = 0; f < 4; ++f) {
                if (!Tetrahedron3.getFaceIndices(f).includes(v)) {
                    opposite = f;
                }
            }
            const expected = mul(outwardFaceNormal(t, opposite), -1);
            const actual = t.computeVertexNormal(v);
            for (let k = 0; k < 3; ++k) {
                expect(actual.values[k]).toBeCloseTo(expected.values[k], 12);
            }
        }

        // The edge normal is the normalized sum of the two adjacent
        // unnormalized (that is, area-weighted) face normals.
        for (let e = 0; e < 6; ++e) {
            const edge = Tetrahedron3.getEdgeIndices(e);
            const adjacent: number[] = [];
            for (let f = 0; f < 4; ++f) {
                const fi = Tetrahedron3.getFaceIndices(f);
                if (fi.includes(edge[0]) && fi.includes(edge[1])) {
                    adjacent.push(f);
                }
            }
            expect(adjacent.length).toBe(2);
            const sum = add(rawFaceNormal(t, adjacent[0]),
                rawFaceNormal(t, adjacent[1]));
            const expected = mul(sum, 1 / length(sum));
            const actual = t.computeEdgeNormal(e);
            for (let k = 0; k < 3; ++k) {
                expect(actual.values[k]).toBeCloseTo(expected.values[k], 12);
            }
        }
    });

    it('computes the centroid as the average of the vertices', () => {
        const t = Tetrahedron3.fromVertices(
            V(1, 0, 0), V(3, 2, -1), V(0, 4, 5), V(-4, 2, 8));
        expect(t.computeCentroid().values).toEqual([0, 2, 3]);
    });

    it('builds outward-pointing face planes for both orientations', () => {
        // A positively oriented tetrahedron and a negatively oriented one
        // (obtained by swapping two vertices).
        const positive = new Tetrahedron3();
        const negative = Tetrahedron3.fromVertices(V(0, 0, 0), V(1, 0, 0),
            V(0, 0, 1), V(0, 1, 0));
        expect(signedVolume(positive)).toBeGreaterThan(0);
        expect(signedVolume(negative)).toBeLessThan(0);

        for (const t of [positive, negative]) {
            const plane = t.getPlanes();
            expect(plane.length).toBe(4);
            for (let f = 0; f < 4; ++f) {
                // Every vertex is in the negative halfspace of the plane, so
                // the normals point outward.
                for (let j = 0; j < 4; ++j) {
                    expect(dot(plane[f].normal, t.v[j]))
                        .toBeLessThanOrEqual(plane[f].constant + 1e-12);
                }
                // The three vertices of the face are on the plane.
                for (const index of Tetrahedron3.getFaceIndices(f)) {
                    expect(dot(plane[f].normal, t.v[index]))
                        .toBeCloseTo(plane[f].constant, 12);
                }
                expect(length(plane[f].normal)).toBeCloseTo(1, 12);
                // The port's origin is a point of the plane.
                expect(dot(plane[f].normal, plane[f].origin))
                    .toBeCloseTo(plane[f].constant, 12);
            }
        }

        // The canonical tetrahedron's planes have known coefficients.
        const plane = positive.getPlanes();
        expect(plane[0].normal.values).toEqual([0, 0, -1]);
        expect(plane[0].constant).toBe(0);
        expect(plane[1].normal.values).toEqual([0, -1, 0]);
        expect(plane[1].constant).toBe(0);
        expect(plane[2].normal.values).toEqual([-1, 0, 0]);
        expect(plane[2].constant).toBe(0);
        expect(plane[3].constant).toBeCloseTo(s3, 15);
    });

    it('agrees with an independent computation on random tetrahedra', () => {
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const R = (): number => 4 * rand() - 2;

        for (let trial = 0; trial < 60; ++trial) {
            const t = Tetrahedron3.fromVertices(V(R(), R(), R()),
                V(R(), R(), R()), V(R(), R(), R()), V(R(), R(), R()));
            if (Math.abs(signedVolume(t)) < 1e-3) {
                continue;
            }

            // The centroid is strictly inside, so it is in the negative
            // halfspace of every face plane.
            const centroid = t.computeCentroid();
            const plane = t.getPlanes();
            for (let f = 0; f < 4; ++f) {
                expect(dot(plane[f].normal, centroid))
                    .toBeLessThan(plane[f].constant);
                const expected = outwardFaceNormal(t, f);
                for (let k = 0; k < 3; ++k) {
                    expect(plane[f].normal.values[k])
                        .toBeCloseTo(expected.values[k], 10);
                }
                // computeFaceNormal uses the raw winding, which is the
                // outward normal exactly when the volume is positive.
                const sign = (signedVolume(t) > 0 ? 1 : -1);
                const faceNormal = t.computeFaceNormal(f);
                for (let k = 0; k < 3; ++k) {
                    expect(sign * faceNormal.values[k])
                        .toBeCloseTo(expected.values[k], 10);
                }
            }

            // The four outward face areas sum to zero (Minkowski's identity).
            let sum = new Vector(3);
            for (let f = 0; f < 4; ++f) {
                const fi = Tetrahedron3.getFaceIndices(f);
                const area = 0.5 * length(cross(
                    sub(t.v[fi[1]], t.v[fi[0]]), sub(t.v[fi[2]], t.v[fi[0]])));
                sum = add(sum, mul(outwardFaceNormal(t, f), area));
            }
            expect(length(sum)).toBeLessThan(1e-10);
        }
    });

    it('supports the lexicographic comparisons', () => {
        const a = new Tetrahedron3();
        const b = new Tetrahedron3();
        expect(a.equals(b)).toBe(true);
        expect(a.notEquals(b)).toBe(false);
        expect(a.lessThan(b)).toBe(false);
        expect(a.lessThanOrEqual(b)).toBe(true);
        expect(a.greaterThanOrEqual(b)).toBe(true);
        expect(a.greaterThan(b)).toBe(false);

        b.v[2].values[0] = 1;
        expect(a.equals(b)).toBe(false);
        expect(a.lessThan(b)).toBe(true);
        expect(b.greaterThan(a)).toBe(true);
        expect(b.lessThanOrEqual(a)).toBe(false);
    });
});
