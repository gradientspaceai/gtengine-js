import { describe, it, expect } from 'vitest';
import { RectanglePatchMesh } from '../src/RectanglePatchMesh';
import { MeshDescription, MeshTopology } from '../src/Mesh';
import { IndexAttribute } from '../src/IndexAttribute';
import { VertexAttribute } from '../src/VertexAttribute';
import { ParametricSurface } from '../src/ParametricSurface';
import { Vector, dot, length, normalize, sub } from '../src/Vector';
import { cross, unitCross } from '../src/Vector3';

function V(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A surface whose partial derivatives are supplied in closed form. The
// upstream RectanglePatchMesh accepts any rectangular ParametricSurface<3,T>.
class FunctionSurface extends ParametricSurface {
    private mF: (u: number, v: number) => number[][];

    constructor(umin: number, umax: number, vmin: number, vmax: number,
        rectangular: boolean, f: (u: number, v: number) => number[][]) {
        super(3, umin, umax, vmin, vmax, rectangular);
        this.mF = f;
        this.mConstructed = true;
    }

    evaluate(u: number, v: number, order: number, jet: Vector[]): void {
        const values = this.mF(u, v);
        const maxOrder = Math.min(order === 0 ? 0 : (order === 1 ? 2 : 5),
            values.length - 1);
        for (let i = 0; i <= maxOrder; ++i) {
            for (let k = 0; k < 3; ++k) {
                jet[i].values[k] = values[i][k];
            }
        }
        for (let i = maxOrder + 1; i < ParametricSurface.SUP_ORDER; ++i) {
            jet[i].makeZero();
        }
    }
}

// The plane X(u,v) = C + u*A0 + v*A1.
function planeSurface(C: number[], A0: number[], A1: number[],
    umin: number, umax: number, vmin: number, vmax: number): FunctionSurface {
    return new FunctionSurface(umin, umax, vmin, vmax, true, (u, v) => [
        [C[0] + u * A0[0] + v * A1[0], C[1] + u * A0[1] + v * A1[1],
            C[2] + u * A0[2] + v * A1[2]],
        A0, A1, [0, 0, 0], [0, 0, 0], [0, 0, 0]
    ]);
}

// The paraboloid X(u,v) = (u, v, s*(u^2+v^2)); 's' is mutable so the dynamic
// update path can be exercised.
class Paraboloid extends ParametricSurface {
    scale: number;

    constructor(scale: number) {
        super(3, -1, 1, -1, 1, true);
        this.scale = scale;
        this.mConstructed = true;
    }

    evaluate(u: number, v: number, order: number, jet: Vector[]): void {
        const s = this.scale;
        jet[0].values[0] = u;
        jet[0].values[1] = v;
        jet[0].values[2] = s * (u * u + v * v);
        if (order >= 1) {
            jet[1].values[0] = 1;
            jet[1].values[1] = 0;
            jet[1].values[2] = 2 * s * u;
            jet[2].values[0] = 0;
            jet[2].values[1] = 1;
            jet[2].values[2] = 2 * s * v;
        }
    }
}

interface Storage {
    description: MeshDescription;
    indices: Uint32Array;
    positions: Float64Array;
    normals: Float64Array | null;
    tcoords: Float64Array | null;
    tangents: Float64Array | null;
    bitangents: Float64Array | null;
    dpdus: Float64Array | null;
    dpdvs: Float64Array | null;
}

function makeStorage(numRows: number, numCols: number, options: {
    normals?: boolean, tcoords?: boolean, frame?: boolean
} = {}): Storage {
    const description = new MeshDescription(MeshTopology.RECTANGLE,
        numRows, numCols);
    const indices = new Uint32Array(3 * description.numTriangles);
    description.indexAttribute = new IndexAttribute(indices, 4);

    const positions = new Float64Array(3 * description.numVertices);
    description.vertexAttributes.push(
        new VertexAttribute('position', positions, 24));

    let normals: Float64Array | null = null;
    if (options.normals !== false) {
        normals = new Float64Array(3 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('normal', normals, 24));
    }

    let tcoords: Float64Array | null = null;
    if (options.tcoords) {
        tcoords = new Float64Array(2 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('tcoord', tcoords, 16));
    }

    let tangents: Float64Array | null = null;
    let bitangents: Float64Array | null = null;
    let dpdus: Float64Array | null = null;
    let dpdvs: Float64Array | null = null;
    if (options.frame) {
        tangents = new Float64Array(3 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('tangent', tangents, 24));
        bitangents = new Float64Array(3 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('bitangent', bitangents, 24));
        dpdus = new Float64Array(3 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('dpdu', dpdus, 24));
        dpdvs = new Float64Array(3 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('dpdv', dpdvs, 24));
        description.wantDynamicTangentSpaceUpdate = true;
    }

    return { description, indices, positions, normals, tcoords, tangents,
        bitangents, dpdus, dpdvs };
}

function get(buffer: Float64Array, i: number): Vector {
    return V(buffer[3 * i], buffer[3 * i + 1], buffer[3 * i + 2]);
}

function P(storage: Storage, i: number): Vector {
    return get(storage.positions, i);
}

describe('RectanglePatchMesh', () => {
    it('samples the parameter domain uniformly into the texture coordinates',
        () => {
            const storage = makeStorage(3, 5, { tcoords: true });
            const surface = planeSurface([1, 2, 3], [1, 0, 0], [0, 1, 0],
                -1, 3, 10, 20);
            const mesh = new RectanglePatchMesh(storage.description, surface);
            expect(mesh.getDescription().numVertices).toBe(15);
            expect(mesh.getSurface()).toBe(surface);

            const tcoords = storage.tcoords as Float64Array;
            for (let r = 0, i = 0; r < 3; ++r) {
                for (let c = 0; c < 5; ++c, ++i) {
                    const u = -1 + (3 - -1) * c / 4;
                    const v = 10 + (20 - 10) * r / 2;
                    expect(tcoords[2 * i]).toBeCloseTo(u, 12);
                    expect(tcoords[2 * i + 1]).toBeCloseTo(v, 12);

                    // The position is the surface evaluated at (u,v).
                    const expected = V(1 + u, 2 + v, 3);
                    const actual = P(storage, i);
                    for (let k = 0; k < 3; ++k) {
                        expect(actual.values[k])
                            .toBeCloseTo(expected.values[k], 12);
                    }
                }
            }
        });

    it('computes normals from the surface derivatives', () => {
        const storage = makeStorage(5, 5, { tcoords: true });
        const surface = new Paraboloid(1);
        const mesh = new RectanglePatchMesh(storage.description, surface);

        const normals = storage.normals as Float64Array;
        const tcoords = storage.tcoords as Float64Array;
        for (let i = 0; i < mesh.getDescription().numVertices; ++i) {
            const u = tcoords[2 * i];
            const v = tcoords[2 * i + 1];
            // The unit normal of z = u^2+v^2 is proportional to (-2u,-2v,1).
            const expected = V(-2 * u, -2 * v, 1);
            normalize(expected);
            const actual = get(normals, i);
            expect(length(actual)).toBeCloseTo(1, 12);
            for (let k = 0; k < 3; ++k) {
                expect(actual.values[k]).toBeCloseTo(expected.values[k], 10);
            }
            // The position is on the paraboloid.
            const p = P(storage, i);
            expect(p.values[2]).toBeCloseTo(u * u + v * v, 12);
        }
    });

    it('builds a tangent-space frame that is orthonormal and consistent', () => {
        const storage = makeStorage(4, 4, { frame: true, tcoords: true });
        const surface = new Paraboloid(0.5);
        const mesh = new RectanglePatchMesh(storage.description, surface);
        expect(mesh.getDescription().allowUpdateFrame).toBe(true);

        const tcoords = storage.tcoords as Float64Array;
        const tangents = storage.tangents as Float64Array;
        const bitangents = storage.bitangents as Float64Array;
        const normals = storage.normals as Float64Array;
        const dpdus = storage.dpdus as Float64Array;
        const dpdvs = storage.dpdvs as Float64Array;
        for (let i = 0; i < mesh.getDescription().numVertices; ++i) {
            const u = tcoords[2 * i];
            const v = tcoords[2 * i + 1];
            const T = get(tangents, i);
            const B = get(bitangents, i);
            const N = get(normals, i);

            // The frame is right-handed and orthonormal.
            expect(length(T)).toBeCloseTo(1, 12);
            expect(length(B)).toBeCloseTo(1, 12);
            expect(length(N)).toBeCloseTo(1, 12);
            expect(dot(T, B)).toBeCloseTo(0, 12);
            expect(dot(T, N)).toBeCloseTo(0, 12);
            expect(dot(B, N)).toBeCloseTo(0, 12);
            const rightHanded = cross(T, B);
            for (let k = 0; k < 3; ++k) {
                expect(rightHanded.values[k]).toBeCloseTo(N.values[k], 12);
            }

            // The dpdu/dpdv channels hold the normalized surface derivatives,
            // recorded before the orthogonalization.
            const dpdu = V(1, 0, 2 * 0.5 * u);
            normalize(dpdu);
            const dpdv = V(0, 1, 2 * 0.5 * v);
            normalize(dpdv);
            for (let k = 0; k < 3; ++k) {
                expect(get(dpdus, i).values[k])
                    .toBeCloseTo(dpdu.values[k], 12);
                expect(get(dpdvs, i).values[k])
                    .toBeCloseTo(dpdv.values[k], 12);
            }

            // The tangent is the normalized dX/du and the normal agrees with
            // the cross product of the derivatives.
            for (let k = 0; k < 3; ++k) {
                expect(T.values[k]).toBeCloseTo(dpdu.values[k], 12);
            }
            const expectedNormal = unitCross(dpdu, dpdv);
            for (let k = 0; k < 3; ++k) {
                expect(N.values[k]).toBeCloseTo(expectedNormal.values[k], 12);
            }
        }
    });

    it('recomputes positions, normals and frame in update()', () => {
        const storage = makeStorage(4, 4, { tcoords: true });
        const surface = new Paraboloid(1);
        const mesh = new RectanglePatchMesh(storage.description, surface);
        expect(P(storage, 0).values[2]).toBeCloseTo(2, 12);

        surface.scale = 3;
        mesh.update();
        const tcoords = storage.tcoords as Float64Array;
        const normals = storage.normals as Float64Array;
        for (let i = 0; i < mesh.getDescription().numVertices; ++i) {
            const u = tcoords[2 * i];
            const v = tcoords[2 * i + 1];
            expect(P(storage, i).values[2]).toBeCloseTo(3 * (u * u + v * v), 12);
            const expected = V(-6 * u, -6 * v, 1);
            normalize(expected);
            const actual = get(normals, i);
            for (let k = 0; k < 3; ++k) {
                expect(actual.values[k]).toBeCloseTo(expected.values[k], 10);
            }
        }

        // The frame path of update().
        const frameStorage = makeStorage(3, 3, { frame: true });
        const frameSurface = new Paraboloid(1);
        const frameMesh = new RectanglePatchMesh(frameStorage.description,
            frameSurface);
        frameSurface.scale = 2;
        frameMesh.update();
        const T = get(frameStorage.tangents as Float64Array, 0);
        const dpdu = V(1, 0, 2 * 2 * -1);
        normalize(dpdu);
        for (let k = 0; k < 3; ++k) {
            expect(T.values[k]).toBeCloseTo(dpdu.values[k], 12);
        }
    });

    it('emits triangles that tile the patch consistently', () => {
        const storage = makeStorage(4, 6);
        const surface = planeSurface([0, 0, 0], [1, 0, 0], [0, 1, 0],
            0, 2, 0, 3);
        const mesh = new RectanglePatchMesh(storage.description, surface);

        expect(storage.indices.length)
            .toBe(3 * mesh.getDescription().numTriangles);
        let area = 0;
        for (let t = 0; t < mesh.getDescription().numTriangles; ++t) {
            const a = P(storage, storage.indices[3 * t]);
            const b = P(storage, storage.indices[3 * t + 1]);
            const c = P(storage, storage.indices[3 * t + 2]);
            const n = cross(sub(b, a), sub(c, a));
            expect(n.values[2]).toBeGreaterThan(0);
            area += 0.5 * length(n);
        }
        expect(area).toBeCloseTo(2 * 3, 10);
    });

    it('requires a nonnull rectangular 3D surface', () => {
        const storage = makeStorage(3, 3);
        const nonRectangular = new FunctionSurface(0, 1, 0, 1, false,
            (u, v) => [[u, v, 0], [1, 0, 0], [0, 1, 0], [0, 0, 0], [0, 0, 0],
                [0, 0, 0]]);
        expect(() => new RectanglePatchMesh(storage.description,
            nonRectangular)).toThrow();
    });
});
