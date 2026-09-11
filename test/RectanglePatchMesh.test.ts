import { describe, it, expect } from 'vitest';
import { RectanglePatchMesh } from '../src/RectanglePatchMesh.js';
import { MeshDescription, MeshTopology } from '../src/Mesh.js';
import { IndexAttribute } from '../src/IndexAttribute.js';
import { VertexAttribute } from '../src/VertexAttribute.js';
import { ParametricSurface } from '../src/ParametricSurface.js';
import { Vector, dot, length, normalize, sub } from '../src/Vector.js';
import { cross, unitCross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// V43 verification: the sampled grid against the surface parametrization, the
// normals and frame against the surface derivatives, and the update paths.
// ---------------------------------------------------------------------------

describe('RectanglePatchMesh verification', () => {
    const gridArb = fc.tuple(fc.integer({ min: 2, max: 6 }),
        fc.integer({ min: 2, max: 6 }));

    // A paraboloid with a random curvature; its partial derivatives are known
    // in closed form, so every channel can be checked independently.
    const scaleArb = finite(-2, 2);

    function analytic(scale: number, u: number, v: number) {
        const position = V(u, v, scale * (u * u + v * v));
        const dpdu = V(1, 0, 2 * scale * u);
        const dpdv = V(0, 1, 2 * scale * v);
        return { position, dpdu, dpdv };
    }

    it('samples the domain uniformly and evaluates the surface there', () => {
        check(fc.tuple(scaleArb, gridArb), ([scale, [nr, nc]]) => {
            const surface = new Paraboloid(scale);
            const storage = makeStorage(nr, nc, { tcoords: true });
            const mesh = new RectanglePatchMesh(storage.description, surface);
            expect(mesh.getSurface()).toBe(surface);
            expect(mesh.getDescription().numVertices).toBe(nr * nc);

            const uMin = surface.getUMin();
            const uDelta = (surface.getUMax() - uMin) / (nc - 1);
            const vMin = surface.getVMin();
            const vDelta = (surface.getVMax() - vMin) / (nr - 1);
            const tcoords = storage.tcoords as Float64Array;
            for (let r = 0, i = 0; r < nr; ++r) {
                for (let c = 0; c < nc; ++c, ++i) {
                    // The columns walk u and the rows walk v.
                    const u = uMin + uDelta * c;
                    const v = vMin + vDelta * r;
                    expectClose(tcoords[2 * i], u, 0, 0);
                    expectClose(tcoords[2 * i + 1], v, 0, 0);
                    expectVectorClose(P(storage, i),
                        analytic(scale, u, v).position, 1e-12, 1e-12);
                }
            }
            // The grid spans the whole parameter rectangle.
            expectClose(tcoords[0], surface.getUMin(), 0, 0);
            expectClose(tcoords[1], surface.getVMin(), 0, 0);
            expectClose(tcoords[2 * (nr * nc - 1)], surface.getUMax(),
                1e-14, 1e-14);
            expectClose(tcoords[2 * (nr * nc - 1) + 1], surface.getVMax(),
                1e-14, 1e-14);
        }, 100);
    });

    it('computes normals as the unit cross of the surface derivatives', () => {
        check(fc.tuple(scaleArb, gridArb), ([scale, [nr, nc]]) => {
            const surface = new Paraboloid(scale);
            const storage = makeStorage(nr, nc);
            const mesh = new RectanglePatchMesh(storage.description, surface);
            // Without a tangent-space request, InitializeNormals runs.
            expect(mesh.getDescription().allowUpdateFrame).toBe(false);

            const uMin = surface.getUMin();
            const uDelta = (surface.getUMax() - uMin) / (nc - 1);
            const vMin = surface.getVMin();
            const vDelta = (surface.getVMax() - vMin) / (nr - 1);
            for (let r = 0, i = 0; r < nr; ++r) {
                for (let c = 0; c < nc; ++c, ++i) {
                    const { dpdu, dpdv } = analytic(scale, uMin + uDelta * c,
                        vMin + vDelta * r);
                    const expected = unitCross(dpdu, dpdv, true);
                    const n = get(storage.normals as Float64Array, i);
                    expectVectorClose(n, expected, 1e-12, 1e-12);
                    expectClose(length(n), 1, 1e-12, 1e-12);
                    // The normal is orthogonal to the tangent plane.
                    expectClose(dot(n, dpdu), 0, 1e-12, 1e-12);
                    expectClose(dot(n, dpdv), 0, 1e-12, 1e-12);
                }
            }
        }, 100);
    });

    it('builds an orthonormal frame from the surface derivatives', () => {
        check(fc.tuple(scaleArb, gridArb), ([scale, [nr, nc]]) => {
            const surface = new Paraboloid(scale);
            const storage = makeStorage(nr, nc, { frame: true });
            const mesh = new RectanglePatchMesh(storage.description, surface);
            expect(mesh.getDescription().allowUpdateFrame).toBe(true);

            const uMin = surface.getUMin();
            const uDelta = (surface.getUMax() - uMin) / (nc - 1);
            const vMin = surface.getVMin();
            const vDelta = (surface.getVMax() - vMin) / (nr - 1);
            for (let r = 0, i = 0; r < nr; ++r) {
                for (let c = 0; c < nc; ++c, ++i) {
                    const { dpdu, dpdv } = analytic(scale, uMin + uDelta * c,
                        vMin + vDelta * r);
                    const unitU = dpdu.clone();
                    normalize(unitU, true);
                    const unitV = dpdv.clone();
                    normalize(unitV, true);

                    // DPDU and DPDV are the normalized partials, written
                    // before ComputeOrthogonalComplement reorthogonalizes.
                    expectVectorClose(get(storage.dpdus as Float64Array, i),
                        unitU, 1e-12, 1e-12);
                    expectVectorClose(get(storage.dpdvs as Float64Array, i),
                        unitV, 1e-12, 1e-12);

                    const t = get(storage.tangents as Float64Array, i);
                    const b = get(storage.bitangents as Float64Array, i);
                    const n = get(storage.normals as Float64Array, i);
                    // {tangent, bitangent, normal} is right handed and
                    // orthonormal.
                    expectClose(length(t), 1, 1e-12, 1e-12);
                    expectClose(length(b), 1, 1e-12, 1e-12);
                    expectClose(length(n), 1, 1e-12, 1e-12);
                    expectClose(dot(t, b), 0, 1e-12, 1e-12);
                    expectClose(dot(t, n), 0, 1e-12, 1e-12);
                    expectClose(dot(b, n), 0, 1e-12, 1e-12);
                    expectVectorClose(cross(t, b), n, 1e-12, 1e-12);
                    // The tangent is the normalized u-partial (Gram-Schmidt
                    // leaves the first vector alone) and the normal spans the
                    // same line as dpdu x dpdv.
                    expectVectorClose(t, unitU, 1e-12, 1e-12);
                    expectVectorClose(n, unitCross(dpdu, dpdv, true),
                        1e-12, 1e-12);
                    // The bitangent stays in the tangent plane.
                    expectClose(dot(b, cross(dpdu, dpdv)), 0, 1e-12, 1e-12);
                }
            }
        }, 60);
    });

    it('recomputes every channel when the surface changes', () => {
        check(fc.tuple(scaleArb, scaleArb, gridArb),
            ([scale0, scale1, [nr, nc]]) => {
                const surface = new Paraboloid(scale0);
                const storage = makeStorage(nr, nc, { frame: true });
                const mesh = new RectanglePatchMesh(storage.description,
                    surface);
                surface.scale = scale1;
                mesh.update();

                const fresh = new Paraboloid(scale1);
                const freshStorage = makeStorage(nr, nc, { frame: true });
                const freshMesh = new RectanglePatchMesh(
                    freshStorage.description, fresh);
                expect(freshMesh.getDescription().numVertices).toBe(nr * nc);
                for (let i = 0; i < nr * nc; ++i) {
                    expectVectorClose(P(storage, i), P(freshStorage, i),
                        1e-12, 1e-12);
                    expectVectorClose(get(storage.normals as Float64Array, i),
                        get(freshStorage.normals as Float64Array, i),
                        1e-12, 1e-12);
                    expectVectorClose(get(storage.tangents as Float64Array, i),
                        get(freshStorage.tangents as Float64Array, i),
                        1e-12, 1e-12);
                }
            }, 40);
    });

    it('emits a consistently wound manifold grid of triangles', () => {
        check(fc.tuple(scaleArb, gridArb), ([scale, [nr, nc]]) => {
            const surface = new Paraboloid(scale);
            const storage = makeStorage(nr, nc);
            const mesh = new RectanglePatchMesh(storage.description, surface);
            const numTriangles = mesh.getDescription().numTriangles;
            expect(numTriangles).toBe(2 * (nr - 1) * (nc - 1));

            const directed = new Set<string>();
            for (let t = 0; t < numTriangles; ++t) {
                const i0 = storage.indices[3 * t];
                const i1 = storage.indices[3 * t + 1];
                const i2 = storage.indices[3 * t + 2];
                expect(i0 !== i1 && i1 !== i2 && i2 !== i0).toBe(true);
                for (const i of [i0, i1, i2]) {
                    expect(i).toBeLessThan(nr * nc);
                }
                // The triangle is nondegenerate.
                const n = cross(sub(P(storage, i1), P(storage, i0)),
                    sub(P(storage, i2), P(storage, i0)));
                expect(length(n)).toBeGreaterThan(0);
                for (const [a, b] of [[i0, i1], [i1, i2], [i2, i0]]) {
                    const key = `${a},${b}`;
                    expect(directed.has(key)).toBe(false);
                    directed.add(key);
                }
            }
            let boundary = 0;
            for (const key of directed) {
                const [a, b] = key.split(',');
                if (!directed.has(`${b},${a}`)) {
                    ++boundary;
                }
            }
            expect(boundary).toBe(2 * (nr - 1) + 2 * (nc - 1));
        }, 100);
    });
});
