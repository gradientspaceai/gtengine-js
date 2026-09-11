import { describe, it, expect } from 'vitest';
import { RectangleMesh } from '../src/RectangleMesh.js';
import { MeshDescription, MeshTopology } from '../src/Mesh.js';
import { IndexAttribute } from '../src/IndexAttribute.js';
import { VertexAttribute } from '../src/VertexAttribute.js';
import { Rectangle } from '../src/Rectangle.js';
import {
    Vector, add, dot, length, mul, normalize, sub
} from '../src/Vector.js';
import { cross, unitCross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    wellScaledVector
} from './helpers/arbitraries.js';

function V(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function unit(x: number, y: number, z: number): Vector {
    const v = V(x, y, z);
    normalize(v);
    return v;
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

function getVec(buffer: Float64Array, i: number): Vector {
    return V(buffer[3 * i], buffer[3 * i + 1], buffer[3 * i + 2]);
}

function P(storage: Storage, i: number): Vector {
    return V(storage.positions[3 * i], storage.positions[3 * i + 1],
        storage.positions[3 * i + 2]);
}

function N(storage: Storage, i: number): Vector {
    const n = storage.normals as Float64Array;
    return V(n[3 * i], n[3 * i + 1], n[3 * i + 2]);
}

// The unit square in the plane z = 0, centered at the origin.
function unitRectangle(): Rectangle {
    return Rectangle.fromCenterAxisExtent(V(0, 0, 0),
        [V(1, 0, 0), V(0, 1, 0)], Vector.fromArray([1, 1]));
}

describe('RectangleMesh', () => {
    it('produces the vertex and triangle counts of a grid', () => {
        const storage = makeStorage(3, 4);
        const mesh = new RectangleMesh(storage.description, unitRectangle());
        const description = mesh.getDescription();
        expect(description.constructed).toBe(true);
        expect(description.numRows).toBe(3);
        expect(description.numCols).toBe(4);
        expect(description.numVertices).toBe(12);
        expect(description.numTriangles).toBe(2 * 2 * 3);
        expect(storage.indices.length).toBe(3 * 12);
        for (const index of storage.indices) {
            expect(index).toBeLessThan(12);
        }

        // The topology is clamped to at least 2 rows and 2 columns.
        const small = makeStorage(1, 1);
        const smallMesh = new RectangleMesh(small.description, unitRectangle());
        expect(smallMesh.getDescription().numRows).toBe(2);
        expect(smallMesh.getDescription().numCols).toBe(2);
        expect(smallMesh.getDescription().numVertices).toBe(4);
        expect(smallMesh.getDescription().numTriangles).toBe(2);
    });

    it('lays out positions on the rectangle in row-major order', () => {
        const storage = makeStorage(3, 5, { tcoords: true });
        const rectangle = Rectangle.fromCenterAxisExtent(V(1, 2, 3),
            [unit(1, 0, 0), unit(0, 1, 0)], Vector.fromArray([2, 4]));
        new RectangleMesh(storage.description, rectangle);

        const numRows = 3;
        const numCols = 5;
        for (let r = 0, i = 0; r < numRows; ++r) {
            for (let c = 0; c < numCols; ++c, ++i) {
                const u = c / (numCols - 1);
                const v = r / (numRows - 1);
                const tcoords = storage.tcoords as Float64Array;
                expect(tcoords[2 * i]).toBeCloseTo(u, 15);
                expect(tcoords[2 * i + 1]).toBeCloseTo(v, 15);

                const expected = V(1 + (2 * u - 1) * 2, 2 + (2 * v - 1) * 4, 3);
                const actual = P(storage, i);
                for (let k = 0; k < 3; ++k) {
                    expect(actual.values[k]).toBeCloseTo(expected.values[k], 12);
                }
            }
        }

        // The four corners are the rectangle's corners.
        expect(P(storage, 0).values).toEqual([-1, -2, 3]);
        expect(P(storage, numCols - 1).values).toEqual([3, -2, 3]);
        expect(P(storage, numRows * numCols - numCols).values)
            .toEqual([-1, 6, 3]);
        expect(P(storage, numRows * numCols - 1).values).toEqual([3, 6, 3]);
    });

    it('assigns the constant rectangle normal to every vertex', () => {
        const storage = makeStorage(4, 4);
        const rectangle = Rectangle.fromCenterAxisExtent(V(0, 0, 0),
            [unit(1, 1, 0), unit(-1, 1, 1)], Vector.fromArray([1, 2]));
        const mesh = new RectangleMesh(storage.description, rectangle);
        const expected = unitCross(rectangle.axis[0], rectangle.axis[1]);
        for (let i = 0; i < mesh.getDescription().numVertices; ++i) {
            const n = N(storage, i);
            expect(length(n)).toBeCloseTo(1, 12);
            for (let k = 0; k < 3; ++k) {
                expect(n.values[k]).toBeCloseTo(expected.values[k], 12);
            }
        }
    });

    it('copies the rectangle and exposes it', () => {
        const storage = makeStorage(2, 2);
        const rectangle = unitRectangle();
        const mesh = new RectangleMesh(storage.description, rectangle);
        rectangle.center.values[0] = 100;
        expect(mesh.getRectangle().center.values).toEqual([0, 0, 0]);

        // A non-3D rectangle is rejected.
        const flat = makeStorage(2, 2);
        expect(() => new RectangleMesh(flat.description, new Rectangle(2)))
            .toThrow();
    });

    it('emits counterclockwise triangles whose areas sum to the rectangle', () => {
        const storage = makeStorage(4, 6);
        const rectangle = Rectangle.fromCenterAxisExtent(V(0, 0, 0),
            [unit(1, 0, 0), unit(0, 1, 0)], Vector.fromArray([1.5, 2.5]));
        const mesh = new RectangleMesh(storage.description, rectangle);
        const normal = unitCross(rectangle.axis[0], rectangle.axis[1]);

        let area = 0;
        for (let t = 0; t < mesh.getDescription().numTriangles; ++t) {
            const a = P(storage, storage.indices[3 * t]);
            const b = P(storage, storage.indices[3 * t + 1]);
            const c = P(storage, storage.indices[3 * t + 2]);
            const n = cross(sub(b, a), sub(c, a));
            // The winding agrees with the rectangle normal.
            expect(dot(n, normal)).toBeGreaterThan(0);
            area += 0.5 * length(n);
        }
        expect(area).toBeCloseTo(4 * 1.5 * 2.5, 10);

        // With wantCCW = false, the winding is reversed.
        const cw = makeStorage(4, 6);
        cw.description.wantCCW = false;
        new RectangleMesh(cw.description, rectangle);
        for (let t = 0; t < cw.description.numTriangles; ++t) {
            const a = P(cw, cw.indices[3 * t]);
            const b = P(cw, cw.indices[3 * t + 1]);
            const c = P(cw, cw.indices[3 * t + 2]);
            expect(dot(cross(sub(b, a), sub(c, a)), normal)).toBeLessThan(0);
        }
    });

    it('initializes the frame when dynamic tangent-space updates are on', () => {
        const storage = makeStorage(3, 3, { frame: true });
        const rectangle = Rectangle.fromCenterAxisExtent(V(0, 0, 0),
            [unit(0, 1, 0), unit(0, 0, 1)], Vector.fromArray([1, 1]));
        const mesh = new RectangleMesh(storage.description, rectangle);
        expect(mesh.getDescription().allowUpdateFrame).toBe(true);

        // Port fix: upstream hardcodes tangent = (1,0,0) and
        // bitangent = (0,1,0); the port uses the rectangle's own axes.
        const normal = unitCross(rectangle.axis[0], rectangle.axis[1]);
        for (let i = 0; i < mesh.getDescription().numVertices; ++i) {
            const n = N(storage, i);
            for (let k = 0; k < 3; ++k) {
                expect(n.values[k]).toBeCloseTo(normal.values[k], 12);
            }
            const t = storage.tangents as Float64Array;
            const b = storage.bitangents as Float64Array;
            expect([t[3 * i], t[3 * i + 1], t[3 * i + 2]]).toEqual([0, 1, 0]);
            expect([b[3 * i], b[3 * i + 1], b[3 * i + 2]]).toEqual([0, 0, 1]);
        }
    });

    it('builds an orthonormal frame for a rectangle rotated out of z = 0', () => {
        // Regression test for the upstream hardcoded tangent/bitangent. The
        // rectangle's plane contains neither (1,0,0) nor (0,1,0), so the
        // upstream values would produce a frame that is neither orthonormal
        // nor tangent to the rectangle.
        const axis0 = unit(1, 1, 1);
        // A second axis orthogonal to axis0 and not axis aligned.
        const axis1 = unit(1, -2, 1);
        expect(dot(axis0, axis1)).toBeCloseTo(0, 15);
        const rectangle = Rectangle.fromCenterAxisExtent(V(-2, 3, 5),
            [axis0, axis1], Vector.fromArray([1.5, 0.75]));
        const storage = makeStorage(4, 3, { frame: true });
        const mesh = new RectangleMesh(storage.description, rectangle);
        expect(mesh.getDescription().allowUpdateFrame).toBe(true);

        const expectedNormal = unitCross(axis0, axis1);
        const tangents = storage.tangents as Float64Array;
        const bitangents = storage.bitangents as Float64Array;
        const dpdus = storage.dpdus as Float64Array;
        const dpdvs = storage.dpdvs as Float64Array;
        for (let i = 0; i < mesh.getDescription().numVertices; ++i) {
            const T = getVec(tangents, i);
            const B = getVec(bitangents, i);
            const Nv = N(storage, i);

            // The frame is orthonormal and right-handed.
            expect(length(T)).toBeCloseTo(1, 12);
            expect(length(B)).toBeCloseTo(1, 12);
            expect(length(Nv)).toBeCloseTo(1, 12);
            expect(dot(T, B)).toBeCloseTo(0, 12);
            expect(dot(T, Nv)).toBeCloseTo(0, 12);
            expect(dot(B, Nv)).toBeCloseTo(0, 12);
            const rightHanded = cross(T, B);
            for (let k = 0; k < 3; ++k) {
                expect(rightHanded.values[k]).toBeCloseTo(Nv.values[k], 12);
            }

            // The frame is the rectangle's own frame.
            for (let k = 0; k < 3; ++k) {
                expect(T.values[k]).toBeCloseTo(axis0.values[k], 12);
                expect(B.values[k]).toBeCloseTo(axis1.values[k], 12);
                expect(Nv.values[k]).toBeCloseTo(expectedNormal.values[k], 12);
                // dpdu/dpdv carry the same tangent and bitangent.
                expect(getVec(dpdus, i).values[k]).toBeCloseTo(axis0.values[k], 12);
                expect(getVec(dpdvs, i).values[k]).toBeCloseTo(axis1.values[k], 12);
            }

            // The vertex lies in the rectangle's plane, so the frame is
            // tangent to the surface it describes.
            const delta = sub(P(storage, i), rectangle.center);
            expect(dot(delta, Nv)).toBeCloseTo(0, 12);
        }
    });

    it('works without a normal channel and with client texture coordinates', () => {
        // No normals: only the positions and texture coordinates are written.
        const storage = makeStorage(3, 3, { normals: false, tcoords: true });
        const mesh = new RectangleMesh(storage.description, unitRectangle());
        expect(mesh.getDescription().allowUpdateFrame).toBe(false);
        expect(P(storage, 0).values).toEqual([-1, -1, 0]);
        expect(P(storage, 8).values).toEqual([1, 1, 0]);
        const tcoords = storage.tcoords as Float64Array;
        expect(tcoords[0]).toBe(0);
        expect(tcoords[1]).toBe(0);
        expect(tcoords[16]).toBe(1);
        expect(tcoords[17]).toBe(1);
    });
});

// ---------------------------------------------------------------------------
// V43 verification: the tessellation against the closed-form parametrization
// of the rectangle, and the frame against the rectangle axes (the #268 fix).
// ---------------------------------------------------------------------------

describe('RectangleMesh verification', () => {
    // A rectangle with a rotated orthonormal axis pair and positive extents.
    const rectangleArb = fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
        fc.array(finite(0.25, 4), { minLength: 2, maxLength: 2 }))
        .map(([center, frame, e]) => Rectangle.fromCenterAxisExtent(center,
            [frame[0], frame[1]], Vector.fromArray(e)));

    const gridArb = fc.tuple(fc.integer({ min: 2, max: 6 }),
        fc.integer({ min: 2, max: 6 }));

    it('places every vertex at the closed-form grid position', () => {
        check(fc.tuple(rectangleArb, gridArb), ([rectangle, [nr, nc]]) => {
            const storage = makeStorage(nr, nc, { tcoords: true });
            const mesh = new RectangleMesh(storage.description, rectangle);
            const description = mesh.getDescription();
            expect(description.numVertices).toBe(nr * nc);
            expect(description.numTriangles).toBe(2 * (nr - 1) * (nc - 1));

            const normal = unitCross(rectangle.axis[0], rectangle.axis[1]);
            const tcoords = storage.tcoords as Float64Array;
            for (let r = 0, i = 0; r < nr; ++r) {
                for (let c = 0; c < nc; ++c, ++i) {
                    // The texture coordinates are the normalized grid
                    // parameters, with the columns along u.
                    expectClose(tcoords[2 * i], c / (nc - 1), 0, 0);
                    expectClose(tcoords[2 * i + 1], r / (nr - 1), 0, 0);

                    const w0 = (2 * (c / (nc - 1)) - 1)
                        * rectangle.extent.get(0);
                    const w1 = (2 * (r / (nr - 1)) - 1)
                        * rectangle.extent.get(1);
                    const expected = add(add(rectangle.center,
                        mul(rectangle.axis[0], w0)),
                        mul(rectangle.axis[1], w1));
                    expectVectorClose(P(storage, i), expected, 1e-12, 1e-12);

                    // The vertex lies in the plane of the rectangle.
                    expectClose(dot(normal, sub(P(storage, i),
                        rectangle.center)), 0, 1e-12, 1e-12);
                    // Every vertex carries the rectangle normal.
                    expectVectorClose(N(storage, i), normal, 1e-12, 1e-12);
                }
            }
        }, 100);
    });

    it('emits a consistently wound manifold grid of triangles', () => {
        check(fc.tuple(rectangleArb, gridArb), ([rectangle, [nr, nc]]) => {
            const storage = makeStorage(nr, nc);
            const mesh = new RectangleMesh(storage.description, rectangle);
            const numTriangles = mesh.getDescription().numTriangles;
            const normal = unitCross(rectangle.axis[0], rectangle.axis[1]);

            let area = 0;
            const edgeCount = new Map<string, number>();
            for (let t = 0; t < numTriangles; ++t) {
                const i0 = storage.indices[3 * t];
                const i1 = storage.indices[3 * t + 1];
                const i2 = storage.indices[3 * t + 2];
                expect(i0 !== i1 && i1 !== i2 && i2 !== i0).toBe(true);
                const p0 = P(storage, i0);
                const p1 = P(storage, i1);
                const p2 = P(storage, i2);
                const n = cross(sub(p1, p0), sub(p2, p0));
                // The triangle is wound so its normal agrees with the
                // rectangle normal (positive area, no flipped triangle).
                expect(dot(n, normal)).toBeGreaterThan(0);
                area += 0.5 * length(n);

                // Every directed edge is used at most once: a manifold grid
                // with consistent orientation.
                for (const [a, b] of [[i0, i1], [i1, i2], [i2, i0]]) {
                    const key = `${a},${b}`;
                    expect(edgeCount.has(key)).toBe(false);
                    edgeCount.set(key, 1);
                }
            }
            // The triangles tile the rectangle exactly.
            expectClose(area,
                4 * rectangle.extent.get(0) * rectangle.extent.get(1),
                1e-12, 1e-12);

            // Each interior edge is shared by exactly two triangles (it
            // appears once in each direction); the boundary edges appear
            // once. A grid has 2*(nr-1)*(nc-1) triangles and
            // 2*(nr-1) + 2*(nc-1) boundary edges.
            let boundary = 0;
            for (const key of edgeCount.keys()) {
                const [a, b] = key.split(',');
                if (!edgeCount.has(`${b},${a}`)) {
                    ++boundary;
                }
            }
            expect(boundary).toBe(2 * (nr - 1) + 2 * (nc - 1));
        }, 100);
    });

    it('builds the tangent frame from the rectangle axes', () => {
        // Upstream hardcodes tangent = (1,0,0) and bitangent = (0,1,0) for
        // every rectangle, which is not in the tangent plane unless the
        // rectangle happens to lie in z = 0 with the standard axes. The port
        // uses the rectangle's own axes (upstream bug, see #268).
        check(fc.tuple(rectangleArb, gridArb), ([rectangle, [nr, nc]]) => {
            const storage = makeStorage(nr, nc, { frame: true });
            const mesh = new RectangleMesh(storage.description, rectangle);
            expect(mesh.getDescription().allowUpdateFrame).toBe(true);

            const normal = unitCross(rectangle.axis[0], rectangle.axis[1]);
            for (let i = 0; i < nr * nc; ++i) {
                const t = getVec(storage.tangents as Float64Array, i);
                const b = getVec(storage.bitangents as Float64Array, i);
                const n = N(storage, i);
                expectVectorClose(t, rectangle.axis[0], 1e-12, 1e-12);
                expectVectorClose(b, rectangle.axis[1], 1e-12, 1e-12);
                expectVectorClose(n, normal, 1e-12, 1e-12);
                // The frame is orthonormal and right handed.
                expectClose(dot(t, b), 0, 1e-12, 1e-12);
                expectClose(dot(t, n), 0, 1e-12, 1e-12);
                expectClose(dot(b, n), 0, 1e-12, 1e-12);
                expectClose(length(t), 1, 1e-12, 1e-12);
                expectClose(length(b), 1, 1e-12, 1e-12);
                expectVectorClose(cross(t, b), n, 1e-12, 1e-12);
                // dpdu and dpdv repeat the tangent and bitangent.
                expectVectorClose(getVec(storage.dpdus as Float64Array, i), t,
                    0, 0);
                expectVectorClose(getVec(storage.dpdvs as Float64Array, i), b,
                    0, 0);
            }
        }, 100);
    });

    it('copies the rectangle passed to the constructor', () => {
        check(rectangleArb, rectangle => {
            const storage = makeStorage(3, 3);
            const mesh = new RectangleMesh(storage.description, rectangle);
            const before = mesh.getRectangle().center.get(0);
            rectangle.center.set(0, rectangle.center.get(0) + 100);
            expect(mesh.getRectangle().center.get(0)).toBe(before);
        });
    });
});
