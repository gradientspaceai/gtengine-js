import { describe, it, expect } from 'vitest';
import { RectangleMesh } from '../src/RectangleMesh';
import { MeshDescription, MeshTopology } from '../src/Mesh';
import { IndexAttribute } from '../src/IndexAttribute';
import { VertexAttribute } from '../src/VertexAttribute';
import { Rectangle } from '../src/Rectangle';
import { Vector, dot, length, normalize, sub } from '../src/Vector';
import { cross, unitCross } from '../src/Vector3';

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
    if (options.frame) {
        tangents = new Float64Array(3 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('tangent', tangents, 24));
        bitangents = new Float64Array(3 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('bitangent', bitangents, 24));
        description.wantDynamicTangentSpaceUpdate = true;
    }

    return { description, indices, positions, normals, tcoords, tangents,
        bitangents };
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

        const normal = unitCross(rectangle.axis[0], rectangle.axis[1]);
        for (let i = 0; i < mesh.getDescription().numVertices; ++i) {
            const n = N(storage, i);
            for (let k = 0; k < 3; ++k) {
                expect(n.values[k]).toBeCloseTo(normal.values[k], 12);
            }
            const t = storage.tangents as Float64Array;
            const b = storage.bitangents as Float64Array;
            // Upstream hardcodes tangent = (1,0,0), bitangent = (0,1,0).
            expect([t[3 * i], t[3 * i + 1], t[3 * i + 2]]).toEqual([1, 0, 0]);
            expect([b[3 * i], b[3 * i + 1], b[3 * i + 2]]).toEqual([0, 1, 0]);
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
