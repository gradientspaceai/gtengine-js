import { describe, it, expect } from 'vitest';
import {
    Mesh, MeshChannel, MeshDescription, MeshTopology
} from '../src/Mesh';
import { IndexAttribute } from '../src/IndexAttribute';
import { VertexAttribute } from '../src/VertexAttribute';
import { Vector, dot, length, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

// A minimal concrete Mesh used to drive the protected base-class algorithms.
// The positions are supplied by a parametric function of the grid indices.
class TestMesh extends Mesh {
    computeIndicesPublic(): void {
        this.computeIndices();
    }

    updateNormalsPublic(): void {
        this.updateNormals();
    }

    updateFramePublic(): void {
        this.updateFrame();
    }

    positionPublic(i: number): Vector {
        return this.position(i);
    }

    setPositionPublic(i: number, v: Vector): void {
        this.setPosition(i, v);
    }

    normalPublic(i: number): Vector {
        return this.normal(i);
    }

    tangentPublic(i: number): Vector {
        return this.tangent(i);
    }

    bitangentPublic(i: number): Vector {
        return this.bitangent(i);
    }

    dpduPublic(i: number): Vector {
        return this.dpdu(i);
    }

    dpdvPublic(i: number): Vector {
        return this.dpdv(i);
    }

    tcoordPublic(i: number): Vector {
        return this.tcoord(i);
    }

    setTCoordPublic(i: number, v: Vector): void {
        this.setTCoord(i, v);
    }
}

// A Mesh whose updatePositions() writes a caller-supplied surface function
// of the vertex index.
class SurfaceMesh extends TestMesh {
    surface: (i: number) => Vector = () => Vector.fromArray([0, 0, 0]);

    protected override updatePositions(): void {
        for (let i = 0; i < this.getDescription().numVertices; ++i) {
            this.setPosition(i, this.surface(i));
        }
    }
}

// Storage for the attributes of one test mesh.
interface Storage {
    description: MeshDescription;
    indices: Uint32Array;
    channels: Map<string, Float64Array>;
}

const SENTINEL = 0xFFFFFFFF;

function makeStorage(topology: MeshTopology, numRows: number, numCols: number,
    semantics: readonly string[], wantCCW: boolean = true,
    wantDynamic: boolean = false): Storage {
    const description = new MeshDescription(topology, numRows, numCols);
    // One extra triangle slot filled with a sentinel so that a test can
    // verify computeIndices() writes exactly numTriangles triangles.
    const indices = new Uint32Array(3 * (description.numTriangles + 1));
    indices.fill(SENTINEL);
    description.indexAttribute = new IndexAttribute(indices, 4);
    description.wantCCW = wantCCW;
    description.wantDynamicTangentSpaceUpdate = wantDynamic;

    const channels = new Map<string, Float64Array>();
    for (const semantic of semantics) {
        const numComponents = (semantic === 'tcoord' ? 2 : 3);
        const data = new Float64Array(description.numVertices * numComponents);
        channels.set(semantic, data);
        description.vertexAttributes.push(
            new VertexAttribute(semantic, data, 8 * numComponents));
    }

    return { description, indices, channels };
}

// The triangles actually written by computeIndices().
function getTriangles(storage: Storage): Array<[number, number, number]> {
    const triangles: Array<[number, number, number]> = [];
    for (let t = 0; t < storage.description.numTriangles; ++t) {
        const { v0, v1, v2 } = storage.description.indexAttribute.getTriangle(t);
        triangles.push([v0, v1, v2]);
    }
    return triangles;
}

function expectWellFormed(storage: Storage): void {
    const numVertices = storage.description.numVertices;
    const triangles = getTriangles(storage);
    for (const [v0, v1, v2] of triangles) {
        for (const v of [v0, v1, v2]) {
            expect(v).toBeGreaterThanOrEqual(0);
            expect(v).toBeLessThan(numVertices);
        }
        expect(v0).not.toBe(v1);
        expect(v1).not.toBe(v2);
        expect(v0).not.toBe(v2);
    }

    // The slot past the last triangle must still hold the sentinel, which
    // verifies that exactly numTriangles triangles were written.
    const t = storage.description.numTriangles;
    expect(storage.indices[3 * t]).toBe(SENTINEL);
    expect(storage.indices[3 * t + 1]).toBe(SENTINEL);
    expect(storage.indices[3 * t + 2]).toBe(SENTINEL);
}

// Every vertex must be used by at least one triangle.
function usedVertices(storage: Storage): Set<number> {
    const used = new Set<number>();
    for (const [v0, v1, v2] of getTriangles(storage)) {
        used.add(v0);
        used.add(v1);
        used.add(v2);
    }
    return used;
}

// A planar grid position for topologies whose rows have 'rowStride' vertices.
function planePosition(r: number, c: number, numRows: number,
    numCols: number): Vector {
    const u = c / (numCols - 1);
    const v = r / (numRows - 1);
    return Vector.fromArray([2 * u, 3 * v, 0]);
}

describe('MeshTopology and MeshDescription', () => {
    it('computes the upstream vertex and triangle counts', () => {
        const rows = 5, cols = 7;

        const rectangle = new MeshDescription(MeshTopology.RECTANGLE, rows, cols);
        expect(rectangle.numRows).toBe(rows);
        expect(rectangle.numCols).toBe(cols);
        expect(rectangle.rMax).toBe(rows - 1);
        expect(rectangle.cMax).toBe(cols - 1);
        expect(rectangle.rIncrement).toBe(cols);
        expect(rectangle.numVertices).toBe(rows * cols);
        expect(rectangle.numTriangles).toBe(2 * (rows - 1) * (cols - 1));

        const cylinder = new MeshDescription(MeshTopology.CYLINDER, rows, cols);
        expect(cylinder.rMax).toBe(rows - 1);
        expect(cylinder.cMax).toBe(cols);
        expect(cylinder.rIncrement).toBe(cols + 1);
        expect(cylinder.numVertices).toBe(rows * (cols + 1));
        expect(cylinder.numTriangles).toBe(2 * (rows - 1) * cols);

        const torus = new MeshDescription(MeshTopology.TORUS, rows, cols);
        expect(torus.rMax).toBe(rows);
        expect(torus.cMax).toBe(cols);
        expect(torus.rIncrement).toBe(cols + 1);
        expect(torus.numVertices).toBe((rows + 1) * (cols + 1));
        expect(torus.numTriangles).toBe(2 * rows * cols);

        const disk = new MeshDescription(MeshTopology.DISK, rows, cols);
        expect(disk.rMax).toBe(rows - 1);
        expect(disk.cMax).toBe(cols);
        expect(disk.rIncrement).toBe(cols + 1);
        expect(disk.numVertices).toBe(rows * (cols + 1) + 1);
        expect(disk.numTriangles).toBe(2 * (rows - 1) * cols + cols);

        const sphere = new MeshDescription(MeshTopology.SPHERE, rows, cols);
        expect(sphere.rMax).toBe(rows - 1);
        expect(sphere.cMax).toBe(cols);
        expect(sphere.rIncrement).toBe(cols + 1);
        expect(sphere.numVertices).toBe(rows * (cols + 1) + 2);
        expect(sphere.numTriangles).toBe(2 * (rows - 1) * cols + 2 * cols);
    });

    it('clamps invalid row and column counts', () => {
        const rectangle = new MeshDescription(MeshTopology.RECTANGLE, 0, 1);
        expect(rectangle.numRows).toBe(2);
        expect(rectangle.numCols).toBe(2);

        for (const topology of [MeshTopology.CYLINDER, MeshTopology.TORUS]) {
            const description = new MeshDescription(topology, 1, 2);
            expect(description.numRows).toBe(2);
            expect(description.numCols).toBe(3);
        }

        for (const topology of [MeshTopology.DISK, MeshTopology.SPHERE]) {
            const description = new MeshDescription(topology, 0, 0);
            expect(description.numRows).toBe(1);
            expect(description.numCols).toBe(3);
        }
    });

    it('supports the ARBITRARY topology constructors', () => {
        const viaCtor = new MeshDescription(MeshTopology.ARBITRARY, 10, 4);
        expect(viaCtor.numVertices).toBe(10);
        expect(viaCtor.numTriangles).toBe(4);
        expect(viaCtor.numRows).toBe(0);
        expect(viaCtor.numCols).toBe(0);
        expect(viaCtor.rMax).toBe(0);
        expect(viaCtor.cMax).toBe(0);
        expect(viaCtor.rIncrement).toBe(0);

        const viaFactory = MeshDescription.arbitrary(10, 4);
        expect(viaFactory.topology).toBe(MeshTopology.ARBITRARY);
        expect(viaFactory.numVertices).toBe(10);
        expect(viaFactory.numTriangles).toBe(4);

        expect(() => MeshDescription.arbitrary(2, 1)).toThrow('Invalid input.');
        expect(() => MeshDescription.arbitrary(3, 0)).toThrow('Invalid input.');
    });

    it('has the upstream defaults and clones all members', () => {
        const description = new MeshDescription(MeshTopology.RECTANGLE, 3, 4);
        expect(description.wantCCW).toBe(true);
        expect(description.wantDynamicTangentSpaceUpdate).toBe(false);
        expect(description.hasTangentSpaceVectors).toBe(false);
        expect(description.allowUpdateFrame).toBe(false);
        expect(description.constructed).toBe(false);
        expect(description.vertexAttributes.length).toBe(0);
        expect(description.indexAttribute.source).toBeNull();

        description.wantCCW = false;
        description.vertexAttributes.push(new VertexAttribute('position',
            new Float64Array(3 * description.numVertices), 24));
        const clone = description.clone();
        expect(clone.topology).toBe(description.topology);
        expect(clone.numVertices).toBe(description.numVertices);
        expect(clone.numTriangles).toBe(description.numTriangles);
        expect(clone.wantCCW).toBe(false);
        expect(clone.rIncrement).toBe(description.rIncrement);
        expect(clone.vertexAttributes).not.toBe(description.vertexAttributes);
        expect(clone.vertexAttributes[0]).toBe(description.vertexAttributes[0]);
        expect(clone.indexAttribute).toBe(description.indexAttribute);
    });
});

describe('Mesh construction', () => {
    it('rejects descriptions without indices or positions', () => {
        const missingIndices = new MeshDescription(MeshTopology.RECTANGLE, 3, 3);
        missingIndices.vertexAttributes.push(new VertexAttribute('position',
            new Float64Array(27), 24));
        expect(() => new TestMesh(missingIndices, [MeshTopology.RECTANGLE]))
            .toThrow('The mesh needs triangles/indices in Mesh constructor.');

        const missingPositions = makeStorage(MeshTopology.RECTANGLE, 3, 3,
            ['normal']).description;
        expect(() => new TestMesh(missingPositions, [MeshTopology.RECTANGLE]))
            .toThrow('The mesh needs positions in Mesh constructor.');

        // An attribute with a nonpositive stride is ignored, as upstream.
        const zeroStride = makeStorage(MeshTopology.RECTANGLE, 3, 3, []).description;
        zeroStride.vertexAttributes.push(new VertexAttribute('position',
            new Float64Array(27), 0));
        expect(() => new TestMesh(zeroStride, [MeshTopology.RECTANGLE]))
            .toThrow('The mesh needs positions in Mesh constructor.');
    });

    it('sets constructed only for a valid topology and does not modify the input',
        () => {
            const storage = makeStorage(MeshTopology.TORUS, 3, 4, ['position']);
            const mesh = new TestMesh(storage.description,
                [MeshTopology.RECTANGLE, MeshTopology.CYLINDER]);
            expect(mesh.getDescription().constructed).toBe(false);
            // The description is copied by value upstream, so the caller's
            // object is unchanged.
            expect(storage.description.constructed).toBe(false);
            expect(() => mesh.update()).toThrow(
                'The Mesh object failed the construction.');

            const good = new TestMesh(storage.description,
                [MeshTopology.RECTANGLE, MeshTopology.TORUS]);
            expect(good.getDescription().constructed).toBe(true);
            expect(storage.description.constructed).toBe(false);
            expect(() => good.update()).not.toThrow();
        });

    it('determines hasTangentSpaceVectors and allowUpdateFrame as upstream', () => {
        // No tangent-space channels: no frame updates.
        let storage = makeStorage(MeshTopology.RECTANGLE, 3, 4,
            ['position', 'normal', 'tcoord'], true, true);
        let mesh = new TestMesh(storage.description, [MeshTopology.RECTANGLE]);
        expect(mesh.getDescription().hasTangentSpaceVectors).toBe(false);
        expect(mesh.getDescription().allowUpdateFrame).toBe(false);

        // Tangent-space channels but no texture coordinates and no normals.
        storage = makeStorage(MeshTopology.RECTANGLE, 3, 4,
            ['position', 'dpdu', 'dpdv'], true, true);
        mesh = new TestMesh(storage.description, [MeshTopology.RECTANGLE]);
        expect(mesh.getDescription().hasTangentSpaceVectors).toBe(true);
        expect(mesh.getDescription().allowUpdateFrame).toBe(false);

        // Tangent-space channels with normals: local coordinates are used.
        storage = makeStorage(MeshTopology.RECTANGLE, 3, 4,
            ['position', 'normal', 'tangent', 'bitangent'], true, true);
        mesh = new TestMesh(storage.description, [MeshTopology.RECTANGLE]);
        expect(mesh.getDescription().hasTangentSpaceVectors).toBe(true);
        expect(mesh.getDescription().allowUpdateFrame).toBe(true);

        // Tangent-space channels with texture coordinates.
        storage = makeStorage(MeshTopology.RECTANGLE, 3, 4,
            ['position', 'tcoord', 'dpdu', 'dpdv'], true, true);
        mesh = new TestMesh(storage.description, [MeshTopology.RECTANGLE]);
        expect(mesh.getDescription().allowUpdateFrame).toBe(true);

        // The client must request dynamic updates.
        storage = makeStorage(MeshTopology.RECTANGLE, 3, 4,
            ['position', 'tcoord', 'dpdu', 'dpdv'], true, false);
        mesh = new TestMesh(storage.description, [MeshTopology.RECTANGLE]);
        expect(mesh.getDescription().hasTangentSpaceVectors).toBe(true);
        expect(mesh.getDescription().allowUpdateFrame).toBe(false);
    });
});

describe('Mesh.computeIndices', () => {
    it('builds a well-formed RECTANGLE index buffer with consistent winding',
        () => {
            const rows = 4, cols = 5;
            const storage = makeStorage(MeshTopology.RECTANGLE, rows, cols,
                ['position']);
            const mesh = new SurfaceMesh(storage.description,
                [MeshTopology.RECTANGLE]);
            mesh.surface = (i: number) => planePosition(Math.floor(i / cols),
                i % cols, rows, cols);
            mesh.computeIndicesPublic();
            mesh.update();

            expectWellFormed(storage);
            expect(storage.description.numTriangles).toBe(2 * (rows - 1) * (cols - 1));
            expect(usedVertices(storage).size).toBe(storage.description.numVertices);

            // Each quad of the grid produces the upstream pair of triangles.
            const triangles = getTriangles(storage);
            let t = 0;
            for (let r = 0; r < rows - 1; ++r) {
                for (let c = 0; c < cols - 1; ++c) {
                    const v0 = r * cols + c, v1 = v0 + 1;
                    const v2 = (r + 1) * cols + c, v3 = v2 + 1;
                    expect(triangles[t++]).toEqual([v0, v1, v2]);
                    expect(triangles[t++]).toEqual([v1, v3, v2]);
                }
            }

            // The plane lies in z = 0 and the winding is counterclockwise, so
            // every triangle normal is +z.
            for (const [v0, v1, v2] of triangles) {
                const n = cross(sub(mesh.positionPublic(v1), mesh.positionPublic(v0)),
                    sub(mesh.positionPublic(v2), mesh.positionPublic(v0)));
                expect(n.values[2]).toBeGreaterThan(0);
            }
        });

    it('reverses the winding when wantCCW is false', () => {
        const rows = 3, cols = 4;
        const ccw = makeStorage(MeshTopology.RECTANGLE, rows, cols, ['position'],
            true);
        const cw = makeStorage(MeshTopology.RECTANGLE, rows, cols, ['position'],
            false);
        new TestMesh(ccw.description, [MeshTopology.RECTANGLE]).computeIndicesPublic();
        new TestMesh(cw.description, [MeshTopology.RECTANGLE]).computeIndicesPublic();

        const ccwTriangles = getTriangles(ccw);
        const cwTriangles = getTriangles(cw);
        for (let t = 0; t < ccwTriangles.length; t += 2) {
            const [v0, v1, v2] = ccwTriangles[t];
            const [w1, w3, w2] = ccwTriangles[t + 1];
            expect(cwTriangles[t]).toEqual([v0, v2, v1]);
            expect(cwTriangles[t + 1]).toEqual([w1, w2, w3]);
        }
    });

    it('closes the CYLINDER seam with the duplicated last column', () => {
        const rows = 3, cols = 5;
        const storage = makeStorage(MeshTopology.CYLINDER, rows, cols, ['position']);
        const mesh = new TestMesh(storage.description, [MeshTopology.CYLINDER]);
        mesh.computeIndicesPublic();
        expectWellFormed(storage);
        expect(usedVertices(storage).size).toBe(storage.description.numVertices);

        const triangles = getTriangles(storage);
        const rIncrement = cols + 1;
        let t = 0;
        for (let r = 0; r < rows - 1; ++r) {
            for (let c = 0; c < cols; ++c) {
                const v0 = r * rIncrement + c, v1 = v0 + 1;
                const v2 = (r + 1) * rIncrement + c, v3 = v2 + 1;
                expect(triangles[t++]).toEqual([v0, v1, v2]);
                expect(triangles[t++]).toEqual([v1, v3, v2]);
            }
        }

        // The last quad of each row references the seam duplicate, which is
        // the last vertex of that row.
        for (let r = 0; r < rows - 1; ++r) {
            const seam = r * rIncrement + cols;
            const wrapTriangle = triangles[2 * (r * cols + cols - 1)];
            expect(wrapTriangle[1]).toBe(seam);
        }
    });

    it('closes the TORUS seams in both directions', () => {
        const rows = 4, cols = 5;
        const storage = makeStorage(MeshTopology.TORUS, rows, cols, ['position']);
        const mesh = new TestMesh(storage.description, [MeshTopology.TORUS]);
        mesh.computeIndicesPublic();
        expectWellFormed(storage);

        const used = usedVertices(storage);
        expect(used.size).toBe(storage.description.numVertices);
        // The last row (r = numRows) duplicates row 0 and the last column
        // duplicates column 0; both are referenced.
        const rIncrement = cols + 1;
        expect(used.has(rows * rIncrement)).toBe(true);
        expect(used.has(rows * rIncrement + cols)).toBe(true);
        expect(Math.max(...used)).toBe(storage.description.numVertices - 1);
    });

    it('builds the DISK fan around the center vertex', () => {
        const rows = 3, cols = 6;
        const storage = makeStorage(MeshTopology.DISK, rows, cols, ['position']);
        const mesh = new TestMesh(storage.description, [MeshTopology.DISK]);
        mesh.computeIndicesPublic();
        expectWellFormed(storage);
        expect(usedVertices(storage).size).toBe(storage.description.numVertices);

        const triangles = getTriangles(storage);
        const center = storage.description.numVertices - 1;
        const fan = triangles.slice(2 * (rows - 1) * cols);
        expect(fan.length).toBe(cols);
        for (let c = 0; c < cols; ++c) {
            expect(fan[c]).toEqual([c, center, c + 1]);
        }
    });

    it('handles a single-row DISK (no grid, fan only)', () => {
        const cols = 4;
        const storage = makeStorage(MeshTopology.DISK, 1, cols, ['position']);
        expect(storage.description.numTriangles).toBe(cols);
        expect(storage.description.numVertices).toBe(cols + 2);
        const mesh = new TestMesh(storage.description, [MeshTopology.DISK]);
        mesh.computeIndicesPublic();
        expectWellFormed(storage);
        const center = storage.description.numVertices - 1;
        const triangles = getTriangles(storage);
        for (let c = 0; c < cols; ++c) {
            expect(triangles[c]).toEqual([c, center, c + 1]);
        }
    });

    it('builds both SPHERE pole fans against the correct rows', () => {
        const rows = 4, cols = 5;
        const storage = makeStorage(MeshTopology.SPHERE, rows, cols, ['position']);
        const mesh = new TestMesh(storage.description, [MeshTopology.SPHERE]);
        mesh.computeIndicesPublic();
        expectWellFormed(storage);
        expect(usedVertices(storage).size).toBe(storage.description.numVertices);

        const triangles = getTriangles(storage);
        const numVertices = storage.description.numVertices;
        const rIncrement = cols + 1;
        const first = 2 * (rows - 1) * cols;

        // The first fan joins row 0 to the pole numVertices - 2.
        for (let c = 0; c < cols; ++c) {
            expect(triangles[first + c]).toEqual([c, numVertices - 2, c + 1]);
        }

        // The second fan joins the last row to the pole numVertices - 1.
        // Upstream computes the first vertex of the last row as
        // (numRows - 1) * numCols; the port uses rMax * rIncrement.
        const lastRow = (rows - 1) * rIncrement;
        for (let c = 0; c < cols; ++c) {
            expect(triangles[first + cols + c]).toEqual(
                [lastRow + c, numVertices - 1, lastRow + c + 1]);
        }
    });

    it('handles a single-row SPHERE where both fans share row 0', () => {
        const cols = 3;
        const storage = makeStorage(MeshTopology.SPHERE, 1, cols, ['position']);
        expect(storage.description.numTriangles).toBe(2 * cols);
        const mesh = new TestMesh(storage.description, [MeshTopology.SPHERE]);
        mesh.computeIndicesPublic();
        expectWellFormed(storage);
        const numVertices = storage.description.numVertices;
        const triangles = getTriangles(storage);
        for (let c = 0; c < cols; ++c) {
            expect(triangles[c]).toEqual([c, numVertices - 2, c + 1]);
            expect(triangles[cols + c]).toEqual([c, numVertices - 1, c + 1]);
        }
    });
});

describe('Mesh.updateNormals', () => {
    it('computes exact unit normals for a planar rectangle mesh', () => {
        const rows = 4, cols = 5;
        const storage = makeStorage(MeshTopology.RECTANGLE, rows, cols,
            ['position', 'normal']);
        const mesh = new SurfaceMesh(storage.description, [MeshTopology.RECTANGLE]);
        mesh.surface = (i: number) => planePosition(Math.floor(i / cols), i % cols,
            rows, cols);
        mesh.computeIndicesPublic();
        mesh.update();

        for (let i = 0; i < storage.description.numVertices; ++i) {
            const n = mesh.normalPublic(i);
            expect(length(n)).toBeCloseTo(1, 12);
            expect(n.values[0]).toBeCloseTo(0, 12);
            expect(n.values[1]).toBeCloseTo(0, 12);
            expect(n.values[2]).toBeCloseTo(1, 12);
        }
    });

    it('flips the normals when the winding is clockwise', () => {
        const rows = 3, cols = 3;
        const storage = makeStorage(MeshTopology.RECTANGLE, rows, cols,
            ['position', 'normal'], false);
        const mesh = new SurfaceMesh(storage.description, [MeshTopology.RECTANGLE]);
        mesh.surface = (i: number) => planePosition(Math.floor(i / cols), i % cols,
            rows, cols);
        mesh.computeIndicesPublic();
        mesh.update();

        for (let i = 0; i < storage.description.numVertices; ++i) {
            expect(mesh.normalPublic(i).values[2]).toBeCloseTo(-1, 12);
        }
    });

    it('matches the analytic radial normals of a cylinder', () => {
        const rows = 4, cols = 8, radius = 2, height = 3;
        const storage = makeStorage(MeshTopology.CYLINDER, rows, cols,
            ['position', 'normal']);
        const mesh = new SurfaceMesh(storage.description, [MeshTopology.CYLINDER]);
        const rIncrement = cols + 1;
        const angleOf = (c: number) => 2 * Math.PI * c / cols;
        mesh.surface = (i: number) => {
            const r = Math.floor(i / rIncrement);
            const c = i - r * rIncrement;
            const angle = angleOf(c);
            return Vector.fromArray([radius * Math.cos(angle),
                radius * Math.sin(angle), height * r / (rows - 1)]);
        };
        mesh.computeIndicesPublic();
        mesh.update();

        // Every facet normal of a prism is horizontal, so every vertex
        // normal is horizontal and unit length. The area-weighted average is
        // not exactly radial because the quads are split into triangles
        // asymmetrically, but it stays within half a facet angle of the
        // analytic radial direction.
        const halfFacet = Math.cos(Math.PI / cols);
        for (let r = 0; r < rows; ++r) {
            for (let c = 0; c <= cols; ++c) {
                const n = mesh.normalPublic(r * rIncrement + c);
                expect(length(n)).toBeCloseTo(1, 12);
                expect(n.values[2]).toBeCloseTo(0, 12);
                const angle = angleOf(c);
                const radial = Vector.fromArray(
                    [Math.cos(angle), Math.sin(angle), 0]);
                expect(dot(n, radial)).toBeGreaterThanOrEqual(halfFacet - 1e-12);
            }
        }
    });

    it('matches an independent accumulation for random positions', () => {
        let seed = 987654321;
        const random = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 20; ++trial) {
            const rows = 3 + (trial % 3), cols = 3 + ((trial + 1) % 4);
            const storage = makeStorage(MeshTopology.RECTANGLE, rows, cols,
                ['position', 'normal']);
            const numVertices = storage.description.numVertices;
            const points: Vector[] = [];
            for (let i = 0; i < numVertices; ++i) {
                points.push(Vector.fromArray([
                    (i % cols) + 0.2 * random(),
                    Math.floor(i / cols) + 0.2 * random(),
                    0.5 * random()
                ]));
            }

            const mesh = new SurfaceMesh(storage.description,
                [MeshTopology.RECTANGLE]);
            mesh.surface = (i: number) => points[i];
            mesh.computeIndicesPublic();
            mesh.update();

            // Independent computation of the area-weighted vertex normals.
            const expected: Vector[] = [];
            for (let i = 0; i < numVertices; ++i) {
                expected.push(Vector.fromArray([0, 0, 0]));
            }
            for (const [v0, v1, v2] of getTriangles(storage)) {
                const n = cross(sub(points[v1], points[v0]),
                    sub(points[v2], points[v0]));
                for (const v of [v0, v1, v2]) {
                    for (let k = 0; k < 3; ++k) {
                        expected[v].values[k] += n.values[k];
                    }
                }
            }
            for (let i = 0; i < numVertices; ++i) {
                const len = length(expected[i]);
                const actual = mesh.normalPublic(i);
                expect(length(actual)).toBeCloseTo(1, 10);
                for (let k = 0; k < 3; ++k) {
                    expect(actual.values[k]).toBeCloseTo(
                        expected[i].values[k] / len, 10);
                }
            }
        }
    });
});

describe('Mesh.updateFrame', () => {
    it('recovers the exact surface derivatives of a plane from tcoords', () => {
        const rows = 4, cols = 5;
        const storage = makeStorage(MeshTopology.RECTANGLE, rows, cols,
            ['position', 'normal', 'tcoord', 'dpdu', 'dpdv', 'tangent',
                'bitangent'], true, true);
        const mesh = new SurfaceMesh(storage.description, [MeshTopology.RECTANGLE]);
        expect(mesh.getDescription().allowUpdateFrame).toBe(true);

        // P(u,v) = (2u, 3v, 0) with (u,v) the grid texture coordinates.
        for (let r = 0, i = 0; r < rows; ++r) {
            for (let c = 0; c < cols; ++c, ++i) {
                mesh.setTCoordPublic(i, Vector.fromArray(
                    [c / (cols - 1), r / (rows - 1)]));
            }
        }
        mesh.surface = (i: number) => planePosition(Math.floor(i / cols), i % cols,
            rows, cols);
        mesh.computeIndicesPublic();
        mesh.update();

        for (let i = 0; i < storage.description.numVertices; ++i) {
            const dpdu = mesh.dpduPublic(i);
            const dpdv = mesh.dpdvPublic(i);
            expect(dpdu.values[0]).toBeCloseTo(2, 10);
            expect(dpdu.values[1]).toBeCloseTo(0, 10);
            expect(dpdu.values[2]).toBeCloseTo(0, 10);
            expect(dpdv.values[0]).toBeCloseTo(0, 10);
            expect(dpdv.values[1]).toBeCloseTo(3, 10);
            expect(dpdv.values[2]).toBeCloseTo(0, 10);

            const tangent = mesh.tangentPublic(i);
            const bitangent = mesh.bitangentPublic(i);
            const normal = mesh.normalPublic(i);
            expect(tangent.values[0]).toBeCloseTo(1, 10);
            expect(bitangent.values[1]).toBeCloseTo(1, 10);
            expect(normal.values[2]).toBeCloseTo(1, 10);
            expect(length(tangent)).toBeCloseTo(1, 10);
            expect(length(bitangent)).toBeCloseTo(1, 10);
            expect(dot(tangent, bitangent)).toBeCloseTo(0, 10);
        }

        // The texture coordinates are unchanged, in [0,1] and monotone along
        // the grid.
        for (let r = 0, i = 0; r < rows; ++r) {
            for (let c = 0; c < cols; ++c, ++i) {
                const tc = mesh.tcoordPublic(i);
                expect(tc.values[0]).toBeGreaterThanOrEqual(0);
                expect(tc.values[0]).toBeLessThanOrEqual(1);
                expect(tc.values[1]).toBeGreaterThanOrEqual(0);
                expect(tc.values[1]).toBeLessThanOrEqual(1);
                if (c > 0) {
                    expect(tc.values[0]).toBeGreaterThan(
                        mesh.tcoordPublic(i - 1).values[0]);
                }
                if (r > 0) {
                    expect(tc.values[1]).toBeGreaterThan(
                        mesh.tcoordPublic(i - cols).values[1]);
                }
            }
        }
    });

    it('uses local coordinates when no texture coordinates are supplied', () => {
        const rows = 4, cols = 4;
        const storage = makeStorage(MeshTopology.RECTANGLE, rows, cols,
            ['position', 'normal', 'dpdu', 'dpdv', 'tangent', 'bitangent'],
            true, true);
        const mesh = new SurfaceMesh(storage.description, [MeshTopology.RECTANGLE]);
        expect(mesh.getDescription().allowUpdateFrame).toBe(true);
        mesh.surface = (i: number) => planePosition(Math.floor(i / cols), i % cols,
            rows, cols);
        mesh.computeIndicesPublic();
        mesh.update();

        for (let i = 0; i < storage.description.numVertices; ++i) {
            // The local frame of a planar patch is an isometry, so the
            // estimated derivatives are orthonormal and span the plane.
            const dpdu = mesh.dpduPublic(i);
            const dpdv = mesh.dpdvPublic(i);
            expect(length(dpdu)).toBeCloseTo(1, 8);
            expect(length(dpdv)).toBeCloseTo(1, 8);
            expect(dot(dpdu, dpdv)).toBeCloseTo(0, 8);
            expect(dpdu.values[2]).toBeCloseTo(0, 8);
            expect(dpdv.values[2]).toBeCloseTo(0, 8);

            const normal = mesh.normalPublic(i);
            expect(Math.abs(normal.values[2])).toBeCloseTo(1, 8);
            expect(dot(normal, mesh.tangentPublic(i))).toBeCloseTo(0, 8);
            expect(dot(normal, mesh.bitangentPublic(i))).toBeCloseTo(0, 8);
        }
    });

    it('produces an orthonormal frame for random cylinder-like surfaces', () => {
        let seed = 24680;
        const random = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 10; ++trial) {
            const rows = 4, cols = 6;
            const storage = makeStorage(MeshTopology.CYLINDER, rows, cols,
                ['position', 'normal', 'tcoord', 'tangent', 'bitangent'],
                true, true);
            const rIncrement = cols + 1;
            const mesh = new SurfaceMesh(storage.description,
                [MeshTopology.CYLINDER]);
            const radius = 1 + random();
            for (let r = 0, i = 0; r < rows; ++r) {
                for (let c = 0; c <= cols; ++c, ++i) {
                    mesh.setTCoordPublic(i, Vector.fromArray(
                        [c / cols, r / (rows - 1)]));
                }
            }
            mesh.surface = (i: number) => {
                const r = Math.floor(i / rIncrement);
                const c = i - r * rIncrement;
                const angle = 2 * Math.PI * c / cols;
                return Vector.fromArray([radius * Math.cos(angle),
                    radius * Math.sin(angle), r / (rows - 1)]);
            };
            mesh.computeIndicesPublic();
            mesh.update();

            for (let i = 0; i < storage.description.numVertices; ++i) {
                const t = mesh.tangentPublic(i);
                const b = mesh.bitangentPublic(i);
                const n = mesh.normalPublic(i);
                expect(length(t)).toBeCloseTo(1, 8);
                expect(length(b)).toBeCloseTo(1, 8);
                expect(length(n)).toBeCloseTo(1, 8);
                expect(dot(t, b)).toBeCloseTo(0, 8);
                expect(dot(t, n)).toBeCloseTo(0, 8);
                expect(dot(b, n)).toBeCloseTo(0, 8);
                // Right-handed: n = t x b.
                const c = cross(t, b);
                for (let k = 0; k < 3; ++k) {
                    expect(c.values[k]).toBeCloseTo(n.values[k], 8);
                }
                // Away from the two boundary rows, where the least-squares
                // fit is one-sided, the estimated normal of a cylinder of
                // revolution is horizontal and within half a facet angle of
                // the radial direction.
                const r = Math.floor(i / rIncrement);
                if (0 < r && r < rows - 1) {
                    expect(n.values[2]).toBeCloseTo(0, 8);
                    const p = mesh.positionPublic(i);
                    expect(dot(n, Vector.fromArray(
                        [p.values[0] / radius, p.values[1] / radius, 0])))
                        .toBeGreaterThan(Math.cos(Math.PI / cols) - 1e-12);
                }
            }
        }
    });
});

describe('MeshChannel', () => {
    it('addresses interleaved float64 vertex records', () => {
        // Six numbers per vertex: position followed by normal.
        const data = new Float64Array(6 * 4);
        const positions = new MeshChannel(data, 48, 3);
        const normals = new MeshChannel(data.subarray(3), 48, 3);
        expect(positions.stride).toBe(48);
        expect(positions.numComponents).toBe(3);

        for (let i = 0; i < 4; ++i) {
            positions.set(i, Vector.fromArray([i, 2 * i, 3 * i]));
            normals.set(i, Vector.fromArray([0, 0, 1]));
        }
        for (let i = 0; i < 4; ++i) {
            expect(positions.get(i).values).toEqual([i, 2 * i, 3 * i]);
            expect(normals.get(i).values).toEqual([0, 0, 1]);
            expect(positions.getComponent(i, 1)).toBe(2 * i);
            expect(data[6 * i]).toBe(i);
            expect(data[6 * i + 5]).toBe(1);
        }

        positions.setComponent(2, 2, -7);
        expect(positions.get(2).values[2]).toBe(-7);
    });

    it('uses float32 elements when the source is a Float32Array', () => {
        const data = new Float32Array(3 * 2);
        const channel = new MeshChannel(data, 12, 3);
        channel.set(0, Vector.fromArray([0.1, 0.2, 0.3]));
        expect(data[0]).toBe(Math.fround(0.1));
        expect(channel.get(0).values[0]).toBe(Math.fround(0.1));
        expect(channel.get(0).values[0]).not.toBe(0.1);
    });

    it('allocates tightly packed storage and rejects invalid strides', () => {
        const channel = MeshChannel.allocate(5, 2);
        expect(channel.stride).toBe(16);
        expect(channel.numComponents).toBe(2);
        channel.set(4, Vector.fromArray([1, 2]));
        expect(channel.get(4).values).toEqual([1, 2]);
        expect(channel.get(0).values).toEqual([0, 0]);

        expect(() => new MeshChannel(new Float64Array(8), 0, 3))
            .toThrow('The stride must be positive.');
        expect(() => new MeshChannel(new Float64Array(8), 24, 0))
            .toThrow('The number of components must be positive.');
        expect(() => new MeshChannel(new Float64Array(8), 20, 3))
            .toThrow('The source offset and stride must be multiples of the element size.');
        expect(() => MeshChannel.fromAttribute(new VertexAttribute('position'), 3))
            .toThrow('The attribute source is null.');
    });

    it('works with an ArrayBuffer source', () => {
        const buffer = new ArrayBuffer(8 * 6);
        const channel = new MeshChannel(buffer, 24, 3);
        channel.set(1, Vector.fromArray([4, 5, 6]));
        expect(new Float64Array(buffer)[3]).toBe(4);
        expect(channel.get(1).values).toEqual([4, 5, 6]);
    });
});

describe('Mesh.update', () => {
    it('does nothing beyond positions when no frame data is present', () => {
        const rows = 3, cols = 3;
        const storage = makeStorage(MeshTopology.RECTANGLE, rows, cols,
            ['position']);
        const mesh = new SurfaceMesh(storage.description, [MeshTopology.RECTANGLE]);
        mesh.surface = (i: number) => Vector.fromArray([i, 0, 0]);
        mesh.computeIndicesPublic();
        mesh.update();
        for (let i = 0; i < storage.description.numVertices; ++i) {
            expect(mesh.positionPublic(i).values).toEqual([i, 0, 0]);
        }
    });

    it('supports the ARBITRARY topology with client-supplied triangles', () => {
        const description = MeshDescription.arbitrary(4, 2);
        const indices = new Uint32Array([0, 1, 2, 0, 2, 3]);
        description.indexAttribute = new IndexAttribute(indices, 4);
        const positions = new Float64Array(12);
        const normals = new Float64Array(12);
        description.vertexAttributes.push(
            new VertexAttribute('position', positions, 24),
            new VertexAttribute('normal', normals, 24));

        const mesh = new SurfaceMesh(description, [MeshTopology.ARBITRARY]);
        expect(mesh.getDescription().constructed).toBe(true);
        const corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
        mesh.surface = (i: number) => Vector.fromArray(corners[i]);
        mesh.update();

        for (let i = 0; i < 4; ++i) {
            expect(mesh.normalPublic(i).values[2]).toBeCloseTo(1, 12);
        }
    });
});
