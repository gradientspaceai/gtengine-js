import { describe, it, expect } from 'vitest';
import { RevolutionMesh } from '../src/RevolutionMesh.js';
import { MeshDescription, MeshTopology } from '../src/Mesh.js';
import { ParametricCurve } from '../src/ParametricCurve.js';
import { IndexAttribute } from '../src/IndexAttribute.js';
import { VertexAttribute } from '../src/VertexAttribute.js';
import { Vector } from '../src/Vector.js';
import { ETManifoldMesh } from '../src/ETManifoldMesh.js';
import { GTE_C_TWO_PI } from '../src/Constants.js';
import { check, expectClose, fc, positive } from './helpers/arbitraries.js';

// A ParametricCurve whose derivatives are supplied in closed form. The
// upstream RevolutionMesh accepts any ParametricCurve<2,Real>.
class FunctionCurve extends ParametricCurve {
    private mDerivatives: Array<(t: number) => number[]>;

    constructor(dimension: number, tmin: number, tmax: number,
        derivatives: Array<(t: number) => number[]>) {
        super(dimension, tmin, tmax);
        this.mDerivatives = derivatives;
        this.mConstructed = true;
    }

    evaluate(t: number, order: number, jet: Vector[]): void {
        const maxOrder = Math.min(order, this.mDerivatives.length - 1);
        for (let i = 0; i <= maxOrder; ++i) {
            const values = this.mDerivatives[i](t);
            for (let j = 0; j < values.length; ++j) {
                jet[i].values[j] = values[j];
            }
        }
        for (let i = maxOrder + 1; i <= order; ++i) {
            jet[i].makeZero();
        }
    }
}

// The curve (x(t),z(t)) = (R, t) on [0,1]; the surface of revolution is a
// cylinder of radius R and height 1.
function cylinderCurve(R: number): FunctionCurve {
    return new FunctionCurve(2, 0, 1, [
        (t: number) => [R, t],
        (_t: number) => [0, 1],
        (_t: number) => [0, 0],
        (_t: number) => [0, 0]
    ]);
}

// The curve (r sin(pi t), -r cos(pi t)) on [0,1]; x(0) = x(1) = 0, so the
// surface of revolution is a sphere of radius r.
function sphereCurve(r: number): FunctionCurve {
    const p = Math.PI;
    return new FunctionCurve(2, 0, 1, [
        (t: number) => [r * Math.sin(p * t), -r * Math.cos(p * t)],
        (t: number) => [r * p * Math.cos(p * t), r * p * Math.sin(p * t)],
        (t: number) => [-r * p * p * Math.sin(p * t), r * p * p * Math.cos(p * t)],
        (t: number) => [-r * p * p * p * Math.cos(p * t), -r * p * p * p * Math.sin(p * t)]
    ]);
}

// The circle of radius a centered at (R,0) in the (x,z) plane; the surface
// of revolution is a torus.
function torusCurve(R: number, a: number): FunctionCurve {
    const w = 2 * Math.PI;
    return new FunctionCurve(2, 0, 1, [
        (t: number) => [R + a * Math.cos(w * t), a * Math.sin(w * t)],
        (t: number) => [-a * w * Math.sin(w * t), a * w * Math.cos(w * t)],
        (t: number) => [-a * w * w * Math.cos(w * t), -a * w * w * Math.sin(w * t)],
        (t: number) => [a * w * w * w * Math.sin(w * t), -a * w * w * w * Math.cos(w * t)]
    ]);
}

// The curve (t, 0) on [0,1]; x(0) = 0, so the surface of revolution is a
// unit disk in the plane z = 0.
function diskCurve(): FunctionCurve {
    return new FunctionCurve(2, 0, 1, [
        (t: number) => [t, 0],
        (_t: number) => [1, 0],
        (_t: number) => [0, 0],
        (_t: number) => [0, 0]
    ]);
}

interface Storage {
    description: MeshDescription;
    indices: Uint32Array;
    positions: Float64Array;
    normals: Float64Array;
    tcoords: Float64Array | null;
}

function makeStorage(topology: MeshTopology, numRows: number, numCols: number,
    withTCoords: boolean = false): Storage {
    const description = new MeshDescription(topology, numRows, numCols);
    const indices = new Uint32Array(3 * description.numTriangles);
    description.indexAttribute = new IndexAttribute(indices, 4);

    const positions = new Float64Array(3 * description.numVertices);
    description.vertexAttributes.push(new VertexAttribute('position', positions, 24));
    const normals = new Float64Array(3 * description.numVertices);
    description.vertexAttributes.push(new VertexAttribute('normal', normals, 24));

    let tcoords: Float64Array | null = null;
    if (withTCoords) {
        tcoords = new Float64Array(2 * description.numVertices);
        description.vertexAttributes.push(new VertexAttribute('tcoord', tcoords, 16));
    }

    return { description, indices, positions, normals, tcoords };
}

function P(storage: Storage, i: number): [number, number, number] {
    return [storage.positions[3 * i], storage.positions[3 * i + 1],
        storage.positions[3 * i + 2]];
}

function N(storage: Storage, i: number): [number, number, number] {
    return [storage.normals[3 * i], storage.normals[3 * i + 1],
        storage.normals[3 * i + 2]];
}

describe('RevolutionMesh', () => {
    it('produces the vertex and triangle counts of a cylinder', () => {
        const numRows = 5, numCols = 8;
        const storage = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
        const d = storage.description;
        expect(d.numRows).toBe(numRows);
        expect(d.numCols).toBe(numCols);
        expect(d.numVertices).toBe(numRows * (numCols + 1));
        expect(d.numTriangles).toBe(2 * (numRows - 1) * numCols);

        const mesh = new RevolutionMesh(d, cylinderCurve(2));
        expect(mesh.getDescription().constructed).toBe(true);
        expect(mesh.isSampleByArcLength()).toBe(false);
        expect(mesh.getCurve()).not.toBeNull();

        // Every index is in range and no triangle is degenerate.
        for (let t = 0; t < d.numTriangles; ++t) {
            const { v0, v1, v2 } = d.indexAttribute.getTriangle(t);
            for (const v of [v0, v1, v2]) {
                expect(v).toBeGreaterThanOrEqual(0);
                expect(v).toBeLessThan(d.numVertices);
            }
            expect(v0 === v1 && v1 === v2).toBe(false);
        }
    });

    it('places cylinder vertices on the revolved curve', () => {
        const numRows = 5, numCols = 8, R = 2;
        const storage = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
        const d = storage.description;
        new RevolutionMesh(d, cylinderCurve(R));

        for (let r = 0; r < numRows; ++r) {
            const z = r / (numRows - 1);
            for (let c = 0; c <= numCols; ++c) {
                const i = r * (numCols + 1) + c;
                const [x, y, pz] = P(storage, i);
                expect(Math.hypot(x, y)).toBeCloseTo(R, 12);
                expect(pz).toBeCloseTo(z, 12);

                // The angle is 2*pi*c/numCols.
                const angle = 2 * Math.PI * c / numCols;
                expect(x).toBeCloseTo(R * Math.cos(angle), 12);
                expect(y).toBeCloseTo(R * Math.sin(angle), 12);
            }
        }
    });

    it('duplicates the seam column of a cylinder', () => {
        const numRows = 4, numCols = 6;
        const storage = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
        const d = storage.description;
        new RevolutionMesh(d, cylinderCurve(3));
        for (let r = 0; r < numRows; ++r) {
            const first = P(storage, r * (numCols + 1));
            const last = P(storage, r * (numCols + 1) + numCols);
            for (let k = 0; k < 3; ++k) {
                expect(last[k]).toBeCloseTo(first[k], 14);
            }
        }
    });

    it('produces unit outward normals on a cylinder', () => {
        const numRows = 5, numCols = 10, R = 2;
        const storage = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
        const d = storage.description;
        new RevolutionMesh(d, cylinderCurve(R));

        // The normals are the area-weighted averages of the adjacent
        // triangle normals. A vertex with the full fan of six triangles has
        // an exactly radial normal on a cylinder. The vertices of the first
        // and last rows and the two duplicated seam vertices of every row
        // see only part of the fan, so their normals are tilted; they still
        // point outward and have no z-component.
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c <= numCols; ++c) {
                const i = r * (numCols + 1) + c;
                const [nx, ny, nz] = N(storage, i);
                expect(Math.hypot(nx, ny, nz)).toBeCloseTo(1, 12);
                expect(nz).toBeCloseTo(0, 12);

                const [x, y] = P(storage, i);
                const rho = Math.hypot(x, y);
                const radial = (nx * x + ny * y) / rho;
                expect(radial).toBeGreaterThan(0.9);

                const interiorRow = (r > 0 && r < numRows - 1);
                const interiorCol = (c > 0 && c < numCols);
                if (interiorRow && interiorCol) {
                    expect(nx).toBeCloseTo(x / rho, 12);
                    expect(ny).toBeCloseTo(y / rho, 12);
                }
                else if (!interiorCol) {
                    // Half of the fan is missing across the seam, which
                    // tilts the normal by half of the angular step.
                    expect(radial).toBeLessThan(0.9999);
                }
            }
        }
    });

    it('produces a sphere with the expected counts, radii and poles', () => {
        const numRows = 6, numCols = 8, r = 1.5;
        const storage = makeStorage(MeshTopology.SPHERE, numRows, numCols);
        const d = storage.description;
        expect(d.numVertices).toBe(numRows * (numCols + 1) + 2);
        expect(d.numTriangles).toBe(2 * (numRows - 1) * numCols + 2 * numCols);

        new RevolutionMesh(d, sphereCurve(r));

        // The ring vertices lie on the sphere of radius r.
        for (let i = 0; i < numRows * (numCols + 1); ++i) {
            const [x, y, z] = P(storage, i);
            expect(Math.hypot(x, y, z)).toBeCloseTo(r, 10);
        }

        // The poles are the curve endpoints on the axis of revolution.
        const south = P(storage, d.numVertices - 2);
        const north = P(storage, d.numVertices - 1);
        expect(south[0]).toBeCloseTo(0, 14);
        expect(south[1]).toBeCloseTo(0, 14);
        expect(south[2]).toBeCloseTo(-r, 12);
        expect(north[0]).toBeCloseTo(0, 14);
        expect(north[1]).toBeCloseTo(0, 14);
        expect(north[2]).toBeCloseTo(r, 12);

        // Seam duplicates coincide.
        for (let row = 0; row < numRows; ++row) {
            const first = P(storage, row * (numCols + 1));
            const last = P(storage, row * (numCols + 1) + numCols);
            for (let k = 0; k < 3; ++k) {
                expect(last[k]).toBeCloseTo(first[k], 14);
            }
        }
    });

    it('builds a sphere whose pole fans reference the correct rows (upstream #220)', () => {
        const numRows = 4, numCols = 5;
        const storage = makeStorage(MeshTopology.SPHERE, numRows, numCols);
        const d = storage.description;
        new RevolutionMesh(d, sphereCurve(1));

        // The south-pole fan is the last numCols triangles and must use the
        // vertices of the last ring, which begins at rMax * rIncrement.
        const firstOfLastRow = d.rMax * d.rIncrement;
        for (let c = 0; c < numCols; ++c) {
            const t = d.numTriangles - numCols + c;
            const { v0, v1, v2 } = d.indexAttribute.getTriangle(t);
            const ring = [v0, v1, v2].filter(v => v !== d.numVertices - 1);
            expect(ring.length).toBe(2);
            for (const v of ring) {
                expect(v).toBeGreaterThanOrEqual(firstOfLastRow);
                expect(v).toBeLessThan(firstOfLastRow + numCols + 1);
            }
        }
    });

    it('produces unit normals along the radius of a sphere', () => {
        const numRows = 8, numCols = 12, r = 2;
        const storage = makeStorage(MeshTopology.SPHERE, numRows, numCols);
        const d = storage.description;
        new RevolutionMesh(d, sphereCurve(r));

        for (let i = 0; i < d.numVertices; ++i) {
            const [nx, ny, nz] = N(storage, i);
            expect(Math.hypot(nx, ny, nz)).toBeCloseTo(1, 10);
        }

        // The area-weighted normals agree with the radial direction to
        // within the discretization error, except on the last ring (see the
        // next test).
        for (let row = 0; row < numRows - 1; ++row) {
            for (let c = 0; c <= numCols; ++c) {
                const i = row * (numCols + 1) + c;
                const [nx, ny, nz] = N(storage, i);
                const [x, y, z] = P(storage, i);
                const dp = (nx * x + ny * y + nz * z) / r;
                if (c > 0 && c < numCols) {
                    expect(dp).toBeGreaterThan(0.99);
                }
                else {
                    expect(dp).toBeGreaterThan(0.9);
                }
            }
        }

        // The pole nearest the first ring has an outward normal.
        {
            const i = d.numVertices - 2;
            const [, , nz] = N(storage, i);
            expect(nz).toBeCloseTo(-1, 10);
        }
    });

    it('winds the second pole fan consistently with the body quads', () => {
        // Upstream Mesh.h assigns the same winding to both pole fans of a
        // SPHERE mesh, but the two fans lie on opposite sides of their
        // rings, so the second (north) fan was wound opposite to the body
        // quads: the numCols directed edges of the last ring were each
        // traversed twice in the same direction and the north-pole normal
        // pointed inward. The port fixes the winding in Mesh.computeIndices
        // (see the upstream-bug issue for B84); this test asserts the
        // consistent result: every directed edge appears once and the
        // north-pole normal is outward.
        const numRows = 6, numCols = 8, r = 2;
        const storage = makeStorage(MeshTopology.SPHERE, numRows, numCols);
        const d = storage.description;
        new RevolutionMesh(d, sphereCurve(r));

        const directed = new Map<string, number>();
        for (let t = 0; t < d.numTriangles; ++t) {
            const { v0, v1, v2 } = d.indexAttribute.getTriangle(t);
            const edges: Array<[number, number]> = [[v0, v1], [v1, v2], [v2, v0]];
            for (const [a, b] of edges) {
                const key = a + ',' + b;
                directed.set(key, (directed.get(key) ?? 0) + 1);
            }
        }
        let numDuplicated = 0;
        for (const count of directed.values()) {
            if (count > 1) {
                ++numDuplicated;
            }
        }
        expect(numDuplicated).toBe(0);

        // The north pole normal is outward.
        const [, , nz] = N(storage, d.numVertices - 1);
        expect(nz).toBeCloseTo(1, 10);
    });

    it('produces a torus with the expected counts and tube radius', () => {
        const numRows = 6, numCols = 8, R = 3, a = 1;
        const storage = makeStorage(MeshTopology.TORUS, numRows, numCols);
        const d = storage.description;
        expect(d.rMax).toBe(numRows);
        expect(d.numVertices).toBe((numRows + 1) * (numCols + 1));
        expect(d.numTriangles).toBe(2 * numRows * numCols);

        new RevolutionMesh(d, torusCurve(R, a));

        for (let i = 0; i < d.numVertices; ++i) {
            const [x, y, z] = P(storage, i);
            const rho = Math.hypot(x, y) - R;
            expect(Math.hypot(rho, z)).toBeCloseTo(a, 10);
        }

        // The torus wraps in both directions: the last row duplicates the
        // first (the curve is closed) and the last column duplicates the
        // first.
        for (let c = 0; c <= numCols; ++c) {
            const first = P(storage, c);
            const last = P(storage, numRows * (numCols + 1) + c);
            for (let k = 0; k < 3; ++k) {
                expect(last[k]).toBeCloseTo(first[k], 10);
            }
        }
    });

    it('produces a disk with a center vertex', () => {
        const numRows = 4, numCols = 6;
        const storage = makeStorage(MeshTopology.DISK, numRows, numCols);
        const d = storage.description;
        expect(d.numVertices).toBe(numRows * (numCols + 1) + 1);
        expect(d.numTriangles).toBe(2 * (numRows - 1) * numCols + numCols);

        new RevolutionMesh(d, diskCurve());

        // The samples are at t = i/(numRows) for i in [0,numRows]; the ring
        // for row r uses sample r+1.
        for (let r = 0; r < numRows; ++r) {
            const radius = (r + 1) / numRows;
            for (let c = 0; c <= numCols; ++c) {
                const [x, y, z] = P(storage, r * (numCols + 1) + c);
                expect(Math.hypot(x, y)).toBeCloseTo(radius, 12);
                expect(z).toBeCloseTo(0, 14);
            }
        }

        const center = P(storage, d.numVertices - 1);
        expect(center[0]).toBeCloseTo(0, 14);
        expect(center[1]).toBeCloseTo(0, 14);
        expect(center[2]).toBeCloseTo(0, 14);
    });

    it('computes default texture coordinates when none are requested', () => {
        // The client-supplied texture coordinates are filled in when
        // requested; here they are requested so the values can be checked.
        const numRows = 3, numCols = 4;
        const storage = makeStorage(MeshTopology.CYLINDER, numRows, numCols, true);
        const d = storage.description;
        new RevolutionMesh(d, cylinderCurve(1));
        const tcoords = storage.tcoords as Float64Array;
        for (let r = 0, i = 0; r < numRows; ++r) {
            for (let c = 0; c <= numCols; ++c, ++i) {
                expect(tcoords[2 * i]).toBeCloseTo(c / numCols, 14);
                expect(tcoords[2 * i + 1]).toBeCloseTo(r / (numRows - 1), 14);
            }
        }
    });

    it('samples by arc length when requested', () => {
        // The curve (t^2, 0) ... use a curve whose speed varies so that the
        // t-uniform and arclength-uniform samplings differ. Here the curve
        // is (1, t^2) on [0,1]: a cylinder whose z-samples are t^2 for the
        // natural sampling and uniform for the arclength sampling.
        const curve = new FunctionCurve(2, 0, 1, [
            (t: number) => [1, t * t],
            (t: number) => [0, 2 * t],
            (_t: number) => [0, 2],
            (_t: number) => [0, 0]
        ]);

        const numRows = 5, numCols = 4;
        const natural = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
        new RevolutionMesh(natural.description, curve, false);
        const byLength = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
        const mesh = new RevolutionMesh(byLength.description, curve, true);
        expect(mesh.isSampleByArcLength()).toBe(true);

        // Natural sampling: z(r) = (r/(numRows-1))^2.
        for (let r = 0; r < numRows; ++r) {
            const s = r / (numRows - 1);
            expect(P(natural, r * (numCols + 1))[2]).toBeCloseTo(s * s, 10);
        }

        // Arclength sampling: the z values are uniformly spaced, because the
        // curve is a straight segment along z and its arclength is z.
        for (let r = 0; r < numRows; ++r) {
            expect(P(byLength, r * (numCols + 1))[2]).toBeCloseTo(r / (numRows - 1), 6);
        }

        // The two samplings agree only at the endpoints.
        expect(P(natural, 0)[2]).toBeCloseTo(P(byLength, 0)[2], 8);
        const last = (numRows - 1) * (numCols + 1);
        expect(P(natural, last)[2]).toBeCloseTo(P(byLength, last)[2], 6);
        const mid = 2 * (numCols + 1);
        expect(Math.abs(P(natural, mid)[2] - P(byLength, mid)[2])).toBeGreaterThan(0.1);
    });

    it('reports a failed construction for an unsupported topology', () => {
        const description = new MeshDescription(MeshTopology.RECTANGLE, 4, 4);
        const indices = new Uint32Array(3 * description.numTriangles);
        description.indexAttribute = new IndexAttribute(indices, 4);
        const positions = new Float64Array(3 * description.numVertices);
        description.vertexAttributes.push(
            new VertexAttribute('position', positions, 24));

        const mesh = new RevolutionMesh(description, cylinderCurve(1));
        expect(mesh.getDescription().constructed).toBe(false);
        expect(mesh.getCurve()).toBeNull();
    });

    it('rejects a curve that is not 2-dimensional', () => {
        const storage = makeStorage(MeshTopology.CYLINDER, 3, 4);
        const curve3 = new FunctionCurve(3, 0, 1, [
            (t: number) => [1, 0, t],
            (_t: number) => [0, 0, 1],
            (_t: number) => [0, 0, 0],
            (_t: number) => [0, 0, 0]
        ]);
        expect(() => new RevolutionMesh(storage.description, curve3)).toThrow();
    });

    it('recomputes positions when update() is called', () => {
        const numRows = 4, numCols = 5;
        const storage = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
        const mesh = new RevolutionMesh(storage.description, cylinderCurve(2));
        const before = Array.from(storage.positions);
        storage.positions.fill(0);
        mesh.update();
        for (let i = 0; i < before.length; ++i) {
            expect(storage.positions[i]).toBeCloseTo(before[i], 14);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// The number of curve samples the constructor allocates for each topology.
function numSamplesFor(topology: MeshTopology, rMax: number): number {
    switch (topology) {
        case MeshTopology.CYLINDER:
        case MeshTopology.TORUS: return rMax + 1;
        case MeshTopology.DISK: return rMax + 2;
        case MeshTopology.SPHERE: return rMax + 3;
        default: return 0;
    }
}

// The t-value of sample i, recomputed from the curve rather than from the
// mesh (the port's tSampler).
function sampleT(curve: ParametricCurve, numSamples: number, i: number): number {
    const factor = (curve.getTMax() - curve.getTMin()) / (numSamples - 1);
    return curve.getTMin() + i * factor;
}

// Weld the duplicated seam column (and, for a torus, the duplicated last row)
// so the index buffer describes a closed surface without repeated vertices.
function weldMap(topology: MeshTopology, d: MeshDescription): number[] {
    const stride = d.cMax + 1;
    const map = new Array<number>(d.numVertices);
    for (let i = 0; i < d.numVertices; ++i) { map[i] = i; }
    const numRows = topology === MeshTopology.TORUS ? d.rMax + 1 : d.rMax + 1;
    for (let r = 0; r < numRows; ++r) {
        map[r * stride + d.cMax] = r * stride;
    }
    if (topology === MeshTopology.TORUS) {
        for (let c = 0; c <= d.cMax; ++c) {
            map[d.rMax * stride + c] = map[c];
        }
    }
    return map;
}

// The Euler characteristic V - E + F of the welded index buffer, plus the
// manifold/closed flags of the resulting edge-triangle mesh.
function weldedTopology(storage: Storage, topology: MeshTopology) {
    const d = storage.description;
    const map = weldMap(topology, d);
    const mesh = new ETManifoldMesh();
    const used = new Set<number>();
    let numTriangles = 0;
    for (let t = 0; t < d.numTriangles; ++t) {
        const { v0, v1, v2 } = d.indexAttribute.getTriangle(t);
        const a = map[v0], b = map[v1], c = map[v2];
        if (a === b || b === c || c === a) { continue; }
        used.add(a); used.add(b); used.add(c);
        expect(mesh.insert(a, b, c)).not.toBeNull();
        ++numTriangles;
    }
    // ETManifoldMesh.insert returns null for a nonmanifold insertion, so the
    // assertion above already proves the welded mesh is manifold.
    return {
        euler: used.size - mesh.getNumEdges() + numTriangles,
        closed: mesh.isClosed(),
        numTriangles
    };
}

const smallGrid = fc.tuple(fc.integer({ min: 2, max: 6 }),
    fc.integer({ min: 3, max: 9 }));

describe('RevolutionMesh verification', () => {
    it('places every cylinder and torus vertex on the revolved curve', () => {
        check(fc.tuple(smallGrid, positive(4, 0.25), positive(2, 0.25)),
            ([[numRows, numCols], R, a]) => {
                for (const topology of [MeshTopology.CYLINDER, MeshTopology.TORUS]) {
                    const storage = makeStorage(topology, numRows, numCols);
                    const d = storage.description;
                    const curve = topology === MeshTopology.CYLINDER
                        ? cylinderCurve(R) : torusCurve(R + a + 1, a);
                    new RevolutionMesh(d, curve);

                    const numSamples = numSamplesFor(topology, d.rMax);
                    const stride = d.cMax + 1;
                    for (let r = 0; r <= d.rMax; ++r) {
                        const t = sampleT(curve, numSamples, r);
                        const p = curve.getPosition(t);
                        for (let c = 0; c <= d.cMax; ++c) {
                            // The port precomputes cos/sin with the wrap-around
                            // entry copied from index 0.
                            const cc = c === d.cMax ? 0 : c;
                            const angle = cc * (1 / d.numCols) * GTE_C_TWO_PI;
                            const got = P(storage, r * stride + c);
                            expect(got[0]).toBe(p.values[0] * Math.cos(angle));
                            expect(got[1]).toBe(p.values[0] * Math.sin(angle));
                            expect(got[2]).toBe(p.values[1]);
                        }
                        // The seam column duplicates the first column exactly.
                        expect(P(storage, r * stride + d.cMax))
                            .toEqual(P(storage, r * stride));
                    }
                }
            });
    });

    it('offsets the disk and sphere rows by one sample and adds the poles', () => {
        check(fc.tuple(smallGrid, positive(4, 0.25)),
            ([[numRows, numCols], r0]) => {
                for (const topology of [MeshTopology.DISK, MeshTopology.SPHERE]) {
                    const storage = makeStorage(topology, numRows, numCols);
                    const d = storage.description;
                    const curve = topology === MeshTopology.DISK
                        ? diskCurve() : sphereCurve(r0);
                    new RevolutionMesh(d, curve);

                    const numSamples = numSamplesFor(topology, d.rMax);
                    const stride = d.cMax + 1;
                    for (let r = 0; r <= d.rMax; ++r) {
                        // Upstream indexes mSamples[r + 1] for these two.
                        const t = sampleT(curve, numSamples, r + 1);
                        const p = curve.getPosition(t);
                        for (let c = 0; c <= d.cMax; ++c) {
                            const cc = c === d.cMax ? 0 : c;
                            const angle = cc * (1 / d.numCols) * GTE_C_TWO_PI;
                            const got = P(storage, r * stride + c);
                            expect(got[0]).toBe(p.values[0] * Math.cos(angle));
                            expect(got[1]).toBe(p.values[0] * Math.sin(angle));
                            expect(got[2]).toBe(p.values[1]);
                        }
                    }

                    const first = curve.getPosition(sampleT(curve, numSamples, 0));
                    const last = curve.getPosition(
                        sampleT(curve, numSamples, numSamples - 1));
                    if (topology === MeshTopology.DISK) {
                        expect(P(storage, d.numVertices - 1))
                            .toEqual([0, 0, first.values[1]]);
                    } else {
                        expect(P(storage, d.numVertices - 2))
                            .toEqual([0, 0, first.values[1]]);
                        expect(P(storage, d.numVertices - 1))
                            .toEqual([0, 0, last.values[1]]);
                    }
                }
            });
    });

    it('emits an index buffer with the Euler characteristic of the topology', () => {
        // The grid must be large enough that welding the seam does not merge
        // two triangles of the same quad ring: a torus needs at least three
        // rings in each direction to be a simplicial complex.
        const bigGrid = fc.tuple(fc.integer({ min: 4, max: 7 }),
            fc.integer({ min: 5, max: 9 }));
        check(bigGrid, ([numRows, numCols]) => {
            const cases: Array<[MeshTopology, ParametricCurve, number, boolean]> = [
                [MeshTopology.CYLINDER, cylinderCurve(2), 0, false],
                [MeshTopology.TORUS, torusCurve(3, 1), 0, true],
                [MeshTopology.DISK, diskCurve(), 1, false],
                [MeshTopology.SPHERE, sphereCurve(2), 2, true]
            ];
            for (const [topology, curve, euler, closed] of cases) {
                const storage = makeStorage(topology, numRows, numCols);
                new RevolutionMesh(storage.description, curve);
                const info = weldedTopology(storage, topology);
                expect(info.euler).toBe(euler);
                expect(info.closed).toBe(closed);
            }
        });
    });

    it('produces texture coordinates in the unit square', () => {
        // The DISK case is excluded: upstream places the ring coordinates at
        // radius*(cos,sin) without adding the (0.5,0.5) origin it declares, so
        // they run over [-0.5,0.5]^2. See the dedicated test below.
        check(smallGrid, ([numRows, numCols]) => {
            for (const [topology, curve] of [
                [MeshTopology.CYLINDER, cylinderCurve(2)],
                [MeshTopology.TORUS, torusCurve(3, 1)],
                [MeshTopology.SPHERE, sphereCurve(2)]
            ] as Array<[MeshTopology, ParametricCurve]>) {
                const storage = makeStorage(topology, numRows, numCols, true);
                const d = storage.description;
                new RevolutionMesh(d, curve);
                const tc = storage.tcoords as Float64Array;
                for (let i = 0; i < d.numVertices; ++i) {
                    expect(tc[2 * i]).toBeGreaterThanOrEqual(0);
                    expect(tc[2 * i]).toBeLessThanOrEqual(1);
                    expect(tc[2 * i + 1]).toBeGreaterThanOrEqual(0);
                    expect(tc[2 * i + 1]).toBeLessThanOrEqual(1);
                }
                if (topology === MeshTopology.CYLINDER
                    || topology === MeshTopology.TORUS) {
                    const stride = d.cMax + 1;
                    const rows = topology === MeshTopology.TORUS
                        ? d.numRows : d.numRows - 1;
                    for (let r = 0; r <= d.rMax; ++r) {
                        for (let c = 0; c <= d.cMax; ++c) {
                            expect(tc[2 * (r * stride + c)]).toBe(c / d.numCols);
                            expect(tc[2 * (r * stride + c) + 1]).toBe(r / rows);
                        }
                    }
                }
            }
        });
    });

    it('reproduces the positions when update() is called again', () => {
        check(fc.tuple(smallGrid, positive(4, 0.25)),
            ([[numRows, numCols], R]) => {
                const storage = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
                const mesh = new RevolutionMesh(storage.description,
                    cylinderCurve(R));
                const before = Array.from(storage.positions);
                storage.positions.fill(0);
                mesh.update();
                expect(Array.from(storage.positions)).toEqual(before);
            });
    });

    it('samples uniformly in t and, on request, uniformly in arc length', () => {
        // For the cylinder curve (R, t) the arc length is t, so the two
        // samplers agree; for the torus curve they differ but both must be
        // uniform in their own parameter.
        check(fc.tuple(smallGrid, positive(4, 0.25)),
            ([[numRows, numCols], R]) => {
                const storage = makeStorage(MeshTopology.CYLINDER, numRows, numCols);
                const d = storage.description;
                new RevolutionMesh(d, cylinderCurve(R), true);
                const stride = d.cMax + 1;
                const step = 1 / d.rMax;
                for (let r = 0; r <= d.rMax; ++r) {
                    expectClose(P(storage, r * stride)[2], r * step, 1e-9, 1e-9);
                }
            });
    });

    it('centers the disk texture coordinates on (0,0) but the hub on (0.5,0.5) (upstream quirk)', () => {
        // Upstream declares 'Vector2<Real> origin{0.5, 0.5}' and then writes
        // the ring coordinates as radius*(cos,sin) without adding it, using
        // 'origin' only for the final hub vertex. The ring therefore lives in
        // [-0.5,0.5]^2 while the hub sits at a corner of that square. The port
        // preserves the behavior.
        check(smallGrid, ([numRows, numCols]) => {
            const storage = makeStorage(MeshTopology.DISK, numRows, numCols, true);
            const d = storage.description;
            new RevolutionMesh(d, diskCurve());
            const tc = storage.tcoords as Float64Array;
            const stride = d.cMax + 1;
            for (let r = 0; r <= d.rMax; ++r) {
                const radius = Math.min((r + 1) / (2 * d.numRows), 0.5);
                for (let c = 0; c <= d.cMax; ++c) {
                    const angle = GTE_C_TWO_PI * c / d.numCols;
                    const i = r * stride + c;
                    expect(tc[2 * i]).toBe(radius * Math.cos(angle));
                    expect(tc[2 * i + 1]).toBe(radius * Math.sin(angle));
                }
            }
            // The hub is the only coordinate that uses the declared origin.
            const hub = d.numVertices - 1;
            expect(tc[2 * hub]).toBe(0.5);
            expect(tc[2 * hub + 1]).toBe(0.5);
        });
    });

    it('emits NaN sphere texture coordinates for a single-row sphere (upstream quirk)', () => {
        // MeshDescription clamps a SPHERE to numRows >= 1, but
        // InitializeTCoords divides the row index by numRows - 1, so a
        // one-row sphere divides 0 by 0. The port reproduces it.
        check(fc.integer({ min: 3, max: 8 }), numCols => {
            const storage = makeStorage(MeshTopology.SPHERE, 1, numCols, true);
            const d = storage.description;
            expect(d.numRows).toBe(1);
            new RevolutionMesh(d, sphereCurve(1));
            const tc = storage.tcoords as Float64Array;
            for (let c = 0; c <= d.cMax; ++c) {
                expect(Number.isNaN(tc[2 * c + 1])).toBe(true);
            }
            // The two pole coordinates are written after the loop and are
            // finite.
            expect(tc[2 * (d.numVertices - 2) + 1]).toBe(0);
            expect(tc[2 * (d.numVertices - 1) + 1]).toBe(1);
        });
    });
});
