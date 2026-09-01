import { describe, it, expect } from 'vitest';
import { TubeMesh } from '../src/TubeMesh';
import { MeshDescription, MeshTopology } from '../src/Mesh';
import { ParametricCurve } from '../src/ParametricCurve';
import { FrenetFrame3 } from '../src/FrenetFrame';
import { IndexAttribute } from '../src/IndexAttribute';
import { VertexAttribute } from '../src/VertexAttribute';
import { Vector, dot, length, sub } from '../src/Vector';

// A ParametricCurve whose derivatives are supplied in closed form.
class FunctionCurve3 extends ParametricCurve {
    private mDerivatives: Array<(t: number) => number[]>;

    constructor(tmin: number, tmax: number,
        derivatives: Array<(t: number) => number[]>, dimension: number = 3) {
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

// p(t) = (0,0,t) on [0,1].
function lineCurve(): FunctionCurve3 {
    return new FunctionCurve3(0, 1, [
        (t: number) => [0, 0, t],
        (_t: number) => [0, 0, 1],
        (_t: number) => [0, 0, 0],
        (_t: number) => [0, 0, 0]
    ]);
}

// p(t) = (0,0,t^2) on [0,1]: a straight segment with nonuniform speed.
function acceleratingLineCurve(): FunctionCurve3 {
    return new FunctionCurve3(0, 1, [
        (t: number) => [0, 0, t * t],
        (t: number) => [0, 0, 2 * t],
        (_t: number) => [0, 0, 2],
        (_t: number) => [0, 0, 0]
    ]);
}

// p(t) = (R cos(2 pi t), R sin(2 pi t), 0) on [0,1]: a closed circle.
function circleCurve(R: number): FunctionCurve3 {
    const w = 2 * Math.PI;
    return new FunctionCurve3(0, 1, [
        (t: number) => [R * Math.cos(w * t), R * Math.sin(w * t), 0],
        (t: number) => [-R * w * Math.sin(w * t), R * w * Math.cos(w * t), 0],
        (t: number) => [-R * w * w * Math.cos(w * t), -R * w * w * Math.sin(w * t), 0],
        (t: number) => [R * w * w * w * Math.sin(w * t), -R * w * w * w * Math.cos(w * t), 0]
    ]);
}

// A helix, used to exercise the Frenet-frame sampler.
function helixCurve(): FunctionCurve3 {
    const w = 2 * Math.PI;
    return new FunctionCurve3(0, 1, [
        (t: number) => [Math.cos(w * t), Math.sin(w * t), t],
        (t: number) => [-w * Math.sin(w * t), w * Math.cos(w * t), 1],
        (t: number) => [-w * w * Math.cos(w * t), -w * w * Math.sin(w * t), 0],
        (t: number) => [w * w * w * Math.sin(w * t), -w * w * w * Math.cos(w * t), 0]
    ]);
}

interface Storage {
    description: MeshDescription;
    indices: Uint32Array;
    positions: Float64Array;
    normals: Float64Array;
    tcoords: Float64Array | null;
}

function makeStorage(numRows: number, numCols: number,
    withTCoords: boolean = false,
    topology: MeshTopology = MeshTopology.CYLINDER): Storage {
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

function P(storage: Storage, i: number): Vector {
    return Vector.fromArray([storage.positions[3 * i],
        storage.positions[3 * i + 1], storage.positions[3 * i + 2]]);
}

const ZERO = Vector.fromArray([0, 0, 0]);
const UPY = Vector.fromArray([0, 1, 0]);
const UPZ = Vector.fromArray([0, 0, 1]);

describe('TubeMesh', () => {
    it('produces the vertex and triangle counts of a cylinder', () => {
        const numRows = 5, numCols = 9;
        const storage = makeStorage(numRows, numCols);
        const d = storage.description;
        expect(d.numVertices).toBe(numRows * (numCols + 1));
        expect(d.numTriangles).toBe(2 * (numRows - 1) * numCols);

        const mesh = new TubeMesh(d, lineCurve(), () => 0.25, false, false, UPY);
        expect(mesh.getDescription().constructed).toBe(true);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.isSampleByArcLength()).toBe(false);
        expect(mesh.getMedial()).not.toBeNull();
        expect(mesh.getRadial()(0.5)).toBe(0.25);
        expect(mesh.getUpVector().values).toEqual([0, 1, 0]);

        for (let t = 0; t < d.numTriangles; ++t) {
            const { v0, v1, v2 } = d.indexAttribute.getTriangle(t);
            for (const v of [v0, v1, v2]) {
                expect(v).toBeGreaterThanOrEqual(0);
                expect(v).toBeLessThan(d.numVertices);
            }
        }
    });

    it('places the vertices at the tube radius around the medial curve', () => {
        const numRows = 6, numCols = 9, radius = 0.3;
        const storage = makeStorage(numRows, numCols);
        const d = storage.description;
        const medial = lineCurve();
        new TubeMesh(d, medial, () => radius, false, false, UPY);

        for (let r = 0; r < numRows; ++r) {
            const t = r / (numRows - 1);
            const center = medial.getPosition(t);
            const tangent = medial.getTangent(t);
            for (let c = 0; c <= numCols; ++c) {
                const p = P(storage, r * (numCols + 1) + c);
                const offset = sub(p, center);
                expect(length(offset)).toBeCloseTo(radius, 12);
                // The ring lies in the plane through the medial point that
                // is perpendicular to the tangent.
                expect(dot(offset, tangent)).toBeCloseTo(0, 12);
            }
        }
    });

    it('samples the ring at numCols-1 distinct angles with two duplicates', () => {
        // The upstream ring has numCols samples, the last of which repeats
        // the first, and updatePositions() appends one more repeat of the
        // first to close the seam.
        const numRows = 3, numCols = 7, radius = 0.5;
        const storage = makeStorage(numRows, numCols);
        const d = storage.description;
        const medial = lineCurve();
        new TubeMesh(d, medial, () => radius, false, false, UPY);

        // The frame for a line with tangent (0,0,1) and up (0,1,0) is
        // binormal = tangent x up = (-1,0,0) and normal = binormal x
        // tangent = (0,1,0).
        for (let r = 0; r < numRows; ++r) {
            const t = r / (numRows - 1);
            for (let c = 0; c < numCols; ++c) {
                const angle = 2 * Math.PI * (c % (numCols - 1)) / (numCols - 1);
                const p = P(storage, r * (numCols + 1) + c);
                expect(p.values[0]).toBeCloseTo(-radius * Math.sin(angle), 12);
                expect(p.values[1]).toBeCloseTo(radius * Math.cos(angle), 12);
                expect(p.values[2]).toBeCloseTo(t, 12);
            }

            // The last two ring vertices repeat the first.
            const first = P(storage, r * (numCols + 1));
            for (const c of [numCols - 1, numCols]) {
                const dup = P(storage, r * (numCols + 1) + c);
                expect(length(sub(dup, first))).toBeCloseTo(0, 12);
            }
        }
    });

    it('honors a varying radial function', () => {
        const numRows = 5, numCols = 9;
        const radial = (t: number) => 0.2 + 0.5 * t;
        const storage = makeStorage(numRows, numCols);
        const medial = lineCurve();
        new TubeMesh(storage.description, medial, radial, false, false, UPY);

        for (let r = 0; r < numRows; ++r) {
            const t = r / (numRows - 1);
            const center = medial.getPosition(t);
            for (let c = 0; c <= numCols; ++c) {
                const offset = sub(P(storage, r * (numCols + 1) + c), center);
                expect(length(offset)).toBeCloseTo(radial(t), 12);
            }
        }
    });

    it('uses the Frenet frame when the up vector is zero', () => {
        const numRows = 7, numCols = 9, radius = 0.2;
        const storage = makeStorage(numRows, numCols);
        const medial = helixCurve();
        new TubeMesh(storage.description, medial, () => radius, false, false, ZERO);

        const frenet = new FrenetFrame3(medial);
        for (let r = 0; r < numRows; ++r) {
            const t = r / (numRows - 1);
            const frame = frenet.compute(t);
            for (let c = 0; c < numCols - 1; ++c) {
                const angle = 2 * Math.PI * c / (numCols - 1);
                const expected = Vector.fromArray([0, 0, 0]);
                for (let k = 0; k < 3; ++k) {
                    expected.values[k] = frame.position.values[k] + radius *
                        (Math.cos(angle) * frame.normal.values[k] +
                            Math.sin(angle) * frame.binormal.values[k]);
                }
                const p = P(storage, r * (numCols + 1) + c);
                expect(length(sub(p, expected))).toBeCloseTo(0, 12);
            }
        }
    });

    it('closes the tube without corrupting the interior rows (upstream stride)', () => {
        // Upstream indexes the last row with a stride of numCols instead of
        // rIncrement = numCols + 1, which writes over interior rows. With
        // the corrected stride, every row of a closed tube is a ring of the
        // requested radius about its own medial point, and the last row
        // repeats the first.
        const numRows = 6, numCols = 7, radius = 0.25, R = 3;
        const storage = makeStorage(numRows, numCols);
        const d = storage.description;
        const medial = circleCurve(R);
        const mesh = new TubeMesh(d, medial, () => radius, true, false, UPZ);
        expect(mesh.isClosed()).toBe(true);

        // For the closed case the parameter step is (tmax-tmin)/numRows.
        for (let r = 0; r < numRows - 1; ++r) {
            const t = r / numRows;
            const center = medial.getPosition(t);
            const tangent = medial.getTangent(t);
            for (let c = 0; c <= numCols; ++c) {
                const offset = sub(P(storage, r * (numCols + 1) + c), center);
                expect(length(offset)).toBeCloseTo(radius, 10);
                expect(dot(offset, tangent)).toBeCloseTo(0, 10);
            }
        }

        // The last row repeats the first, seam vertex included.
        for (let c = 0; c <= numCols; ++c) {
            const first = P(storage, c);
            const last = P(storage, (numRows - 1) * (numCols + 1) + c);
            expect(length(sub(last, first))).toBeCloseTo(0, 14);
        }
    });

    it('leaves the last row distinct when the tube is open', () => {
        const numRows = 5, numCols = 7, radius = 0.25, R = 3;
        const storage = makeStorage(numRows, numCols);
        const d = storage.description;
        new TubeMesh(d, circleCurve(R), () => radius, false, false, UPZ);

        // With closed = false the parameter step is (tmax-tmin)/(numRows-1),
        // so the last row is at t = tmax. The circle is closed, so the last
        // ring coincides with the first up to the frame at t = 1 versus
        // t = 0, which are the same for a circle.
        for (let c = 0; c <= numCols; ++c) {
            const first = P(storage, c);
            const last = P(storage, (numRows - 1) * (numCols + 1) + c);
            expect(length(sub(last, first))).toBeCloseTo(0, 8);
        }

        // The interior rows are at distinct medial points.
        const p0 = P(storage, 0);
        const p1 = P(storage, numCols + 1);
        expect(length(sub(p1, p0))).toBeGreaterThan(0.5);
    });

    it('samples the medial curve by arc length when requested', () => {
        const numRows = 5, numCols = 7, radius = 0.1;
        const medial = acceleratingLineCurve();

        const natural = makeStorage(numRows, numCols);
        new TubeMesh(natural.description, medial, () => radius, false, false, UPY);
        const byLength = makeStorage(numRows, numCols);
        const mesh = new TubeMesh(byLength.description, medial, () => radius,
            false, true, UPY);
        expect(mesh.isSampleByArcLength()).toBe(true);

        // The medial curve is (0,0,t^2), so the natural sampling gives
        // z = (r/(numRows-1))^2 and the arclength sampling gives uniformly
        // spaced z (the arclength of the segment is z).
        for (let r = 0; r < numRows; ++r) {
            const s = r / (numRows - 1);
            expect(P(natural, r * (numCols + 1)).values[2]).toBeCloseTo(s * s, 10);
            expect(P(byLength, r * (numCols + 1)).values[2]).toBeCloseTo(s, 5);
        }

        // The two samplings differ away from the endpoints.
        const mid = 2 * (numCols + 1);
        expect(Math.abs(P(natural, mid).values[2] -
            P(byLength, mid).values[2])).toBeGreaterThan(0.1);
    });

    it('computes texture coordinates over the row/column grid', () => {
        const numRows = 4, numCols = 5;
        const storage = makeStorage(numRows, numCols, true);
        new TubeMesh(storage.description, lineCurve(), () => 0.5, false, false, UPY);
        const tcoords = storage.tcoords as Float64Array;
        for (let r = 0, i = 0; r < numRows; ++r) {
            for (let c = 0; c <= numCols; ++c, ++i) {
                expect(tcoords[2 * i]).toBeCloseTo(c / numCols, 14);
                expect(tcoords[2 * i + 1]).toBeCloseTo(r / (numRows - 1), 14);
            }
        }
    });

    it('produces unit normals', () => {
        const numRows = 6, numCols = 9, radius = 0.4;
        const storage = makeStorage(numRows, numCols);
        const d = storage.description;
        new TubeMesh(d, lineCurve(), () => radius, false, false, UPY);
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c <= numCols; ++c) {
                const i = r * (numCols + 1) + c;
                const n = Vector.fromArray([storage.normals[3 * i],
                    storage.normals[3 * i + 1], storage.normals[3 * i + 2]]);
                if (c === numCols) {
                    // The ring has only numCols-1 distinct angles: column
                    // numCols-1 repeats column 0 and updatePositions()
                    // appends one more repeat. The quad between columns
                    // numCols-1 and numCols is therefore degenerate, so the
                    // area-weighted normal of the trailing column is zero.
                    // This is upstream behavior.
                    expect(length(n)).toBeCloseTo(0, 14);
                    continue;
                }
                expect(length(n)).toBeCloseTo(1, 12);
                // The tube around a straight line has normals perpendicular
                // to the line, pointing away from it.
                expect(n.values[2]).toBeCloseTo(0, 12);
                const p = P(storage, i);
                const radial = Vector.fromArray([p.values[0], p.values[1], 0]);
                expect(dot(n, radial) / radius).toBeGreaterThan(0.9);
            }
        }
    });

    it('reports a failed construction for an unsupported topology', () => {
        const storage = makeStorage(4, 5, false, MeshTopology.TORUS);
        const mesh = new TubeMesh(storage.description, lineCurve(), () => 1,
            false, false, UPY);
        expect(mesh.getDescription().constructed).toBe(false);
        expect(mesh.getMedial()).toBeNull();
    });

    it('rejects a medial curve that is not 3-dimensional', () => {
        const storage = makeStorage(4, 5);
        const curve2 = new FunctionCurve3(0, 1, [
            (t: number) => [1, t],
            (_t: number) => [0, 1],
            (_t: number) => [0, 0],
            (_t: number) => [0, 0]
        ], 2);
        expect(() => new TubeMesh(storage.description, curve2, () => 1,
            false, false, UPY)).toThrow();
    });

    it('recomputes positions when update() is called', () => {
        const numRows = 4, numCols = 6;
        const storage = makeStorage(numRows, numCols);
        const mesh = new TubeMesh(storage.description, lineCurve(), () => 0.3,
            false, false, UPY);
        const before = Array.from(storage.positions);
        storage.positions.fill(0);
        mesh.update();
        for (let i = 0; i < before.length; ++i) {
            expect(storage.positions[i]).toBeCloseTo(before[i], 14);
        }
    });
});
