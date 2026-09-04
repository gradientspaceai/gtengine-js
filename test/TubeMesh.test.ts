import { describe, it, expect } from 'vitest';
import { TubeMesh } from '../src/TubeMesh.js';
import { MeshDescription, MeshTopology } from '../src/Mesh.js';
import { ParametricCurve } from '../src/ParametricCurve.js';
import { FrenetFrame3 } from '../src/FrenetFrame.js';
import { IndexAttribute } from '../src/IndexAttribute.js';
import { VertexAttribute } from '../src/VertexAttribute.js';
import { Vector, dot, length, sub } from '../src/Vector.js';
import { ETManifoldMesh } from '../src/ETManifoldMesh.js';
import { check, expectClose, fc, positive } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// The t-value of row 'row', recomputed from the curve rather than from the
// mesh (the port's tSampler).
function tubeSampleT(curve: ParametricCurve, numRows: number, closed: boolean,
    row: number): number {
    const invDenom = closed ? 1 / numRows : 1 / (numRows - 1);
    const factor = (curve.getTMax() - curve.getTMin()) * invDenom;
    return curve.getTMin() + row * factor;
}

// Weld the duplicated ring columns (the last two both repeat column 0) and,
// for a closed tube, the duplicated last row.
function tubeWeldMap(d: MeshDescription, closed: boolean): number[] {
    const stride = d.cMax + 1;
    const map = new Array<number>(d.numVertices);
    for (let i = 0; i < d.numVertices; ++i) { map[i] = i; }
    for (let r = 0; r < d.numRows; ++r) {
        map[r * stride + d.numCols - 1] = r * stride;
        map[r * stride + d.numCols] = r * stride;
    }
    if (closed) {
        for (let c = 0; c <= d.cMax; ++c) {
            map[(d.numRows - 1) * stride + c] = map[c];
        }
    }
    return map;
}

function tubeWeldedEuler(storage: Storage, closed: boolean): number {
    const d = storage.description;
    const map = tubeWeldMap(d, closed);
    const mesh = new ETManifoldMesh();
    const used = new Set<number>();
    let numTriangles = 0;
    for (let t = 0; t < d.numTriangles; ++t) {
        const { v0, v1, v2 } = d.indexAttribute.getTriangle(t);
        const a = map[v0], b = map[v1], c = map[v2];
        if (a === b || b === c || c === a) { continue; }
        used.add(a); used.add(b); used.add(c);
        // ETManifoldMesh.insert returns null for a nonmanifold insertion.
        expect(mesh.insert(a, b, c)).not.toBeNull();
        ++numTriangles;
    }
    return used.size - mesh.getNumEdges() + numTriangles;
}

const tubeGrid = fc.tuple(fc.integer({ min: 3, max: 7 }),
    fc.integer({ min: 4, max: 9 }));

describe('TubeMesh verification', () => {
    it('puts every vertex at the radial distance from the medial curve, orthogonal to the tangent', () => {
        check(fc.tuple(tubeGrid, positive(3, 0.25), positive(2, 0.25)),
            ([[numRows, numCols], R, r0]) => {
                const radial = (t: number) => r0 * (1 + 0.5 * t);
                const storage = makeStorage(numRows, numCols);
                const d = storage.description;
                const medial = circleCurve(R + r0 + 1);
                new TubeMesh(d, medial, radial, false, false, UPZ);

                const stride = d.cMax + 1;
                for (let row = 0; row < d.numRows; ++row) {
                    const t = tubeSampleT(medial, d.numRows, false, row);
                    const center = medial.getPosition(t);
                    const tangent = medial.getTangent(t);
                    const expectedRadius = radial(t);
                    for (let col = 0; col <= d.numCols; ++col) {
                        const diff = sub(P(storage, row * stride + col), center);
                        expectClose(length(diff), expectedRadius, 1e-12, 1e-12);
                        expectClose(dot(diff, tangent), 0, 1e-12, 1e-12);
                    }
                    // The last two columns both repeat the first one, because
                    // mCosAngle[numCols-1] is copied from mCosAngle[0] and the
                    // extra seam vertex is an explicit copy.
                    expect(P(storage, row * stride + d.numCols - 1).values)
                        .toEqual(P(storage, row * stride).values);
                    expect(P(storage, row * stride + d.numCols).values)
                        .toEqual(P(storage, row * stride).values);
                }
            });
    });

    it('closes a closed tube by copying the whole first row, seam vertex included', () => {
        // Upstream indexes the last row as 'col + numCols * (numRows - 1)',
        // but the CYLINDER row stride is numCols + 1, so it overwrites an
        // interior row and leaves the real last row stale. The port uses the
        // description's rIncrement and copies every column of the ring.
        check(fc.tuple(tubeGrid, positive(3, 0.5)), ([[numRows, numCols], R]) => {
            const storage = makeStorage(numRows, numCols);
            const d = storage.description;
            const medial = circleCurve(R + 2);
            new TubeMesh(d, medial, () => 0.5, true, false, UPZ);

            const stride = d.cMax + 1;
            const last = (d.numRows - 1) * stride;
            for (let col = 0; col <= d.cMax; ++col) {
                expect(P(storage, last + col).values)
                    .toEqual(P(storage, col).values);
            }
            // No interior row was clobbered: every interior row still lies on
            // its own medial sample.
            for (let row = 1; row + 1 < d.numRows; ++row) {
                const t = tubeSampleT(medial, d.numRows, true, row);
                const center = medial.getPosition(t);
                for (let col = 0; col <= d.cMax; ++col) {
                    expectClose(length(sub(P(storage, row * stride + col),
                        center)), 0.5, 1e-12, 1e-12);
                }
            }
        });
    });

    it('leaves the last row on its own sample when the tube is open', () => {
        check(fc.tuple(tubeGrid, positive(2, 0.25)),
            ([[numRows, numCols], r0]) => {
                const storage = makeStorage(numRows, numCols);
                const d = storage.description;
                const medial = lineCurve();
                new TubeMesh(d, medial, () => r0, false, false, UPY);
                const stride = d.cMax + 1;
                for (let row = 0; row < d.numRows; ++row) {
                    const t = tubeSampleT(medial, d.numRows, false, row);
                    for (let col = 0; col <= d.cMax; ++col) {
                        // p(t) = (0,0,t), so the ring for row lies in z = t.
                        expectClose(P(storage, row * stride + col).values[2], t,
                            1e-12, 1e-12);
                    }
                }
                expect(P(storage, (d.numRows - 1) * stride).values[2])
                    .not.toBe(P(storage, 0).values[2]);
            });
    });

    it('emits an index buffer with the Euler characteristic of the tube', () => {
        // A welded open tube is an annulus (chi = 0) and a welded closed tube
        // is a torus (chi = 0); both must be manifold, which insert() checks.
        const bigTube = fc.tuple(fc.integer({ min: 5, max: 8 }),
            fc.integer({ min: 5, max: 9 }));
        check(bigTube, ([numRows, numCols]) => {
            for (const closed of [false, true]) {
                const storage = makeStorage(numRows, numCols);
                const medial = closed ? circleCurve(4) : lineCurve();
                new TubeMesh(storage.description, medial, () => 0.4, closed,
                    false, UPZ);
                expect(tubeWeldedEuler(storage, closed)).toBe(0);
            }
        });
    });

    it('scales the tube with the radial function', () => {
        check(fc.tuple(tubeGrid, positive(2, 0.25), positive(4, 0.5)),
            ([[numRows, numCols], r0, s]) => {
                const a = makeStorage(numRows, numCols);
                new TubeMesh(a.description, lineCurve(), () => r0, false,
                    false, UPY);
                const b = makeStorage(numRows, numCols);
                new TubeMesh(b.description, lineCurve(), () => s * r0, false,
                    false, UPY);
                const d = a.description;
                const stride = d.cMax + 1;
                for (let row = 0; row < d.numRows; ++row) {
                    const t = tubeSampleT(lineCurve(), d.numRows, false, row);
                    const center = Vector.fromArray([0, 0, t]);
                    for (let col = 0; col <= d.cMax; ++col) {
                        const i = row * stride + col;
                        const da = sub(P(a, i), center);
                        const db = sub(P(b, i), center);
                        for (let k = 0; k < 3; ++k) {
                            expectClose(db.values[k], s * da.values[k],
                                1e-12, 1e-12);
                        }
                    }
                }
            });
    });

    it('reproduces the positions when update() is called again', () => {
        check(fc.tuple(tubeGrid, positive(2, 0.25)),
            ([[numRows, numCols], r0]) => {
                const storage = makeStorage(numRows, numCols);
                const mesh = new TubeMesh(storage.description, helixCurve(),
                    () => r0, false, false, ZERO);
                const before = Array.from(storage.positions);
                storage.positions.fill(0);
                mesh.update();
                expect(Array.from(storage.positions)).toEqual(before);
            });
    });

    it('produces texture coordinates over the row/column grid', () => {
        check(tubeGrid, ([numRows, numCols]) => {
            const storage = makeStorage(numRows, numCols, true);
            const d = storage.description;
            new TubeMesh(d, lineCurve(), () => 0.5, false, false, UPY);
            const tc = storage.tcoords as Float64Array;
            const stride = d.cMax + 1;
            for (let r = 0; r < d.numRows; ++r) {
                for (let c = 0; c <= d.numCols; ++c) {
                    const i = r * stride + c;
                    expect(tc[2 * i]).toBe(c / d.numCols);
                    expect(tc[2 * i + 1]).toBe(r / d.rMax);
                }
            }
        });
    });

    it('agrees with the Frenet frame when the up vector is zero', () => {
        // With a zero up vector upstream builds a FrenetFrame3 of the medial
        // curve; the ring must then lie in the plane spanned by the Frenet
        // normal and binormal.
        check(fc.tuple(tubeGrid, positive(2, 0.25)),
            ([[numRows, numCols], r0]) => {
                const storage = makeStorage(numRows, numCols);
                const d = storage.description;
                const medial = helixCurve();
                new TubeMesh(d, medial, () => r0, false, false, ZERO);
                const frenet = new FrenetFrame3(medial);
                const stride = d.cMax + 1;
                for (let row = 0; row < d.numRows; ++row) {
                    const t = tubeSampleT(medial, d.numRows, false, row);
                    const frame = frenet.compute(t);
                    for (let col = 0; col <= d.cMax; ++col) {
                        const diff = sub(P(storage, row * stride + col),
                            frame.position);
                        expectClose(length(diff), r0, 1e-12, 1e-12);
                        expectClose(dot(diff, frame.tangent), 0, 1e-12, 1e-12);
                    }
                }
            });
    });
});
