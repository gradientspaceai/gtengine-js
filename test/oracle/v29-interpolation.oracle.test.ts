// Replays oracle/cpp/cases/v29-interpolation.cpp. Keep the two files in the
// same order. See the C++ file comment for the mesh adapters: the port's
// Delaunay2Mesh / Delaunay3Mesh already number the simplices in sorted
// TriangleKey / TetrahedronKey order and start searches at simplex 0, which
// is what the C++ SortedMesh2 / SortedMesh3 reproduce over upstream's
// triangulation, so both the "real mesh" and the "sortedMesh" cases replay
// through the port's own mesh classes.
import { describe } from 'vitest';
import { BSRational } from '../../src/BSRational.js';
import { Delaunay2 } from '../../src/Delaunay2.js';
import { Delaunay2Mesh } from '../../src/Delaunay2Mesh.js';
import { Delaunay3 } from '../../src/Delaunay3.js';
import { Delaunay3Mesh } from '../../src/Delaunay3Mesh.js';
import { IntpAkimaNonuniform1 } from '../../src/IntpAkimaNonuniform1.js';
import { IntpAkimaUniform1 } from '../../src/IntpAkimaUniform1.js';
import {
    IntpBSplineUniform, IntpBSplineUniform1, IntpBSplineUniform2, IntpBSplineUniform3,
    type IntpBSplineUniformControls
} from '../../src/IntpBSplineUniform.js';
import {
    IntpLinearNonuniform2, type IntpLinearNonuniform2TriangleMesh
} from '../../src/IntpLinearNonuniform2.js';
import {
    IntpLinearNonuniform3, type IntpLinearNonuniform3TetrahedronMesh
} from '../../src/IntpLinearNonuniform3.js';
import {
    IntpQuadraticNonuniform2, type IntpQuadraticNonuniform2TriangleMesh
} from '../../src/IntpQuadraticNonuniform2.js';
import { IntpSphere2 } from '../../src/IntpSphere2.js';
import { IntpThinPlateSpline2 } from '../../src/IntpThinPlateSpline2.js';
import { IntpThinPlateSpline3 } from '../../src/IntpThinPlateSpline3.js';
import { IntpTricubic3 } from '../../src/IntpTricubic3.js';
import { IntpTrilinear3 } from '../../src/IntpTrilinear3.js';
import { IntpVectorField2 } from '../../src/IntpVectorField2.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- uniform 3D grids -----------------------------------------------------
interface Grid3 {
    bound: number[];
    min: number[];
    spacing: number[];
    F: number[];
}

function grid3(io: OracleIO): Grid3 {
    const bound = [io.integer(), io.integer(), io.integer()];
    const min: number[] = [];
    const spacing: number[] = [];
    for (let d = 0; d < 3; ++d) {
        min.push(io.real());
        spacing.push(io.real());
    }
    const F = io.reals(bound[0] * bound[1] * bound[2]);
    return { bound, min, spacing, F };
}

interface GridIntp {
    getXBound(): number; getYBound(): number; getZBound(): number;
    getQuantity(): number;
    getXMin(): number; getXMax(): number; getXSpacing(): number;
    getYMin(): number; getYMax(): number; getYSpacing(): number;
    getZMin(): number; getZMax(): number; getZSpacing(): number;
}

function emitGridAccessors(io: OracleIO, intp: GridIntp): void {
    io.outInt(intp.getXBound());
    io.outInt(intp.getYBound());
    io.outInt(intp.getZBound());
    io.outInt(intp.getQuantity());
    io.outReal(intp.getXMin());
    io.outReal(intp.getXMax());
    io.outReal(intp.getXSpacing());
    io.outReal(intp.getYMin());
    io.outReal(intp.getYMax());
    io.outReal(intp.getYSpacing());
    io.outReal(intp.getZMin());
    io.outReal(intp.getZMax());
    io.outReal(intp.getZSpacing());
}

// ---- Akima -------------------------------------------------------------------
// The port's subclasses expose the protected polynomials like the C++
// AkimaExposed template.
class AkimaUniform1Exposed extends IntpAkimaUniform1 {
    coeff(i: number, k: number): number { return this.mPoly[i].getCoeff(k); }
}

class AkimaNonuniform1Exposed extends IntpAkimaNonuniform1 {
    coeff(i: number, k: number): number { return this.mPoly[i].getCoeff(k); }
}

function emitAkima(io: OracleIO, intp: AkimaUniform1Exposed | AkimaNonuniform1Exposed): void {
    const q = intp.getQuantity();
    io.outInt(q);
    io.outReal(intp.getXMin());
    io.outReal(intp.getXMax());
    for (let i = 0; i + 1 < q; ++i) {
        for (let k = 0; k < 4; ++k) {
            io.outReal(intp.coeff(i, k));
        }
    }
    for (let j = 0; j < 5; ++j) {
        const x = io.real();
        io.outReal(intp.evaluate(x));
        for (let order = -1; order <= 4; ++order) {
            io.outReal(intp.evaluate(order, x));
        }
    }
}

// ---- B-splines -------------------------------------------------------------
// The C++ GridControls: a lexicographic array, dimension 0 fastest.
class GridControls implements IntpBSplineUniformControls<number> {
    constructor(readonly size: number[], readonly data: number[]) {}

    getSize(dimension: number): number { return this.size[dimension]; }

    get(indices: readonly number[]): number {
        const n = this.size.length;
        let index = indices[n - 1];
        for (let d = n - 1; d-- > 0;) {
            index = index * this.size[d] + indices[d];
        }
        return this.data[index];
    }

    add(c0: number, c1: number): number { return c0 + c1; }

    mul(c0: number, s: number): number { return c0 * s; }
}

function controls(io: OracleIO, size: number[]): GridControls {
    const q = size.reduce((a, b) => a * b, 1);
    return new GridControls(size, io.reals(q));
}

describe('oracle: v29-interpolation', () => {
    const family = new OracleFamily('v29-interpolation');

    family.case('IntpTrilinear3.evaluate', (io) => {
        const g = grid3(io);
        const intp = new IntpTrilinear3(g.bound[0], g.bound[1], g.bound[2], g.min[0],
            g.spacing[0], g.min[1], g.spacing[1], g.min[2], g.spacing[2], g.F);
        emitGridAccessors(io, intp);
        for (let q = 0; q < 3; ++q) {
            const x = io.real(), y = io.real(), z = io.real();
            io.outReal(intp.evaluate(x, y, z));
            for (let zo = -1; zo <= 2; ++zo) {
                for (let yo = -1; yo <= 2; ++yo) {
                    for (let xo = -1; xo <= 2; ++xo) {
                        io.outReal(intp.evaluate(xo, yo, zo, x, y, z));
                    }
                }
            }
        }
    }, { exact: true });

    family.case('IntpTricubic3.evaluate', (io) => {
        const g = grid3(io);
        const catmullRom = io.boolean();
        const intp = new IntpTricubic3(g.bound[0], g.bound[1], g.bound[2], g.min[0],
            g.spacing[0], g.min[1], g.spacing[1], g.min[2], g.spacing[2], g.F, catmullRom);
        emitGridAccessors(io, intp);
        for (let q = 0; q < 2; ++q) {
            const x = io.real(), y = io.real(), z = io.real();
            io.outReal(intp.evaluate(x, y, z));
            for (let zo = 0; zo <= 3; ++zo) {
                for (let yo = 0; yo <= 3; ++yo) {
                    for (let xo = 0; xo <= 3; ++xo) {
                        io.outReal(intp.evaluate(xo, yo, zo, x, y, z));
                    }
                }
            }
            io.outReal(intp.evaluate(4, 0, 0, x, y, z));
            io.outReal(intp.evaluate(0, -1, 0, x, y, z));
            io.outReal(intp.evaluate(1, 2, 4, x, y, z));
        }
    }, { exact: true });

    family.case('IntpTrilinear3.construct.invalid', (io) => {
        const g = grid3(io);
        const intp = new IntpTrilinear3(g.bound[0], g.bound[1], g.bound[2], g.min[0],
            g.spacing[0], g.min[1], g.spacing[1], g.min[2], g.spacing[2], g.F);
        io.outInt(intp.getQuantity());
        io.outReal(intp.getXMax());
        io.outReal(intp.getYMax());
        io.outReal(intp.getZMax());
        io.outReal(intp.evaluate(g.min[0] + 0.25, g.min[1] + 0.5, g.min[2] + 0.75));
    }, { exact: true });

    family.case('IntpTricubic3.construct.invalid', (io) => {
        const g = grid3(io);
        const catmullRom = io.boolean();
        const intp = new IntpTricubic3(g.bound[0], g.bound[1], g.bound[2], g.min[0],
            g.spacing[0], g.min[1], g.spacing[1], g.min[2], g.spacing[2], g.F, catmullRom);
        io.outInt(intp.getQuantity());
        io.outReal(intp.getXMax());
        io.outReal(intp.getYMax());
        io.outReal(intp.getZMax());
        io.outReal(intp.evaluate(g.min[0] + 0.25, g.min[1] + 0.5, g.min[2] + 0.75));
    }, { exact: true });

    family.case('IntpAkimaUniform1.evaluate', (io) => {
        const q = io.integer();
        const xMin = io.real();
        const spacing = io.real();
        const F = io.reals(q);
        const intp = new AkimaUniform1Exposed(q, xMin, spacing, F);
        io.outReal(intp.getXSpacing());
        emitAkima(io, intp);
    }, { exact: true });

    family.case('IntpAkimaNonuniform1.evaluate', (io) => {
        const q = io.integer();
        const X = io.reals(q);
        const F = io.reals(q);
        emitAkima(io, new AkimaNonuniform1Exposed(q, X, F));
    }, { exact: true });

    family.case('IntpAkimaUniform1.construct.invalid', (io) => {
        const q = io.integer();
        const xMin = io.real();
        const spacing = io.real();
        const F = io.reals(q);
        const intp = new IntpAkimaUniform1(q, xMin, spacing, F);
        io.outReal(intp.getXMax());
        io.outReal(intp.evaluate(xMin + 0.75));
    }, { exact: true });

    family.case('IntpAkimaNonuniform1.construct.invalid', (io) => {
        const q = io.integer();
        const X = io.reals(q);
        const F = io.reals(q);
        const intp = new IntpAkimaNonuniform1(q, X, F);
        io.outReal(intp.getXMax());
        io.outReal(intp.evaluate(X[0] + 0.5));
    }, { exact: true });

    family.case('IntpBSplineUniform1.evaluate', (io) => {
        const degree = io.integer();
        const n = io.integer();
        const cacheMode = io.integer();
        const c = controls(io, [n]);
        const ctZero = io.real();
        const intp = new IntpBSplineUniform1<number>(degree, c, ctZero, cacheMode);
        io.outInt(intp.getDegree(0));
        io.outInt(intp.getNumControls(0));
        io.outReal(intp.getTMin(0));
        io.outReal(intp.getTMax(0));
        io.outInt(intp.getCacheMode());
        for (let j = 0; j < 6; ++j) {
            const t = io.real();
            for (let order = -1; order <= degree + 1; ++order) {
                io.outReal(intp.evaluate([order], [t]));
            }
        }
    }, { exact: true });

    family.case('IntpBSplineUniform2.evaluate', (io) => {
        const degrees = [io.integer(), io.integer()];
        const sizes = [io.integer(), io.integer()];
        const cacheMode = io.integer();
        const c = controls(io, sizes);
        const ctZero = io.real();
        const intp = new IntpBSplineUniform2<number>(degrees, c, ctZero, cacheMode);
        for (let d = 0; d < 2; ++d) {
            io.outInt(intp.getDegree(d));
            io.outInt(intp.getNumControls(d));
            io.outReal(intp.getTMin(d));
            io.outReal(intp.getTMax(d));
        }
        io.outInt(intp.getCacheMode());
        for (let j = 0; j < 4; ++j) {
            const t = [io.real(), io.real()];
            for (let o1 = -1; o1 <= degrees[1] + 1; ++o1) {
                for (let o0 = -1; o0 <= degrees[0] + 1; ++o0) {
                    io.outReal(intp.evaluate([o0, o1], t));
                }
            }
        }
    }, { exact: true });

    family.case('IntpBSplineUniform3.evaluate', (io) => {
        const degrees = [io.integer(), io.integer(), io.integer()];
        const sizes = [io.integer(), io.integer(), io.integer()];
        const cacheMode = io.integer();
        const c = controls(io, sizes);
        const ctZero = io.real();
        const intp = new IntpBSplineUniform3<number>(degrees, c, ctZero, cacheMode);
        for (let d = 0; d < 3; ++d) {
            io.outInt(intp.getDegree(d));
            io.outInt(intp.getNumControls(d));
            io.outReal(intp.getTMin(d));
            io.outReal(intp.getTMax(d));
        }
        io.outInt(intp.getCacheMode());
        for (let j = 0; j < 3; ++j) {
            const t = io.reals(3);
            for (let o2 = 0; o2 <= degrees[2] + 1; ++o2) {
                for (let o1 = 0; o1 <= degrees[1] + 1; ++o1) {
                    for (let o0 = 0; o0 <= degrees[0] + 1; ++o0) {
                        io.outReal(intp.evaluate([o0, o1, o2], t));
                    }
                }
            }
            io.outReal(intp.evaluate([-1, 0, 0], t));
        }
    }, { exact: true });

    family.case('IntpBSplineUniform.evaluate.general', (io) => {
        const N = 1 + io.index % 4;
        const degrees: number[] = [];
        for (let d = 0; d < N; ++d) { degrees.push(io.integer()); }
        const sizes: number[] = [];
        for (let d = 0; d < N; ++d) { sizes.push(io.integer()); }
        const cacheMode = io.integer();
        const c = controls(io, sizes);
        const ctZero = io.real();
        const intp = new IntpBSplineUniform<number>(degrees, c, ctZero, cacheMode);
        for (let d = 0; d < N; ++d) {
            io.outInt(intp.getDegree(d));
            io.outInt(intp.getNumControls(d));
            io.outReal(intp.getTMin(d));
            io.outReal(intp.getTMax(d));
        }
        io.outInt(intp.getCacheMode());
        for (let j = 0; j < 3; ++j) {
            const t = io.reals(N);
            const order = new Array<number>(N).fill(0);
            for (;;) {
                io.outReal(intp.evaluate(order, t));
                let d = 0;
                for (; d < N; ++d) {
                    if (++order[d] <= degrees[d] + 1) { break; }
                    order[d] = 0;
                }
                if (d === N) { break; }
            }
            io.outReal(intp.evaluate(new Array<number>(N - 1).fill(0), t));
            io.outReal(intp.evaluate(new Array<number>(N).fill(0),
                new Array<number>(N - 1).fill(0)));
            const negative = new Array<number>(N).fill(0);
            negative[N - 1] = -1;
            io.outReal(intp.evaluate(negative, t));
        }
    }, { exact: true });

    family.case('IntpBSplineUniform.evaluate.fixedN', (io) => {
        const degrees = [io.integer(), io.integer(), io.integer(), io.integer()];
        const sizes = [io.integer(), io.integer(), io.integer(), io.integer()];
        const cacheMode = io.integer();
        const c = controls(io, sizes);
        const ctZero = io.real();
        const intp = new IntpBSplineUniform<number>(degrees, c, ctZero, cacheMode);
        for (let j = 0; j < 3; ++j) {
            const t = io.reals(4);
            const order = [io.integer(), io.integer(), io.integer(), io.integer()];
            io.outReal(intp.evaluate(order, t));
            io.outReal(intp.evaluate([0, 0, 0, 0], t));
        }
    }, { exact: true });

    family.case('IntpBSplineUniform.construct.invalid', (io) => {
        const which = io.integer();
        const degree = io.integer();
        const n = io.integer();
        const n2 = io.integer();
        const cacheMode = io.integer();
        const sizes = (which === 0 || which === 3) ? [n] : (which === 1 ? [n2, n] : [n2, n2, n]);
        const c = controls(io, sizes);
        let value: number;
        if (which === 0) {
            value = new IntpBSplineUniform1<number>(degree, c, 0, cacheMode).evaluate([0], [0.25]);
        } else if (which === 1) {
            value = new IntpBSplineUniform2<number>([degree, degree], c, 0, cacheMode)
                .evaluate([0, 0], [0.25, 0.25]);
        } else if (which === 2) {
            value = new IntpBSplineUniform3<number>([degree, degree, degree], c, 0, cacheMode)
                .evaluate([0, 0, 0], [0.25, 0.25, 0.25]);
        } else {
            value = new IntpBSplineUniform<number>([degree], c, 0, cacheMode).evaluate([0], [0.25]);
        }
        io.outReal(value);
    }, { exact: true });

    // ---- nonuniform linear interpolation -------------------------------------
    function points(io: OracleIO, dim: number): Vector[] {
        const n = io.integer();
        const pts: Vector[] = [];
        for (let i = 0; i < n; ++i) { pts.push(io.vec(dim)); }
        return pts;
    }

    function mesh2(pts: Vector[]): Delaunay2Mesh {
        const del = new Delaunay2();
        del.compute(pts);
        return new Delaunay2Mesh(del);
    }

    function mesh3(pts: Vector[]): Delaunay3Mesh {
        const del = new Delaunay3();
        del.compute(pts);
        return new Delaunay3Mesh(del);
    }

    function linear2(io: OracleIO, mesh: IntpLinearNonuniform2TriangleMesh, F: number[],
        numQueries: number): void {
        const intp = new IntpLinearNonuniform2(mesh, F);
        for (let q = 0; q < numQueries; ++q) {
            const r = intp.evaluate(io.vec(2));
            io.outBool(r.valid);
            if (r.valid) { io.outReal(r.F); }
        }
    }

    function linear3(io: OracleIO, mesh: IntpLinearNonuniform3TetrahedronMesh, F: number[],
        numQueries: number): void {
        const intp = new IntpLinearNonuniform3(mesh, F);
        for (let q = 0; q < numQueries; ++q) {
            const r = intp.evaluate(io.vec(3));
            io.outBool(r.valid);
            if (r.valid) { io.outReal(r.F); }
        }
    }

    family.case('IntpLinearNonuniform2.evaluate', (io) => {
        const pts = points(io, 2);
        const F = io.reals(pts.length);
        linear2(io, mesh2(pts), F, 6);
    }, { exact: true });

    family.case('IntpLinearNonuniform2.evaluate.sortedMesh', (io) => {
        const pts = points(io, 2);
        const F = io.reals(pts.length);
        linear2(io, mesh2(pts), F, 6);
    }, { exact: true });

    family.case('IntpLinearNonuniform2.deviation.getIndices', (io) => {
        const pts = points(io, 2);
        const F = io.reals(pts.length);
        const mesh = mesh2(pts);
        linear2(io, {
            getContainingTriangle: (P) => mesh.getContainingTriangle(P),
            getBarycentrics: (t, P) => mesh.getBarycentrics(t, P),
            getTriangleIndices: () => null
        }, F, 4);
    }, { deviation: '#135 (IntpLinearNonuniform2 discards the GetIndices failure flag)' });

    family.case('IntpLinearNonuniform3.evaluate', (io) => {
        const pts = points(io, 3);
        const F = io.reals(pts.length);
        linear3(io, mesh3(pts), F, 6);
    }, { exact: true });

    family.case('IntpLinearNonuniform3.evaluate.sortedMesh', (io) => {
        const pts = points(io, 3);
        const F = io.reals(pts.length);
        linear3(io, mesh3(pts), F, 6);
    }, { exact: true });

    family.case('IntpLinearNonuniform3.deviation.getIndices', (io) => {
        const pts = points(io, 3);
        const F = io.reals(pts.length);
        const mesh = mesh3(pts);
        linear3(io, {
            getContainingTetrahedron: (P) => mesh.getContainingTetrahedron(P),
            getBarycentrics: (t, P) => mesh.getBarycentrics(t, P),
            getTetrahedronIndices: () => null
        }, F, 4);
    }, { deviation: '#135 (IntpLinearNonuniform3 discards the GetIndices failure flag)' });

    // ---- IntpQuadraticNonuniform2 ----------------------------------------------
    function gradients(io: OracleIO, n: number): { FX: number[], FY: number[] } {
        const FX: number[] = [], FY: number[] = [];
        for (let i = 0; i < n; ++i) {
            FX.push(io.real());
            FY.push(io.real());
        }
        return { FX, FY };
    }

    function emitQuadratic(io: OracleIO, intp: IntpQuadraticNonuniform2, P: Vector): void {
        const r = intp.evaluate(P);
        io.outBool(r.valid);
        if (r.valid) {
            io.outReal(r.F);
            io.outReal(r.FX);
            io.outReal(r.FY);
        }
    }

    function quadraticFromDerivatives(io: OracleIO): void {
        const pts = points(io, 2);
        const F = io.reals(pts.length);
        const { FX, FY } = gradients(io, pts.length);
        const intp = IntpQuadraticNonuniform2.fromDerivatives(mesh2(pts), F, FX, FY);
        for (let q = 0; q < 5; ++q) { emitQuadratic(io, intp, io.vec(2)); }
    }

    function quadraticFromSpatialDelta(io: OracleIO): void {
        const pts = points(io, 2);
        const F = io.reals(pts.length);
        const spatialDelta = io.real();
        const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh2(pts), F, spatialDelta);
        for (let q = 0; q < 5; ++q) { emitQuadratic(io, intp, io.vec(2)); }
    }

    family.case('IntpQuadraticNonuniform2.fromDerivatives', quadraticFromDerivatives,
        { exact: true });

    family.case('IntpQuadraticNonuniform2.fromDerivatives.sortedMesh', quadraticFromDerivatives,
        { exact: true });

    family.case('IntpQuadraticNonuniform2.fromSpatialDelta.sortedMesh', quadraticFromSpatialDelta,
        { exact: true });


    // The C++ FailingMesh2: one per-triangle accessor fails for every
    // triangle (mode 0 vertices, 1 indices, 2 adjacencies, 3 barycentrics).
    function failingMesh(mesh: Delaunay2Mesh, mode: number): IntpQuadraticNonuniform2TriangleMesh {
        return {
            getNumVertices: () => mesh.getNumVertices(),
            getNumTriangles: () => mesh.getNumTriangles(),
            getVertices: () => mesh.getVertices(),
            getIndices: () => mesh.getIndices(),
            getTriangleVertices: (t) => (mode === 0 ? null : mesh.getTriangleVertices(t)),
            getTriangleIndices: (t) => (mode === 1 ? null : mesh.getTriangleIndices(t)),
            getTriangleAdjacencies: (t) => (mode === 2 ? null : mesh.getTriangleAdjacencies(t)),
            getBarycentrics: (t, P) => (mode === 3 ? null : mesh.getBarycentrics(t, P)),
            getContainingTriangle: (P) => mesh.getContainingTriangle(P),
            getInvalidIndex: () => -1
        };
    }

    family.case('IntpQuadraticNonuniform2.deviation.meshFlags', (io) => {
        const mode = io.integer();
        const pts = points(io, 2);
        const F = io.reals(pts.length);
        const { FX, FY } = gradients(io, pts.length);
        const intp = IntpQuadraticNonuniform2.fromDerivatives(failingMesh(mesh2(pts), mode),
            F, FX, FY);
        for (let q = 0; q < 4; ++q) { emitQuadratic(io, intp, io.vec(2)); }
    }, { deviation: '#337 (IntpQuadraticNonuniform2 discards every mesh-accessor failure flag)' });

    // The C++ ListMesh2: explicit triangles and adjacencies, containment by
    // a scan for the first triangle with nonnegative exact barycentrics.
    function listMesh(V: Vector[], indices: number[], adjacencies: number[]):
        IntpQuadraticNonuniform2TriangleMesh {
        const tri = (t: number) => [indices[3 * t], indices[3 * t + 1], indices[3 * t + 2]];
        const bary = (t: number, P: Vector): number[] | null => {
            const r = (x: number) => BSRational.fromNumber(x);
            const [a, b, c] = tri(t).map((i) => V[i].values);
            const d0 = [r(a[0]).sub(r(c[0])), r(a[1]).sub(r(c[1]))];
            const d1 = [r(b[0]).sub(r(c[0])), r(b[1]).sub(r(c[1]))];
            const d2 = [r(P.values[0]).sub(r(c[0])), r(P.values[1]).sub(r(c[1]))];
            const det = d0[0].mul(d1[1]).sub(d0[1].mul(d1[0]));
            if (det.getSign() === 0) { return null; }
            const b0 = d2[0].mul(d1[1]).sub(d2[1].mul(d1[0])).div(det);
            const b1 = d0[0].mul(d2[1]).sub(d0[1].mul(d2[0])).div(det);
            const b2 = r(1).sub(b0).sub(b1);
            return [b0.toNumber(), b1.toNumber(), b2.toNumber()];
        };
        const numTriangles = indices.length / 3;
        return {
            getNumVertices: () => V.length,
            getNumTriangles: () => numTriangles,
            getVertices: () => V,
            getIndices: () => indices,
            getTriangleVertices: (t) => tri(t).map((i) => V[i]),
            getTriangleIndices: (t) => tri(t),
            getTriangleAdjacencies: (t) => adjacencies.slice(3 * t, 3 * t + 3),
            getBarycentrics: bary,
            getContainingTriangle: (P) => {
                for (let t = 0; t < numTriangles; ++t) {
                    const b = bary(t, P);
                    if (b !== null && b[0] >= 0 && b[1] >= 0 && b[2] >= 0) { return t; }
                }
                return -1;
            },
            getInvalidIndex: () => -1
        };
    }

    family.case('IntpQuadraticNonuniform2.degenerateTriangle', (io) => {
        const V0 = io.vec(2);
        const dir = io.vec(2);
        const V1 = Vector.fromArray([V0.get(0) + dir.get(0), V0.get(1) + dir.get(1)]);
        const k = io.integer() * 0.25;
        const V2 = Vector.fromArray([V0.get(0) + k * (V1.get(0) - V0.get(0)),
            V0.get(1) + k * (V1.get(1) - V0.get(1))]);
        io.real();  // h: V3 is recorded below
        const V3 = io.vec(2);
        const F = io.reals(4);
        const useSpatialDelta = io.boolean();
        const { FX, FY } = gradients(io, 4);
        const mesh = listMesh([V0, V1, V2, V3], [0, 2, 3, 2, 1, 3, 0, 1, 2],
            [2, 1, -1, 2, -1, 0, -1, 1, 0]);
        const intp = useSpatialDelta ? IntpQuadraticNonuniform2.fromSpatialDelta(mesh, F, 1)
            : IntpQuadraticNonuniform2.fromDerivatives(mesh, F, FX, FY);
        for (let q = 0; q < 4; ++q) { emitQuadratic(io, intp, io.vec(2)); }
    }, { exact: true });

    // ---- thin-plate splines -----------------------------------------------------
    // std::log in the 2D kernel (MSVC log vs Math.log, 1 ulp on ~3.6% of the
    // arguments), amplified by the two inverted systems. The generator keeps
    // cond1(A) * cond1(Q) <= 1e4, where the measured scaled error is below
    // 1e-12; compared at 1e-11. IsInitialized compares exactly.
    family.case('IntpThinPlateSpline2.evaluate', (io) => {
        const n = io.integer();
        const X = io.reals(n);
        const Y = io.reals(n);
        const F = io.reals(n);
        const smooth = io.real();
        const transform = io.boolean();
        const tps = new IntpThinPlateSpline2(n, X, Y, F, smooth, transform);
        io.outBool(tps.isInitialized());
        for (let i = 0; i < n; ++i) { io.outReal(tps.evaluate(X[i], Y[i])); }
        for (let q = 0; q < 3; ++q) {
            const x = io.real();
            const y = io.real();
            io.outReal(tps.evaluate(x, y));
        }
        io.outReal(tps.computeFunctional());
    }, { tol: 1e-11 });

    family.case('IntpThinPlateSpline3.evaluate', (io) => {
        const n = io.integer();
        const X = io.reals(n);
        const Y = io.reals(n);
        const Z = io.reals(n);
        const F = io.reals(n);
        const smooth = io.real();
        const transform = io.boolean();
        const tps = new IntpThinPlateSpline3(n, X, Y, Z, F, smooth, transform);
        io.outBool(tps.isInitialized());
        for (let i = 0; i < n; ++i) { io.outReal(tps.evaluate(X[i], Y[i], Z[i])); }
        for (let q = 0; q < 3; ++q) {
            const x = io.real();
            const y = io.real();
            const z = io.real();
            io.outReal(tps.evaluate(x, y, z));
        }
        io.outReal(tps.computeFunctional());
    }, { exact: true });

    family.case('IntpThinPlateSpline2.construct.invalid', (io) => {
        const n = io.integer();
        const X: number[] = [], Y: number[] = [], F: number[] = [];
        for (let i = 0; i < n; ++i) {
            X.push(io.real());
            Y.push(io.real());
            F.push(io.real());
        }
        const smooth = io.real();
        io.outBool(new IntpThinPlateSpline2(n, X, Y, F, smooth, false).isInitialized());
    }, { exact: true });

    family.case('IntpThinPlateSpline3.construct.invalid', (io) => {
        const n = io.integer();
        const X: number[] = [], Y: number[] = [], Z: number[] = [], F: number[] = [];
        for (let i = 0; i < n; ++i) {
            X.push(io.real());
            Y.push(io.real());
            Z.push(io.real());
            F.push(io.real());
        }
        const smooth = io.real();
        io.outBool(new IntpThinPlateSpline3(n, X, Y, Z, F, smooth, false).isInitialized());
    }, { exact: true });

    // ---- IntpSphere2 ---------------------------------------------------------------
    function sphereSamples(io: OracleIO): { theta: number[], phi: number[], F: number[] } {
        const n = io.integer();
        const theta = io.reals(n);
        const phi = io.reals(n);
        const F = io.reals(n);
        return { theta, phi, F };
    }

    // std::atan2 and std::acos (MSVC) against Math.atan2 and Math.acos.
    family.case('IntpSphere2.getSphericalCoordinates', (io) => {
        const v = io.reals(3);
        const r = IntpSphere2.getSphericalCoordinates(v[0], v[1], v[2]);
        io.outReal(r.theta);
        io.outReal(r.phi);
    });

    family.case('IntpSphere2.deviation.constructorThrows', (io) => {
        const s = sphereSamples(io);
        const queries: Vector[] = [];
        for (let q = 0; q < 5; ++q) { queries.push(io.vec(2)); }
        const intp = new IntpSphere2(s.theta, s.phi, s.F);
        for (const a of queries) {
            const r = intp.evaluate(a.get(0), a.get(1));
            io.outBool(r.valid);
            if (r.valid) { io.outReal(r.F); }
        }
    }, { deviation: 'v29 finding: IntpSphere2<T> constructs Delaunay2Mesh<T> before the triangulation' });

    family.case('IntpSphere2.evaluate.sortedMesh', (io) => {
        const s = sphereSamples(io);
        const intp = new IntpSphere2(s.theta, s.phi, s.F);
        for (let q = 0; q < 5; ++q) {
            const a = io.vec(2);
            const r = intp.evaluate(a.get(0), a.get(1));
            io.outBool(r.valid);
            if (r.valid) { io.outReal(r.F); }
        }
    }, { exact: true });

    // ---- IntpVectorField2 -------------------------------------------------------------
    function vectorField(io: OracleIO): IntpVectorField2 {
        const domain = points(io, 2);
        const range: Vector[] = [];
        for (let i = 0; i < domain.length; ++i) { range.push(io.vec(2)); }
        return new IntpVectorField2(domain, range);
    }

    family.case('IntpVectorField2.deviation.constructorThrows', (io) => {
        const domain = points(io, 2);
        const range: Vector[] = [];
        for (let i = 0; i < domain.length; ++i) { range.push(io.vec(2)); }
        const queries: Vector[] = [];
        for (let q = 0; q < 5; ++q) { queries.push(io.vec(2)); }
        const intp = new IntpVectorField2(domain, range);
        for (const P of queries) {
            const r = intp.evaluate(P);
            io.outBool(r.valid);
            if (r.valid) { io.outVec(r.output); }
        }
    }, { deviation: 'v29 finding: IntpVectorField2<T> constructs Delaunay2Mesh<T> before the triangulation' });

    family.case('IntpVectorField2.evaluate.sortedMesh', (io) => {
        const intp = vectorField(io);
        for (let q = 0; q < 5; ++q) {
            const r = intp.evaluate(io.vec(2));
            io.outBool(r.valid);
            if (r.valid) { io.outVec(r.output); }
        }
    }, { exact: true });

    family.case('IntpVectorField2.deviation.staleOutput', (io) => {
        const intp = vectorField(io);
        for (let q = 0; q < 3; ++q) {
            const r = intp.evaluate(io.vec(2));
            io.outBool(r.valid);
            io.outVec(r.output);
        }
    }, { deviation: '#337 (IntpVectorField2 leaves the caller output stale on failure)' });

    family.finish();
});
