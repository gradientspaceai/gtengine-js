// Replays oracle/cpp/cases/v28-interpolation.cpp. Keep the two files in the
// same order.
//
// Every case is exact: the only operations on these paths are + - * / and
// the truncating float-to-int conversion of the cell index (Math.trunc). The
// '.maxBoundary' cases of IntpAkimaUniform2/3 are deliberate deviations: the
// port fixes the sign of the max-boundary mixed-derivative stencils (#58).
import { describe } from 'vitest';
import { HermiteBicubic, HermiteBicubicSample } from '../../src/HermiteBicubic.js';
import { HermiteBiquintic, HermiteBiquinticSample } from '../../src/HermiteBiquintic.js';
import { HermiteCubic, HermiteCubicSample } from '../../src/HermiteCubic.js';
import { HermiteQuintic, HermiteQuinticSample } from '../../src/HermiteQuintic.js';
import { HermiteTricubic, HermiteTricubicSample } from '../../src/HermiteTricubic.js';
import { HermiteTriquintic, HermiteTriquinticSample } from '../../src/HermiteTriquintic.js';
import { IntpAkimaUniform1 } from '../../src/IntpAkimaUniform1.js';
import { IntpAkimaUniform2 } from '../../src/IntpAkimaUniform2.js';
import { IntpAkimaUniform3 } from '../../src/IntpAkimaUniform3.js';
import { IntpBicubic2 } from '../../src/IntpBicubic2.js';
import { IntpBilinear2 } from '../../src/IntpBilinear2.js';
import { OracleFamily, type OracleIO } from './harness.js';

const exact = { exact: true };

type Pair<T> = readonly [T, T];

function cubicBlocks(io: OracleIO): Pair<HermiteCubicSample> {
    const b0 = new HermiteCubicSample(io.real(), io.real());
    const b1 = new HermiteCubicSample(io.real(), io.real());
    return [b0, b1];
}

function quinticBlocks(io: OracleIO): Pair<HermiteQuinticSample> {
    const b0 = new HermiteQuinticSample(io.real(), io.real(), io.real());
    const b1 = new HermiteQuinticSample(io.real(), io.real(), io.real());
    return [b0, b1];
}

function bicubicBlocks(io: OracleIO): Pair<Pair<HermiteBicubicSample>> {
    const s = (): HermiteBicubicSample =>
        new HermiteBicubicSample(io.real(), io.real(), io.real(), io.real());
    const b00 = s(), b01 = s(), b10 = s(), b11 = s();
    return [[b00, b01], [b10, b11]];
}

function biquinticBlocks(io: OracleIO): Pair<Pair<HermiteBiquinticSample>> {
    const s = (): HermiteBiquinticSample => {
        const v = io.reals(9);
        return new HermiteBiquinticSample(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
    };
    const b00 = s(), b01 = s(), b10 = s(), b11 = s();
    return [[b00, b01], [b10, b11]];
}

function tricubicBlocks(io: OracleIO): Pair<Pair<Pair<HermiteTricubicSample>>> {
    const s = (): HermiteTricubicSample => {
        const v = io.reals(8);
        return new HermiteTricubicSample(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
    };
    const b000 = s(), b001 = s(), b010 = s(), b011 = s();
    const b100 = s(), b101 = s(), b110 = s(), b111 = s();
    return [[[b000, b001], [b010, b011]], [[b100, b101], [b110, b111]]];
}

function triquinticBlocks(io: OracleIO): Pair<Pair<Pair<HermiteTriquinticSample>>> {
    // The 27 fields in declaration order, which is also the constructor's
    // parameter order.
    const s = (): HermiteTriquinticSample => {
        const v = io.reals(27);
        return new HermiteTriquinticSample(...(v as [number, number, number, number,
            number, number, number, number, number, number, number, number, number,
            number, number, number, number, number, number, number, number, number,
            number, number, number, number, number]));
    };
    const b000 = s(), b001 = s(), b010 = s(), b011 = s();
    const b100 = s(), b101 = s(), b110 = s(), b111 = s();
    return [[[b000, b001], [b010, b011]], [[b100, b101], [b110, b111]]];
}

function outCoeff2(io: OracleIO, c: number[][]): void {
    for (const row of c) { io.outReals(row); }
}

function outCoeff3(io: OracleIO, c: number[][][]): void {
    for (const plane of c) { for (const row of plane) { io.outReals(row); } }
}

function outOrders2(io: OracleIO, h: { evaluate(xo: number, yo: number, x: number, y: number): number },
    maxOrder: number, x: number, y: number): void {
    for (let xo = 0; xo <= maxOrder; ++xo) {
        for (let yo = 0; yo <= maxOrder; ++yo) { io.outReal(h.evaluate(xo, yo, x, y)); }
    }
}

interface Evaluate3 {
    evaluate(xo: number, yo: number, zo: number, x: number, y: number, z: number): number;
}

function outDrawnOrders3(io: OracleIO, h: Evaluate3, count: number,
    x: number, y: number, z: number): void {
    for (let i = 0; i < count; ++i) {
        const xo = io.integer();
        const yo = io.integer();
        const zo = io.integer();
        io.outReal(h.evaluate(xo, yo, zo, x, y, z));
    }
}

describe('oracle: v28-interpolation', () => {
    const family = new OracleFamily('v28-interpolation');

    // ------------------------------------------------------ Hermite 1D

    family.case('HermiteCubic.evaluate', (io) => {
        const hermite = new HermiteCubic(cubicBlocks(io));
        for (let q = 0; q < 3; ++q) {
            const x = io.real();
            for (let order = 0; order <= 4; ++order) { io.outReal(hermite.evaluate(order, x)); }
        }
    }, exact);

    family.case('HermiteCubic.generate', (io) => {
        const hermite = new HermiteCubic();
        const x = io.real();
        io.outReal(hermite.evaluate(0, x));
        hermite.generate(cubicBlocks(io));
        for (let order = 0; order <= 3; ++order) { io.outReal(hermite.evaluate(order, x)); }
        hermite.generate(cubicBlocks(io));
        for (let order = 0; order <= 3; ++order) { io.outReal(hermite.evaluate(order, x)); }
    }, exact);

    family.case('HermiteCubic.P', (io) => {
        const t = io.real();
        for (let select = 0; select < 4; ++select) {
            for (let order = 0; order <= 4; ++order) { io.outReal(HermiteCubic.p(select, order, t)); }
        }
    }, exact);

    family.case('HermiteQuintic.evaluate', (io) => {
        const hermite = new HermiteQuintic(quinticBlocks(io));
        for (let q = 0; q < 3; ++q) {
            const x = io.real();
            for (let order = 0; order <= 6; ++order) { io.outReal(hermite.evaluate(order, x)); }
        }
    }, exact);

    family.case('HermiteQuintic.generate', (io) => {
        const hermite = new HermiteQuintic();
        const x = io.real();
        io.outReal(hermite.evaluate(0, x));
        hermite.generate(quinticBlocks(io));
        for (let order = 0; order <= 5; ++order) { io.outReal(hermite.evaluate(order, x)); }
        hermite.generate(quinticBlocks(io));
        for (let order = 0; order <= 5; ++order) { io.outReal(hermite.evaluate(order, x)); }
    }, exact);

    family.case('HermiteQuintic.P', (io) => {
        const t = io.real();
        for (let select = 0; select < 6; ++select) {
            for (let order = 0; order <= 6; ++order) { io.outReal(HermiteQuintic.p(select, order, t)); }
        }
    }, exact);

    // ------------------------------------------------------ Hermite 2D/3D

    family.case('HermiteBicubic.evaluate', (io) => {
        const hermite = new HermiteBicubic(bicubicBlocks(io));
        outCoeff2(io, hermite.c);
        for (let q = 0; q < 2; ++q) {
            const x = io.real();
            const y = io.real();
            outOrders2(io, hermite, 4, x, y);
        }
    }, exact);

    family.case('HermiteBicubic.manual', (io) => {
        const hermite = new HermiteBicubic();
        const x0 = io.real();
        const y0 = io.real();
        io.outReal(hermite.evaluate(0, 0, x0, y0));
        for (let i = 0; i < 4; ++i) {
            for (let j = 0; j < 4; ++j) { hermite.c[i][j] = io.real(); }
        }
        const x = io.real();
        const y = io.real();
        outOrders2(io, hermite, 4, x, y);
        hermite.generate(bicubicBlocks(io));
        outCoeff2(io, hermite.c);
    }, exact);

    family.case('HermiteBiquintic.evaluate', (io) => {
        const hermite = new HermiteBiquintic(biquinticBlocks(io));
        outCoeff2(io, hermite.c);
        for (let q = 0; q < 2; ++q) {
            const x = io.real();
            const y = io.real();
            outOrders2(io, hermite, 6, x, y);
        }
    }, exact);

    family.case('HermiteBiquintic.manual', (io) => {
        const hermite = new HermiteBiquintic();
        const x0 = io.real();
        const y0 = io.real();
        io.outReal(hermite.evaluate(0, 0, x0, y0));
        for (let i = 0; i < 6; ++i) {
            for (let j = 0; j < 6; ++j) { hermite.c[i][j] = io.real(); }
        }
        const x = io.real();
        const y = io.real();
        outOrders2(io, hermite, 6, x, y);
        hermite.generate(biquinticBlocks(io));
        outCoeff2(io, hermite.c);
    }, exact);

    family.case('HermiteTricubic.evaluate', (io) => {
        const hermite = new HermiteTricubic(tricubicBlocks(io));
        outCoeff3(io, hermite.c);
        const x = io.real();
        const y = io.real();
        const z = io.real();
        for (let xo = 0; xo <= 3; ++xo) {
            for (let yo = 0; yo <= 3; ++yo) {
                for (let zo = 0; zo <= 3; ++zo) { io.outReal(hermite.evaluate(xo, yo, zo, x, y, z)); }
            }
        }
        outDrawnOrders3(io, hermite, 4, x, y, z);
    }, exact);

    family.case('HermiteTricubic.manual', (io) => {
        const hermite = new HermiteTricubic();
        const x0 = io.real();
        const y0 = io.real();
        const z0 = io.real();
        io.outReal(hermite.evaluate(0, 0, 0, x0, y0, z0));
        for (let i = 0; i < 4; ++i) {
            for (let j = 0; j < 4; ++j) {
                for (let k = 0; k < 4; ++k) { hermite.c[i][j][k] = io.real(); }
            }
        }
        const x = io.real();
        const y = io.real();
        const z = io.real();
        outDrawnOrders3(io, hermite, 16, x, y, z);
    }, exact);

    family.case('HermiteTriquintic.evaluate', (io) => {
        const hermite = new HermiteTriquintic(triquinticBlocks(io));
        outCoeff3(io, hermite.c);
        const x = io.real();
        const y = io.real();
        const z = io.real();
        io.outReal(hermite.evaluate(0, 0, 0, x, y, z));
        outDrawnOrders3(io, hermite, 20, x, y, z);
    }, exact);

    family.case('HermiteTriquintic.manual', (io) => {
        const hermite = new HermiteTriquintic();
        const x0 = io.real();
        const y0 = io.real();
        const z0 = io.real();
        io.outReal(hermite.evaluate(0, 0, 0, x0, y0, z0));
        for (let i = 0; i < 6; ++i) {
            for (let j = 0; j < 6; ++j) {
                for (let k = 0; k < 6; ++k) { hermite.c[i][j][k] = io.real(); }
            }
        }
        const x = io.real();
        const y = io.real();
        const z = io.real();
        io.outReal(hermite.evaluate(0, 0, 0, x, y, z));
        outDrawnOrders3(io, hermite, 12, x, y, z);
    }, exact);

    // ------------------------------------------------------ IntpAkima1

    family.case('IntpAkima1.evaluate.uniform1', (io) => {
        const quantity = io.integer();
        const xMin = io.real();
        const xSpacing = io.real();
        const F = io.reals(quantity);
        const interp = new IntpAkimaUniform1(quantity, xMin, xSpacing, F);
        io.outInt(interp.getQuantity());
        io.outReal(interp.getXMin());
        io.outReal(interp.getXMax());
        for (let q = 0; q < 5; ++q) {
            const x = io.real();
            io.outReal(interp.evaluate(x));
            for (let order = -1; order <= 4; ++order) { io.outReal(interp.evaluate(order, x)); }
        }
    }, exact);

    family.case('IntpAkima1.construct.validation', (io) => {
        const quantity = io.integer();
        const xSpacing = io.real();
        const F = io.reals(4);
        const interp = new IntpAkimaUniform1(quantity, 0, xSpacing, F);
        io.outReal(interp.evaluate(0.25));
        io.outReal(interp.evaluate(1, 0.75));
    }, exact);

    // ------------------------------------------------------ IntpAkimaUniform2

    function akima2(io: OracleIO): IntpAkimaUniform2 {
        const nx = io.integer();
        const ny = io.integer();
        const xMin = io.real();
        const xSpacing = io.real();
        const yMin = io.real();
        const ySpacing = io.real();
        const F = io.reals(nx * ny);
        const interp = new IntpAkimaUniform2(nx, ny, xMin, xSpacing, yMin, ySpacing, F);
        io.outInt(interp.getQuantity());
        io.outReal(interp.getXMax());
        io.outReal(interp.getYMax());
        return interp;
    }

    function outAkima2(io: OracleIO, interp: IntpAkimaUniform2, queries: number): void {
        for (let q = 0; q < queries; ++q) {
            const x = io.real();
            const y = io.real();
            io.outReal(interp.evaluate(x, y));
            for (let xo = 0; xo <= 3; ++xo) {
                for (let yo = 0; yo <= 3; ++yo) { io.outReal(interp.evaluate(xo, yo, x, y)); }
            }
            io.outReal(interp.evaluate(4, 0, x, y));
            io.outReal(interp.evaluate(0, 4, x, y));
            io.outReal(interp.evaluate(-1, 1, x, y));
        }
    }

    family.case('IntpAkimaUniform2.evaluate', (io) => {
        outAkima2(io, akima2(io), 3);
    }, exact);

    family.case('IntpAkimaUniform2.evaluate.separable', (io) => {
        outAkima2(io, akima2(io), 3);
    }, exact);

    // The port negates upstream's max-boundary mixed-derivative stencils
    // (docs/UPSTREAM-FINDINGS.md, IntpAkimaUniform2.h).
    family.case('IntpAkimaUniform2.evaluate.maxBoundary', (io) => {
        outAkima2(io, akima2(io), 2);
    }, { exact: true, deviation: '#58: max-boundary FXY stencils have the wrong sign' });

    family.case('IntpAkimaUniform2.construct.validation', (io) => {
        const xBound = io.integer();
        const yBound = io.integer();
        const xSpacing = io.real();
        const ySpacing = io.real();
        const F = io.reals(xBound * yBound);
        const interp = new IntpAkimaUniform2(xBound, yBound, 0, xSpacing, 0, ySpacing, F);
        io.outReal(interp.evaluate(0.25, 0.125));
    }, exact);

    // ------------------------------------------------------ IntpAkimaUniform3

    function akima3(io: OracleIO): IntpAkimaUniform3 {
        const nx = io.integer();
        const ny = io.integer();
        const nz = io.integer();
        const xMin = io.real();
        const xSpacing = io.real();
        const yMin = io.real();
        const ySpacing = io.real();
        const zMin = io.real();
        const zSpacing = io.real();
        const F = io.reals(nx * ny * nz);
        const interp = new IntpAkimaUniform3(nx, ny, nz, xMin, xSpacing, yMin, ySpacing,
            zMin, zSpacing, F);
        io.outInt(interp.getQuantity());
        io.outReal(interp.getXMax());
        io.outReal(interp.getYMax());
        io.outReal(interp.getZMax());
        return interp;
    }

    function outAkima3(io: OracleIO, interp: IntpAkimaUniform3): void {
        for (let q = 0; q < 2; ++q) {
            const x = io.real();
            const y = io.real();
            const z = io.real();
            io.outReal(interp.evaluate(x, y, z));
            for (let xo = 0; xo <= 3; ++xo) {
                for (let yo = 0; yo <= 3; ++yo) {
                    for (let zo = 0; zo <= 3; ++zo) {
                        io.outReal(interp.evaluate(xo, yo, zo, x, y, z));
                    }
                }
            }
            io.outReal(interp.evaluate(4, 0, 0, x, y, z));
            io.outReal(interp.evaluate(0, 0, 4, x, y, z));
            io.outReal(interp.evaluate(1, -1, 0, x, y, z));
        }
    }

    family.case('IntpAkimaUniform3.evaluate', (io) => {
        outAkima3(io, akima3(io));
    }, exact);

    family.case('IntpAkimaUniform3.evaluate.separable', (io) => {
        outAkima3(io, akima3(io));
    }, exact);

    // The port negates upstream's max-boundary FXY, FXZ, FYZ and FXYZ
    // stencils (docs/UPSTREAM-FINDINGS.md, IntpAkimaUniform2.h).
    family.case('IntpAkimaUniform3.evaluate.maxBoundary', (io) => {
        outAkima3(io, akima3(io));
    }, { exact: true, deviation: '#58: max-boundary mixed-derivative stencils have the wrong sign' });

    family.case('IntpAkimaUniform3.construct.validation', (io) => {
        const xBound = io.integer();
        const yBound = io.integer();
        const zBound = io.integer();
        const xSpacing = io.real();
        const ySpacing = io.real();
        const zSpacing = io.real();
        const F = io.reals(xBound * yBound * zBound);
        const interp = new IntpAkimaUniform3(xBound, yBound, zBound, 0, xSpacing, 0, ySpacing,
            0, zSpacing, F);
        io.outReal(interp.evaluate(0.25, 0.125, 0.0625));
    }, exact);

    // ------------------------------------------------------ IntpBicubic2, IntpBilinear2

    family.case('IntpBicubic2.evaluate', (io) => {
        const xBound = io.integer();
        const yBound = io.integer();
        const xMin = io.real();
        const xSpacing = io.real();
        const yMin = io.real();
        const ySpacing = io.real();
        const catmullRom = io.boolean();
        const F = io.reals(xBound * yBound);
        const interp = new IntpBicubic2(xBound, yBound, xMin, xSpacing, yMin, ySpacing, F, catmullRom);
        io.outInt(interp.getQuantity());
        io.outReal(interp.getXMax());
        io.outReal(interp.getYMax());
        for (let q = 0; q < 3; ++q) {
            const x = io.real();
            const y = io.real();
            io.outReal(interp.evaluate(x, y));
            outOrders2(io, interp, 4, x, y);
            io.outReal(interp.evaluate(-1, 0, x, y));
        }
    }, exact);

    family.case('IntpBicubic2.construct.validation', (io) => {
        const xBound = io.integer();
        const yBound = io.integer();
        const xSpacing = io.real();
        const ySpacing = io.real();
        const catmullRom = io.boolean();
        const F = io.reals(xBound * yBound);
        const interp = new IntpBicubic2(xBound, yBound, 0, xSpacing, 0, ySpacing, F, catmullRom);
        io.outReal(interp.evaluate(0.25, 0.5));
    }, exact);

    family.case('IntpBilinear2.evaluate', (io) => {
        const xBound = io.integer();
        const yBound = io.integer();
        const xMin = io.real();
        const xSpacing = io.real();
        const yMin = io.real();
        const ySpacing = io.real();
        const F = io.reals(xBound * yBound);
        const interp = new IntpBilinear2(xBound, yBound, xMin, xSpacing, yMin, ySpacing, F);
        io.outInt(interp.getQuantity());
        io.outReal(interp.getXMax());
        io.outReal(interp.getYMax());
        for (let q = 0; q < 4; ++q) {
            const x = io.real();
            const y = io.real();
            io.outReal(interp.evaluate(x, y));
            outOrders2(io, interp, 2, x, y);
            io.outReal(interp.evaluate(0, -1, x, y));
        }
    }, exact);

    family.case('IntpBilinear2.construct.validation', (io) => {
        const xBound = io.integer();
        const yBound = io.integer();
        const xSpacing = io.real();
        const ySpacing = io.real();
        const F = io.reals(xBound * yBound);
        const interp = new IntpBilinear2(xBound, yBound, 0, xSpacing, 0, ySpacing, F);
        io.outReal(interp.evaluate(0.25, 0.5));
    }, exact);

    family.finish();
});
