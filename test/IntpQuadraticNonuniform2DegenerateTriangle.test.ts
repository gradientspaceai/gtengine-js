// Regression test for IntpQuadraticNonuniform2 on a mesh with a collinear
// triangle (found by the C++ oracle, group 29 / v29-interpolation).
//
// Upstream ProcessTriangles passes a value-initialized Circle2 to Inscribe
// and ignores the returned bool. Inscribe (ContScribeCircle2.h) stores
//   circle.center = len21 * v0 + len20 * v1 + len10 * v2
// (the edge lengths normalized by the perimeter) *before* it tests the
// radius, so for a collinear triangle with a positive perimeter the center
// is that weighted point and only a triangle with three coincident vertices
// keeps (0,0). The port used inscribeCircle2, which returns null for every
// degenerate triangle, and left (0,0) for the collinear case too. The center
// of a triangle enters its neighbours' cross-edge intersections and Bezier
// coefficients, so the neighbours' values and gradients were wrong (scaled
// errors up to 0.67 in the oracle's committed records).
//
// The mesh is the oracle's: A = V0, B = V1, D = V2 on segment AB, C = V3,
// t0 = <A,D,C>, t1 = <D,B,C> and the collinear t2 = <A,B,D>, which shares
// <D,A> with t0 and <B,D> with t1. The expected values are the MSVC build's
// (oracle/golden/v29-interpolation.txt, IntpQuadraticNonuniform2.
// degenerateTriangle, record 7); PRE_FIX holds what the port produced before
// the fix and is asserted to differ.
import { describe, expect, it } from 'vitest';
import { BSRational } from '../src/BSRational.js';
import {
    IntpQuadraticNonuniform2, type IntpQuadraticNonuniform2TriangleMesh
} from '../src/IntpQuadraticNonuniform2.js';
import { Vector } from '../src/Vector.js';

function fromHex(hex: string): number {
    const view = new DataView(new ArrayBuffer(8));
    view.setBigUint64(0, BigInt('0x' + hex));
    return view.getFloat64(0);
}

// An explicit triangle list with exact barycentrics; containment is the
// first triangle whose barycentrics are all nonnegative.
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
        return [b0.toNumber(), b1.toNumber(), r(1).sub(b0).sub(b1).toNumber()];
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

describe('IntpQuadraticNonuniform2 with a collinear triangle', () => {
    const A = Vector.fromArray([-1, -3]);
    const B = Vector.fromArray([1, 0]);
    const D = Vector.fromArray([-1 + 0.75 * 2, -3 + 0.75 * 3]);
    const C = Vector.fromArray([fromHex('c0144103346cfc27'), fromHex('3ffe02b33677f5be')]);
    const F = ['bffeb060d0e61c3c', '3fd7612115f60dd0', '3ff2623192871c98', '3fc0f2def87c4ca0']
        .map(fromHex);
    const FX = [0.25, -0.5, 1, 0];
    const FY = [0.5, 0, -1, 0.75];
    const mesh = listMesh([A, B, D, C], [0, 2, 3, 2, 1, 3, 0, 1, 2],
        [2, 1, -1, 2, -1, 0, -1, 1, 0]);

    it('matches the MSVC build of upstream bit for bit', () => {
        // Record 7 uses the spatialDelta constructor with spatialDelta = 1.
        const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, F, 1);
        const queries = [['3fbdf7e65c981ec8', '3fa80accd9dfd6f8'], ['c006c103346cfc27', 'bfd1fa9993101484'],
            ['c002b0c26751bd1d', 'bfec7df9972607b2'], ['bfe84103346cfc27', 'bff43fa999310148']];
        const expected = [
            ['3fec2d84e58ba03f', 'bfcbb2abd57fc08f', 'bfb926f7af5a753b'],
            ['bfcbfb284841eff4', '3fdebca70e87822f', '3fe83eaf3558df08'],
            ['bfdd1b3f048f8520', '3fde213dc169d9a6', '3fea2eaecea798e9'],
            ['3fc2789cfa619cde', '3fe1ccd91cefb872', '3ff236e2f21f70ba']];
        queries.forEach((q, k) => {
            const r = intp.evaluate(Vector.fromArray(q.map(fromHex)));
            expect(r.valid).toBe(true);
            expect([r.F, r.FX, r.FY]).toEqual(expected[k].map(fromHex));
        });
        // The pre-fix port, which gave the collinear triangle the center (0,0).
        const PRE_FIX = [0.12473321430698642, 0.43932153506347726, 1.132150443416905];
        const r = intp.evaluate(Vector.fromArray(queries[3].map(fromHex)));
        expect([r.F, r.FX, r.FY]).not.toEqual(PRE_FIX);
    });

    it('gives the collinear triangle upstream\'s weighted center, (0,0) only when coincident', () => {
        const intp = IntpQuadraticNonuniform2.fromDerivatives(mesh, F, FX, FY);
        const tData = (intp as unknown as { mTData: Array<{ center: Vector }> }).mTData;
        // |AB| = sqrt(13), |AD| = 0.75 sqrt(13), |DB| = 0.25 sqrt(13); the
        // weights are the opposite edge lengths over the perimeter 2 sqrt(13).
        const center = tData[2].center.values;
        expect(center[0]).toBeCloseTo(0.125 * -1 + 0.375 * 1 + 0.5 * D.values[0], 14);
        expect(center[1]).toBeCloseTo(0.125 * -3 + 0.375 * 0 + 0.5 * D.values[1], 14);

        const P = Vector.fromArray([2, 2]);
        const coincident = listMesh([A, B, C, P, P], [0, 1, 2, 3, 4, 3], [-1, -1, -1, -1, -1, -1]);
        const intp2 = IntpQuadraticNonuniform2.fromDerivatives(coincident, [1, 2, 3, 4, 5],
            [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]);
        const t2 = (intp2 as unknown as { mTData: Array<{ center: Vector }> }).mTData;
        expect(t2[1].center.values).toEqual([0, 0]);
    });
});
