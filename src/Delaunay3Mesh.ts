// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Delaunay3Mesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A lightweight mesh interface to a Delaunay3 tetrahedralization. The class
// stores a reference to the tetrahedralization; it must outlive the mesh
// object and must not be recomputed while the mesh object is in use.
//
// Port notes:
// * Upstream Delaunay3Mesh.h declares a variadic class template with two
//   specializations: the deprecated
//   Delaunay3Mesh<InputType, ComputeType, RationalType> and the replacement
//   Delaunay3Mesh<T>. Only the replacement is ported, matching the Delaunay3
//   port.
// * The upstream overloads GetVertices(), GetVertices(t, out), and so on,
//   cannot be distinguished by name in TypeScript, so the per-tetrahedron
//   forms are named getTetrahedronVertices(t), getTetrahedronIndices(t) and
//   getTetrahedronAdjacencies(t), matching the naming already used by
//   Delaunay3. The bool/out-parameter forms return the value or null.
// * The upstream barycentric computation uses BSRational<UIntegerAP32> for
//   exactness; the port does the same with BSRational, then converts the
//   result to double.
// * The deprecated upstream specialization has a GetInvalidIndex() only in
//   the replacement specialization; the port provides it, returning -1.

import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';
import { BSRational } from './BSRational.js';
import { Delaunay3, Delaunay3SearchInfo } from './Delaunay3.js';

// The port of Vector3<Rational>.
type RationalPoint3 = [BSRational, BSRational, BSRational];

function rationalSub3(v0: RationalPoint3, v1: RationalPoint3): RationalPoint3 {
    return [v0[0].sub(v1[0]), v0[1].sub(v1[1]), v0[2].sub(v1[2])];
}

// The port of DotCross(v0, v1, v2) = Dot(v0, Cross(v1, v2)) for rational
// vectors.
function rationalDotCross3(v0: RationalPoint3, v1: RationalPoint3,
    v2: RationalPoint3): BSRational {
    const cross: RationalPoint3 = [
        v1[1].mul(v2[2]).sub(v1[2].mul(v2[1])),
        v1[2].mul(v2[0]).sub(v1[0].mul(v2[2])),
        v1[0].mul(v2[1]).sub(v1[1].mul(v2[0]))
    ];
    return v0[0].mul(cross[0]).add(v0[1].mul(cross[1])).add(v0[2].mul(cross[2]));
}

export class Delaunay3Mesh {
    private mDelaunay: Delaunay3;

    // Construction. The Delaunay tetrahedralization must have dimension 3.
    constructor(delaunay: Delaunay3) {
        logAssert(delaunay.getDimension() === 3, 'Invalid Delaunay dimension.');
        this.mDelaunay = delaunay;
    }

    // Mesh information.
    getNumVertices(): number {
        return this.mDelaunay.getNumVertices();
    }

    getNumTetrahedra(): number {
        return this.mDelaunay.getNumTetrahedra();
    }

    getVertices(): readonly Vector[] {
        return this.mDelaunay.getVertices();
    }

    getIndices(): readonly number[] {
        return this.mDelaunay.getIndices();
    }

    getAdjacencies(): readonly number[] {
        return this.mDelaunay.getAdjacencies();
    }

    // The port of Delaunay3::negOne, the "no such tetrahedron" sentinel.
    getInvalidIndex(): number {
        return -1;
    }

    // Containment queries. The returned value is the index of the tetrahedron
    // containing P, or getInvalidIndex() when there is no such tetrahedron.
    getContainingTetrahedron(p: Vector): number {
        const info = new Delaunay3SearchInfo();
        return this.mDelaunay.getContainingTetrahedron(p, info);
    }

    // The vertices of tetrahedron t, or null when t is not a valid
    // tetrahedron index.
    getTetrahedronVertices(t: number): [Vector, Vector, Vector, Vector] | null {
        const indices = this.mDelaunay.getTetrahedronIndices(t);
        if (indices === null) {
            return null;
        }
        const vertices = this.mDelaunay.getVertices();
        return [vertices[indices[0]].clone(), vertices[indices[1]].clone(),
            vertices[indices[2]].clone(), vertices[indices[3]].clone()];
    }

    // The vertex indices of tetrahedron t, or null when t is not a valid
    // tetrahedron index.
    getTetrahedronIndices(t: number): [number, number, number, number] | null {
        return this.mDelaunay.getTetrahedronIndices(t);
    }

    // The adjacent tetrahedra of tetrahedron t, or null when t is not a valid
    // tetrahedron index. An adjacency of -1 indicates a boundary face.
    getTetrahedronAdjacencies(t: number): [number, number, number, number] | null {
        return this.mDelaunay.getTetrahedronAdjacencies(t);
    }

    // The barycentric coordinates of P with respect to tetrahedron t. The
    // result is null when t is not a valid tetrahedron index or when the
    // tetrahedron is degenerate. The computation is exact (rational
    // arithmetic) and the result is rounded to double at the end.
    getBarycentrics(t: number, p: Vector): [number, number, number, number] | null {
        logAssert(p.size === 3, 'Delaunay3Mesh requires 3D points.');
        const indices = this.mDelaunay.getTetrahedronIndices(t);
        if (indices === null) {
            return null;
        }

        const vertices = this.mDelaunay.getVertices();
        const rtV: RationalPoint3[] = [];
        for (let i = 0; i < 4; ++i) {
            const v = vertices[indices[i]].values;
            rtV.push([BSRational.fromNumber(v[0]), BSRational.fromNumber(v[1]),
                BSRational.fromNumber(v[2])]);
        }
        const rtP: RationalPoint3 = [BSRational.fromNumber(p.values[0]),
            BSRational.fromNumber(p.values[1]), BSRational.fromNumber(p.values[2])];

        // The port of ComputeBarycentrics(P, V0, V1, V2, V3, bary) with the
        // default epsilon of zero.
        const diff0 = rationalSub3(rtV[0], rtV[3]);
        const diff1 = rationalSub3(rtV[1], rtV[3]);
        const diff2 = rationalSub3(rtV[2], rtV[3]);
        const diff3 = rationalSub3(rtP, rtV[3]);
        const det = rationalDotCross3(diff0, diff1, diff2);
        if (det.getSign() === 0) {
            return null;
        }

        const b0 = rationalDotCross3(diff3, diff1, diff2).div(det);
        const b1 = rationalDotCross3(diff3, diff2, diff0).div(det);
        const b2 = rationalDotCross3(diff3, diff0, diff1).div(det);
        const b3 = BSRational.fromNumber(1).sub(b0).sub(b1).sub(b2);
        return [b0.toNumber(), b1.toNumber(), b2.toNumber(), b3.toNumber()];
    }
}
