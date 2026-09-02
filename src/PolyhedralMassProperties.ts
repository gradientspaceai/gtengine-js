// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PolyhedralMassProperties.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Fast computation of centroid and inertia tensor are described in
// https://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
//
// Port notes: upstream 'Vector3<Real>' is the runtime-sized Vector with size
// 3 and 'Matrix3x3<Real>' is the runtime-sized Matrix with size 3-by-3. The
// output reference parameters mass, center and inertia become the fields of
// the returned ComputeMassPropertiesResult.

import { Matrix } from './Matrix';
import { Vector, div } from './Vector';
import { cross } from './Vector3';

export interface ComputeMassPropertiesResult {
    // The volume of the polyhedron (the mass for a body of density 1).
    mass: number;

    // The center of mass.
    center: Vector;

    // The inertia tensor, relative to the world origin when bodyCoords is
    // false and relative to the center of mass when bodyCoords is true.
    inertia: Matrix;
}

// The input triangle mesh must represent a polyhedron. The triangles are
// represented as triples of indices <V0,V1,V2> into the vertex array. The
// index array has numTriangles such triples. The Boolean value 'bodyCoords'
// is 'true' if you want the inertia tensor to be relative to body coordinates
// but 'false' if you want it to be relative to world coordinates.
//
// The code assumes the rigid body has a constant density of 1. If your
// application assigns a constant density of 'd', then you must multiply the
// output 'mass' by 'd' and the output 'inertia' by 'd'.
export function computeMassProperties(vertices: readonly Vector[],
    numTriangles: number, indices: readonly number[],
    bodyCoords: boolean): ComputeMassPropertiesResult {
    const oneDiv6 = 1 / 6;
    const oneDiv24 = 1 / 24;
    const oneDiv60 = 1 / 60;
    const oneDiv120 = 1 / 120;

    // order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
    const integral = new Array<number>(10).fill(0);

    let index = 0;
    for (let i = 0; i < numTriangles; ++i) {
        // Get vertices of triangle i.
        const v0 = vertices[indices[index++]].values;
        const v1 = vertices[indices[index++]].values;
        const v2 = vertices[indices[index++]].values;

        // Get cross product of edges and normal vector.
        const V1mV0 = Vector.fromArray([v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]]);
        const V2mV0 = Vector.fromArray([v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]]);
        const N = cross(V1mV0, V2mV0).values;

        // Compute integral terms.
        let tmp0: number;
        let tmp1: number;
        let tmp2: number;

        tmp0 = v0[0] + v1[0];
        const f1x = tmp0 + v2[0];
        tmp1 = v0[0] * v0[0];
        tmp2 = tmp1 + v1[0] * tmp0;
        const f2x = tmp2 + v2[0] * f1x;
        const f3x = v0[0] * tmp1 + v1[0] * tmp2 + v2[0] * f2x;
        const g0x = f2x + v0[0] * (f1x + v0[0]);
        const g1x = f2x + v1[0] * (f1x + v1[0]);
        const g2x = f2x + v2[0] * (f1x + v2[0]);

        tmp0 = v0[1] + v1[1];
        const f1y = tmp0 + v2[1];
        tmp1 = v0[1] * v0[1];
        tmp2 = tmp1 + v1[1] * tmp0;
        const f2y = tmp2 + v2[1] * f1y;
        const f3y = v0[1] * tmp1 + v1[1] * tmp2 + v2[1] * f2y;
        const g0y = f2y + v0[1] * (f1y + v0[1]);
        const g1y = f2y + v1[1] * (f1y + v1[1]);
        const g2y = f2y + v2[1] * (f1y + v2[1]);

        tmp0 = v0[2] + v1[2];
        const f1z = tmp0 + v2[2];
        tmp1 = v0[2] * v0[2];
        tmp2 = tmp1 + v1[2] * tmp0;
        const f2z = tmp2 + v2[2] * f1z;
        const f3z = v0[2] * tmp1 + v1[2] * tmp2 + v2[2] * f2z;
        const g0z = f2z + v0[2] * (f1z + v0[2]);
        const g1z = f2z + v1[2] * (f1z + v1[2]);
        const g2z = f2z + v2[2] * (f1z + v2[2]);

        // Update integrals.
        integral[0] += N[0] * f1x;
        integral[1] += N[0] * f2x;
        integral[2] += N[1] * f2y;
        integral[3] += N[2] * f2z;
        integral[4] += N[0] * f3x;
        integral[5] += N[1] * f3y;
        integral[6] += N[2] * f3z;
        integral[7] += N[0] * (v0[1] * g0x + v1[1] * g1x + v2[1] * g2x);
        integral[8] += N[1] * (v0[2] * g0y + v1[2] * g1y + v2[2] * g2y);
        integral[9] += N[2] * (v0[0] * g0z + v1[0] * g1z + v2[0] * g2z);
    }

    integral[0] *= oneDiv6;
    integral[1] *= oneDiv24;
    integral[2] *= oneDiv24;
    integral[3] *= oneDiv24;
    integral[4] *= oneDiv60;
    integral[5] *= oneDiv60;
    integral[6] *= oneDiv60;
    integral[7] *= oneDiv120;
    integral[8] *= oneDiv120;
    integral[9] *= oneDiv120;

    // mass
    const mass = integral[0];

    // center of mass
    const center = div(Vector.fromArray([integral[1], integral[2], integral[3]]), mass);

    // inertia relative to world origin
    const inertia = new Matrix(3, 3);
    inertia.set(0, 0, integral[5] + integral[6]);
    inertia.set(0, 1, -integral[7]);
    inertia.set(0, 2, -integral[9]);
    inertia.set(1, 0, inertia.get(0, 1));
    inertia.set(1, 1, integral[4] + integral[6]);
    inertia.set(1, 2, -integral[8]);
    inertia.set(2, 0, inertia.get(0, 2));
    inertia.set(2, 1, inertia.get(1, 2));
    inertia.set(2, 2, integral[4] + integral[5]);

    // inertia relative to center of mass
    if (bodyCoords) {
        const c = center.values;
        inertia.set(0, 0, inertia.get(0, 0) - mass * (c[1] * c[1] + c[2] * c[2]));
        inertia.set(0, 1, inertia.get(0, 1) + mass * c[0] * c[1]);
        inertia.set(0, 2, inertia.get(0, 2) + mass * c[2] * c[0]);
        inertia.set(1, 0, inertia.get(0, 1));
        inertia.set(1, 1, inertia.get(1, 1) - mass * (c[2] * c[2] + c[0] * c[0]));
        inertia.set(1, 2, inertia.get(1, 2) + mass * c[1] * c[2]);
        inertia.set(2, 0, inertia.get(0, 2));
        inertia.set(2, 1, inertia.get(1, 2));
        inertia.set(2, 2, inertia.get(2, 2) - mass * (c[0] * c[0] + c[1] * c[1]));
    }

    return { mass, center, inertia };
}
