// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// ContScribeCircle3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Circumscribed and inscribed circles (for triangles) and spheres (for
// tetrahedra) in 3D.
//
// Port notes: upstream overloads Circumscribe/Inscribe on the container type
// written through the output reference. Following the Cont* naming precedent,
// the port suffixes each function with that container type
// (circumscribeCircle3, circumscribeSphere3, inscribeCircle3,
// inscribeSphere3). Upstream returns 'true' when the circle/sphere was
// constructed and writes it through an output parameter; the port returns the
// circle/sphere, or null when the input points are linearly dependent.

import { Circle3 } from './Circle3';
import { Hypersphere, type Sphere3 } from './Hypersphere';
import { LinearSystem } from './LinearSystem';
import { logAssert } from './Logger';
import { Matrix } from './Matrix';
import { Vector, add, dot, length, mul, normalize, sub } from './Vector';
import { cross, unitCross } from './Vector3';

function assert3D(name: string, ...v: Vector[]): void {
    for (const u of v) {
        logAssert(u.size === 3, name + ': points must be 3D.');
    }
}

// Circle circumscribing a triangle in 3D. Returns null when the input points
// are linearly dependent (a degenerate triangle).
export function circumscribeCircle3(v0: Vector, v1: Vector,
    v2: Vector): Circle3 | null {
    assert3D('circumscribeCircle3', v0, v1, v2);

    const E02 = sub(v0, v2);
    const E12 = sub(v1, v2);
    const e02e02 = dot(E02, E02);
    const e02e12 = dot(E02, E12);
    const e12e12 = dot(E12, E12);
    const det = e02e02 * e12e12 - e02e12 * e02e12;
    if (det !== 0) {
        const halfInvDet = 0.5 / det;
        const u0 = halfInvDet * e12e12 * (e02e02 - e02e12);
        const u1 = halfInvDet * e02e02 * (e12e12 - e02e12);
        const tmp = add(mul(u0, E02), mul(u1, E12));
        const circle = new Circle3();
        circle.center = add(v2, tmp);
        circle.normal = unitCross(E02, E12);
        circle.radius = length(tmp);
        return circle;
    }
    return null;
}

// Sphere circumscribing a tetrahedron. Returns null when the input points are
// linearly dependent (a degenerate tetrahedron).
export function circumscribeSphere3(v0: Vector, v1: Vector, v2: Vector,
    v3: Vector): Sphere3 | null {
    assert3D('circumscribeSphere3', v0, v1, v2, v3);

    const E10 = sub(v1, v0);
    const E20 = sub(v2, v0);
    const E30 = sub(v3, v0);

    const A = Matrix.zero(3, 3);
    A.setRow(0, E10);
    A.setRow(1, E20);
    A.setRow(2, E30);

    const B = Vector.fromArray([
        0.5 * dot(E10, E10),
        0.5 * dot(E20, E20),
        0.5 * dot(E30, E30)
    ]);

    const { X, invertible } = LinearSystem.solve3x3(A, B);
    if (invertible) {
        const sphere = new Hypersphere(3);
        sphere.center = add(v0, X);
        sphere.radius = length(X);
        return sphere;
    }
    return null;
}

// Circle inscribing a triangle in 3D. Returns null when the triangle is
// degenerate.
export function inscribeCircle3(v0: Vector, v1: Vector,
    v2: Vector): Circle3 | null {
    assert3D('inscribeCircle3', v0, v1, v2);

    const circle = new Circle3();

    // Edges.
    const E0 = sub(v1, v0);
    const E1 = sub(v2, v1);
    const E2 = sub(v0, v2);

    // Plane normal.
    circle.normal = cross(E1, E0);

    // Edge normals within the plane.
    const N0 = unitCross(circle.normal, E0);
    const N1 = unitCross(circle.normal, E1);
    const N2 = unitCross(circle.normal, E2);

    const a0 = dot(N1, E0);
    if (a0 === 0) {
        return null;
    }

    const a1 = dot(N2, E1);
    if (a1 === 0) {
        return null;
    }

    const a2 = dot(N0, E2);
    if (a2 === 0) {
        return null;
    }

    const invA0 = 1 / a0;
    const invA1 = 1 / a1;
    const invA2 = 1 / a2;

    circle.radius = 1 / (invA0 + invA1 + invA2);
    circle.center = mul(circle.radius,
        add(add(mul(invA0, v0), mul(invA1, v1)), mul(invA2, v2)));
    normalize(circle.normal);
    return circle;
}

// Sphere inscribing a tetrahedron. Returns null when a face normal is
// degenerate or when the linear system is not invertible.
export function inscribeSphere3(v0: Vector, v1: Vector, v2: Vector,
    v3: Vector): Sphere3 | null {
    assert3D('inscribeSphere3', v0, v1, v2, v3);

    // Edges.
    const E10 = sub(v1, v0);
    const E20 = sub(v2, v0);
    const E30 = sub(v3, v0);
    const E21 = sub(v2, v1);
    const E31 = sub(v3, v1);

    // Normals.
    const N0 = cross(E31, E21);
    const N1 = cross(E20, E30);
    const N2 = cross(E30, E10);
    const N3 = cross(E10, E20);

    // Normalize the normals.
    if (normalize(N0) === 0) {
        return null;
    }
    if (normalize(N1) === 0) {
        return null;
    }
    if (normalize(N2) === 0) {
        return null;
    }
    if (normalize(N3) === 0) {
        return null;
    }

    // The insphere center C = v3 + S satisfies equal signed distances to the
    // four faces. Faces 1, 2 and 3 contain v0 and face 0 contains v1, and
    // N1, N2 are perpendicular to E30 while N0 is perpendicular to E31, so
    // the distance differences reduce to the system below.
    const A = Matrix.zero(3, 3);
    A.setRow(0, sub(N1, N0));
    A.setRow(1, sub(N2, N0));
    A.setRow(2, sub(N3, N0));
    const B = Vector.fromArray([0, 0, -dot(N3, E30)]);
    const { X, invertible } = LinearSystem.solve3x3(A, B);
    if (invertible) {
        const sphere = new Hypersphere(3);
        sphere.center = add(v3, X);
        sphere.radius = Math.abs(dot(N0, X));
        return sphere;
    }
    return null;
}
