// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContTetrahedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment query for a point in a tetrahedron.
//
// Port notes: the Cont* naming precedent suffixes the query with the
// container type, so InContainer becomes inContainerTetrahedron3.

import { logAssert } from './Logger';
import { Tetrahedron3 } from './Tetrahedron3';
import { Vector, sub } from './Vector';
import { dotCross } from './Vector3';

// Test for containment of a point by a tetrahedron. The test assumes the
// vertices are ordered so that the triple scalar products below are
// nonpositive for interior points; that is, the tetrahedron has the ordering
// used by Tetrahedron3.
export function inContainerTetrahedron3(point: Vector,
    tetra: Tetrahedron3): boolean {
    logAssert(point.size === 3, 'inContainerTetrahedron3: point must be 3D.');

    const zero = 0;

    // A loop over the faces is not used in order to avoid redundant
    // computations of edge directions. The difference vector is the same for
    // the first 3 faces but differs for the last face.

    // face <0,2,1>
    const edge20 = sub(tetra.v[2], tetra.v[0]);
    const edge10 = sub(tetra.v[1], tetra.v[0]);
    const diffP0 = sub(point, tetra.v[0]);
    if (dotCross(edge20, edge10, diffP0) > zero) {
        return false;
    }

    // face <0,1,3>
    const edge30 = sub(tetra.v[3], tetra.v[0]);
    if (dotCross(edge10, edge30, diffP0) > zero) {
        return false;
    }

    // face <0,3,2>
    if (dotCross(edge30, edge20, diffP0) > zero) {
        return false;
    }

    // face <1,2,3>
    const edge21 = sub(tetra.v[2], tetra.v[1]);
    const edge31 = sub(tetra.v[3], tetra.v[1]);
    const diffP1 = sub(point, tetra.v[1]);
    if (dotCross(edge21, edge31, diffP1) > zero) {
        return false;
    }

    return true;
}
