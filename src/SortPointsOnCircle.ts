// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SortPointsOnCircle.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The sorting algorithms are described in
// https://www.geometrictools.com/Documentation/SortPointsOnCircle.pdf
// The input P[] are points and C is the center point about which the points
// are to be sorted. The reference ray (angle 0) is C+t*D for initial point
// C and nonzero direction D. The direction does not have to be unit length.
// If sortCCW is 'true', the angles counterclockwise from the reference ray
// are positive and in [0,pi]. The angles clockwise from the reference ray
// are negative and in (-pi,0]. If sortCCW is 'false', the angles clockwise
// from the reference ray are positive and in [0,pi]. The angles
// counterclockwise from the reference ray are negative and in (-pi,0]. The
// output 'indices[]' provides an indirect sorting. The sorted points are
// P[indices[0]], P[indices[1]], ..., P[indices[P.length-1]].
//
// Port notes: points are [x, y] tuples (the port of std::array<T, 2> with
// T = number). The output index array is returned rather than written to an
// out parameter. Array.prototype.sort with a three-way comparator derived
// from the upstream strict-weak-order predicates replaces std::sort; JS sort
// is stable, so points that compare equivalent keep their input order (an
// allowed ordering for std::sort as well).

interface SortObject {
    W: [number, number];
    index: number;
}

function makeSortObjects(
    P: readonly (readonly [number, number])[],
    C: readonly [number, number],
    D: readonly [number, number],
    sortCCW: boolean): SortObject[] {
    const Dperp: [number, number] = (sortCCW ? [-D[1], D[0]] : [D[1], -D[0]]);
    const objects: SortObject[] = new Array(P.length);
    for (let i = 0; i < P.length; ++i) {
        const V: [number, number] = [P[i][0] - C[0], P[i][1] - C[1]];
        objects[i] = {
            W: [D[0] * V[0] + D[1] * V[1], Dperp[0] * V[0] + Dperp[1] * V[1]],
            index: i
        };
    }
    return objects;
}

// The port of the upstream strict less-than predicates as three-way
// comparators (negative when object0 sorts before object1).
function lessThanByAngle(object0: SortObject, object1: SortObject): number {
    const x0 = object0.W[0], y0 = object0.W[1];
    const x1 = object1.W[0], y1 = object1.W[1];

    const angle0 = Math.atan2(y0, x0);
    const angle1 = Math.atan2(y1, x1);
    if (angle0 < angle1) {
        return -1;
    }
    if (angle0 > angle1) {
        return 1;
    }

    // Equal angles; sort by squared distance from the center.
    const s0 = (x0 - x1) * (x0 + x1);
    const s1 = (y1 - y0) * (y1 + y0);
    return s0 < s1 ? -1 : (s1 < s0 ? 1 : 0);
}

function lessThanByGeometry(object0: SortObject, object1: SortObject): number {
    // Derive the three-way result from the upstream strict predicate
    // less(a, b): negative if less(0, 1), positive if less(1, 0), else 0.
    if (strictLessThanByGeometry(object0, object1)) {
        return -1;
    }
    if (strictLessThanByGeometry(object1, object0)) {
        return 1;
    }
    return 0;
}

function strictLessThanByGeometry(object0: SortObject, object1: SortObject): boolean {
    const x0 = object0.W[0], y0 = object0.W[1];
    const x1 = object1.W[0], y1 = object1.W[1];
    const zero = 0;

    if (y0 < zero && y1 >= zero) {
        return true;
    }

    if (y1 < zero && y0 >= zero) {
        return false;
    }

    if (y0 > zero && y1 === zero) {
        return x1 < zero;
    }

    if (y1 > zero && y0 === zero) {
        return x0 > zero;
    }

    if (y0 === zero && y1 === zero) {
        return (x1 < zero && x1 < x0) || (x0 > zero && x1 > x0);
    }

    const c = x0 * y1 - x1 * y0;
    if (c > zero) {
        return true;
    }

    if (c < zero) {
        return false;
    }

    // c == 0, compare squared distances s0 < s1.
    return (x0 - x1) * (x0 + x1) < (y1 - y0) * (y1 + y0);
}

export class SortPointsOnCircle {
    // The sorting algorithm uses Math.atan2 and contains arithmetic
    // operations, all subject to floating-point rounding errors. An exact
    // rational type does not fix the problem because atan2 has mathematical
    // errors: the function cannot be implemented to produce exact angles
    // using only arithmetic operations.
    static byAngle(
        P: readonly (readonly [number, number])[],
        C: readonly [number, number],
        D: readonly [number, number],
        sortCCW: boolean): number[] {
        const objects = makeSortObjects(P, C, D, sortCCW);
        objects.sort(lessThanByAngle);
        return objects.map((object) => object.index);
    }

    // The sorting algorithm uses only arithmetic operations. With 'number'
    // inputs the correctness is not guaranteed to be theoretically correct
    // because of rounding errors. (Upstream: if T is an exact rational type,
    // the output is theoretically correct.)
    static byGeometry(
        P: readonly (readonly [number, number])[],
        C: readonly [number, number],
        D: readonly [number, number],
        sortCCW: boolean): number[] {
        const objects = makeSortObjects(P, C, D, sortCCW);
        objects.sort(lessThanByGeometry);
        return objects.map((object) => object.index);
    }
}
