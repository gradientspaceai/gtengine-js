// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BoxManager.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Sort-and-sweep management of a collection of axis-aligned boxes in 3D. The
// initial set of overlapping pairs is computed by a full sort of the box
// endpoints along each coordinate axis followed by a sweep of the x-endpoints.
// After the boxes are moved by small amounts, the endpoint arrays are nearly
// sorted, so an insertion sort restores the ordering quickly and the overlap
// set is updated incrementally as endpoints are transposed.
//
// Port notes:
// - Upstream stores 'std::vector<AlignedBox3<Real>>&', a reference to the
//   caller's array. The port stores the array itself, which has reference
//   semantics in TypeScript, so setBox() writes through to the caller's array
//   exactly as upstream does. The AlignedBox elements are cloned on the way in
//   and out (setBox/getBox) because C++ assigns box values.
// - The overlap container 'std::set<EdgeKey<false>>' becomes a Map from
//   EdgeKey.mapKey() to the EdgeKey, because JavaScript Set compares object
//   keys by identity. getOverlap() returns the keys sorted lexicographically
//   by (V[0], V[1]), which is the order std::set iterates them.
// - std::sort is not stable while Array.prototype.sort is. Endpoints that
//   compare equal (same value and same type) may therefore end up in a
//   different relative order than in the C++ build. The set of overlapping
//   pairs is unaffected.

import { AlignedBox } from './AlignedBox';
import { EdgeKey } from './EdgeKey';
import { IntrAlignedBox3AlignedBox3TI } from './IntrAlignedBox3AlignedBox3';

// An endpoint of a box interval along one coordinate axis.
interface Endpoint {
    // The endpoint value.
    value: number;

    // 0 if an interval minimum, 1 if an interval maximum.
    type: number;

    // The index of the interval (box) containing this endpoint.
    index: number;
}

// The port of the upstream Endpoint::operator<: order by value, breaking
// ties by type so that a 'begin' endpoint precedes an 'end' endpoint.
function endpointLess(e0: Endpoint, e1: Endpoint): boolean {
    if (e0.value < e1.value) {
        return true;
    }
    if (e0.value > e1.value) {
        return false;
    }
    return e0.type < e1.type;
}

function endpointCompare(e0: Endpoint, e1: Endpoint): number {
    if (endpointLess(e0, e1)) {
        return -1;
    }
    if (endpointLess(e1, e0)) {
        return 1;
    }
    return 0;
}

// The port of 'EdgeKey<false>(v0, v1)', which stores the pair with
// V[0] < V[1].
function makeOverlapKey(v0: number, v1: number): EdgeKey {
    return new EdgeKey(false, v0, v1);
}

export class BoxManager {
    private mBoxes: AlignedBox[];
    private mXEndpoints: Endpoint[];
    private mYEndpoints: Endpoint[];
    private mZEndpoints: Endpoint[];
    private mOverlap: Map<string, EdgeKey>;

    // The intervals are indexed 0 <= i < n. The endpoint array has 2*n
    // entries. The original 2*n interval values are ordered as
    //   b[0], e[0], b[1], e[1], ..., b[n-1], e[n-1]
    // When the endpoint array is sorted, the mapping between interval values
    // and endpoints is lost. In order to modify interval values that are
    // stored in the endpoint array, the mapping must be maintained. This is
    // done by the following lookup tables of 2*n entries. The value
    // mLookup[2*i] is the index of b[i] in the endpoint array. The value
    // mLookup[2*i+1] is the index of e[i] in the endpoint array.
    private mXLookup: number[];
    private mYLookup: number[];
    private mZLookup: number[];

    // Construction. Upstream has no default, copy or assignment
    // construction.
    constructor(boxes: AlignedBox[]) {
        this.mBoxes = boxes;
        this.mXEndpoints = [];
        this.mYEndpoints = [];
        this.mZEndpoints = [];
        this.mOverlap = new Map<string, EdgeKey>();
        this.mXLookup = [];
        this.mYLookup = [];
        this.mZLookup = [];
        this.initialize();
    }

    // This function is called by the constructor and does the sort-and-sweep
    // to initialize the update system. However, if you add or remove items
    // from the array of boxes after the constructor call, you will need to
    // call this function once before you start the multiple calls of the
    // update function.
    initialize(): void {
        // Get the box endpoints.
        const intrSize = this.mBoxes.length;
        const endpSize = 2 * intrSize;
        this.mXEndpoints = new Array<Endpoint>(endpSize);
        this.mYEndpoints = new Array<Endpoint>(endpSize);
        this.mZEndpoints = new Array<Endpoint>(endpSize);
        for (let i = 0, j = 0; i < intrSize; ++i) {
            const box = this.mBoxes[i];

            this.mXEndpoints[j] = { type: 0, value: box.min.values[0], index: i };
            this.mYEndpoints[j] = { type: 0, value: box.min.values[1], index: i };
            this.mZEndpoints[j] = { type: 0, value: box.min.values[2], index: i };
            ++j;

            this.mXEndpoints[j] = { type: 1, value: box.max.values[0], index: i };
            this.mYEndpoints[j] = { type: 1, value: box.max.values[1], index: i };
            this.mZEndpoints[j] = { type: 1, value: box.max.values[2], index: i };
            ++j;
        }

        // Sort the box endpoints.
        this.mXEndpoints.sort(endpointCompare);
        this.mYEndpoints.sort(endpointCompare);
        this.mZEndpoints.sort(endpointCompare);

        // Create the interval-to-endpoint lookup tables.
        this.mXLookup = new Array<number>(endpSize).fill(0);
        this.mYLookup = new Array<number>(endpSize).fill(0);
        this.mZLookup = new Array<number>(endpSize).fill(0);
        for (let j = 0; j < endpSize; ++j) {
            this.mXLookup[2 * this.mXEndpoints[j].index + this.mXEndpoints[j].type] = j;
            this.mYLookup[2 * this.mYEndpoints[j].index + this.mYEndpoints[j].type] = j;
            this.mZLookup[2 * this.mZEndpoints[j].index + this.mZEndpoints[j].type] = j;
        }

        // Active set of boxes (stored by index in array). Upstream uses a
        // std::set, which iterates in ascending order; the sweep result does
        // not depend on the iteration order, but the order is replicated so
        // the port matches upstream step for step.
        const active = new Set<number>();

        // Set of overlapping boxes (stored by pairs of indices in array).
        this.mOverlap.clear();

        // Sweep through the endpoints to determine overlapping x-intervals.
        for (let i = 0; i < endpSize; ++i) {
            const endpoint = this.mXEndpoints[i];
            const index = endpoint.index;
            if (endpoint.type === 0) {
                // An interval 'begin' value. In the 1D problem, the current
                // interval overlaps with all the active intervals. In 3D we
                // also need to check for y-overlap and z-overlap.
                const sorted = Array.from(active).sort((a, b) => a - b);
                for (const activeIndex of sorted) {
                    // Boxes activeIndex and index overlap in the x-dimension.
                    // Test for overlap in the y-dimension and z-dimension.
                    const b0 = this.mBoxes[activeIndex];
                    const b1 = this.mBoxes[index];
                    if (b0.max.values[1] >= b1.min.values[1]
                        && b0.min.values[1] <= b1.max.values[1]
                        && b0.max.values[2] >= b1.min.values[2]
                        && b0.min.values[2] <= b1.max.values[2]) {
                        const key = (activeIndex < index
                            ? makeOverlapKey(activeIndex, index)
                            : makeOverlapKey(index, activeIndex));
                        this.mOverlap.set(key.mapKey(), key);
                    }
                }
                active.add(index);
            } else {
                // An interval 'end' value.
                active.delete(index);
            }
        }
    }

    // After the system is initialized, you can move the boxes using this
    // function. It is not enough to modify the input array of boxes because
    // the endpoint values stored internally by this class must also change.
    setBox(i: number, box: AlignedBox): void {
        this.mBoxes[i] = box.clone();
        const twoI = 2 * i;
        this.mXEndpoints[this.mXLookup[twoI]].value = box.min.values[0];
        this.mXEndpoints[this.mXLookup[twoI + 1]].value = box.max.values[0];
        this.mYEndpoints[this.mYLookup[twoI]].value = box.min.values[1];
        this.mYEndpoints[this.mYLookup[twoI + 1]].value = box.max.values[1];
        this.mZEndpoints[this.mZLookup[twoI]].value = box.min.values[2];
        this.mZEndpoints[this.mZLookup[twoI + 1]].value = box.max.values[2];
    }

    // Retrieve the current box information. Upstream copies into an output
    // parameter; the port returns a copy.
    getBox(i: number): AlignedBox {
        return this.mBoxes[i].clone();
    }

    // When you are finished moving boxes, call this function to determine the
    // overlapping boxes. An incremental update is applied to determine the
    // new set of overlapping boxes.
    update(): void {
        this.insertionSort(this.mXEndpoints, this.mXLookup);
        this.insertionSort(this.mYEndpoints, this.mYLookup);
        this.insertionSort(this.mZEndpoints, this.mZLookup);
    }

    // If (i,j) is in the overlap set, then box i and box j are overlapping.
    // The indices are those for the input array. The keys (i,j) are stored so
    // that i < j. The array is sorted lexicographically by (V[0], V[1]), the
    // order in which the upstream std::set iterates.
    getOverlap(): EdgeKey[] {
        const keys = Array.from(this.mOverlap.values());
        keys.sort((e0, e1) => {
            if (e0.V[0] !== e1.V[0]) {
                return e0.V[0] - e1.V[0];
            }
            return e0.V[1] - e1.V[1];
        });
        return keys;
    }

    private insertionSort(endpoint: Endpoint[], lookup: number[]): void {
        // Apply an insertion sort. Under the assumption that the boxes have
        // not changed much since the last call, the endpoints are nearly
        // sorted. The insertion sort should be very fast in this case.

        const query = new IntrAlignedBox3AlignedBox3TI();
        const endpSize = endpoint.length;
        for (let j = 1; j < endpSize; ++j) {
            const key = endpoint[j];
            let i = j - 1;
            while (i >= 0 && endpointLess(key, endpoint[i])) {
                const e0 = endpoint[i];
                const e1 = endpoint[i + 1];

                // Update the overlap status.
                if (e0.type === 0) {
                    if (e1.type === 1) {
                        // The 'b' of interval e0.index was smaller than the
                        // 'e' of interval e1.index, and the intervals *might
                        // have been* overlapping. Now 'b' and 'e' are
                        // swapped, and the intervals cannot overlap. Remove
                        // the pair from the overlap set.
                        this.mOverlap.delete(makeOverlapKey(e0.index, e1.index).mapKey());
                    }
                } else {
                    if (e1.type === 0) {
                        // The 'b' of interval e1.index was larger than the
                        // 'e' of interval e0.index, and the intervals were
                        // not overlapping. Now 'b' and 'e' are swapped, and
                        // the intervals *might be* overlapping. Determine if
                        // they are overlapping and then insert.
                        if (query.test(this.mBoxes[e0.index], this.mBoxes[e1.index]).intersect) {
                            const okey = makeOverlapKey(e0.index, e1.index);
                            this.mOverlap.set(okey.mapKey(), okey);
                        }
                    }
                }

                // Reorder the items to maintain the sorted list.
                endpoint[i] = e1;
                endpoint[i + 1] = e0;
                lookup[2 * e1.index + e1.type] = i;
                lookup[2 * e0.index + e0.type] = i + 1;
                --i;
            }
            endpoint[i + 1] = key;
            lookup[2 * key.index + key.type] = i + 1;
        }
    }
}
