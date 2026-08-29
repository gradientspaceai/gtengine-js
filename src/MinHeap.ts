// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinHeap.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A min-heap is a binary tree whose nodes have weights and with the
// constraint that the weight of a parent node is less than or equal to the
// weights of its children. This data structure may be used as a priority
// queue. Unlike a plain priority queue, this min-heap supports updating the
// weight of an element that is already in the heap without removing and
// reinserting it: Insert returns a Record whose weight can later be changed
// through Update. For example, when decimating a polyline, removing the
// minimum-weight vertex changes the weights of its two neighbors; those
// neighbors are already in the min-heap and their records can be updated in
// place. Usage for (open) polyline decimation:
//
//    interface Vertex { previous: number, current: number, next: number }
//    const numVertices = <number of polyline vertices>;
//    const positions: Vector[] = <the polyline vertices>;
//    const minHeap = new MinHeap<Vertex>(numVertices);
//    const records: MinHeapRecord<Vertex, number>[] = [];
//    for (let i = 0; i < numVertices; ++i) {
//        const vertex = {
//            previous: (i + numVertices - 1) % numVertices,
//            current: i,
//            next: (i + 1) % numVertices
//        };
//        records[i] = minHeap.insert(vertex, weight(positions, vertex))!;
//    }
//
//    while (minHeap.getNumElements() >= 2) {
//        const { key: vertex } = minHeap.remove()!;
//        // <consume 'vertex' according to your application's needs>
//
//        // Remove 'vertex' from the doubly linked list.
//        const vp = records[vertex.previous].key;
//        const vc = records[vertex.current].key;
//        const vn = records[vertex.next].key;
//        vp.next = vc.next;
//        vn.previous = vc.previous;
//
//        // Update the neighbors' weights in the min-heap.
//        minHeap.update(records[vertex.previous], weight(positions, vp));
//        minHeap.update(records[vertex.next], weight(positions, vn));
//    }
//
// Port notes: upstream ValueType requires operators "<" and "<=". The port
// defaults ValueType to number and compares with "<"; a custom lessThan
// comparator may be passed for other value types, with "a <= b" derived as
// "!(b < a)" (equivalent for any strict total order). Upstream out-parameter
// functions return object literals or null instead: GetMinimum/Remove return
// { key, value } or null, and Insert returns the Record or null when the
// heap is full. Records are preallocated by reset(); their key and value are
// undefined until the record is used by insert, matching the upstream
// warning that record values are uninitialized.

// The port of MinHeap<KeyType, ValueType>::Record. The 'key' member may be
// read and written by users (see the polyline example). The 'value' member
// must be modified only through MinHeap.update. The 'index' member is
// internal bookkeeping and must not be modified.
export interface MinHeapRecord<KeyType, ValueType> {
    key: KeyType;
    value: ValueType;
    index: number;
}

export class MinHeap<KeyType, ValueType = number> {
    // A 2-level storage system is used. The records are unique to each
    // inserted value in order to support the update() capability of the
    // min-heap. mPointers is the permutation of references to the records
    // that stores the heap ordering.
    private mNumElements: number;
    private mRecords: MinHeapRecord<KeyType, ValueType>[];
    private mPointers: MinHeapRecord<KeyType, ValueType>[];
    private readonly lt: (a: ValueType, b: ValueType) => boolean;

    constructor(maxElements: number = 0,
        lessThan?: (a: ValueType, b: ValueType) => boolean) {
        this.lt = lessThan ?? ((a: ValueType, b: ValueType) =>
            (a as unknown as number) < (b as unknown as number));
        this.mNumElements = 0;
        this.mRecords = [];
        this.mPointers = [];
        this.reset(maxElements);
    }

    // The port of "a <= b", derived from the strict order.
    private le(a: ValueType, b: ValueType): boolean {
        return !this.lt(b, a);
    }

    // Clear the min-heap so that it has the specified max elements, the
    // number of elements is zero, and the pointers are set to the natural
    // ordering of the records.
    reset(maxElements: number): void {
        this.mNumElements = 0;
        this.mRecords = [];
        this.mPointers = [];
        for (let i = 0; i < maxElements; ++i) {
            const record: MinHeapRecord<KeyType, ValueType> = {
                key: undefined as unknown as KeyType,
                value: undefined as unknown as ValueType,
                index: i
            };
            this.mRecords.push(record);
            this.mPointers.push(record);
        }
    }

    // Get the remaining number of elements in the min-heap. This number is
    // in the range {0..maxElements}.
    getNumElements(): number {
        return this.mNumElements;
    }

    // Get the root of the min-heap without removing it. The return value is
    // null whenever the min-heap is empty.
    getMinimum(): { key: KeyType, value: ValueType } | null {
        if (this.mNumElements > 0) {
            return { key: this.mPointers[0].key, value: this.mPointers[0].value };
        }
        return null;
    }

    // Insert into the min-heap the 'value' that corresponds to the 'key'.
    // The return value is the heap record that stores the value, and the
    // record identity is constant for the life of the min-heap. If you must
    // update a member of the min-heap, say, as illustrated in the polyline
    // decimation example, pass the record to update():
    //    const valueRecord = minHeap.insert(key, value);
    //    // <do whatever>
    //    minHeap.update(valueRecord, newValue);
    // The return value is null when the heap is full.
    insert(key: KeyType, value: ValueType): MinHeapRecord<KeyType, ValueType> | null {
        // Return immediately when the heap is full.
        if (this.mNumElements === this.mRecords.length) {
            return null;
        }

        // Store the input information in the last heap record, which is the
        // last leaf in the tree.
        let child = this.mNumElements++;
        const record = this.mPointers[child];
        record.key = key;
        record.value = value;

        // Propagate the information toward the root of the tree until it
        // reaches its correct position, thus restoring the tree to a valid
        // heap.
        while (child > 0) {
            const parent = (child - 1) >> 1;
            if (this.le(this.mPointers[parent].value, value)) {
                // The parent has a value smaller than or equal to the
                // child's value, so we now have a valid heap.
                break;
            }

            // The parent has a larger value than the child's value. Swap the
            // parent and child:

            // Move the parent into the child's slot.
            this.mPointers[child] = this.mPointers[parent];
            this.mPointers[child].index = child;

            // Move the child into the parent's slot.
            this.mPointers[parent] = record;
            this.mPointers[parent].index = parent;

            child = parent;
        }

        return this.mPointers[child];
    }

    // Remove the root of the heap and return its key and value. The root
    // contains the minimum value of all heap elements. The return value is
    // null whenever the min-heap is empty.
    remove(): { key: KeyType, value: ValueType } | null {
        // Return immediately when the heap is empty.
        if (this.mNumElements === 0) {
            return null;
        }

        // Get the information from the root of the heap.
        const root = this.mPointers[0];
        const result = { key: root.key, value: root.value };

        // Restore the tree to a heap. Abstractly, record is the new root of
        // the heap. It is moved down the tree via parent-child swaps until
        // it is in a location that restores the tree to a heap.
        const last = --this.mNumElements;
        const record = this.mPointers[last];
        let parent = 0;
        let child = 1;
        while (child <= last) {
            if (child < last) {
                // Select the child with smallest value to be the one that is
                // swapped with the parent, if necessary.
                const childP1 = child + 1;
                if (this.lt(this.mPointers[childP1].value, this.mPointers[child].value)) {
                    child = childP1;
                }
            }

            if (this.le(record.value, this.mPointers[child].value)) {
                // The tree is now a heap.
                break;
            }

            // Move the child into the parent's slot.
            this.mPointers[parent] = this.mPointers[child];
            this.mPointers[parent].index = parent;

            parent = child;
            child = 2 * child + 1;
        }

        // The previous 'last' record was moved to the root and propagated
        // down the tree to its final resting place, restoring the tree to a
        // heap. The slot mPointers[parent] is that resting place.
        this.mPointers[parent] = record;
        this.mPointers[parent].index = parent;

        // The old root record must not be lost. Attach it to the slot that
        // contained the old last record.
        this.mPointers[last] = root;
        this.mPointers[last].index = last;
        return result;
    }

    // The value of a heap record must be modified through this function
    // call. The side effect is that the heap is updated accordingly to
    // restore the data structure to a min-heap. The input 'record' should be
    // a record returned by insert(key, value); see the comments for the
    // insert() function.
    update(record: MinHeapRecord<KeyType, ValueType> | null, value: ValueType): void {
        // Return immediately on invalid record.
        if (!record) {
            return;
        }

        if (this.lt(record.value, value)) {
            record.value = value;

            // The new value is larger than the old value. Propagate it
            // toward the leaves.
            let parent = record.index;
            let child = 2 * parent + 1;
            while (child < this.mNumElements) {
                // At least one child exists. Locate the one of maximum
                // value.
                let maxChild: number;
                const childP1 = child + 1;
                if (childP1 < this.mNumElements) {
                    // Two children exist.
                    if (this.le(this.mPointers[child].value, this.mPointers[childP1].value)) {
                        maxChild = child;
                    } else {
                        maxChild = childP1;
                    }
                } else {
                    // One child exists.
                    maxChild = child;
                }

                if (this.le(value, this.mPointers[maxChild].value)) {
                    // The new value is in the correct place to restore the
                    // tree to a heap.
                    break;
                }

                // The child has a larger value than the parent's value. Swap
                // the parent and child:

                // Move the child into the parent's slot.
                this.mPointers[parent] = this.mPointers[maxChild];
                this.mPointers[parent].index = parent;

                // Move the parent into the child's slot.
                this.mPointers[maxChild] = record;
                this.mPointers[maxChild].index = maxChild;

                parent = maxChild;
                child = 2 * parent + 1;
            }
        } else if (this.lt(value, record.value)) {
            record.value = value;

            // The new weight is smaller than the old weight. Propagate it
            // toward the root.
            let child = record.index;
            while (child > 0) {
                // A parent exists.
                const parent = (child - 1) >> 1;

                if (this.le(this.mPointers[parent].value, value)) {
                    // The new value is in the correct place to restore the
                    // tree to a heap.
                    break;
                }

                // The parent has a larger value than the child's value. Swap
                // the child and parent:

                // Move the parent into the child's slot.
                this.mPointers[child] = this.mPointers[parent];
                this.mPointers[child].index = child;

                // Move the child into the parent's slot.
                this.mPointers[parent] = record;
                this.mPointers[parent].index = parent;

                child = parent;
            }
        }
    }

    // Support for debugging. The function tests whether the data structure
    // is a valid min-heap.
    isValid(): boolean {
        for (let child = 0; child < this.mNumElements; ++child) {
            const parent = (child - 1) >> 1;
            if (child > 0) {
                if (this.lt(this.mPointers[child].value, this.mPointers[parent].value)) {
                    return false;
                }

                if (this.mPointers[parent].index !== parent) {
                    return false;
                }
            }
        }

        return true;
    }
}
