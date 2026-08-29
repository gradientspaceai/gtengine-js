// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DisjointIntervals.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute Boolean operations of disjoint sets of half-open intervals of the
// form [xmin,xmax) with xmin < xmax.
//
// Port notes: the upstream template parameter Scalar becomes 'number'. The
// friend operators |, &, -, ^ become the static methods 'union',
// 'intersection', 'difference' and 'exclusiveOr'. The out-parameter function
// GetInterval returns { xmin, xmax } or null instead of a bool with
// reference outputs.

export class DisjointIntervals {
    // The array of endpoints has an even number of elements. The i-th
    // interval is [mEndpoints[2*i], mEndpoints[2*i+1]).
    private mEndpoints: number[];

    // Construction. The two-argument form requires xmin < xmax; otherwise
    // the set is empty.
    constructor(xmin?: number, xmax?: number) {
        this.mEndpoints = [];
        if (xmin !== undefined && xmax !== undefined && xmin < xmax) {
            this.mEndpoints = [xmin, xmax];
        }
    }

    // The port of the C++ copy constructor/assignment (TS objects alias, so
    // copying must be explicit).
    clone(): DisjointIntervals {
        const copy = new DisjointIntervals();
        copy.mEndpoints = this.mEndpoints.slice();
        return copy;
    }

    // The number of intervals in the set.
    getNumIntervals(): number {
        return Math.trunc(this.mEndpoints.length / 2);
    }

    // The i-th interval is [xmin,xmax). The values are valid (non-null) only
    // when 0 <= i < getNumIntervals().
    getInterval(i: number): { xmin: number, xmax: number } | null {
        const index = 2 * i;
        if (0 <= index && index < this.mEndpoints.length) {
            return { xmin: this.mEndpoints[index], xmax: this.mEndpoints[index + 1] };
        }
        return null;
    }

    // Make this set empty.
    clear(): void {
        this.mEndpoints = [];
    }

    // Insert [xmin,xmax) into the set. This is a Boolean 'union' operation.
    // The operation is successful only when xmin < xmax.
    insert(xmin: number, xmax: number): boolean {
        if (xmin < xmax) {
            const input = new DisjointIntervals(xmin, xmax);
            const output = DisjointIntervals.union(this, input);
            this.mEndpoints = output.mEndpoints;
            return true;
        }
        return false;
    }

    // Remove [xmin,xmax) from the set. This is a Boolean 'difference'
    // operation. The operation is successful only when xmin < xmax.
    remove(xmin: number, xmax: number): boolean {
        if (xmin < xmax) {
            const input = new DisjointIntervals(xmin, xmax);
            const output = DisjointIntervals.difference(this, input);
            this.mEndpoints = output.mEndpoints;
            return true;
        }
        return false;
    }

    // Get the union of the interval sets, input0 union input1. The port of
    // the upstream friend operator|.
    static union(input0: DisjointIntervals, input1: DisjointIntervals): DisjointIntervals {
        const output = new DisjointIntervals();

        const numEndpoints0 = input0.mEndpoints.length;
        const numEndpoints1 = input1.mEndpoints.length;
        let i0 = 0, i1 = 0;
        let parity0 = 0, parity1 = 0;
        while (i0 < numEndpoints0 && i1 < numEndpoints1) {
            const value0 = input0.mEndpoints[i0];
            const value1 = input1.mEndpoints[i1];

            if (value0 < value1) {
                if (parity0 === 0) {
                    parity0 = 1;
                    if (parity1 === 0) {
                        output.mEndpoints.push(value0);
                    }
                } else {
                    if (parity1 === 0) {
                        output.mEndpoints.push(value0);
                    }
                    parity0 = 0;
                }
                ++i0;
            } else if (value1 < value0) {
                if (parity1 === 0) {
                    parity1 = 1;
                    if (parity0 === 0) {
                        output.mEndpoints.push(value1);
                    }
                } else {
                    if (parity0 === 0) {
                        output.mEndpoints.push(value1);
                    }
                    parity1 = 0;
                }
                ++i1;
            } else {
                // value0 === value1
                if (parity0 === parity1) {
                    output.mEndpoints.push(value0);
                }
                parity0 ^= 1;
                parity1 ^= 1;
                ++i0;
                ++i1;
            }
        }

        while (i0 < numEndpoints0) {
            output.mEndpoints.push(input0.mEndpoints[i0]);
            ++i0;
        }

        while (i1 < numEndpoints1) {
            output.mEndpoints.push(input1.mEndpoints[i1]);
            ++i1;
        }

        return output;
    }

    // Get the intersection of the interval sets, input0 intersect input1.
    // The port of the upstream friend operator&.
    static intersection(input0: DisjointIntervals, input1: DisjointIntervals): DisjointIntervals {
        const output = new DisjointIntervals();

        const numEndpoints0 = input0.mEndpoints.length;
        const numEndpoints1 = input1.mEndpoints.length;
        let i0 = 0, i1 = 0;
        let parity0 = 0, parity1 = 0;
        while (i0 < numEndpoints0 && i1 < numEndpoints1) {
            const value0 = input0.mEndpoints[i0];
            const value1 = input1.mEndpoints[i1];

            if (value0 < value1) {
                if (parity0 === 0) {
                    parity0 = 1;
                    if (parity1 === 1) {
                        output.mEndpoints.push(value0);
                    }
                } else {
                    if (parity1 === 1) {
                        output.mEndpoints.push(value0);
                    }
                    parity0 = 0;
                }
                ++i0;
            } else if (value1 < value0) {
                if (parity1 === 0) {
                    parity1 = 1;
                    if (parity0 === 1) {
                        output.mEndpoints.push(value1);
                    }
                } else {
                    if (parity0 === 1) {
                        output.mEndpoints.push(value1);
                    }
                    parity1 = 0;
                }
                ++i1;
            } else {
                // value0 === value1
                if (parity0 === parity1) {
                    output.mEndpoints.push(value0);
                }
                parity0 ^= 1;
                parity1 ^= 1;
                ++i0;
                ++i1;
            }
        }

        return output;
    }

    // Get the difference of the interval sets, input0 minus input1. The port
    // of the upstream friend operator-.
    static difference(input0: DisjointIntervals, input1: DisjointIntervals): DisjointIntervals {
        const output = new DisjointIntervals();

        const numEndpoints0 = input0.mEndpoints.length;
        const numEndpoints1 = input1.mEndpoints.length;
        let i0 = 0, i1 = 0;
        let parity0 = 0, parity1 = 1;
        while (i0 < numEndpoints0 && i1 < numEndpoints1) {
            const value0 = input0.mEndpoints[i0];
            const value1 = input1.mEndpoints[i1];

            if (value0 < value1) {
                if (parity0 === 0) {
                    parity0 = 1;
                    if (parity1 === 1) {
                        output.mEndpoints.push(value0);
                    }
                } else {
                    if (parity1 === 1) {
                        output.mEndpoints.push(value0);
                    }
                    parity0 = 0;
                }
                ++i0;
            } else if (value1 < value0) {
                if (parity1 === 0) {
                    parity1 = 1;
                    if (parity0 === 1) {
                        output.mEndpoints.push(value1);
                    }
                } else {
                    if (parity0 === 1) {
                        output.mEndpoints.push(value1);
                    }
                    parity1 = 0;
                }
                ++i1;
            } else {
                // value0 === value1
                if (parity0 === parity1) {
                    output.mEndpoints.push(value0);
                }
                parity0 ^= 1;
                parity1 ^= 1;
                ++i0;
                ++i1;
            }
        }

        while (i0 < numEndpoints0) {
            output.mEndpoints.push(input0.mEndpoints[i0]);
            ++i0;
        }

        return output;
    }

    // Get the exclusive or of the interval sets, input0 xor input1 =
    // (input0 minus input1) or (input1 minus input0). The port of the
    // upstream friend operator^.
    static exclusiveOr(input0: DisjointIntervals, input1: DisjointIntervals): DisjointIntervals {
        const output = new DisjointIntervals();

        const numEndpoints0 = input0.mEndpoints.length;
        const numEndpoints1 = input1.mEndpoints.length;
        let i0 = 0, i1 = 0;
        while (i0 < numEndpoints0 && i1 < numEndpoints1) {
            const value0 = input0.mEndpoints[i0];
            const value1 = input1.mEndpoints[i1];

            if (value0 < value1) {
                output.mEndpoints.push(value0);
                ++i0;
            } else if (value1 < value0) {
                output.mEndpoints.push(value1);
                ++i1;
            } else {
                // value0 === value1
                ++i0;
                ++i1;
            }
        }

        while (i0 < numEndpoints0) {
            output.mEndpoints.push(input0.mEndpoints[i0]);
            ++i0;
        }

        while (i1 < numEndpoints1) {
            output.mEndpoints.push(input1.mEndpoints[i1]);
            ++i1;
        }

        return output;
    }
}
