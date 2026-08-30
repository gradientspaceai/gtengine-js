// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DisjointRectangles.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute Boolean operations of disjoint sets of half-open rectangles of the
// form [xmin,xmax)x[ymin,ymax) with xmin < xmax and ymin < ymax.
//
// Port notes: the upstream template parameter Scalar becomes 'number'. The
// friend operators |, &, -, ^ become the static methods 'union',
// 'intersection', 'difference' and 'exclusiveOr' (matching the
// DisjointIntervals port). The out-parameter functions GetRectangle and
// GetStrip return object literals or null instead of a bool with reference
// outputs. The nested Strip class becomes the DisjointRectanglesStrip
// interface; getStrip returns a copy (upstream copies via assignment).

import { DisjointIntervals } from './DisjointIntervals';

// The rectangle set consists of y-strips of x-interval sets. This is the
// port of the upstream nested class DisjointRectangles::Strip.
export interface DisjointRectanglesStrip {
    ymin: number;
    ymax: number;
    intervalSet: DisjointIntervals;
}

export class DisjointRectangles {
    // The number of rectangles in the set.
    private mNumRectangles: number;

    // The y-strips of the set, each containing an x-interval set.
    private mStrips: DisjointRectanglesStrip[];

    // Construction. The four-argument form requires xmin < xmax and
    // ymin < ymax; otherwise the set is empty.
    constructor(xmin?: number, xmax?: number, ymin?: number, ymax?: number) {
        this.mStrips = [];
        if (xmin !== undefined && xmax !== undefined &&
            ymin !== undefined && ymax !== undefined &&
            xmin < xmax && ymin < ymax) {
            this.mNumRectangles = 1;
            this.mStrips.push({ ymin, ymax, intervalSet: new DisjointIntervals(xmin, xmax) });
        } else {
            this.mNumRectangles = 0;
        }
    }

    // The port of the C++ copy constructor/assignment (TS objects alias, so
    // copying must be explicit).
    clone(): DisjointRectangles {
        const copy = new DisjointRectangles();
        copy.mNumRectangles = this.mNumRectangles;
        copy.mStrips = this.mStrips.map((strip) => ({
            ymin: strip.ymin,
            ymax: strip.ymax,
            intervalSet: strip.intervalSet.clone()
        }));
        return copy;
    }

    // The number of rectangles in the set.
    getNumRectangles(): number {
        return this.mNumRectangles;
    }

    // The i-th rectangle is [xmin,xmax)x[ymin,ymax). The values are valid
    // (non-null) only when 0 <= i < getNumRectangles().
    getRectangle(i: number): { xmin: number, xmax: number, ymin: number, ymax: number } | null {
        let totalQuantity = 0;
        for (const strip of this.mStrips) {
            const intervalSet = strip.intervalSet;
            const xQuantity = intervalSet.getNumIntervals();
            const nextTotalQuantity = totalQuantity + xQuantity;
            if (i < nextTotalQuantity) {
                const interval = intervalSet.getInterval(i - totalQuantity);
                if (interval === null) {
                    return null;
                }
                return {
                    xmin: interval.xmin,
                    xmax: interval.xmax,
                    ymin: strip.ymin,
                    ymax: strip.ymax
                };
            }
            totalQuantity = nextTotalQuantity;
        }
        return null;
    }

    // Make this set empty.
    clear(): void {
        this.mNumRectangles = 0;
        this.mStrips = [];
    }

    // The number of y-strips in the set.
    getNumStrips(): number {
        return this.mStrips.length;
    }

    // The i-th strip. The returned value is valid (non-null) only when
    // 0 <= i < getNumStrips(). The interval set is a copy, matching the C++
    // out-parameter assignment.
    getStrip(i: number): DisjointRectanglesStrip | null {
        if (0 <= i && i < this.getNumStrips()) {
            const strip = this.mStrips[i];
            return {
                ymin: strip.ymin,
                ymax: strip.ymax,
                intervalSet: strip.intervalSet.clone()
            };
        }
        return null;
    }

    // Insert [xmin,xmax)x[ymin,ymax) into the set. This is a Boolean 'union'
    // operation. The operation is successful only when xmin < xmax and
    // ymin < ymax.
    insert(xmin: number, xmax: number, ymin: number, ymax: number): boolean {
        if (xmin < xmax && ymin < ymax) {
            const input = new DisjointRectangles(xmin, xmax, ymin, ymax);
            const output = DisjointRectangles.union(this, input);
            this.mNumRectangles = output.mNumRectangles;
            this.mStrips = output.mStrips;
            return true;
        }
        return false;
    }

    // Remove [xmin,xmax)x[ymin,ymax) from the set. This is a Boolean
    // 'difference' operation. The operation is successful only when
    // xmin < xmax and ymin < ymax.
    remove(xmin: number, xmax: number, ymin: number, ymax: number): boolean {
        if (xmin < xmax && ymin < ymax) {
            const input = new DisjointRectangles(xmin, xmax, ymin, ymax);
            const output = DisjointRectangles.difference(this, input);
            this.mNumRectangles = output.mNumRectangles;
            this.mStrips = output.mStrips;
            return true;
        }
        return false;
    }

    // Get the union of the rectangle sets, input0 union input1. The port of
    // the upstream friend operator|.
    static union(input0: DisjointRectangles, input1: DisjointRectangles): DisjointRectangles {
        return DisjointRectangles.execute(
            (i0, i1) => DisjointIntervals.union(i0, i1),
            true, true, input0, input1);
    }

    // Get the intersection of the rectangle sets, input0 intersect input1.
    // The port of the upstream friend operator&.
    static intersection(input0: DisjointRectangles, input1: DisjointRectangles): DisjointRectangles {
        return DisjointRectangles.execute(
            (i0, i1) => DisjointIntervals.intersection(i0, i1),
            false, false, input0, input1);
    }

    // Get the difference of the rectangle sets, input0 minus input1. The
    // port of the upstream friend operator-.
    static difference(input0: DisjointRectangles, input1: DisjointRectangles): DisjointRectangles {
        return DisjointRectangles.execute(
            (i0, i1) => DisjointIntervals.difference(i0, i1),
            false, true, input0, input1);
    }

    // Get the exclusive or of the rectangle sets, input0 xor input1 =
    // (input0 minus input1) or (input1 minus input0). The port of the
    // upstream friend operator^.
    static exclusiveOr(input0: DisjointRectangles, input1: DisjointRectangles): DisjointRectangles {
        return DisjointRectangles.execute(
            (i0, i1) => DisjointIntervals.exclusiveOr(i0, i1),
            true, true, input0, input1);
    }

    private static execute(
        operation: (i0: DisjointIntervals, i1: DisjointIntervals) => DisjointIntervals,
        unionExclusiveOr: boolean, unionExclusiveOrDifference: boolean,
        input0: DisjointRectangles, input1: DisjointRectangles): DisjointRectangles {
        const output = new DisjointRectangles();

        const numStrips0 = input0.getNumStrips();
        const numStrips1 = input1.getNumStrips();
        let i0 = 0, i1 = 0;
        let getOriginal0 = true, getOriginal1 = true;
        let ymin0 = 0;
        let ymax0 = 0;
        let ymin1 = 0;
        let ymax1 = 0;

        while (i0 < numStrips0 && i1 < numStrips1) {
            const intr0 = input0.mStrips[i0].intervalSet;
            if (getOriginal0) {
                ymin0 = input0.mStrips[i0].ymin;
                ymax0 = input0.mStrips[i0].ymax;
            }

            const intr1 = input1.mStrips[i1].intervalSet;
            if (getOriginal1) {
                ymin1 = input1.mStrips[i1].ymin;
                ymax1 = input1.mStrips[i1].ymax;
            }

            // Case 1.
            if (ymax1 <= ymin0) {
                // operation(empty,strip1)
                if (unionExclusiveOr) {
                    output.mStrips.push({ ymin: ymin1, ymax: ymax1, intervalSet: intr1.clone() });
                }

                ++i1;
                getOriginal0 = false;
                getOriginal1 = true;
                continue;  // using next ymin1/ymax1
            }

            // Case 11.
            if (ymin1 >= ymax0) {
                // operation(strip0,empty)
                if (unionExclusiveOrDifference) {
                    output.mStrips.push({ ymin: ymin0, ymax: ymax0, intervalSet: intr0.clone() });
                }

                ++i0;
                getOriginal0 = true;
                getOriginal1 = false;
                continue;  // using next ymin0/ymax0
            }

            // Reduce cases 2, 3, 4 to cases 5, 6, 7.
            if (ymin1 < ymin0) {
                // operation(empty,[ymin1,ymin0))
                if (unionExclusiveOr) {
                    output.mStrips.push({ ymin: ymin1, ymax: ymin0, intervalSet: intr1.clone() });
                }

                ymin1 = ymin0;
                getOriginal1 = false;
            }

            // Reduce cases 8, 9, 10 to cases 5, 6, 7.
            if (ymin1 > ymin0) {
                // operation([ymin0,ymin1),empty)
                if (unionExclusiveOrDifference) {
                    output.mStrips.push({ ymin: ymin0, ymax: ymin1, intervalSet: intr0.clone() });
                }

                ymin0 = ymin1;
                getOriginal0 = false;
            }

            // Case 5.
            if (ymax1 < ymax0) {
                // operation(strip0,[ymin1,ymax1))
                const result = operation(intr0, intr1);
                output.mStrips.push({ ymin: ymin1, ymax: ymax1, intervalSet: result });

                ymin0 = ymax1;
                ++i1;
                getOriginal0 = false;
                getOriginal1 = true;
                continue;  // using next ymin1/ymax1
            }

            // Case 6.
            if (ymax1 === ymax0) {
                // operation(strip0,[ymin1,ymax1))
                const result = operation(intr0, intr1);
                output.mStrips.push({ ymin: ymin1, ymax: ymax1, intervalSet: result });

                ++i0;
                ++i1;
                getOriginal0 = true;
                getOriginal1 = true;
                continue;  // using next ymin0/ymax0 and ymin1/ymax1
            }

            // Case 7.
            if (ymax1 > ymax0) {
                // operation(strip0,[ymin1,ymax0))
                const result = operation(intr0, intr1);
                output.mStrips.push({ ymin: ymin1, ymax: ymax0, intervalSet: result });

                ymin1 = ymax0;
                ++i0;
                getOriginal0 = true;
                getOriginal1 = false;
                // continue;  using current ymin1/ymax1
            }
        }

        if (unionExclusiveOrDifference) {
            while (i0 < numStrips0) {
                if (getOriginal0) {
                    ymin0 = input0.mStrips[i0].ymin;
                    ymax0 = input0.mStrips[i0].ymax;
                } else {
                    getOriginal0 = true;
                }

                // operation(strip0,empty)
                output.mStrips.push({
                    ymin: ymin0, ymax: ymax0,
                    intervalSet: input0.mStrips[i0].intervalSet.clone()
                });

                ++i0;
            }
        }

        if (unionExclusiveOr) {
            while (i1 < numStrips1) {
                if (getOriginal1) {
                    ymin1 = input1.mStrips[i1].ymin;
                    ymax1 = input1.mStrips[i1].ymax;
                } else {
                    getOriginal1 = true;
                }

                // operation(empty,strip1)
                output.mStrips.push({
                    ymin: ymin1, ymax: ymax1,
                    intervalSet: input1.mStrips[i1].intervalSet.clone()
                });

                ++i1;
            }
        }

        output.computeRectangleQuantity();
        return output;
    }

    private computeRectangleQuantity(): void {
        this.mNumRectangles = 0;
        for (const strip of this.mStrips) {
            this.mNumRectangles += strip.intervalSet.getNumIntervals();
        }
    }
}
