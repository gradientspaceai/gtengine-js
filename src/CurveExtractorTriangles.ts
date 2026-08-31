// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CurveExtractorTriangles.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The level set extraction algorithm implemented here is described in
// Section 2 of the document
// https://www.geometrictools.com/Documentation/ExtractLevelCurves.pdf
//
// Each image square is split into two triangles, with the diagonal chosen
// according to the parity of the square, and the level curve is extracted
// per triangle. Because the function is linear on each triangle, the level
// curve within a triangle is a line segment (or a vertex or the entire
// triangle boundary in the degenerate cases).
//
// Port notes: see src/CurveExtractor.ts and src/CurveExtractorSquares.ts.

import {
    CurveExtractor, CurveExtractorEdge, CurveExtractorVertex
} from './CurveExtractor';

export class CurveExtractorTriangles extends CurveExtractor {
    // The input is a 2D image with lexicographically ordered pixels (x,y)
    // stored in a linear array. Pixel (x,y) is stored in the array at
    // location index = x + xBound * y. The inputs xBound and yBound must
    // each be 2 or larger so that there is at least one image square to
    // process. The inputPixels must contain at least xBound * yBound
    // elements.
    constructor(xBound: number, yBound: number, inputPixels: ArrayLike<number>) {
        super(xBound, yBound, inputPixels);
    }

    // Extract level curves and return rational vertices. Use extractReal of
    // the base class if you want real-valued vertices.
    extract(level: number): {
        vertices: CurveExtractorVertex[];
        edges: CurveExtractorEdge[];
    } {
        // Adjust the image so that the level set is F(x,y) = 0.
        for (let i = 0; i < this.mPixels.length; ++i) {
            this.mPixels[i] = this.mInputPixels[i] - level;
        }

        const vertices: CurveExtractorVertex[] = [];
        const edges: CurveExtractorEdge[] = [];
        for (let y = 0, yp = 1; yp < this.mYBound; ++y, ++yp) {
            const yParity = y & 1;

            for (let x = 0, xp = 1; xp < this.mXBound; ++x, ++xp) {
                const xParity = x & 1;

                // Get the image values at the corners of the square.
                const i00 = x + this.mXBound * y;
                const i10 = i00 + 1;
                const i01 = i00 + this.mXBound;
                const i11 = i10 + this.mXBound;
                const f00 = this.mPixels[i00];
                const f10 = this.mPixels[i10];
                const f01 = this.mPixels[i01];
                const f11 = this.mPixels[i11];

                // Construct the vertices and edges of the level curve in the
                // square. The diagonal of the split alternates with the
                // parity of the square so that the extraction is symmetric.
                if (xParity === yParity) {
                    this.processTriangle(vertices, edges, x, y, f00, x, yp, f01, xp, y, f10);
                    this.processTriangle(vertices, edges, xp, yp, f11, xp, y, f10, x, yp, f01);
                } else {
                    this.processTriangle(vertices, edges, x, yp, f01, xp, yp, f11, x, y, f00);
                    this.processTriangle(vertices, edges, xp, y, f10, x, y, f00, xp, yp, f11);
                }
            }
        }
        return { vertices, edges };
    }

    protected processTriangle(vertices: CurveExtractorVertex[], edges: CurveExtractorEdge[],
        x0: number, y0: number, f0: number,
        x1: number, y1: number, f1: number,
        x2: number, y2: number, f2: number): void {
        let xn0: number, yn0: number, xn1: number, yn1: number, d0: number, d1: number;

        if (f0 !== 0) {
            // Convert to case "+**".
            if (f0 < 0) {
                f0 = 0 - f0;
                f1 = 0 - f1;
                f2 = 0 - f2;
            }

            if (f1 > 0) {
                if (f2 > 0) {
                    // +++
                    return;
                } else if (f2 < 0) {
                    // ++-
                    d0 = f0 - f2;
                    xn0 = f0 * x2 - f2 * x0;
                    yn0 = f0 * y2 - f2 * y0;
                    d1 = f1 - f2;
                    xn1 = f1 * x2 - f2 * x1;
                    yn1 = f1 * y2 - f2 * y1;
                    this.addEdge(vertices, edges, xn0, d0, yn0, d0, xn1, d1, yn1, d1);
                } else {
                    // ++0
                    this.addVertex(vertices, x2, 1, y2, 1);
                }
            } else if (f1 < 0) {
                d0 = f0 - f1;
                xn0 = f0 * x1 - f1 * x0;
                yn0 = f0 * y1 - f1 * y0;

                if (f2 > 0) {
                    // +-+
                    d1 = f2 - f1;
                    xn1 = f2 * x1 - f1 * x2;
                    yn1 = f2 * y1 - f1 * y2;
                    this.addEdge(vertices, edges, xn0, d0, yn0, d0, xn1, d1, yn1, d1);
                } else if (f2 < 0) {
                    // +--
                    d1 = f2 - f0;
                    xn1 = f2 * x0 - f0 * x2;
                    yn1 = f2 * y0 - f0 * y2;
                    this.addEdge(vertices, edges, xn0, d0, yn0, d0, xn1, d1, yn1, d1);
                } else {
                    // +-0
                    this.addEdge(vertices, edges, x2, 1, y2, 1, xn0, d0, yn0, d0);
                }
            } else {
                if (f2 > 0) {
                    // +0+
                    this.addVertex(vertices, x1, 1, y1, 1);
                } else if (f2 < 0) {
                    // +0-
                    d0 = f2 - f0;
                    xn0 = f2 * x0 - f0 * x2;
                    yn0 = f2 * y0 - f0 * y2;
                    this.addEdge(vertices, edges, x1, 1, y1, 1, xn0, d0, yn0, d0);
                } else {
                    // +00
                    this.addEdge(vertices, edges, x1, 1, y1, 1, x2, 1, y2, 1);
                }
            }
        } else if (f1 !== 0) {
            // Convert to case "0+*".
            if (f1 < 0) {
                f1 = 0 - f1;
                f2 = 0 - f2;
            }

            if (f2 > 0) {
                // 0++
                this.addVertex(vertices, x0, 1, y0, 1);
            } else if (f2 < 0) {
                // 0+-
                d0 = f1 - f2;
                xn0 = f1 * x2 - f2 * x1;
                yn0 = f1 * y2 - f2 * y1;
                this.addEdge(vertices, edges, x0, 1, y0, 1, xn0, d0, yn0, d0);
            } else {
                // 0+0
                this.addEdge(vertices, edges, x0, 1, y0, 1, x2, 1, y2, 1);
            }
        } else if (f2 !== 0) {
            // Cases 00+ or 00-.
            this.addEdge(vertices, edges, x0, 1, y0, 1, x1, 1, y1, 1);
        } else {
            // Case 000.
            this.addEdge(vertices, edges, x0, 1, y0, 1, x1, 1, y1, 1);
            this.addEdge(vertices, edges, x1, 1, y1, 1, x2, 1, y2, 1);
            this.addEdge(vertices, edges, x2, 1, y2, 1, x0, 1, y0, 1);
        }
    }
}
