// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CurveExtractorSquares.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The level set extraction algorithm implemented here is described in
// Section 3 of the document
// https://www.geometrictools.com/Documentation/ExtractLevelCurves.pdf
//
// The image is partitioned into squares whose corners are the pixel
// centers. Each square is processed by classifying the signs of the four
// corner values of F(x,y) = image(x,y) - level and emitting the level-curve
// vertices and edges for that sign pattern. The vertices are rational
// numbers, so the extraction is exact for integer-valued images.
//
// Port notes: see src/CurveExtractor.ts. The abstract Extract that returns
// rational vertices is named extract(level) and returns the vertices and
// edges rather than filling reference parameters. The upstream int64_t
// arithmetic becomes number arithmetic; the products f * x remain exact for
// images whose pixel magnitudes and dimensions are far below 2^26.

import {
    CurveExtractor, CurveExtractorEdge, CurveExtractorVertex
} from './CurveExtractor.js';

export class CurveExtractorSquares extends CurveExtractor {
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
            for (let x = 0, xp = 1; xp < this.mXBound; ++x, ++xp) {
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
                // square.
                this.processSquare(vertices, edges, x, xp, y, yp, f00, f10, f11, f01);
            }
        }
        return { vertices, edges };
    }

    protected processSquare(vertices: CurveExtractorVertex[], edges: CurveExtractorEdge[],
        x: number, xp: number, y: number, yp: number,
        f00: number, f10: number, f11: number, f01: number): void {
        let xn0: number, yn0: number, xn1: number, yn1: number;
        let d0: number, d1: number, d2: number, d3: number, det: number;

        if (f00 !== 0) {
            // Convert to case "+***".
            if (f00 < 0) {
                f00 = 0 - f00;
                f10 = 0 - f10;
                f11 = 0 - f11;
                f01 = 0 - f01;
            }

            if (f10 > 0) {
                if (f11 > 0) {
                    if (f01 > 0) {
                        // ++++
                        return;
                    } else if (f01 < 0) {
                        // +++-
                        d0 = f11 - f01;
                        xn0 = f11 * x - f01 * xp;
                        d1 = f00 - f01;
                        yn1 = f00 * yp - f01 * y;
                        this.addEdge(vertices, edges, xn0, d0, yp, 1, x, 1, yn1, d1);
                    } else {
                        // +++0
                        this.addVertex(vertices, x, 1, yp, 1);
                    }
                } else if (f11 < 0) {
                    d0 = f10 - f11;
                    yn0 = f10 * yp - f11 * y;

                    if (f01 > 0) {
                        // ++-+
                        d1 = f01 - f11;
                        xn1 = f01 * xp - f11 * x;
                        this.addEdge(vertices, edges, xp, 1, yn0, d0, xn1, d1, yp, 1);
                    } else if (f01 < 0) {
                        // ++--
                        d1 = f01 - f00;
                        yn1 = f01 * y - f00 * yp;
                        this.addEdge(vertices, edges, x, 1, yn1, d1, xp, 1, yn0, d0);
                    } else {
                        // ++-0
                        this.addEdge(vertices, edges, x, 1, yp, 1, xp, 1, yn0, d0);
                    }
                } else {
                    if (f01 > 0) {
                        // ++0+
                        this.addVertex(vertices, xp, 1, yp, 1);
                    } else if (f01 < 0) {
                        // ++0-
                        d0 = f01 - f00;
                        yn0 = f01 * y - f00 * yp;
                        this.addEdge(vertices, edges, xp, 1, yp, 1, x, 1, yn0, d0);
                    } else {
                        // ++00
                        this.addEdge(vertices, edges, xp, 1, yp, 1, x, 1, yp, 1);
                    }
                }
            } else if (f10 < 0) {
                d0 = f00 - f10;
                xn0 = f00 * xp - f10 * x;

                if (f11 > 0) {
                    d1 = f11 - f10;
                    yn1 = f11 * y - f10 * yp;

                    if (f01 > 0) {
                        // +-++
                        this.addEdge(vertices, edges, xn0, d0, y, 1, xp, 1, yn1, d1);
                    } else if (f01 < 0) {
                        // +-+-
                        d3 = f11 - f01;
                        xn1 = f11 * x - f01 * xp;
                        d2 = f01 - f00;
                        yn0 = f01 * y - f00 * yp;

                        if (d0 * d3 > 0) {
                            det = xn1 * d0 - xn0 * d3;
                        } else {
                            det = xn0 * d3 - xn1 * d0;
                        }

                        if (det > 0) {
                            this.addEdge(vertices, edges, xn1, d3, yp, 1, xp, 1, yn1, d1);
                            this.addEdge(vertices, edges, xn0, d0, y, 1, x, 1, yn0, d2);
                        } else if (det < 0) {
                            this.addEdge(vertices, edges, xn1, d3, yp, 1, x, 1, yn0, d2);
                            this.addEdge(vertices, edges, xn0, d0, y, 1, xp, 1, yn1, d1);
                        } else {
                            this.addEdge(vertices, edges, xn0, d0, yn0, d2, xn0, d0, y, 1);
                            this.addEdge(vertices, edges, xn0, d0, yn0, d2, xn0, d0, yp, 1);
                            this.addEdge(vertices, edges, xn0, d0, yn0, d2, x, 1, yn0, d2);
                            this.addEdge(vertices, edges, xn0, d0, yn0, d2, xp, 1, yn0, d2);
                        }
                    } else {
                        // +-+0
                        this.addEdge(vertices, edges, xn0, d0, y, 1, xp, 1, yn1, d1);
                        this.addVertex(vertices, x, 1, yp, 1);
                    }
                } else if (f11 < 0) {
                    if (f01 > 0) {
                        // +--+
                        d1 = f11 - f01;
                        xn1 = f11 * x - f01 * xp;
                        this.addEdge(vertices, edges, xn0, d0, y, 1, xn1, d1, yp, 1);
                    } else if (f01 < 0) {
                        // +---
                        d1 = f01 - f00;
                        yn1 = f01 * y - f00 * yp;
                        this.addEdge(vertices, edges, x, 1, yn1, d1, xn0, d0, y, 1);
                    } else {
                        // +--0
                        this.addEdge(vertices, edges, x, 1, yp, 1, xn0, d0, y, 1);
                    }
                } else {
                    if (f01 > 0) {
                        // +-0+
                        this.addEdge(vertices, edges, xp, 1, yp, 1, xn0, d0, y, 1);
                    } else if (f01 < 0) {
                        // +-0-
                        d1 = f01 - f00;
                        yn1 = f01 * y - f00 * yp;
                        this.addEdge(vertices, edges, x, 1, yn1, d1, xn0, d0, y, 1);
                        this.addVertex(vertices, xp, 1, yp, 1);
                    } else {
                        // +-00
                        this.addEdge(vertices, edges, xp, 1, yp, 1, xn0, d0, yp, 1);
                        this.addEdge(vertices, edges, xn0, d0, yp, 1, x, 1, yp, 1);
                        this.addEdge(vertices, edges, xn0, d0, yp, 1, xn0, d0, y, 1);
                    }
                }
            } else {
                if (f11 > 0) {
                    if (f01 > 0) {
                        // +0++
                        this.addVertex(vertices, xp, 1, y, 1);
                    } else if (f01 < 0) {
                        // +0+-
                        d0 = f11 - f01;
                        xn0 = f11 * x - f01 * xp;
                        d1 = f00 - f01;
                        yn1 = f00 * yp - f01 * y;
                        this.addEdge(vertices, edges, xn0, d0, yp, 1, x, 1, yn1, d1);
                        this.addVertex(vertices, xp, 1, y, 1);
                    } else {
                        // +0+0
                        this.addVertex(vertices, xp, 1, y, 1);
                        this.addVertex(vertices, x, 1, yp, 1);
                    }
                } else if (f11 < 0) {
                    if (f01 > 0) {
                        // +0-+
                        d0 = f11 - f01;
                        xn0 = f11 * x - f01 * xp;
                        this.addEdge(vertices, edges, xp, 1, y, 1, xn0, d0, yp, 1);
                    } else if (f01 < 0) {
                        // +0--
                        d0 = f01 - f00;
                        yn0 = f01 * y - f00 * yp;
                        this.addEdge(vertices, edges, xp, 1, y, 1, x, 1, yn0, d0);
                    } else {
                        // +0-0
                        this.addEdge(vertices, edges, xp, 1, y, 1, x, 1, yp, 1);
                    }
                } else {
                    if (f01 > 0) {
                        // +00+
                        this.addEdge(vertices, edges, xp, 1, y, 1, xp, 1, yp, 1);
                    } else if (f01 < 0) {
                        // +00-
                        d0 = f00 - f01;
                        yn0 = f00 * yp - f01 * y;
                        this.addEdge(vertices, edges, xp, 1, y, 1, xp, 1, yn0, d0);
                        this.addEdge(vertices, edges, xp, 1, yn0, d0, xp, 1, yp, 1);
                        this.addEdge(vertices, edges, xp, 1, yn0, d0, x, 1, yn0, d0);
                    } else {
                        // +000
                        this.addEdge(vertices, edges, x, 1, yp, 1, x, 1, y, 1);
                        this.addEdge(vertices, edges, x, 1, y, 1, xp, 1, y, 1);
                    }
                }
            }
        } else if (f10 !== 0) {
            // Convert to case "0+**".
            if (f10 < 0) {
                f10 = 0 - f10;
                f11 = 0 - f11;
                f01 = 0 - f01;
            }

            if (f11 > 0) {
                if (f01 > 0) {
                    // 0+++
                    this.addVertex(vertices, x, 1, y, 1);
                } else if (f01 < 0) {
                    // 0++-
                    d0 = f11 - f01;
                    xn0 = f11 * x - f01 * xp;
                    this.addEdge(vertices, edges, x, 1, y, 1, xn0, d0, yp, 1);
                } else {
                    // 0++0
                    this.addEdge(vertices, edges, x, 1, yp, 1, x, 1, y, 1);
                }
            } else if (f11 < 0) {
                if (f01 > 0) {
                    // 0+-+
                    d0 = f10 - f11;
                    yn0 = f10 * yp - f11 * y;
                    d1 = f01 - f11;
                    xn1 = f01 * xp - f11 * x;
                    this.addEdge(vertices, edges, xp, 1, yn0, d0, xn1, d1, yp, 1);
                    this.addVertex(vertices, x, 1, y, 1);
                } else if (f01 < 0) {
                    // 0+--
                    d0 = f10 - f11;
                    yn0 = f10 * yp - f11 * y;
                    this.addEdge(vertices, edges, x, 1, y, 1, xp, 1, yn0, d0);
                } else {
                    // 0+-0
                    d0 = f10 - f11;
                    yn0 = f10 * yp - f11 * y;
                    this.addEdge(vertices, edges, x, 1, y, 1, x, 1, yn0, d0);
                    this.addEdge(vertices, edges, x, 1, yn0, d0, x, 1, yp, 1);
                    this.addEdge(vertices, edges, x, 1, yn0, d0, xp, 1, yn0, d0);
                }
            } else {
                if (f01 > 0) {
                    // 0+0+
                    this.addVertex(vertices, x, 1, y, 1);
                    this.addVertex(vertices, xp, 1, yp, 1);
                } else if (f01 < 0) {
                    // 0+0-
                    this.addEdge(vertices, edges, x, 1, y, 1, xp, 1, yp, 1);
                } else {
                    // 0+00
                    this.addEdge(vertices, edges, xp, 1, yp, 1, x, 1, yp, 1);
                    this.addEdge(vertices, edges, x, 1, yp, 1, x, 1, y, 1);
                }
            }
        } else if (f11 !== 0) {
            // Convert to case "00+*".
            if (f11 < 0) {
                f11 = 0 - f11;
                f01 = 0 - f01;
            }

            if (f01 > 0) {
                // 00++
                this.addEdge(vertices, edges, x, 1, y, 1, xp, 1, y, 1);
            } else if (f01 < 0) {
                // 00+-
                d0 = f01 - f11;
                xn0 = f01 * xp - f11 * x;
                this.addEdge(vertices, edges, x, 1, y, 1, xn0, d0, y, 1);
                this.addEdge(vertices, edges, xn0, d0, y, 1, xp, 1, y, 1);
                this.addEdge(vertices, edges, xn0, d0, y, 1, xn0, d0, yp, 1);
            } else {
                // 00+0
                this.addEdge(vertices, edges, xp, 1, y, 1, xp, 1, yp, 1);
                this.addEdge(vertices, edges, xp, 1, yp, 1, x, 1, yp, 1);
            }
        } else if (f01 !== 0) {
            // Cases 000+ or 000-.
            this.addEdge(vertices, edges, x, 1, y, 1, xp, 1, y, 1);
            this.addEdge(vertices, edges, xp, 1, y, 1, xp, 1, yp, 1);
        } else {
            // Case 0000.
            this.addEdge(vertices, edges, x, 1, y, 1, xp, 1, y, 1);
            this.addEdge(vertices, edges, xp, 1, y, 1, xp, 1, yp, 1);
            this.addEdge(vertices, edges, xp, 1, yp, 1, x, 1, yp, 1);
            this.addEdge(vertices, edges, x, 1, yp, 1, x, 1, y, 1);
        }
    }
}
