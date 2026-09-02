// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrConvexPolygonHyperplane.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The intersection queries are based on the document
// https://www.geometrictools.com/Documentation/ClipConvexPolygonByHyperplane.pdf
//
// Port notes: the queries take the convex polygon as an array of vertices
// (upstream 'std::vector<Vector<N,Real>>'). Upstream declares two identical
// 'Configuration' enums, one in the TIQuery and one in the FIQuery; the port
// exports a single shared enum IntrConvexPolygonHyperplaneConfiguration
// (library-wide export uniqueness, PORTING.md). The FI 'SplitPolygon' helper
// is a module-private function. The vertices stored in the result arrays are
// clones, matching C++ value semantics.

import { Hyperplane } from './Hyperplane';
import { Vector, add, dot, mul } from './Vector';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';

// The relationship of the convex polygon to the hyperplane.
export enum IntrConvexPolygonHyperplaneConfiguration {
    // The polygon has vertices strictly on both sides of the hyperplane.
    SPLIT,

    // The polygon is on the positive side of the hyperplane and touches it
    // in exactly one vertex.
    POSITIVE_SIDE_VERTEX,

    // The polygon is on the positive side of the hyperplane and touches it
    // in an edge (two or more vertices have zero height).
    POSITIVE_SIDE_EDGE,

    // The polygon is strictly on the positive side of the hyperplane.
    POSITIVE_SIDE_STRICT,

    // The polygon is on the negative side of the hyperplane and touches it
    // in exactly one vertex.
    NEGATIVE_SIDE_VERTEX,

    // The polygon is on the negative side of the hyperplane and touches it
    // in an edge (two or more vertices have zero height).
    NEGATIVE_SIDE_EDGE,

    // The polygon is strictly on the negative side of the hyperplane.
    NEGATIVE_SIDE_STRICT,

    // All polygon vertices are on the hyperplane.
    CONTAINED,

    // The input has fewer than 3 vertices.
    INVALID_POLYGON
}

// The result of IntrConvexPolygonHyperplaneTI.test.
export interface IntrConvexPolygonHyperplaneTIResult {
    intersect: boolean;
    configuration: IntrConvexPolygonHyperplaneConfiguration;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrConvexPolygonHyperplaneTIResult():
    IntrConvexPolygonHyperplaneTIResult {
    return {
        intersect: false,
        configuration: IntrConvexPolygonHyperplaneConfiguration.INVALID_POLYGON
    };
}

// The result of IntrConvexPolygonHyperplaneFI.find.
export interface IntrConvexPolygonHyperplaneFIResult {
    configuration: IntrConvexPolygonHyperplaneConfiguration;

    // The intersection is either empty, a single vertex, a single edge or
    // the polygon is contained by the hyperplane.
    intersection: Vector[];

    // If 'configuration' is POSITIVE_* or SPLIT, this polygon is the portion
    // of the query input 'polygon' on the positive side of the hyperplane
    // with possibly a vertex or edge on the hyperplane.
    positivePolygon: Vector[];

    // If 'configuration' is NEGATIVE_* or SPLIT, this polygon is the portion
    // of the query input 'polygon' on the negative side of the hyperplane
    // with possibly a vertex or edge on the hyperplane.
    negativePolygon: Vector[];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrConvexPolygonHyperplaneFIResult():
    IntrConvexPolygonHyperplaneFIResult {
    return {
        configuration: IntrConvexPolygonHyperplaneConfiguration.INVALID_POLYGON,
        intersection: [],
        positivePolygon: [],
        negativePolygon: []
    };
}

// Test-intersection query for a convex polygon and a hyperplane.
export class IntrConvexPolygonHyperplaneTI implements
    TIQuery<Vector[], Hyperplane, IntrConvexPolygonHyperplaneTIResult> {

    test(polygon: readonly Vector[], hyperplane: Hyperplane):
        IntrConvexPolygonHyperplaneTIResult {
        const result = defaultIntrConvexPolygonHyperplaneTIResult();
        const C = IntrConvexPolygonHyperplaneConfiguration;

        const numVertices = polygon.length;
        if (numVertices < 3) {
            // The convex polygon must have at least 3 vertices.
            result.intersect = false;
            result.configuration = C.INVALID_POLYGON;
            return result;
        }

        // Determine on which side of the hyperplane each vertex lies.
        let numPositive = 0, numNegative = 0, numZero = 0;
        for (let i = 0; i < numVertices; ++i) {
            const h = dot(hyperplane.normal, polygon[i]) - hyperplane.constant;
            if (h > 0) {
                ++numPositive;
            }
            else if (h < 0) {
                ++numNegative;
            }
            else {
                ++numZero;
            }
        }

        if (numPositive > 0) {
            if (numNegative > 0) {
                result.intersect = true;
                result.configuration = C.SPLIT;
            }
            else if (numZero === 0) {
                result.intersect = false;
                result.configuration = C.POSITIVE_SIDE_STRICT;
            }
            else if (numZero === 1) {
                result.intersect = true;
                result.configuration = C.POSITIVE_SIDE_VERTEX;
            }
            else {  // numZero > 1
                result.intersect = true;
                result.configuration = C.POSITIVE_SIDE_EDGE;
            }
        }
        else if (numNegative > 0) {
            if (numZero === 0) {
                result.intersect = false;
                result.configuration = C.NEGATIVE_SIDE_STRICT;
            }
            else if (numZero === 1) {
                // The polygon touches the plane in a vertex or an edge.
                result.intersect = true;
                result.configuration = C.NEGATIVE_SIDE_VERTEX;
            }
            else {  // numZero > 1
                result.intersect = true;
                result.configuration = C.NEGATIVE_SIDE_EDGE;
            }
        }
        else {  // numZero == numVertices
            result.intersect = true;
            result.configuration = C.CONTAINED;
        }

        return result;
    }
}

// Find-intersection query for a convex polygon and a hyperplane. The query
// also clips the polygon, producing the sub-polygons on the positive and
// negative sides of the hyperplane.
export class IntrConvexPolygonHyperplaneFI implements
    FIQuery<Vector[], Hyperplane, IntrConvexPolygonHyperplaneFIResult> {

    find(polygon: readonly Vector[], hyperplane: Hyperplane):
        IntrConvexPolygonHyperplaneFIResult {
        const result = defaultIntrConvexPolygonHyperplaneFIResult();
        const Cfg = IntrConvexPolygonHyperplaneConfiguration;

        const numVertices = polygon.length;
        if (numVertices < 3) {
            // The convex polygon must have at least 3 vertices.
            result.configuration = Cfg.INVALID_POLYGON;
            return result;
        }

        // Determine on which side of the hyperplane the vertices live. The
        // index maxPosIndex stores the index of the vertex on the positive
        // side of the hyperplane that is farthest from the hyperplane. The
        // index maxNegIndex stores the index of the vertex on the negative
        // side of the hyperplane that is farthest from the hyperplane. If one
        // or the other such vertex does not exist, the corresponding index
        // remains its initial invalid value.
        const height = new Array<number>(numVertices);
        const zeroHeightIndices: number[] = [];
        let numPositive = 0, numNegative = 0;
        let maxPosHeight = -Number.MAX_VALUE;
        let maxNegHeight = Number.MAX_VALUE;
        let maxPosIndex = -1;
        let maxNegIndex = -1;
        for (let i = 0; i < numVertices; ++i) {
            height[i] = dot(hyperplane.normal, polygon[i]) - hyperplane.constant;
            if (height[i] > 0) {
                ++numPositive;
                if (height[i] > maxPosHeight) {
                    maxPosHeight = height[i];
                    maxPosIndex = i;
                }
            }
            else if (height[i] < 0) {
                ++numNegative;
                if (height[i] < maxNegHeight) {
                    maxNegHeight = height[i];
                    maxNegIndex = i;
                }
            }
            else {
                zeroHeightIndices.push(i);
            }
        }

        if (numPositive > 0) {
            if (numNegative > 0) {
                result.configuration = Cfg.SPLIT;

                const doSwap = (maxPosHeight < -maxNegHeight);
                if (doSwap) {
                    for (let i = 0; i < height.length; ++i) {
                        height[i] = -height[i];
                    }
                    const save = maxPosIndex;
                    maxPosIndex = maxNegIndex;
                    maxNegIndex = save;
                }

                splitPolygon(polygon, height, maxPosIndex, result);

                if (doSwap) {
                    const save = result.positivePolygon;
                    result.positivePolygon = result.negativePolygon;
                    result.negativePolygon = save;
                }
            }
            else {
                const numZero = zeroHeightIndices.length;
                if (numZero === 0) {
                    result.configuration = Cfg.POSITIVE_SIDE_STRICT;
                }
                else if (numZero === 1) {
                    result.configuration = Cfg.POSITIVE_SIDE_VERTEX;
                    result.intersection = [polygon[zeroHeightIndices[0]].clone()];
                }
                else {  // numZero > 1
                    result.configuration = Cfg.POSITIVE_SIDE_EDGE;
                    result.intersection = [
                        polygon[zeroHeightIndices[0]].clone(),
                        polygon[zeroHeightIndices[1]].clone()
                    ];
                }
                result.positivePolygon = polygon.map(p => p.clone());
            }
        }
        else if (numNegative > 0) {
            const numZero = zeroHeightIndices.length;
            if (numZero === 0) {
                result.configuration = Cfg.NEGATIVE_SIDE_STRICT;
            }
            else if (numZero === 1) {
                result.configuration = Cfg.NEGATIVE_SIDE_VERTEX;
                result.intersection = [polygon[zeroHeightIndices[0]].clone()];
            }
            else {  // numZero > 1
                result.configuration = Cfg.NEGATIVE_SIDE_EDGE;
                result.intersection = [
                    polygon[zeroHeightIndices[0]].clone(),
                    polygon[zeroHeightIndices[1]].clone()
                ];
            }
            result.negativePolygon = polygon.map(p => p.clone());
        }
        else {  // numZero == numVertices
            result.configuration = Cfg.CONTAINED;
            result.intersection = polygon.map(p => p.clone());
        }

        return result;
    }
}

// Clip the polygon by the hyperplane, where the polygon has vertices strictly
// on both sides of the hyperplane. The vertex of maximum positive height is
// polygon[maxPosIndex].
function splitPolygon(polygon: readonly Vector[], height: readonly number[],
    maxPosIndex: number, result: IntrConvexPolygonHyperplaneFIResult): void {
    // Find the largest contiguous subset of indices for which height[i] >= 0.
    const numVertices = polygon.length;
    const positiveList: Vector[] = [polygon[maxPosIndex].clone()];
    let end0 = maxPosIndex;
    let end0prev = -1;
    for (let i = 0; i < numVertices; ++i) {
        end0prev = (end0 + numVertices - 1) % numVertices;
        if (height[end0prev] >= 0) {
            positiveList.unshift(polygon[end0prev].clone());
            end0 = end0prev;
        }
        else {
            break;
        }
    }

    let end1 = maxPosIndex;
    let end1next = -1;
    for (let i = 0; i < numVertices; ++i) {
        end1next = (end1 + 1) % numVertices;
        if (height[end1next] >= 0) {
            positiveList.push(polygon[end1next].clone());
            end1 = end1next;
        }
        else {
            break;
        }
    }

    let index = end1next;
    const negativeList: Vector[] = [];
    for (let i = 0; i < numVertices; ++i) {
        negativeList.push(polygon[index].clone());
        index = (index + 1) % numVertices;
        if (index === end0) {
            break;
        }
    }

    // Clip the polygon.
    if (height[end0] > 0) {
        const t = -height[end0prev] / (height[end0] - height[end0prev]);
        const omt = 1 - t;
        const V = add(mul(polygon[end0prev], omt), mul(polygon[end0], t));
        positiveList.unshift(V.clone());
        negativeList.push(V.clone());
        result.intersection.push(V.clone());
    }
    else {
        negativeList.push(polygon[end0].clone());
        result.intersection.push(polygon[end0].clone());
    }

    if (height[end1] > 0) {
        const t = -height[end1next] / (height[end1] - height[end1next]);
        const omt = 1 - t;
        const V = add(mul(polygon[end1next], omt), mul(polygon[end1], t));
        positiveList.push(V.clone());
        negativeList.unshift(V.clone());
        result.intersection.push(V.clone());
    }
    else {
        negativeList.unshift(polygon[end1].clone());
        result.intersection.push(polygon[end1].clone());
    }

    result.positivePolygon = positiveList;
    result.negativePolygon = negativeList;
}
