// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinimumWidthPoints2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The width for a set of 2D points is the minimum distance between pairs of
// parallel lines, each pair bounding the points. The width for a set of 2D
// points is equal to the width for the set of vertices of the convex hull of
// the 2D points. It can be computed using the rotating calipers algorithm.
// For details about the rotating calipers algorithm and computing the width
// of a set of 2D points, see
// http://www-cgrl.cs.mcgill.ca/~godfried/research/calipers.html
// https://web.archive.org/web/20150330010154/http://cgm.cs.mcgill.ca/~orm/rotcal.html
//
// Port notes:
// * The four operator() overloads become compute(points, useRotatingCalipers)
//   and computeIndexed(points, indices, useRotatingCalipers). As upstream, an
//   empty 'indices' forwards to compute(points): both paths compute the
//   convex hull of the input, so the "the points already form a convex
//   polygon" overloads are really just a convenience for a subset selection.
// * The exact width comparison uses BSRational (upstream
//   BSRational<UIntegerAP32>), which is bigint-backed in the port.
// * The std::function GetVertex indirection becomes a plain closure.

import { logAssert } from './Logger.js';
import { BSRational } from './BSRational.js';
import { ConvexHull2 } from './ConvexHull2.js';
import { OrientedBox, type OrientedBox2 } from './OrientedBox.js';
import { RotatingCalipers, type RotatingCalipersAntipode } from './RotatingCalipers.js';
import { Vector, dot, normalize, sub } from './Vector.js';
import { dotPerp } from './Vector2.js';

// The exact squared distance from the antipodal vertex to the line through
// the antipodal edge.
function computeSqrWidth(vertices: readonly Vector[],
    antipode: RotatingCalipersAntipode): BSRational {
    const V = vertices[antipode.vertex];
    const E0 = vertices[antipode.edge[0]];
    const E1 = vertices[antipode.edge[1]];
    const rU: [BSRational, BSRational] = [
        BSRational.fromNumber(E1.get(0)).sub(BSRational.fromNumber(E0.get(0))),
        BSRational.fromNumber(E1.get(1)).sub(BSRational.fromNumber(E0.get(1)))
    ];
    const rDiff: [BSRational, BSRational] = [
        BSRational.fromNumber(V.get(0)).sub(BSRational.fromNumber(E0.get(0))),
        BSRational.fromNumber(V.get(1)).sub(BSRational.fromNumber(E0.get(1)))
    ];
    const rDotPerp = rU[1].mul(rDiff[0]).sub(rU[0].mul(rDiff[1]));
    const rSqrLenU = rU[0].mul(rU[0]).add(rU[1].mul(rU[1]));
    return rDotPerp.mul(rDotPerp).div(rSqrLenU);
}

// Compute the origin, the unit-length edge direction, and the extreme
// heights of the polygon along that direction.
function computeExtremeHeights(getVertex: (j: number) => Vector,
    numElements: number, i0Min: number, i1Min: number):
    { origin: Vector, U: Vector, minHeight: number, maxHeight: number } {
    const origin = getVertex(i0Min);
    const U = sub(getVertex(i1Min), origin);
    normalize(U);

    let minHeight = 0;
    let maxHeight = 0;
    for (let j = 0; j < numElements; ++j) {
        const height = dot(U, sub(getVertex(j), origin));
        if (height < minHeight) {
            minHeight = height;
        }
        else if (height > maxHeight) {
            maxHeight = height;
        }
    }

    return { origin, U, minHeight, maxHeight };
}

// Remove duplicate and collinear vertices from a counterclockwise convex
// polygon, returning the retained vertex indices.
//
// Upstream bug (MinimumWidthPoints2.h, ComputeMinWidth, the brute-force
// branch): the "arriving" edge of the DotPerp test is the immediately
// preceding edge of the input array, which may be the zero-length edge
// produced by a duplicate point. DotPerp is then zero and a genuine polygon
// corner is discarded along with the duplicate; this is the same defect
// found in RotatingCalipers::CreatePolygon (gtengine-js issue #286). It is
// unreachable here because the input is always a ConvexHull2 hull, which has
// no duplicate points, but the port applies the same fix for safety. The
// behavior is identical to upstream whenever the input has no duplicates.
function removeDuplicateAndCollinearVertices(
    vertices: readonly Vector[]): number[] {
    const numVertices = vertices.length;

    const edges: Vector[] = [];
    const isZeroEdge: boolean[] = [];
    for (let i = 0; i < numVertices; ++i) {
        const edge = sub(vertices[(i + 1) % numVertices], vertices[i]);
        edges.push(edge);
        isZeroEdge.push(edge.get(0) === 0 && edge.get(1) === 0);
    }

    let prevNonzero = -1;
    for (let i = numVertices - 1; i >= 0; --i) {
        if (!isZeroEdge[i]) {
            prevNonzero = i;
            break;
        }
    }

    const indices: number[] = [];
    if (prevNonzero < 0) {
        return indices;
    }

    for (let i0 = 0; i0 < numVertices; ++i0) {
        if (isZeroEdge[i0]) {
            continue;
        }
        if (dotPerp(edges[prevNonzero], edges[i0]) !== 0) {
            indices.push(i0);
        }
        prevNonzero = i0;
    }

    return indices;
}

export class MinimumWidthPoints2 {
    // The return value is an oriented box in 2D. The width of the point set
    // is in the direction box.axis[0]; the width is 2*box.extent[0]. The
    // corresponding height is in the direction box.axis[1] =
    // -Perp(box.axis[0]); the height is 2*box.extent[1].
    //
    // The points are arbitrary, so the convex hull must be computed from
    // them to obtain the convex polygon whose minimum width is the desired
    // output.
    compute(points: readonly Vector[],
        useRotatingCalipers: boolean = true): OrientedBox2 {
        logAssert(points.length >= 3, 'Invalid input.');

        const box = new OrientedBox(2);

        // Get the convex hull of the points.
        const ch2 = new ConvexHull2();
        ch2.compute(points);
        const dimension = ch2.getDimension();

        if (dimension === 0) {
            box.center = points[0].clone();
            box.axis[0] = Vector.unit(2, 0);
            box.axis[1] = Vector.unit(2, 1);
            box.extent = Vector.zero(2);
            return box;
        }

        if (dimension === 1) {
            // The points lie on a line. Determine the extreme t-values for
            // the points represented as P = origin + t*direction. The
            // 'origin' is an input vertex, so both t-extremes start at zero.
            const line = ch2.getLine();
            let tmin = 0;
            let tmax = 0;
            for (let i = 0; i < points.length; ++i) {
                const t = dot(sub(points[i], line.origin), line.direction);
                if (t > tmax) {
                    tmax = t;
                }
                else if (t < tmin) {
                    tmin = t;
                }
            }

            const half = 0.5;
            box.center = Vector.fromArray([
                line.origin.get(0) + half * (tmin + tmax) * line.direction.get(0),
                line.origin.get(1) + half * (tmin + tmax) * line.direction.get(1)
            ]);
            box.extent = Vector.fromArray([0, half * (tmax - tmin)]);
            // NOTE: upstream uses Perp(direction) here but -Perp(U) in
            // computeMinWidth, so the two branches produce boxes of opposite
            // handedness. The quirk is harmless (extent[0] is zero, so the
            // box is a segment either way) and is preserved.
            box.axis[0] = Vector.fromArray([
                line.direction.get(1), -line.direction.get(0)
            ]);
            box.axis[1] = line.direction.clone();
            return box;
        }

        // Get the indexed convex hull.
        const hull = ch2.getHull();
        const hullPoints = ch2.getPoints();
        const vertices = hull.map(h => hullPoints[h]);

        this.computeMinWidth(vertices, useRotatingCalipers, box);
        return box;
    }

    // The points already form a counterclockwise, nondegenerate convex
    // polygon. If the points directly are the convex polygon, pass an empty
    // 'indices' array. If the polygon vertices are a subset of the incoming
    // points, that subset is identified by an 'indices' array of at least 3
    // elements.
    computeIndexed(points: readonly Vector[], indices: readonly number[],
        useRotatingCalipers: boolean = true): OrientedBox2 {
        logAssert(points.length >= 3 &&
            (indices.length === 0 || indices.length >= 3), 'Invalid input.');

        if (indices.length > 0) {
            const compactPoints = indices.map(i => points[i]);
            return this.compute(compactPoints, useRotatingCalipers);
        }
        return this.compute(points, useRotatingCalipers);
    }

    private computeMinWidth(vertices: readonly Vector[],
        useRotatingCalipers: boolean, box: OrientedBox2): void {
        let getVertex: (j: number) => Vector;
        let numElements: number;
        let i0Min: number;
        let i1Min: number;
        let minWidth: number;

        if (useRotatingCalipers) {
            const antipodes = RotatingCalipers.computeAntipodes(vertices);
            logAssert(antipodes.length > 0, 'Antipodes must exist.');

            let minSqrWidth = computeSqrWidth(vertices, antipodes[0]);
            let minAntipode = 0;
            for (let i = 1; i < antipodes.length; ++i) {
                const sqrWidth = computeSqrWidth(vertices, antipodes[i]);
                if (sqrWidth.lessThan(minSqrWidth)) {
                    minSqrWidth = sqrWidth;
                    minAntipode = i;
                }
            }
            minWidth = Math.sqrt(minSqrWidth.toNumber());

            getVertex = (j: number) => vertices[j];
            numElements = vertices.length;
            i0Min = antipodes[minAntipode].edge[0];
            i1Min = antipodes[minAntipode].edge[1];
        }
        else {
            // Remove duplicate and collinear vertices.
            const indices = removeDuplicateAndCollinearVertices(vertices);
            const numIndices = indices.length;
            logAssert(numIndices >= 3,
                'The convex polygon must have at least 3 noncollinear vertices.');

            // Iterate over the polygon edges to search for the edge that
            // leads to the minimum width.
            minWidth = Number.MAX_VALUE;
            i0Min = numIndices - 1;
            i1Min = 0;
            for (let i0 = numIndices - 1, i1 = 0; i1 < numIndices; i0 = i1++) {
                const origin = vertices[indices[i0]];
                const U = sub(vertices[indices[i1]], origin);
                normalize(U);

                let maxWidth = 0;
                for (let j = 0; j < numIndices; ++j) {
                    const diff = sub(vertices[indices[j]], origin);
                    const width = U.get(0) * diff.get(1) - U.get(1) * diff.get(0);
                    if (width > maxWidth) {
                        maxWidth = width;
                    }
                }

                if (maxWidth < minWidth) {
                    minWidth = maxWidth;
                    i0Min = i0;
                    i1Min = i1;
                }
            }

            getVertex = (j: number) => vertices[indices[j]];
            numElements = numIndices;
        }

        const { origin, U, minHeight, maxHeight } = computeExtremeHeights(
            getVertex, numElements, i0Min, i1Min);

        const half = 0.5;
        box.extent = Vector.fromArray([
            half * minWidth, half * (maxHeight - minHeight)
        ]);
        // axis[0] = -Perp(U) = (-U1,U0), axis[1] = U.
        box.axis[0] = Vector.fromArray([-U.get(1), U.get(0)]);
        box.axis[1] = U;
        const s0 = box.extent.get(0);
        const s1 = half * (maxHeight + minHeight);
        box.center = Vector.fromArray([
            origin.get(0) + s0 * box.axis[0].get(0) + s1 * box.axis[1].get(0),
            origin.get(1) + s0 * box.axis[0].get(1) + s1 * box.axis[1].get(1)
        ]);
    }
}
