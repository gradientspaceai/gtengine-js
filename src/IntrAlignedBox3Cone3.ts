// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrAlignedBox3Cone3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test for intersection of a box and a cone. The cone can be infinite
//   0 <= minHeight < maxHeight = +infinity
// or finite (cone frustum)
//   0 <= minHeight < maxHeight < +infinity.
// The algorithm is described in
//   https://www.geometrictools.com/Documentation/IntersectionBoxCone.pdf
// and reports an intersection only when the intersection set has positive
// volume. For example, let the box be outside the cone. If the box is below
// the minHeight plane at the cone vertex and just touches the cone vertex, no
// intersection is reported. If the box is above the maxHeight plane and just
// touches the disk capping the cone, either at a single point, a line segment
// of points or a polygon of points, no intersection is reported.
//
// Upstream TODO (preserved): these queries were designed when an infinite
// cone was defined by choosing maxHeight of std::numeric_limits<T>::max().
// The Cone class has been redesigned not to use std::numeric_limits so that
// arithmetic systems without infinities (BSNumber, BSRational) are supported.
// The intersection queries have not been rewritten for the new class design.
// As upstream does, this port substitutes the largest finite value
// (Number.MAX_VALUE, the analogue of std::numeric_limits<double>::max()) for
// the maximum height of an infinite cone.
//
// Port notes: the upstream member-function-pointer table mConfiguration[81]
// becomes an array of arrow functions built in the constructor, indexed
// exactly as upstream. The 81 configuration handlers keep their upstream
// names in comments (NNNN_0, NNNZ_1, ...). The nested Result struct becomes
// the exported interface IntrAlignedBox3Cone3TIResult and operator() becomes
// test(), per the Intr* precedents. The protected static helper
// HasPointInsideCone is exported as a module function so that it can be
// tested and reused.

import type { AlignedBox } from './AlignedBox.js';
import type { Cone } from './Cone.js';
import { IntrRay3AlignedBox3TI } from './IntrRay3AlignedBox3.js';
import { IntrSegment3AlignedBox3TI } from './IntrSegment3AlignedBox3.js';
import { Ray } from './Ray.js';
import { Segment } from './Segment.js';
import type { TIQuery } from './TIQuery.js';
import { Vector, add, dot, length, mul, sub } from './Vector.js';
import { cross } from './Vector3.js';

// The result of the test-intersection query.
export interface IntrAlignedBox3Cone3TIResult {
    intersect: boolean;
}

// The port of the upstream Result default constructor.
export function defaultIntrAlignedBox3Cone3TIResult():
    IntrAlignedBox3Cone3TIResult {
    return { intersect: false };
}

// A box face is { v: [v0, v1, v2, v3], e: [e0, e1, e2, e3] }, where the face
// corner vertices are { v0, v1, v2, v3 } (indices into mVertices) and the
// indices { e0, e1, e2, e3 } are into mEdges. The index e0 refers to edge
// { v0, v1 }, e1 to { v1, v2 }, e2 to { v2, v3 } and e3 to { v3, v0 }. The
// ordering of vertices for the faces is counterclockwise when viewed from
// outside the box. The initial vertex has minimum index for all vertices of
// that face.
interface Face {
    v: [number, number, number, number];
    e: [number, number, number, number];
}

const NUM_BOX_VERTICES = 8;
const NUM_BOX_EDGES = 12;
const NUM_BOX_FACES = 6;
const MAX_VERTICES = 32;
const VERTEX_MIN_BASE = 8;
const VERTEX_MAX_BASE = 20;
const MAX_CANDIDATE_EDGES = 496;

// The stand-in for std::numeric_limits<T>::max() used for the maximum height
// of an infinite cone. See the comments at the beginning of the file.
const INFINITE_CONE_MAX_HEIGHT = Number.MAX_VALUE;

function coneMaxHeightOf(cone: Cone): number {
    return cone.isFinite() ? cone.getMaxHeight() : INFINITE_CONE_MAX_HEIGHT;
}

// Compute the interval of heights of the box relative to the cone axis.
export function intrAlignedBox3Cone3ComputeBoxHeightInterval(box: AlignedBox,
    cone: Cone): { boxMinHeight: number, boxMaxHeight: number } {
    const { center: C, extent: e } = box.getCenteredForm();
    const V = cone.ray.origin;
    const U = cone.ray.direction;
    const CmV = sub(C, V);
    const DdCmV = dot(U, CmV);
    const radius =
        e.values[0] * Math.abs(U.values[0]) +
        e.values[1] * Math.abs(U.values[1]) +
        e.values[2] * Math.abs(U.values[2]);
    return { boxMinHeight: DdCmV - radius, boxMaxHeight: DdCmV + radius };
}

// Determine whether the cone axis (the portion of it that is in the cone
// height range) intersects the box.
export function intrAlignedBox3Cone3ConeAxisIntersectsBox(box: AlignedBox,
    cone: Cone): boolean {
    if (cone.isFinite()) {
        const segment = new Segment(3);
        segment.p[0] = add(cone.ray.origin,
            mul(cone.getMinHeight(), cone.ray.direction));
        segment.p[1] = add(cone.ray.origin,
            mul(cone.getMaxHeight(), cone.ray.direction));
        const sbResult = new IntrSegment3AlignedBox3TI().test(segment, box);
        if (sbResult.intersect) {
            return true;
        }
    }
    else {
        const ray = new Ray(3);
        ray.origin = add(cone.ray.origin,
            mul(cone.getMinHeight(), cone.ray.direction));
        ray.direction = cone.ray.direction.clone();
        const rbResult = new IntrRay3AlignedBox3TI().test(ray, box);
        if (rbResult.intersect) {
            return true;
        }
    }
    return false;
}

// Determine whether the segment <P0,P1> has a point strictly inside the cone.
// The points P0 and P1 are relative to the cone vertex V; that is, the
// original points are X0 = P0 + V and X1 = P1 + V.
//
// Define F(X) = Dot(U,X - V)/|X - V|, where U is the unit-length cone axis
// direction and V is the cone vertex. The segment <P0,P1> and cone intersect
// when a segment point X is inside the cone; that is, when F(X) > cosAngle.
// The comparison is converted to an equivalent one that does not involve
// divisions in order to avoid a division by zero if a vertex or edge contains
// (0,0,0). The function is G(X) = Dot(U,X-V) - cosAngle*Length(X-V).
export function intrAlignedBox3Cone3HasPointInsideCone(P0: Vector, P1: Vector,
    cone: Cone): boolean {
    const U = cone.ray.direction;

    // Test whether P0 or P1 is inside the cone.
    let g = dot(U, P0) - cone.cosAngle * length(P0);
    if (g > 0) {
        // X0 = P0 + V is inside the cone.
        return true;
    }

    g = dot(U, P1) - cone.cosAngle * length(P1);
    if (g > 0) {
        // X1 = P1 + V is inside the cone.
        return true;
    }

    // Test whether an interior segment point is inside the cone.
    const E = sub(P1, P0);
    const crossP0U = cross(P0, U);
    const crossP0E = cross(P0, E);
    const dphi0 = dot(crossP0E, crossP0U);
    if (dphi0 > 0) {
        const crossP1U = cross(P1, U);
        const dphi1 = dot(crossP0E, crossP1U);
        if (dphi1 < 0) {
            const t = dphi0 / (dphi0 - dphi1);
            const PMax = add(P0, mul(t, E));
            g = dot(U, PMax) - cone.cosAngle * length(PMax);
            if (g > 0) {
                // The edge point XMax = PMax + V is inside the cone.
                return true;
            }
        }
    }

    return false;
}

// Test-intersection query for a solid aligned box and a solid cone in 3D.
export class IntrAlignedBox3Cone3TI implements
    TIQuery<AlignedBox, Cone, IntrAlignedBox3Cone3TIResult> {

    // The box topology is that of a cube whose vertices have components in
    // {0,1}. The cube vertices are indexed by
    //   0: (0,0,0), 1: (1,0,0), 2: (1,1,0), 3: (0,1,0)
    //   4: (0,0,1), 5: (1,0,1), 6: (1,1,1), 7: (0,1,1)

    // The first 8 vertices are the box corners, the next 12 vertices are
    // reserved for hmin-edge points and the final 12 vertices are reserved
    // for the hmax-edge points. The conservative upper bound of the number of
    // vertices is 8 + 12 + 12 = 32.
    private readonly mVertices: Vector[];

    // The box has 12 edges. An edge is mEdges[i] = [v0, v1], where the
    // indices v0 and v1 are relative to mVertices with v0 < v1.
    private readonly mEdges: [number, number][];

    // The box has 6 faces. The faces are listed as -x face, +x face, -y face,
    // +y face, -z face and +z face.
    private readonly mFaces: Face[];

    // Signed distances from the minimum and maximum height planes for the
    // cone to the projection of the box vertices onto the cone axis.
    private readonly mProjectionMin: number[];
    private readonly mProjectionMax: number[];

    // The mCandidateEdges array stores the edges of the clipped box that are
    // candidates for containing the optimizing point. The maximum number of
    // candidate edges is 1 + 2 + ... + 31 = 496, which is a conservative
    // bound because not all pairs of vertices form edges on box faces. The
    // candidate edges are stored as (v0,v1) with v0 < v1. The implementation
    // is designed so that during a single query, the number of candidate
    // edges can only grow.
    private mNumCandidateEdges: number;
    private readonly mCandidateEdges: [number, number][];

    // The mAdjacencyMatrix is a simple representation of edges in the graph
    // G = (V,E) that represents the (wireframe) clipped box. An edge (r,c)
    // does not exist when mAdjacencyMatrix[r][c] = 0. If an edge (r,c) does
    // exist, it is appended to mCandidateEdges at index k and the adjacency
    // matrix is set to mAdjacencyMatrix[r][c] = 1. This allows for a fast
    // edge-in-graph test and a fast and efficient clear of mCandidateEdges
    // and mAdjacencyMatrix.
    private readonly mAdjacencyMatrix: Uint8Array[];

    // The 81 possible configurations for a box face. The N stands for '-',
    // the Z stands for '0' and the P stands for '+'. These are listed in the
    // order that maps to the array mConfiguration. Thus, NNNN maps to
    // mConfiguration[0], NNNZ maps to mConfiguration[1], and so on.
    private readonly mConfiguration: ((base: number, face: Face) => void)[];

    constructor() {
        this.mVertices = [];
        for (let i = 0; i < MAX_VERTICES; ++i) {
            this.mVertices.push(Vector.zero(3));
        }

        // An edge is [v0, v1], where v0 and v1 are relative to mVertices with
        // v0 < v1.
        this.mEdges = [
            [0, 1], [1, 3], [2, 3], [0, 2],
            [4, 5], [5, 7], [6, 7], [4, 6],
            [0, 4], [1, 5], [3, 7], [2, 6]
        ];

        // For example, mFaces[0] has vertices { 0, 4, 6, 2 }. The edge {0,4}
        // is mEdges[8], the edge {4,6} is mEdges[7], the edge {6,2} is
        // mEdges[11] and the edge {2,0} is mEdges[3]; thus, the edge indices
        // are { 8, 7, 11, 3 }.
        this.mFaces = [
            { v: [0, 4, 6, 2], e: [8, 7, 11, 3] },
            { v: [1, 3, 7, 5], e: [1, 10, 5, 9] },
            { v: [0, 1, 5, 4], e: [0, 9, 4, 8] },
            { v: [2, 6, 7, 3], e: [11, 6, 10, 2] },
            { v: [0, 2, 3, 1], e: [3, 2, 1, 0] },
            { v: [4, 5, 7, 6], e: [4, 5, 6, 7] }
        ];

        this.mProjectionMin = new Array<number>(NUM_BOX_VERTICES).fill(0);
        this.mProjectionMax = new Array<number>(NUM_BOX_VERTICES).fill(0);

        // Clear the edges.
        this.mNumCandidateEdges = 0;
        this.mCandidateEdges = [];
        for (let r = 0; r < MAX_CANDIDATE_EDGES; ++r) {
            this.mCandidateEdges.push([0, 0]);
        }
        this.mAdjacencyMatrix = [];
        for (let r = 0; r < MAX_VERTICES; ++r) {
            this.mAdjacencyMatrix.push(new Uint8Array(MAX_VERTICES));
        }

        this.mConfiguration = [
            // NNNN_0
            () => { },
            // NNNZ_1
            () => { },
            // NNNP_2
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[2], base + face.e[3]);
            },
            // NNZN_3
            () => { },
            // NNZZ_4
            () => { },
            // NNZP_5
            (base: number, face: Face) => {
                this.insertEdge(face.v[2], base + face.e[3]);
            },
            // NNPN_6
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[1], base + face.e[2]);
            },
            // NNPZ_7
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[1], face.v[3]);
            },
            // NNPP_8
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[1], base + face.e[3]);
            },
            // NZNN_9
            () => { },
            // NZNZ_10
            () => { },
            // NZNP_11
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[2], face.v[3]);
                this.insertEdge(base + face.e[3], face.v[3]);
            },
            // NZZN_12
            () => { },
            // NZZZ_13
            () => { },
            // NZZP_14
            (base: number, face: Face) => {
                this.insertEdge(face.v[2], face.v[3]);
                this.insertEdge(base + face.e[3], face.v[3]);
            },
            // NZPN_15
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[2], face.v[1]);
            },
            // NZPZ_16
            (base: number, face: Face) => {
                this.insertEdge(face.v[1], face.v[3]);
            },
            // NZPP_17
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], face.v[1]);
            },
            // NPNN_18
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], base + face.e[1]);
            },
            // NPNZ_19
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], face.v[1]);
                this.insertEdge(base + face.e[1], face.v[1]);
            },
            // NPNP_20
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], face.v[1]);
                this.insertEdge(base + face.e[1], face.v[1]);
                this.insertEdge(base + face.e[2], face.v[3]);
                this.insertEdge(base + face.e[3], face.v[3]);
            },
            // NPZN_21
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], face.v[2]);
            },
            // NPZZ_22
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], face.v[1]);
                this.insertEdge(face.v[1], face.v[2]);
            },
            // NPZP_23
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], face.v[1]);
                this.insertEdge(face.v[1], face.v[2]);
                this.insertEdge(base + face.e[3], face.v[2]);
                this.insertEdge(face.v[2], face.v[3]);
            },
            // NPPN_24
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], base + face.e[2]);
            },
            // NPPZ_25
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], face.v[3]);
            },
            // NPPP_26
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], base + face.e[3]);
            },
            // ZNNN_27
            () => { },
            // ZNNZ_28
            () => { },
            // ZNNP_29
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[2], face.v[0]);
            },
            // ZNZN_30
            () => { },
            // ZNZZ_31
            () => { },
            // ZNZP_32
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], face.v[2]);
            },
            // ZNPN_33
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[1], face.v[2]);
                this.insertEdge(base + face.e[2], face.v[2]);
            },
            // ZNPZ_34
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[1], face.v[2]);
                this.insertEdge(face.v[2], face.v[3]);
            },
            // ZNPP_35
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], base + face.e[1]);
            },
            // ZZNN_36
            () => { },
            // ZZNZ_37
            () => { },
            // ZZNP_38
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], face.v[3]);
                this.insertEdge(face.v[3], base + face.e[2]);
            },
            // ZZZN_39
            () => { },
            // ZZZZ_40
            () => { },
            // ZZZP_41
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], face.v[3]);
                this.insertEdge(face.v[3], face.v[2]);
            },
            // ZZPN_42
            (base: number, face: Face) => {
                this.insertEdge(face.v[1], face.v[2]);
                this.insertEdge(face.v[2], base + face.e[2]);
            },
            // ZZPZ_43
            (base: number, face: Face) => {
                this.insertEdge(face.v[1], face.v[2]);
                this.insertEdge(face.v[2], face.v[3]);
            },
            // ZZPP_44
            () => { },
            // ZPNN_45
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], base + face.e[1]);
            },
            // ZPNZ_46
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], face.v[1]);
                this.insertEdge(face.v[1], base + face.e[1]);
            },
            // ZPNP_47
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], face.v[1]);
                this.insertEdge(face.v[1], base + face.e[1]);
                this.insertEdge(base + face.e[2], face.v[3]);
                this.insertEdge(face.v[3], face.v[0]);
            },
            // ZPZN_48
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], face.v[2]);
            },
            // ZPZZ_49
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], face.v[1]);
                this.insertEdge(face.v[1], face.v[2]);
            },
            // ZPZP_50
            () => { },
            // ZPPN_51
            (base: number, face: Face) => {
                this.insertEdge(face.v[0], base + face.e[2]);
            },
            // ZPPZ_52
            () => { },
            // ZPPP_53
            () => { },
            // PNNN_54
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], base + face.e[0]);
            },
            // PNNZ_55
            (base: number, face: Face) => {
                this.insertEdge(face.v[3], base + face.e[0]);
            },
            // PNNP_56
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[2], base + face.e[0]);
            },
            // PNZN_57
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], face.v[0]);
                this.insertEdge(face.v[0], base + face.e[0]);
            },
            // PNZZ_58
            (base: number, face: Face) => {
                this.insertEdge(face.v[3], face.v[0]);
                this.insertEdge(face.v[0], base + face.e[0]);
            },
            // PNZP_59
            (base: number, face: Face) => {
                this.insertEdge(face.v[2], base + face.e[0]);
            },
            // PNPN_60
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], face.v[0]);
                this.insertEdge(face.v[0], base + face.e[0]);
                this.insertEdge(base + face.e[1], face.v[2]);
                this.insertEdge(face.v[2], base + face.e[2]);
            },
            // PNPZ_61
            (base: number, face: Face) => {
                this.insertEdge(face.v[3], face.v[0]);
                this.insertEdge(face.v[0], base + face.e[0]);
                this.insertEdge(base + face.e[1], face.v[2]);
                this.insertEdge(face.v[2], face.v[3]);
            },
            // PNPP_62
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[0], base + face.e[1]);
            },
            // PZNN_63
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], face.v[1]);
            },
            // PZNZ_64
            (base: number, face: Face) => {
                this.insertEdge(face.v[3], face.v[1]);
            },
            // PZNP_65
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[2], face.v[1]);
            },
            // PZZN_66
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], face.v[0]);
                this.insertEdge(face.v[0], face.v[1]);
            },
            // PZZZ_67
            () => { },
            // PZZP_68
            () => { },
            // PZPN_69
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], face.v[0]);
                this.insertEdge(face.v[0], face.v[1]);
                this.insertEdge(face.v[1], face.v[2]);
                this.insertEdge(face.v[2], base + face.e[2]);
            },
            // PZPZ_70
            () => { },
            // PZPP_71
            () => { },
            // PPNN_72
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], base + face.e[1]);
            },
            // PPNZ_73
            (base: number, face: Face) => {
                this.insertEdge(face.v[3], base + face.e[1]);
            },
            // PPNP_74
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[2], base + face.e[1]);
            },
            // PPZN_75
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[2], face.v[2]);
            },
            // PPZZ_76
            () => { },
            // PPZP_77
            () => { },
            // PPPN_78
            (base: number, face: Face) => {
                this.insertEdge(base + face.e[3], base + face.e[2]);
            },
            // PPPZ_79
            () => { },
            // PPPP_80
            () => { }
        ];
    }

    test(box: AlignedBox, cone: Cone): IntrAlignedBox3Cone3TIResult {
        const result = defaultIntrAlignedBox3Cone3TIResult();

        // Quick-rejectance test. Determine whether the box is outside the
        // slab bounded by the minimum and maximum height planes. When outside
        // the slab, the box vertices are not required by the cone-box
        // intersection query, so the vertices are not yet computed.
        const { boxMinHeight, boxMaxHeight } =
            intrAlignedBox3Cone3ComputeBoxHeightInterval(box, cone);
        const coneMaxHeight = coneMaxHeightOf(cone);
        if (boxMaxHeight <= cone.getMinHeight() || boxMinHeight >= coneMaxHeight) {
            // There is no volumetric overlap of the box and the cone. The box
            // is clipped entirely.
            result.intersect = false;
            return result;
        }

        // Quick-acceptance test. Determine whether the cone axis intersects
        // the box.
        if (intrAlignedBox3Cone3ConeAxisIntersectsBox(box, cone)) {
            result.intersect = true;
            return result;
        }

        // Test for box fully inside the slab. When inside the slab, the box
        // vertices are required by the cone-box intersection query, so they
        // are computed here; they are also required in the remaining cases.
        // Also when inside the slab, the box edges are the candidates, so
        // they are copied to mCandidateEdges.
        if (this.boxFullyInConeSlab(box, boxMinHeight, boxMaxHeight, cone)) {
            result.intersect = this.candidatesHavePointInsideCone(cone);
            return result;
        }

        // Clear the candidates array and adjacency matrix.
        this.clearCandidates();

        // The box intersects at least one plane. Compute the box-plane
        // edge-interior intersection points. Insert the box subedges into the
        // array of candidate edges.
        this.computeCandidatesOnBoxEdges(cone);

        // Insert any relevant box face-interior clipped edges into the array
        // of candidate edges.
        this.computeCandidatesOnBoxFaces();

        result.intersect = this.candidatesHavePointInsideCone(cone);
        return result;
    }

    private boxFullyInConeSlab(box: AlignedBox, boxMinHeight: number,
        boxMaxHeight: number, cone: Cone): boolean {
        // Compute the box vertices relative to cone vertex as origin.
        const bMin = box.min.values;
        const bMax = box.max.values;
        const corners: [number, number, number][] = [
            [bMin[0], bMin[1], bMin[2]],
            [bMax[0], bMin[1], bMin[2]],
            [bMin[0], bMax[1], bMin[2]],
            [bMax[0], bMax[1], bMin[2]],
            [bMin[0], bMin[1], bMax[2]],
            [bMax[0], bMin[1], bMax[2]],
            [bMin[0], bMax[1], bMax[2]],
            [bMax[0], bMax[1], bMax[2]]
        ];
        const origin = cone.ray.origin.values;
        for (let i = 0; i < NUM_BOX_VERTICES; ++i) {
            const values = this.mVertices[i].values;
            values[0] = corners[i][0] - origin[0];
            values[1] = corners[i][1] - origin[1];
            values[2] = corners[i][2] - origin[2];
        }

        const coneMaxHeight = coneMaxHeightOf(cone);
        if (cone.getMinHeight() <= boxMinHeight && boxMaxHeight <= coneMaxHeight) {
            // The box is fully inside, so no clipping is necessary.
            //
            // Port fix for an upstream bug. Upstream copies mEdges over the
            // first 12 entries of mCandidateEdges and sets mNumCandidateEdges
            // to 12 without clearing the adjacency matrix. When the query
            // object is reused, the entries that a previous query set in
            // mAdjacencyMatrix beyond index 11 are never cleared by a later
            // ClearCandidates (which visits only the first
            // mNumCandidateEdges entries of mCandidateEdges). A subsequent
            // clipping query then silently drops those candidate edges,
            // because InsertEdge treats them as already present, and can
            // report a false negative. The port clears the candidates first
            // and inserts the box edges through insertEdge, which keeps
            // mCandidateEdges and mAdjacencyMatrix consistent. Within a
            // single query the candidate list is identical to upstream's.
            this.clearCandidates();
            for (let i = 0; i < NUM_BOX_EDGES; ++i) {
                this.insertEdge(this.mEdges[i][0], this.mEdges[i][1]);
            }
            return true;
        }
        return false;
    }

    private candidatesHavePointInsideCone(cone: Cone): boolean {
        for (let i = 0; i < this.mNumCandidateEdges; ++i) {
            const edge = this.mCandidateEdges[i];
            const P0 = this.mVertices[edge[0]];
            const P1 = this.mVertices[edge[1]];
            if (intrAlignedBox3Cone3HasPointInsideCone(P0, P1, cone)) {
                return true;
            }
        }
        return false;
    }

    private computeCandidatesOnBoxEdges(cone: Cone): void {
        const coneMaxHeight = coneMaxHeightOf(cone);
        for (let i = 0; i < NUM_BOX_VERTICES; ++i) {
            const h = dot(cone.ray.direction, this.mVertices[i]);
            this.mProjectionMin[i] = cone.getMinHeight() - h;
            this.mProjectionMax[i] = h - coneMaxHeight;
        }

        let v0 = VERTEX_MIN_BASE;
        let v1 = VERTEX_MAX_BASE;
        for (let i = 0; i < NUM_BOX_EDGES; ++i, ++v0, ++v1) {
            const edge = this.mEdges[i];

            // In the next blocks, the sign comparisons can be expressed
            // instead as "s0 * s1 < 0". The multiplication could lead to
            // floating-point underflow where the product becomes 0, so that
            // approach is avoided.

            // Process the hmin-plane.
            const p0Min = this.mProjectionMin[edge[0]];
            const p1Min = this.mProjectionMin[edge[1]];
            const clipMin = (p0Min < 0 && p1Min > 0) || (p0Min > 0 && p1Min < 0);
            if (clipMin) {
                const E0 = this.mVertices[edge[0]].values;
                const E1 = this.mVertices[edge[1]].values;
                const denom = p1Min - p0Min;
                const values = this.mVertices[v0].values;
                for (let d = 0; d < 3; ++d) {
                    values[d] = (p1Min * E0[d] - p0Min * E1[d]) / denom;
                }
            }

            // Process the hmax-plane.
            const p0Max = this.mProjectionMax[edge[0]];
            const p1Max = this.mProjectionMax[edge[1]];
            const clipMax = (p0Max < 0 && p1Max > 0) || (p0Max > 0 && p1Max < 0);
            if (clipMax) {
                const E0 = this.mVertices[edge[0]].values;
                const E1 = this.mVertices[edge[1]].values;
                const denom = p1Max - p0Max;
                const values = this.mVertices[v1].values;
                for (let d = 0; d < 3; ++d) {
                    values[d] = (p1Max * E0[d] - p0Max * E1[d]) / denom;
                }
            }

            // Get the candidate edges that are contained by the box edges.
            if (clipMin) {
                if (clipMax) {
                    this.insertEdge(v0, v1);
                }
                else {
                    if (p0Min < 0) {
                        this.insertEdge(edge[0], v0);
                    }
                    else {  // p1Min < 0
                        this.insertEdge(edge[1], v0);
                    }
                }
            }
            else if (clipMax) {
                if (p0Max < 0) {
                    this.insertEdge(edge[0], v1);
                }
                else {  // p1Max < 0
                    this.insertEdge(edge[1], v1);
                }
            }
            else {
                // No clipping has occurred. If the edge is inside the box, it
                // is a candidate edge. To be inside the box, the p*min and
                // p*max values must be nonpositive.
                if (p0Min <= 0 && p1Min <= 0 && p0Max <= 0 && p1Max <= 0) {
                    this.insertEdge(edge[0], edge[1]);
                }
            }
        }
    }

    private computeCandidatesOnBoxFaces(): void {
        for (let i = 0; i < NUM_BOX_FACES; ++i) {
            const face = this.mFaces[i];

            // Process the hmin-plane.
            let p0 = this.mProjectionMin[face.v[0]];
            let p1 = this.mProjectionMin[face.v[1]];
            let p2 = this.mProjectionMin[face.v[2]];
            let p3 = this.mProjectionMin[face.v[3]];
            let b0 = (p0 < 0 ? 0 : (p0 > 0 ? 2 : 1));
            let b1 = (p1 < 0 ? 0 : (p1 > 0 ? 2 : 1));
            let b2 = (p2 < 0 ? 0 : (p2 > 0 ? 2 : 1));
            let b3 = (p3 < 0 ? 0 : (p3 > 0 ? 2 : 1));
            let index = b3 + 3 * (b2 + 3 * (b1 + 3 * b0));
            this.mConfiguration[index](VERTEX_MIN_BASE, face);

            // Process the hmax-plane.
            p0 = this.mProjectionMax[face.v[0]];
            p1 = this.mProjectionMax[face.v[1]];
            p2 = this.mProjectionMax[face.v[2]];
            p3 = this.mProjectionMax[face.v[3]];
            b0 = (p0 < 0 ? 0 : (p0 > 0 ? 2 : 1));
            b1 = (p1 < 0 ? 0 : (p1 > 0 ? 2 : 1));
            b2 = (p2 < 0 ? 0 : (p2 > 0 ? 2 : 1));
            b3 = (p3 < 0 ? 0 : (p3 > 0 ? 2 : 1));
            index = b3 + 3 * (b2 + 3 * (b1 + 3 * b0));
            this.mConfiguration[index](VERTEX_MAX_BASE, face);
        }
    }

    private clearCandidates(): void {
        for (let i = 0; i < this.mNumCandidateEdges; ++i) {
            const edge = this.mCandidateEdges[i];
            this.mAdjacencyMatrix[edge[0]][edge[1]] = 0;
            this.mAdjacencyMatrix[edge[1]][edge[0]] = 0;
        }
        this.mNumCandidateEdges = 0;
    }

    private insertEdge(v0: number, v1: number): void {
        if (this.mAdjacencyMatrix[v0][v1] === 0) {
            this.mAdjacencyMatrix[v0][v1] = 1;
            this.mAdjacencyMatrix[v1][v0] = 1;
            this.mCandidateEdges[this.mNumCandidateEdges][0] = v0;
            this.mCandidateEdges[this.mNumCandidateEdges][1] = v1;
            ++this.mNumCandidateEdges;
        }
    }
}
