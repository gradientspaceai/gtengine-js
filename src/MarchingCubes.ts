// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MarchingCubes.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Create the lookup table for the Marching Cubes algorithm that is used to
// extract a triangular mesh that represents a level surface of a 3D image
// sampled on a regular lattice. The assumption is that no data sample is
// zero, which allows us to have a table with 256 entries: 2 signs per
// sample, 8 samples per volume element (voxel). Each entry corresponds to
// the pattern of 8 signs at the corners of a voxel. The signs are stored as
// bits (b7,b6,b5,b4,b3,b2,b1,b0). The bit assignments to voxel corners is
//   b0 = (x,y,z),   b1 = (x+1,y,z),   b2 = (x,y+1,z),   b3 = (x+1,y+1,z)
//   b4 = (x,y,z+1), b5 = (x+1,y,z+1), b6 = (x,y+1,z+1), b7 = (x+1,y+1,z+1)
// If a bit is zero, then the voxel value at the corresponding corner is
// positive; otherwise, the bit is one and the value is negative. The
// triangles are counterclockwise ordered according to an observer viewing
// the triangle from the negative side of the level surface.
//
// Port notes: the template parameter IndexType is number. The nested struct
// Topology is exported as MarchingCubesTopology. The topology table is built
// per instance by the constructor (upstream builds a function-local static),
// and the upstream member-function pointers in the configuration table are
// replaced by a type-indexed dispatch (the function is uniquely determined by
// the configuration type). The overload GetTable() returning the flat
// 256 x 41 storage is ported as getFlatTable(). The configuration and
// prebuilt tables below are mechanical transcriptions of the upstream
// tables; do not edit them by hand.

// The topology for one voxel sign configuration.
// mTable[i][0] = numVertices
// mTable[i][1] = numTriangles
// mTable[i][2..25] = pairs of corner indices (maximum of 12 pairs)
// mTable[i][26..40] = triples of indices (maximum of 5 triples)
export interface MarchingCubesTopology {
    numVertices: number;
    numTriangles: number;
    // maxVertices (12) pairs of voxel corner indices; each vertex of the
    // level surface lies on the voxel edge joining the two corners.
    vpair: number[][];
    // maxTriangles (5) triples of indices into vpair.
    itriple: number[][];
}

// Configuration types, in the order of the upstream CTBits* constants (the
// value indexes CONFIGURATION_STRING).
const ConfigurationType = {
    Bits0: 0,
    Bits1: 1,
    Bits7: 2,
    Bits2Edge: 3,
    Bits6Edge: 4,
    Bits2FaceDiag: 5,
    Bits6FaceDiag: 6,
    Bits2BoxDiag: 7,
    Bits6BoxDiag: 8,
    Bits3SameFace: 9,
    Bits5SameFace: 10,
    Bits3EdgeFaceDiag: 11,
    Bits5EdgeFaceDiag: 12,
    Bits3FaceDiagFaceDiag: 13,
    Bits5FaceDiagFaceDiag: 14,
    Bits4SameFace: 15,
    Bits4FaceEdge: 16,
    Bits4FaceFaceDiagL: 17,
    Bits4FaceFaceDiagR: 18,
    Bits4FaceBoxDiag: 19,
    Bits4EdgeEdgePara: 20,
    Bits4EdgeEdgePerp: 21
} as const;

const CONFIGURATION_STRING: readonly string[] = [
    'Bits0',
    'Bits1',
    'Bits7',
    'Bits2Edge',
    'Bits6Edge',
    'Bits2FaceDiag',
    'Bits6FaceDiag',
    'Bits2BoxDiag',
    'Bits6BoxDiag',
    'Bits3SameFace',
    'Bits5SameFace',
    'Bits3EdgeFaceDiag',
    'Bits5EdgeFaceDiag',
    'Bits3FaceDiagFaceDiag',
    'Bits5FaceDiagFaceDiag',
    'Bits4SameFace',
    'Bits4FaceEdge',
    'Bits4FaceFaceDiagL',
    'Bits4FaceFaceDiagR',
    'Bits4FaceBoxDiag',
    'Bits4EdgeEdgePara',
    'Bits4EdgeEdgePerp'
];

interface Configuration {
    type: number;
    index: readonly number[];
}

// The precomputed information about the 256 sign configurations for voxels.
// Each entry selects a configuration function (by type) and the voxel corner
// permutation it is applied to. Transcribed verbatim from upstream; the
// comment on each row is the sign configuration (b7...b0).
const CONFIGURATION_TABLE: readonly Configuration[] = [
    /*00000000*/ { type: ConfigurationType.Bits0, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00000001*/ { type: ConfigurationType.Bits1, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00000010*/ { type: ConfigurationType.Bits1, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*00000011*/ { type: ConfigurationType.Bits2Edge, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00000100*/ { type: ConfigurationType.Bits1, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*00000101*/ { type: ConfigurationType.Bits2Edge, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*00000110*/ { type: ConfigurationType.Bits2FaceDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*00000111*/ { type: ConfigurationType.Bits3SameFace, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00001000*/ { type: ConfigurationType.Bits1, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*00001001*/ { type: ConfigurationType.Bits2FaceDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00001010*/ { type: ConfigurationType.Bits2Edge, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*00001011*/ { type: ConfigurationType.Bits3SameFace, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*00001100*/ { type: ConfigurationType.Bits2Edge, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*00001101*/ { type: ConfigurationType.Bits3SameFace, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*00001110*/ { type: ConfigurationType.Bits3SameFace, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*00001111*/ { type: ConfigurationType.Bits4SameFace, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00010000*/ { type: ConfigurationType.Bits1, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*00010001*/ { type: ConfigurationType.Bits2Edge, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*00010010*/ { type: ConfigurationType.Bits2FaceDiag, index: [1, 0, 5, 4, 3, 2, 7, 6] },
    /*00010011*/ { type: ConfigurationType.Bits3SameFace, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*00010100*/ { type: ConfigurationType.Bits2FaceDiag, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*00010101*/ { type: ConfigurationType.Bits3SameFace, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*00010110*/ { type: ConfigurationType.Bits3FaceDiagFaceDiag, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*00010111*/ { type: ConfigurationType.Bits4FaceEdge, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00011000*/ { type: ConfigurationType.Bits2BoxDiag, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*00011001*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00011010*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [1, 0, 5, 4, 3, 2, 7, 6] },
    /*00011011*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*00011100*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [2, 6, 0, 4, 3, 7, 1, 5] },
    /*00011101*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*00011110*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*00011111*/ { type: ConfigurationType.Bits5SameFace, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*00100000*/ { type: ConfigurationType.Bits1, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*00100001*/ { type: ConfigurationType.Bits2FaceDiag, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*00100010*/ { type: ConfigurationType.Bits2Edge, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*00100011*/ { type: ConfigurationType.Bits3SameFace, index: [1, 0, 5, 4, 3, 2, 7, 6] },
    /*00100100*/ { type: ConfigurationType.Bits2BoxDiag, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*00100101*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*00100110*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*00100111*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00101000*/ { type: ConfigurationType.Bits2FaceDiag, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*00101001*/ { type: ConfigurationType.Bits3FaceDiagFaceDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*00101010*/ { type: ConfigurationType.Bits3SameFace, index: [1, 5, 3, 7, 0, 4, 2, 6] },
    /*00101011*/ { type: ConfigurationType.Bits4FaceEdge, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*00101100*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [3, 1, 7, 5, 2, 0, 6, 4] },
    /*00101101*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*00101110*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*00101111*/ { type: ConfigurationType.Bits5SameFace, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*00110000*/ { type: ConfigurationType.Bits2Edge, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*00110001*/ { type: ConfigurationType.Bits3SameFace, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*00110010*/ { type: ConfigurationType.Bits3SameFace, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*00110011*/ { type: ConfigurationType.Bits4SameFace, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*00110100*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*00110101*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*00110110*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*00110111*/ { type: ConfigurationType.Bits5SameFace, index: [7, 6, 3, 2, 5, 4, 1, 0] },
    /*00111000*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*00111001*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*00111010*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*00111011*/ { type: ConfigurationType.Bits5SameFace, index: [6, 2, 7, 3, 4, 0, 5, 1] },
    /*00111100*/ { type: ConfigurationType.Bits4EdgeEdgePara, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*00111101*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [7, 3, 5, 1, 6, 2, 4, 0] },
    /*00111110*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*00111111*/ { type: ConfigurationType.Bits6Edge, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*01000000*/ { type: ConfigurationType.Bits1, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*01000001*/ { type: ConfigurationType.Bits2FaceDiag, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*01000010*/ { type: ConfigurationType.Bits2BoxDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*01000011*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*01000100*/ { type: ConfigurationType.Bits2Edge, index: [6, 2, 7, 3, 4, 0, 5, 1] },
    /*01000101*/ { type: ConfigurationType.Bits3SameFace, index: [2, 6, 0, 4, 3, 7, 1, 5] },
    /*01000110*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*01000111*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*01001000*/ { type: ConfigurationType.Bits2FaceDiag, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*01001001*/ { type: ConfigurationType.Bits3FaceDiagFaceDiag, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*01001010*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*01001011*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*01001100*/ { type: ConfigurationType.Bits3SameFace, index: [2, 3, 6, 7, 0, 1, 4, 5] },
    /*01001101*/ { type: ConfigurationType.Bits4FaceEdge, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*01001110*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*01001111*/ { type: ConfigurationType.Bits5SameFace, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*01010000*/ { type: ConfigurationType.Bits2Edge, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*01010001*/ { type: ConfigurationType.Bits3SameFace, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*01010010*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*01010011*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*01010100*/ { type: ConfigurationType.Bits3SameFace, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*01010101*/ { type: ConfigurationType.Bits4SameFace, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*01010110*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*01010111*/ { type: ConfigurationType.Bits5SameFace, index: [7, 3, 5, 1, 6, 2, 4, 0] },
    /*01011000*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [6, 2, 7, 3, 4, 0, 5, 1] },
    /*01011001*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*01011010*/ { type: ConfigurationType.Bits4EdgeEdgePara, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*01011011*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [7, 6, 3, 2, 5, 4, 1, 0] },
    /*01011100*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*01011101*/ { type: ConfigurationType.Bits5SameFace, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*01011110*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*01011111*/ { type: ConfigurationType.Bits6Edge, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*01100000*/ { type: ConfigurationType.Bits2FaceDiag, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*01100001*/ { type: ConfigurationType.Bits3FaceDiagFaceDiag, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*01100010*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*01100011*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [1, 0, 5, 4, 3, 2, 7, 6] },
    /*01100100*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*01100101*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [2, 6, 0, 4, 3, 7, 1, 5] },
    /*01100110*/ { type: ConfigurationType.Bits4EdgeEdgePara, index: [6, 2, 7, 3, 4, 0, 5, 1] },
    /*01100111*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*01101000*/ { type: ConfigurationType.Bits3FaceDiagFaceDiag, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*01101001*/ { type: ConfigurationType.Bits4EdgeEdgePerp, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*01101010*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [1, 5, 3, 7, 0, 4, 2, 6] },
    /*01101011*/ { type: ConfigurationType.Bits5FaceDiagFaceDiag, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*01101100*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [2, 3, 6, 7, 0, 1, 4, 5] },
    /*01101101*/ { type: ConfigurationType.Bits5FaceDiagFaceDiag, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*01101110*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*01101111*/ { type: ConfigurationType.Bits6FaceDiag, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*01110000*/ { type: ConfigurationType.Bits3SameFace, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*01110001*/ { type: ConfigurationType.Bits4FaceEdge, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*01110010*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*01110011*/ { type: ConfigurationType.Bits5SameFace, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*01110100*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*01110101*/ { type: ConfigurationType.Bits5SameFace, index: [3, 1, 7, 5, 2, 0, 6, 4] },
    /*01110110*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*01110111*/ { type: ConfigurationType.Bits6Edge, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*01111000*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*01111001*/ { type: ConfigurationType.Bits5FaceDiagFaceDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*01111010*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [2, 3, 6, 7, 0, 1, 4, 5] },
    /*01111011*/ { type: ConfigurationType.Bits6FaceDiag, index: [2, 3, 6, 7, 0, 1, 4, 5] },
    /*01111100*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [1, 5, 3, 7, 0, 4, 2, 6] },
    /*01111101*/ { type: ConfigurationType.Bits6FaceDiag, index: [1, 5, 3, 7, 0, 4, 2, 6] },
    /*01111110*/ { type: ConfigurationType.Bits6BoxDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*01111111*/ { type: ConfigurationType.Bits7, index: [7, 3, 5, 1, 6, 2, 4, 0] },
    /*10000000*/ { type: ConfigurationType.Bits1, index: [7, 3, 5, 1, 6, 2, 4, 0] },
    /*10000001*/ { type: ConfigurationType.Bits2BoxDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*10000010*/ { type: ConfigurationType.Bits2FaceDiag, index: [1, 5, 3, 7, 0, 4, 2, 6] },
    /*10000011*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [1, 5, 3, 7, 0, 4, 2, 6] },
    /*10000100*/ { type: ConfigurationType.Bits2FaceDiag, index: [2, 3, 6, 7, 0, 1, 4, 5] },
    /*10000101*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [2, 3, 6, 7, 0, 1, 4, 5] },
    /*10000110*/ { type: ConfigurationType.Bits3FaceDiagFaceDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*10000111*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*10001000*/ { type: ConfigurationType.Bits2Edge, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*10001001*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*10001010*/ { type: ConfigurationType.Bits3SameFace, index: [3, 1, 7, 5, 2, 0, 6, 4] },
    /*10001011*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*10001100*/ { type: ConfigurationType.Bits3SameFace, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*10001101*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*10001110*/ { type: ConfigurationType.Bits4FaceEdge, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*10001111*/ { type: ConfigurationType.Bits5SameFace, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*10010000*/ { type: ConfigurationType.Bits2FaceDiag, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*10010001*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*10010010*/ { type: ConfigurationType.Bits3FaceDiagFaceDiag, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*10010011*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*10010100*/ { type: ConfigurationType.Bits3FaceDiagFaceDiag, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*10010101*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*10010110*/ { type: ConfigurationType.Bits4EdgeEdgePerp, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*10010111*/ { type: ConfigurationType.Bits5FaceDiagFaceDiag, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*10011000*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*10011001*/ { type: ConfigurationType.Bits4EdgeEdgePara, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*10011010*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [3, 1, 7, 5, 2, 0, 6, 4] },
    /*10011011*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*10011100*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*10011101*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*10011110*/ { type: ConfigurationType.Bits5FaceDiagFaceDiag, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*10011111*/ { type: ConfigurationType.Bits6FaceDiag, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*10100000*/ { type: ConfigurationType.Bits2Edge, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*10100001*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*10100010*/ { type: ConfigurationType.Bits3SameFace, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*10100011*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [1, 0, 5, 4, 3, 2, 7, 6] },
    /*10100100*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [7, 6, 3, 2, 5, 4, 1, 0] },
    /*10100101*/ { type: ConfigurationType.Bits4EdgeEdgePara, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*10100110*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*10100111*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [6, 2, 7, 3, 4, 0, 5, 1] },
    /*10101000*/ { type: ConfigurationType.Bits3SameFace, index: [7, 3, 5, 1, 6, 2, 4, 0] },
    /*10101001*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [7, 3, 5, 1, 6, 2, 4, 0] },
    /*10101010*/ { type: ConfigurationType.Bits4SameFace, index: [1, 5, 3, 7, 0, 4, 2, 6] },
    /*10101011*/ { type: ConfigurationType.Bits5SameFace, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*10101100*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*10101101*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*10101110*/ { type: ConfigurationType.Bits5SameFace, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*10101111*/ { type: ConfigurationType.Bits6Edge, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*10110000*/ { type: ConfigurationType.Bits3SameFace, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*10110001*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*10110010*/ { type: ConfigurationType.Bits4FaceEdge, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*10110011*/ { type: ConfigurationType.Bits5SameFace, index: [2, 3, 6, 7, 0, 1, 4, 5] },
    /*10110100*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [5, 4, 7, 6, 1, 0, 3, 2] },
    /*10110101*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*10110110*/ { type: ConfigurationType.Bits5FaceDiagFaceDiag, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*10110111*/ { type: ConfigurationType.Bits6FaceDiag, index: [3, 7, 2, 6, 1, 5, 0, 4] },
    /*10111000*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [7, 3, 5, 1, 6, 2, 4, 0] },
    /*10111001*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*10111010*/ { type: ConfigurationType.Bits5SameFace, index: [2, 6, 0, 4, 3, 7, 1, 5] },
    /*10111011*/ { type: ConfigurationType.Bits6Edge, index: [6, 2, 7, 3, 4, 0, 5, 1] },
    /*10111100*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*10111101*/ { type: ConfigurationType.Bits6BoxDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*10111110*/ { type: ConfigurationType.Bits6FaceDiag, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*10111111*/ { type: ConfigurationType.Bits7, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*11000000*/ { type: ConfigurationType.Bits2Edge, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*11000001*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*11000010*/ { type: ConfigurationType.Bits3EdgeFaceDiag, index: [7, 3, 5, 1, 6, 2, 4, 0] },
    /*11000011*/ { type: ConfigurationType.Bits4EdgeEdgePara, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*11000100*/ { type: ConfigurationType.Bits3SameFace, index: [6, 2, 7, 3, 4, 0, 5, 1] },
    /*11000101*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [2, 6, 0, 4, 3, 7, 1, 5] },
    /*11000110*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [6, 2, 7, 3, 4, 0, 5, 1] },
    /*11000111*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*11001000*/ { type: ConfigurationType.Bits3SameFace, index: [7, 6, 3, 2, 5, 4, 1, 0] },
    /*11001001*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [7, 6, 3, 2, 5, 4, 1, 0] },
    /*11001010*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [7, 6, 3, 2, 5, 4, 1, 0] },
    /*11001011*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*11001100*/ { type: ConfigurationType.Bits4SameFace, index: [2, 3, 6, 7, 0, 1, 4, 5] },
    /*11001101*/ { type: ConfigurationType.Bits5SameFace, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*11001110*/ { type: ConfigurationType.Bits5SameFace, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*11001111*/ { type: ConfigurationType.Bits6Edge, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*11010000*/ { type: ConfigurationType.Bits3SameFace, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*11010001*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*11010010*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*11010011*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [3, 1, 7, 5, 2, 0, 6, 4] },
    /*11010100*/ { type: ConfigurationType.Bits4FaceEdge, index: [6, 4, 2, 0, 7, 5, 3, 1] },
    /*11010101*/ { type: ConfigurationType.Bits5SameFace, index: [1, 5, 3, 7, 0, 4, 2, 6] },
    /*11010110*/ { type: ConfigurationType.Bits5FaceDiagFaceDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*11010111*/ { type: ConfigurationType.Bits6FaceDiag, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*11011000*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [6, 7, 4, 5, 2, 3, 0, 1] },
    /*11011001*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*11011010*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*11011011*/ { type: ConfigurationType.Bits6BoxDiag, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*11011100*/ { type: ConfigurationType.Bits5SameFace, index: [1, 0, 5, 4, 3, 2, 7, 6] },
    /*11011101*/ { type: ConfigurationType.Bits6Edge, index: [5, 1, 4, 0, 7, 3, 6, 2] },
    /*11011110*/ { type: ConfigurationType.Bits6FaceDiag, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*11011111*/ { type: ConfigurationType.Bits7, index: [5, 7, 1, 3, 4, 6, 0, 2] },
    /*11100000*/ { type: ConfigurationType.Bits3SameFace, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*11100001*/ { type: ConfigurationType.Bits4FaceBoxDiag, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*11100010*/ { type: ConfigurationType.Bits4FaceFaceDiagL, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*11100011*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [2, 6, 0, 4, 3, 7, 1, 5] },
    /*11100100*/ { type: ConfigurationType.Bits4FaceFaceDiagR, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*11100101*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [1, 0, 5, 4, 3, 2, 7, 6] },
    /*11100110*/ { type: ConfigurationType.Bits5EdgeFaceDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*11100111*/ { type: ConfigurationType.Bits6BoxDiag, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*11101000*/ { type: ConfigurationType.Bits4FaceEdge, index: [7, 5, 6, 4, 3, 1, 2, 0] },
    /*11101001*/ { type: ConfigurationType.Bits5FaceDiagFaceDiag, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*11101010*/ { type: ConfigurationType.Bits5SameFace, index: [0, 2, 4, 6, 1, 3, 5, 7] },
    /*11101011*/ { type: ConfigurationType.Bits6FaceDiag, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*11101100*/ { type: ConfigurationType.Bits5SameFace, index: [0, 4, 1, 5, 2, 6, 3, 7] },
    /*11101101*/ { type: ConfigurationType.Bits6FaceDiag, index: [1, 0, 5, 4, 3, 2, 7, 6] },
    /*11101110*/ { type: ConfigurationType.Bits6Edge, index: [4, 0, 6, 2, 5, 1, 7, 3] },
    /*11101111*/ { type: ConfigurationType.Bits7, index: [4, 5, 0, 1, 6, 7, 2, 3] },
    /*11110000*/ { type: ConfigurationType.Bits4SameFace, index: [4, 6, 5, 7, 0, 2, 1, 3] },
    /*11110001*/ { type: ConfigurationType.Bits5SameFace, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*11110010*/ { type: ConfigurationType.Bits5SameFace, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*11110011*/ { type: ConfigurationType.Bits6Edge, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*11110100*/ { type: ConfigurationType.Bits5SameFace, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*11110101*/ { type: ConfigurationType.Bits6Edge, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*11110110*/ { type: ConfigurationType.Bits6FaceDiag, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*11110111*/ { type: ConfigurationType.Bits7, index: [3, 2, 1, 0, 7, 6, 5, 4] },
    /*11111000*/ { type: ConfigurationType.Bits5SameFace, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*11111001*/ { type: ConfigurationType.Bits6FaceDiag, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*11111010*/ { type: ConfigurationType.Bits6Edge, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*11111011*/ { type: ConfigurationType.Bits7, index: [2, 0, 3, 1, 6, 4, 7, 5] },
    /*11111100*/ { type: ConfigurationType.Bits6Edge, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*11111101*/ { type: ConfigurationType.Bits7, index: [1, 3, 0, 2, 5, 7, 4, 6] },
    /*11111110*/ { type: ConfigurationType.Bits7, index: [0, 1, 2, 3, 4, 5, 6, 7] },
    /*11111111*/ { type: ConfigurationType.Bits0, index: [0, 1, 2, 3, 4, 5, 6, 7] }
];

// The pre-built topology table, 256 rows of 41 integers. This is for
// reference in case you want to have a GPU-based implementation where you
// load the table as a GPU resource. Transcribed verbatim from upstream.
const PREBUILT_TABLE: readonly (readonly number[])[] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 0, 1, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 1, 3, 1, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 2, 0, 4, 0, 2, 1, 3, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 0, 2, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 2, 2, 6, 2, 3, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 1, 3, 1, 5, 0, 1, 0, 2, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 0, 4, 2, 6, 2, 3, 1, 3, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [3, 1, 2, 3, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 0, 1, 0, 4, 0, 2, 2, 3, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 2, 1, 5, 0, 1, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 1, 5, 0, 4, 0, 2, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [4, 2, 3, 7, 1, 3, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 2, 6, 3, 7, 1, 3, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [5, 3, 3, 7, 1, 5, 0, 1, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [4, 2, 0, 4, 2, 6, 3, 7, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 4, 5, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 2, 4, 5, 4, 6, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 0, 1, 1, 3, 1, 5, 4, 5, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 0, 2, 1, 3, 1, 5, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 2, 0, 4, 4, 5, 4, 6, 2, 6, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 0, 1, 4, 5, 4, 6, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [9, 3, 0, 2, 2, 6, 2, 3, 1, 3, 1, 5, 0, 1, 0, 4, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    [6, 4, 4, 5, 4, 6, 2, 6, 2, 3, 1, 3, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [6, 2, 2, 3, 3, 7, 1, 3, 0, 4, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 0, 1, 4, 5, 4, 6, 0, 2, 2, 3, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [7, 3, 0, 1, 2, 3, 3, 7, 1, 5, 4, 5, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [6, 4, 4, 5, 4, 6, 0, 2, 2, 3, 3, 7, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 3, 2, 6, 3, 7, 1, 3, 0, 2, 0, 4, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [6, 4, 4, 6, 2, 6, 3, 7, 1, 3, 0, 1, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [8, 4, 3, 7, 1, 5, 0, 1, 0, 2, 2, 6, 4, 5, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [5, 3, 3, 7, 2, 6, 4, 6, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [3, 1, 5, 7, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 0, 4, 0, 2, 0, 1, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 2, 5, 7, 4, 5, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 1, 3, 5, 7, 4, 5, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 2, 0, 2, 2, 6, 2, 3, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 0, 4, 2, 6, 2, 3, 0, 1, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [7, 3, 1, 3, 5, 7, 4, 5, 0, 1, 0, 2, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [6, 4, 4, 5, 0, 4, 2, 6, 2, 3, 1, 3, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [6, 2, 5, 7, 4, 5, 1, 5, 1, 3, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 3, 0, 1, 0, 4, 0, 2, 2, 3, 3, 7, 1, 3, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    [5, 3, 0, 1, 2, 3, 3, 7, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 5, 7, 4, 5, 0, 4, 0, 2, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 3, 1, 3, 0, 2, 2, 6, 3, 7, 5, 7, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 2, 6, 3, 7, 1, 3, 0, 1, 0, 4, 5, 7, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [6, 4, 5, 7, 4, 5, 0, 1, 0, 2, 2, 6, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 2, 6, 0, 4, 4, 5, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [4, 2, 4, 6, 0, 4, 1, 5, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 4, 6, 0, 2, 0, 1, 1, 5, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [5, 3, 5, 7, 4, 6, 0, 4, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [4, 2, 0, 2, 1, 3, 5, 7, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 0, 4, 1, 5, 5, 7, 4, 6, 2, 6, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [6, 4, 1, 5, 5, 7, 4, 6, 2, 6, 2, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [8, 4, 5, 7, 4, 6, 0, 4, 0, 1, 1, 3, 2, 6, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [5, 3, 5, 7, 1, 3, 2, 3, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [7, 3, 5, 7, 4, 6, 0, 4, 1, 5, 1, 3, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 4, 6, 0, 2, 0, 1, 1, 5, 5, 7, 2, 3, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [6, 4, 3, 7, 5, 7, 4, 6, 0, 4, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 4, 6, 5, 7, 3, 7, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [8, 4, 3, 7, 1, 3, 0, 2, 2, 6, 1, 5, 5, 7, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 0, 0],
    [7, 5, 3, 7, 2, 6, 4, 6, 5, 7, 1, 5, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [7, 5, 4, 6, 5, 7, 3, 7, 2, 6, 0, 2, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [4, 2, 2, 6, 4, 6, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 6, 7, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 0, 2, 0, 1, 0, 4, 4, 6, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 1, 3, 1, 5, 0, 1, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 0, 2, 1, 3, 1, 5, 0, 4, 4, 6, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [4, 2, 4, 6, 6, 7, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 2, 3, 0, 1, 0, 4, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [7, 3, 0, 2, 4, 6, 6, 7, 2, 3, 1, 3, 1, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [6, 4, 4, 6, 6, 7, 2, 3, 1, 3, 1, 5, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [6, 2, 3, 7, 1, 3, 2, 3, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 3, 2, 3, 3, 7, 1, 3, 0, 1, 0, 4, 0, 2, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    [7, 3, 3, 7, 1, 5, 0, 1, 2, 3, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 1, 5, 0, 4, 0, 2, 2, 3, 3, 7, 4, 6, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [5, 3, 0, 2, 4, 6, 6, 7, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 4, 6, 6, 7, 3, 7, 1, 3, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [6, 4, 6, 7, 3, 7, 1, 5, 0, 1, 0, 2, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 1, 5, 3, 7, 6, 7, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [4, 2, 6, 7, 2, 6, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 4, 5, 6, 7, 2, 6, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [7, 3, 4, 5, 6, 7, 2, 6, 0, 4, 0, 1, 1, 3, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [6, 4, 2, 6, 0, 2, 1, 3, 1, 5, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 6, 7, 2, 3, 0, 2, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [4, 2, 0, 1, 4, 5, 6, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 4, 6, 7, 2, 3, 0, 2, 0, 4, 4, 5, 1, 3, 1, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [5, 3, 6, 7, 4, 5, 1, 5, 1, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [7, 3, 2, 6, 0, 4, 4, 5, 6, 7, 3, 7, 1, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 4, 5, 6, 7, 2, 6, 0, 2, 0, 1, 3, 7, 1, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [8, 4, 1, 5, 0, 1, 2, 3, 3, 7, 0, 4, 4, 5, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 0, 0],
    [7, 5, 6, 7, 4, 5, 1, 5, 3, 7, 2, 3, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 4, 3, 7, 1, 3, 0, 2, 0, 4, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 4, 5, 0, 1, 1, 3, 3, 7, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [7, 5, 1, 5, 3, 7, 6, 7, 4, 5, 0, 4, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [4, 2, 4, 5, 1, 5, 3, 7, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 4, 5, 1, 5, 5, 7, 6, 7, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 3, 4, 5, 1, 5, 5, 7, 6, 7, 2, 6, 4, 6, 0, 4, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    [7, 3, 4, 5, 0, 1, 1, 3, 5, 7, 6, 7, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 1, 3, 5, 7, 4, 5, 0, 4, 0, 2, 6, 7, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [7, 3, 6, 7, 2, 3, 0, 2, 4, 6, 4, 5, 1, 5, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 2, 3, 0, 1, 0, 4, 4, 6, 6, 7, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [8, 4, 4, 6, 6, 7, 2, 3, 0, 2, 5, 7, 4, 5, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 0, 0],
    [7, 5, 5, 7, 1, 3, 2, 3, 6, 7, 4, 6, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [9, 3, 6, 7, 2, 6, 4, 6, 4, 5, 1, 5, 5, 7, 3, 7, 1, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    [12, 4, 0, 1, 0, 4, 0, 2, 2, 6, 4, 6, 6, 7, 2, 3, 3, 7, 1, 3, 1, 5, 5, 7, 4, 5, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 0, 0],
    [8, 4, 0, 1, 2, 3, 3, 7, 5, 7, 4, 5, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [9, 5, 4, 6, 0, 4, 4, 5, 5, 7, 3, 7, 6, 7, 2, 6, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 4, 3, 1, 7, 4, 1, 8, 7, 0, 5, 6],
    [8, 4, 0, 2, 4, 6, 6, 7, 3, 7, 1, 3, 4, 5, 1, 5, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [9, 5, 5, 7, 3, 7, 6, 7, 4, 6, 0, 4, 4, 5, 1, 5, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 4, 3, 1, 7, 4, 1, 8, 7, 0, 5, 6],
    [7, 5, 4, 6, 0, 2, 0, 1, 4, 5, 5, 7, 3, 7, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 4, 5, 7, 3, 7, 6, 7, 4, 6, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [5, 3, 0, 4, 1, 5, 5, 7, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 0, 2, 0, 1, 1, 5, 5, 7, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [6, 4, 6, 7, 2, 6, 0, 4, 0, 1, 1, 3, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 1, 3, 0, 2, 2, 6, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [6, 4, 0, 2, 0, 4, 1, 5, 5, 7, 6, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 2, 3, 6, 7, 5, 7, 1, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [7, 5, 2, 3, 6, 7, 5, 7, 1, 3, 0, 1, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [4, 2, 1, 3, 2, 3, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 4, 0, 4, 1, 5, 5, 7, 6, 7, 2, 6, 1, 3, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [9, 5, 1, 3, 1, 5, 0, 1, 0, 2, 2, 6, 2, 3, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 4, 3, 1, 7, 4, 1, 8, 7, 0, 5, 6],
    [7, 5, 2, 3, 0, 1, 0, 4, 2, 6, 6, 7, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 4, 2, 3, 0, 2, 2, 6, 6, 7, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [7, 5, 1, 5, 0, 4, 0, 2, 1, 3, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 4, 1, 5, 0, 1, 1, 3, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [6, 2, 0, 1, 0, 4, 0, 2, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 0, 1, 0, 4, 0, 2, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 2, 1, 5, 0, 1, 1, 3, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 1, 5, 0, 4, 0, 2, 1, 3, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [6, 2, 2, 3, 0, 2, 2, 6, 6, 7, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 2, 3, 0, 1, 0, 4, 2, 6, 6, 7, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [9, 3, 1, 3, 1, 5, 0, 1, 0, 2, 2, 6, 2, 3, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    [8, 4, 0, 4, 2, 6, 2, 3, 1, 3, 1, 5, 6, 7, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [4, 2, 1, 3, 2, 3, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 2, 3, 6, 7, 5, 7, 1, 3, 0, 1, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [5, 3, 2, 3, 6, 7, 5, 7, 1, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 5, 7, 1, 5, 0, 4, 0, 2, 2, 3, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 1, 3, 0, 2, 2, 6, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 6, 7, 5, 7, 1, 3, 0, 1, 0, 4, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [6, 4, 6, 7, 5, 7, 1, 5, 0, 1, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 0, 4, 1, 5, 5, 7, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [6, 2, 5, 7, 3, 7, 6, 7, 4, 6, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 4, 6, 0, 2, 0, 1, 4, 5, 5, 7, 3, 7, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [9, 3, 5, 7, 3, 7, 6, 7, 4, 6, 0, 4, 4, 5, 1, 5, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    [8, 4, 0, 2, 1, 3, 1, 5, 4, 5, 4, 6, 3, 7, 6, 7, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [9, 3, 4, 6, 0, 4, 4, 5, 5, 7, 3, 7, 6, 7, 2, 6, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0],
    [8, 4, 0, 1, 4, 5, 4, 6, 2, 6, 2, 3, 5, 7, 3, 7, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [12, 4, 1, 3, 1, 5, 0, 1, 0, 4, 4, 5, 4, 6, 0, 2, 2, 6, 2, 3, 3, 7, 6, 7, 5, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 0, 0],
    [9, 5, 6, 7, 2, 6, 4, 6, 4, 5, 1, 5, 5, 7, 3, 7, 1, 3, 2, 3, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 4, 3, 1, 7, 4, 1, 8, 7, 0, 5, 6],
    [7, 3, 5, 7, 1, 3, 2, 3, 6, 7, 4, 6, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 4, 5, 4, 6, 0, 2, 0, 1, 6, 7, 5, 7, 1, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 0, 0],
    [8, 4, 2, 3, 6, 7, 5, 7, 1, 5, 0, 1, 4, 6, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [7, 5, 6, 7, 2, 3, 0, 2, 4, 6, 4, 5, 1, 5, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [8, 4, 1, 3, 0, 2, 2, 6, 6, 7, 5, 7, 0, 4, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [7, 5, 4, 5, 0, 1, 1, 3, 5, 7, 6, 7, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [9, 5, 4, 5, 1, 5, 5, 7, 6, 7, 2, 6, 4, 6, 0, 4, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 4, 3, 1, 7, 4, 1, 8, 7, 0, 5, 6],
    [6, 4, 4, 5, 1, 5, 5, 7, 6, 7, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [4, 2, 4, 5, 1, 5, 3, 7, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 1, 5, 3, 7, 6, 7, 4, 5, 0, 4, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [5, 3, 4, 5, 0, 1, 1, 3, 3, 7, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 3, 7, 6, 7, 4, 5, 0, 4, 0, 2, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 3, 6, 7, 4, 5, 1, 5, 3, 7, 2, 3, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 2, 6, 2, 3, 0, 1, 0, 4, 3, 7, 6, 7, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 0, 0],
    [8, 4, 4, 5, 0, 1, 1, 3, 3, 7, 6, 7, 0, 2, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [7, 5, 2, 6, 0, 4, 4, 5, 6, 7, 3, 7, 1, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [5, 3, 6, 7, 4, 5, 1, 5, 1, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [8, 4, 6, 7, 4, 5, 1, 5, 1, 3, 2, 3, 0, 4, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [4, 2, 0, 1, 2, 3, 6, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 6, 7, 2, 3, 0, 2, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [6, 4, 1, 5, 1, 3, 0, 2, 2, 6, 6, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 5, 4, 5, 6, 7, 2, 6, 0, 4, 0, 1, 1, 3, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [5, 3, 4, 5, 6, 7, 2, 6, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [4, 2, 6, 7, 2, 6, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 1, 5, 3, 7, 6, 7, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 6, 7, 4, 6, 0, 2, 0, 1, 1, 5, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [6, 4, 3, 7, 6, 7, 4, 6, 0, 4, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 0, 2, 4, 6, 6, 7, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [8, 4, 1, 5, 3, 7, 6, 7, 4, 6, 0, 4, 2, 3, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [7, 5, 3, 7, 1, 5, 0, 1, 2, 3, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [9, 5, 2, 3, 3, 7, 1, 3, 0, 1, 0, 4, 0, 2, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 4, 3, 1, 7, 4, 1, 8, 7, 0, 5, 6],
    [6, 4, 3, 7, 1, 3, 2, 3, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [6, 4, 4, 6, 0, 4, 1, 5, 1, 3, 2, 3, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 5, 0, 2, 4, 6, 6, 7, 2, 3, 1, 3, 1, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [5, 3, 2, 3, 0, 1, 0, 4, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [4, 2, 4, 6, 6, 7, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 5, 0, 2, 1, 3, 1, 5, 0, 4, 4, 6, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 2, 1, 3, 1, 5, 0, 1, 2, 6, 4, 6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 4, 0, 2, 0, 1, 0, 4, 4, 6, 6, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [3, 1, 6, 7, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 2, 2, 6, 4, 6, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [7, 3, 4, 6, 5, 7, 3, 7, 2, 6, 0, 2, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [7, 3, 3, 7, 2, 6, 4, 6, 5, 7, 1, 5, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
    [8, 4, 0, 4, 0, 2, 1, 3, 1, 5, 2, 6, 4, 6, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 0, 0],
    [5, 3, 4, 6, 5, 7, 3, 7, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 3, 7, 2, 3, 0, 1, 0, 4, 4, 6, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [8, 4, 4, 6, 5, 7, 3, 7, 2, 3, 0, 2, 1, 5, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [7, 5, 5, 7, 4, 6, 0, 4, 1, 5, 1, 3, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [5, 3, 5, 7, 1, 3, 2, 3, 2, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [8, 4, 5, 7, 1, 3, 2, 3, 2, 6, 4, 6, 0, 1, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [6, 4, 1, 5, 0, 1, 2, 3, 2, 6, 4, 6, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 5, 0, 4, 1, 5, 5, 7, 4, 6, 2, 6, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [4, 2, 0, 2, 4, 6, 5, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 5, 7, 4, 6, 0, 4, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [5, 3, 4, 6, 0, 2, 0, 1, 1, 5, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [4, 2, 4, 6, 0, 4, 1, 5, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 2, 6, 0, 4, 4, 5, 5, 7, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [6, 4, 5, 7, 3, 7, 2, 6, 0, 2, 0, 1, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [8, 4, 2, 6, 0, 4, 4, 5, 5, 7, 3, 7, 0, 1, 1, 3, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [7, 5, 1, 3, 0, 2, 2, 6, 3, 7, 5, 7, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 4, 5, 7, 3, 7, 2, 3, 0, 2, 0, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [5, 3, 0, 1, 2, 3, 3, 7, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [9, 5, 0, 1, 0, 4, 0, 2, 2, 3, 3, 7, 1, 3, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 4, 3, 1, 7, 4, 1, 8, 7, 0, 5, 6],
    [6, 4, 5, 7, 4, 5, 1, 5, 1, 3, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [6, 4, 2, 3, 2, 6, 0, 4, 4, 5, 5, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 5, 1, 3, 5, 7, 4, 5, 0, 1, 0, 2, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [7, 5, 0, 4, 2, 6, 2, 3, 0, 1, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 2, 0, 2, 2, 6, 2, 3, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 1, 3, 5, 7, 4, 5, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [4, 2, 5, 7, 4, 5, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 4, 0, 4, 0, 2, 0, 1, 1, 5, 5, 7, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [3, 1, 5, 7, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 3, 7, 2, 6, 4, 6, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 0, 0, 0, 0, 0],
    [8, 4, 3, 7, 2, 6, 4, 6, 4, 5, 1, 5, 0, 2, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 5, 6, 7, 0, 0, 0],
    [6, 4, 1, 3, 3, 7, 2, 6, 4, 6, 4, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 5, 2, 6, 3, 7, 1, 3, 0, 2, 0, 4, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 4, 2, 3, 0, 2, 4, 6, 4, 5, 1, 5, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [7, 5, 0, 1, 2, 3, 3, 7, 1, 5, 4, 5, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [7, 5, 0, 1, 4, 5, 4, 6, 0, 2, 2, 3, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 6, 5, 1, 0, 5, 2, 1, 5, 3, 2, 5, 4, 3],
    [6, 2, 2, 3, 3, 7, 1, 3, 0, 4, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 4, 1, 3, 2, 3, 2, 6, 4, 6, 4, 5, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 0, 0],
    [9, 5, 0, 2, 2, 6, 2, 3, 1, 3, 1, 5, 0, 1, 0, 4, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 4, 3, 1, 7, 4, 1, 8, 7, 0, 5, 6],
    [5, 3, 0, 1, 4, 5, 4, 6, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [6, 4, 0, 4, 4, 5, 4, 6, 2, 6, 2, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [5, 3, 0, 2, 1, 3, 1, 5, 4, 5, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [6, 4, 0, 1, 1, 3, 1, 5, 4, 5, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [4, 2, 4, 5, 4, 6, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 4, 5, 4, 6, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 2, 0, 4, 1, 5, 3, 7, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 3, 7, 1, 5, 0, 1, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [5, 3, 2, 6, 3, 7, 1, 3, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [4, 2, 3, 7, 1, 3, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 1, 5, 0, 4, 0, 2, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [4, 2, 1, 5, 0, 1, 2, 3, 3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 4, 0, 1, 0, 4, 0, 2, 2, 3, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [3, 1, 2, 3, 3, 7, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 3, 0, 4, 2, 6, 2, 3, 1, 3, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 0, 0, 0, 0, 0],
    [6, 4, 1, 3, 1, 5, 0, 1, 0, 2, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 1, 5, 4, 1, 4, 3, 1, 3, 2, 0, 0, 0],
    [4, 2, 2, 6, 2, 3, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 0, 2, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 2, 0, 4, 0, 2, 1, 3, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 1, 3, 1, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 1, 0, 1, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
];

export class MarchingCubes {
    static readonly maxVertices = 12;
    static readonly maxTriangles = 5;

    // The 256-entry topology table generated by the constructor.
    protected mTable: MarchingCubesTopology[];

    // The constructor iterates mEntry from 0 to 255 and calls configuration
    // functions, each calling setTable(...). The mEntry value is the table
    // index to be used.
    protected mEntry: number;

    constructor() {
        // All members of each Topology are set to zeros.
        this.mTable = [];
        for (let i = 0; i < 256; ++i) {
            this.mTable.push({
                numVertices: 0,
                numTriangles: 0,
                vpair: Array.from({ length: MarchingCubes.maxVertices }, () => [0, 0]),
                itriple: Array.from({ length: MarchingCubes.maxTriangles }, () => [0, 0, 0])
            });
        }

        // Create the lookup table. The dispatch array is indexed by the
        // configuration type (upstream stores member-function pointers in
        // the configuration table).
        const functions: readonly ((index: readonly number[]) => void)[] = [
            (index) => this.bits0(index),
            (index) => this.bits1(index),
            (index) => this.bits7(index),
            (index) => this.bits2Edge(index),
            (index) => this.bits6Edge(index),
            (index) => this.bits2FaceDiag(index),
            (index) => this.bits6FaceDiag(index),
            (index) => this.bits2BoxDiag(index),
            (index) => this.bits6BoxDiag(index),
            (index) => this.bits3SameFace(index),
            (index) => this.bits5SameFace(index),
            (index) => this.bits3EdgeFaceDiag(index),
            (index) => this.bits5EdgeFaceDiag(index),
            (index) => this.bits3FaceDiagFaceDiag(index),
            (index) => this.bits5FaceDiagFaceDiag(index),
            (index) => this.bits4SameFace(index),
            (index) => this.bits4FaceEdge(index),
            (index) => this.bits4FaceFaceDiagL(index),
            (index) => this.bits4FaceFaceDiagR(index),
            (index) => this.bits4FaceBoxDiag(index),
            (index) => this.bits4EdgeEdgePara(index),
            (index) => this.bits4EdgeEdgePerp(index)
        ];

        for (this.mEntry = 0; this.mEntry < 256; ++this.mEntry) {
            const configuration = CONFIGURATION_TABLE[this.mEntry];
            functions[configuration.type](configuration.index);
        }
    }

    // The entry must be in {0..255}.
    getTable(entry: number): MarchingCubesTopology {
        return this.mTable[entry];
    }

    // The table has 256 entries, each 41 integers, returned flattened in
    // row-major order (upstream GetTable() overload returning the raw
    // pointer to the table storage).
    getFlatTable(): number[] {
        const flat: number[] = [];
        for (const topology of this.mTable) {
            flat.push(topology.numVertices, topology.numTriangles);
            for (const pair of topology.vpair) {
                flat.push(pair[0], pair[1]);
            }
            for (const triple of topology.itriple) {
                flat.push(triple[0], triple[1], triple[2]);
            }
        }
        return flat;
    }

    // The pre-built topology table that matches what the constructor
    // generates, stored as table[256][41].
    static getPrebuiltTable(): readonly (readonly number[])[] {
        return PREBUILT_TABLE;
    }

    // Get the configuration type for the voxel, which is one of the string
    // names of the 'bits*(index)' functions.
    static getConfigurationType(entry: number): string {
        if (0 <= entry && entry < 256) {
            return CONFIGURATION_STRING[CONFIGURATION_TABLE[entry].type];
        }
        return '';
    }

    protected setTable(numV: number, vpair: readonly number[],
        numT: number, itriple: readonly number[]): void {
        // The item is already zeroed in the constructor.
        const topology = this.mTable[this.mEntry];
        topology.numVertices = numV;
        topology.numTriangles = numT;

        // Store vertex pairs with minimum index occurring first.
        for (let i = 0, j = 0; i < numV; ++i, j += 2) {
            topology.vpair[i][0] = Math.min(vpair[j], vpair[j + 1]);
            topology.vpair[i][1] = Math.max(vpair[j], vpair[j + 1]);
        }

        // Store triangle triples as is.
        for (let i = 0, j = 0; i < numT; ++i, j += 3) {
            topology.itriple[i][0] = itriple[j];
            topology.itriple[i][1] = itriple[j + 1];
            topology.itriple[i][2] = itriple[j + 2];
        }
    }

    // The configuration functions. Each receives the permutation of voxel
    // corner indices from the configuration table.
    protected bits0(_index: readonly number[]): void {
        this.setTable(0, [], 0, []);
    }

    protected bits1(index: readonly number[]): void {
        const vpair = [
            index[1], index[0],
            index[4], index[0],
            index[2], index[0]
        ];
        const itriple = [
            0, 1, 2
        ];
        this.setTable(3, vpair, 1, itriple);
    }

    protected bits7(index: readonly number[]): void {
        const vpair = [
            index[1], index[0],
            index[4], index[0],
            index[2], index[0]
        ];
        const itriple = [
            0, 2, 1
        ];
        this.setTable(3, vpair, 1, itriple);
    }

    protected bits2Edge(index: readonly number[]): void {
        const vpair = [
            index[4], index[0],
            index[2], index[0],
            index[3], index[1],
            index[5], index[1]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3
        ];
        this.setTable(4, vpair, 2, itriple);
    }

    protected bits6Edge(index: readonly number[]): void {
        const vpair = [
            index[4], index[0],
            index[2], index[0],
            index[3], index[1],
            index[5], index[1]
        ];
        const itriple = [
            0, 2, 1,
            0, 3, 2
        ];
        this.setTable(4, vpair, 2, itriple);
    }

    protected bits2FaceDiag(index: readonly number[]): void {
        const vpair = [
            index[1], index[0],
            index[4], index[0],
            index[2], index[0],
            index[2], index[3],
            index[7], index[3],
            index[1], index[3]
        ];
        const itriple = [
            0, 1, 2,
            3, 4, 5
        ];
        this.setTable(6, vpair, 2, itriple);
    }

    protected bits6FaceDiag(index: readonly number[]): void {
        const vpair = [
            index[1], index[0],
            index[4], index[0],
            index[2], index[0],
            index[2], index[3],
            index[7], index[3],
            index[1], index[3]
        ];
        // Not the reverse ordering from bits2FaceDiag due to ambiguous face
        // handling.
        const itriple = [
            1, 0, 5,
            1, 5, 4,
            1, 4, 3,
            1, 3, 2
        ];
        this.setTable(6, vpair, 4, itriple);
    }

    protected bits2BoxDiag(index: readonly number[]): void {
        const vpair = [
            index[1], index[0],
            index[4], index[0],
            index[2], index[0],
            index[3], index[7],
            index[6], index[7],
            index[5], index[7]
        ];
        const itriple = [
            0, 1, 2,
            3, 4, 5
        ];
        this.setTable(6, vpair, 2, itriple);
    }

    protected bits6BoxDiag(index: readonly number[]): void {
        const vpair = [
            index[1], index[0],
            index[4], index[0],
            index[2], index[0],
            index[3], index[7],
            index[6], index[7],
            index[5], index[7]
        ];
        const itriple = [
            0, 2, 1,
            3, 5, 4
        ];
        this.setTable(6, vpair, 2, itriple);
    }

    protected bits3SameFace(index: readonly number[]): void {
        const vpair = [
            index[4], index[0],
            index[2], index[6],
            index[2], index[3],
            index[1], index[3],
            index[1], index[5]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3,
            0, 3, 4
        ];
        this.setTable(5, vpair, 3, itriple);
    }

    protected bits5SameFace(index: readonly number[]): void {
        const vpair = [
            index[4], index[0],
            index[2], index[6],
            index[2], index[3],
            index[1], index[3],
            index[1], index[5]
        ];
        const itriple = [
            0, 2, 1,
            0, 3, 2,
            0, 4, 3
        ];
        this.setTable(5, vpair, 3, itriple);
    }

    protected bits3EdgeFaceDiag(index: readonly number[]): void {
        const vpair = [
            index[0], index[1],
            index[4], index[5],
            index[4], index[6],
            index[0], index[2],
            index[2], index[3],
            index[3], index[7],
            index[1], index[3]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3,
            4, 5, 6
        ];
        this.setTable(7, vpair, 3, itriple);
    }

    protected bits5EdgeFaceDiag(index: readonly number[]): void {
        const vpair = [
            index[0], index[1],
            index[4], index[5],
            index[4], index[6],
            index[0], index[2],
            index[2], index[3],
            index[3], index[7],
            index[1], index[3]
        ];
        // Not the reverse ordering from bits3EdgeFaceDiag due to ambiguous
        // face handling.
        const itriple = [
            5, 0, 6,
            5, 1, 0,
            5, 2, 1,
            5, 3, 2,
            5, 4, 3
        ];
        this.setTable(7, vpair, 5, itriple);
    }

    protected bits3FaceDiagFaceDiag(index: readonly number[]): void {
        const vpair = [
            index[0], index[1],
            index[0], index[4],
            index[0], index[2],
            index[2], index[3],
            index[3], index[7],
            index[1], index[3],
            index[1], index[5],
            index[5], index[7],
            index[4], index[5]
        ];
        const itriple = [
            0, 1, 2,
            3, 4, 5,
            6, 7, 8
        ];
        this.setTable(9, vpair, 3, itriple);
    }

    protected bits5FaceDiagFaceDiag(index: readonly number[]): void {
        const vpair = [
            index[0], index[1],
            index[0], index[4],
            index[0], index[2],
            index[2], index[3],
            index[3], index[7],
            index[1], index[3],
            index[1], index[5],
            index[5], index[7],
            index[4], index[5]
        ];
        // Not the reverse ordering from bits3FaceDiagFaceDiag due to
        // ambiguous face handling.
        const itriple = [
            1, 3, 2,
            1, 4, 3,
            1, 7, 4,
            1, 8, 7,
            0, 5, 6
        ];
        this.setTable(9, vpair, 5, itriple);
    }

    protected bits4SameFace(index: readonly number[]): void {
        const vpair = [
            index[0], index[4],
            index[2], index[6],
            index[3], index[7],
            index[1], index[5]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3
        ];
        this.setTable(4, vpair, 2, itriple);
    }

    protected bits4FaceEdge(index: readonly number[]): void {
        const vpair = [
            index[4], index[5],
            index[4], index[6],
            index[2], index[6],
            index[2], index[3],
            index[1], index[3],
            index[1], index[5]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3,
            0, 3, 4,
            0, 4, 5
        ];
        this.setTable(6, vpair, 4, itriple);
    }

    protected bits4FaceFaceDiagL(index: readonly number[]): void {
        const vpair = [
            index[4], index[5],
            index[0], index[4],
            index[2], index[6],
            index[2], index[3],
            index[1], index[3],
            index[5], index[7]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3,
            0, 3, 4,
            0, 4, 5
        ];
        this.setTable(6, vpair, 4, itriple);
    }

    protected bits4FaceFaceDiagR(index: readonly number[]): void {
        const vpair = [
            index[4], index[6],
            index[6], index[7],
            index[2], index[3],
            index[1], index[3],
            index[1], index[5],
            index[0], index[4]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3,
            0, 3, 4,
            0, 4, 5
        ];
        this.setTable(6, vpair, 4, itriple);
    }

    protected bits4FaceBoxDiag(index: readonly number[]): void {
        const vpair = [
            index[0], index[4],
            index[2], index[6],
            index[2], index[3],
            index[1], index[3],
            index[1], index[5],
            index[6], index[7],
            index[5], index[7],
            index[3], index[7]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3,
            0, 3, 4,
            5, 6, 7
        ];
        this.setTable(8, vpair, 4, itriple);
    }

    protected bits4EdgeEdgePara(index: readonly number[]): void {
        const vpair = [
            index[0], index[4],
            index[0], index[2],
            index[1], index[3],
            index[1], index[5],
            index[2], index[6],
            index[4], index[6],
            index[5], index[7],
            index[3], index[7]
        ];
        const itriple = [
            0, 1, 2,
            0, 2, 3,
            4, 5, 6,
            4, 6, 7
        ];
        this.setTable(8, vpair, 4, itriple);
    }

    protected bits4EdgeEdgePerp(index: readonly number[]): void {
        const vpair = [
            index[0], index[1],
            index[0], index[4],
            index[0], index[2],
            index[2], index[6],
            index[4], index[6],
            index[6], index[7],
            index[2], index[3],
            index[3], index[7],
            index[1], index[3],
            index[1], index[5],
            index[5], index[7],
            index[4], index[5]
        ];
        const itriple = [
            0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11
        ];
        this.setTable(12, vpair, 4, itriple);
    }
}
