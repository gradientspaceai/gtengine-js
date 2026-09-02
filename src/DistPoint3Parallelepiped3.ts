// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint3Parallelepiped3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Implementation of a point-parallelepiped distance and closest-point query.
// The details are described in
//   https://www.geometrictools.com/Documentation/DistancePointParallelpiped.pdf
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector3<T>, Parallelepiped3<T>>'
// becomes the class DistPoint3Parallelepiped3 with the result type
// DistPoint3Parallelepiped3Result. The public member function 'GetMinimizer'
// stays a public method; the private region helpers become module-private
// functions.
//
// The six upstream face helpers GetClosestRzzm, GetClosestRzzp,
// GetClosestRzmz, GetClosestRzpz, GetClosestRmzz and GetClosestRpzz have
// identical bodies up to a cyclic permutation of the coordinate indices and
// the sign of the fixed coordinate, so they are expressed here by the single
// module-private faceMinimizer(A3, Z3, ia, ib, ic, value). The permutations
// (ia,ib,ic) used upstream are (0,1,2) for a fixed z, (2,0,1) for a fixed y
// and (1,2,0) for a fixed x; 'value' is the fixed coordinate, -1 or +1.

import type { DCPQuery } from './DCPQuery';
import { DistPoint2Parallelogram2 } from './DistPoint2Parallelogram2';
import { Matrix, multiplyATB, mulMatrix } from './Matrix';
import { inverse2x2 } from './Matrix2x2';
import { inverse3x3 } from './Matrix3x3';
import type { Parallelepiped3 } from './Parallelepiped3';
import { Vector, add, dot, mul, sub } from './Vector';

export interface DistPoint3Parallelepiped3Result {
    // The point closest[0] is the query point. The point closest[1] is the
    // parallelepiped point closest to the query point. The two points are the
    // same when the query point is contained by the parallelepiped.
    distance: number;
    sqrDistance: number;
    closest: [Vector, Vector];
}

// Minimize q on the face of the cube [-1,1]^3 on which coordinate ic is
// fixed at 'value'. The restriction of q to that face is a 2D quadratic
// function whose minimizer is computed by the point-parallelogram query.
function faceMinimizer(A3: Matrix, Z3: Vector, ia: number, ib: number,
    ic: number, value: number): Vector {
    const Z2 = Vector.fromArray([Z3.values[ia], Z3.values[ib]]);
    const u = value - Z3.values[ic];
    const A2 = Matrix.fromArray(2, 2, [
        A3.get(ia, ia), A3.get(ia, ib),
        A3.get(ia, ib), A3.get(ib, ib)
    ]);
    const V2 = Vector.fromArray([A3.get(ic, ia), A3.get(ic, ib)]);
    const Zeta2 = sub(Z2, mul(u, mulMatrix(inverse2x2(A2).inverse, V2)));
    const K2 = new DistPoint2Parallelogram2().getMinimizer(A2, Zeta2);
    const K3 = new Vector(3);
    K3.values[ia] = K2.values[0];
    K3.values[ib] = K2.values[1];
    K3.values[ic] = value;
    return K3;
}

// The six face regions.
function getClosestRzzm(A3: Matrix, Z3: Vector): Vector {
    return faceMinimizer(A3, Z3, 0, 1, 2, -1);
}

function getClosestRzzp(A3: Matrix, Z3: Vector): Vector {
    return faceMinimizer(A3, Z3, 0, 1, 2, +1);
}

function getClosestRzmz(A3: Matrix, Z3: Vector): Vector {
    return faceMinimizer(A3, Z3, 2, 0, 1, -1);
}

function getClosestRzpz(A3: Matrix, Z3: Vector): Vector {
    return faceMinimizer(A3, Z3, 2, 0, 1, +1);
}

function getClosestRmzz(A3: Matrix, Z3: Vector): Vector {
    return faceMinimizer(A3, Z3, 1, 2, 0, -1);
}

function getClosestRpzz(A3: Matrix, Z3: Vector): Vector {
    return faceMinimizer(A3, Z3, 1, 2, 0, +1);
}

// The twelve edge regions.
function getClosestRzmm(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzm(A, Z);
    if (K.values[1] === -1) {
        K = getClosestRzmz(A, Z);
    }
    return K;
}

function getClosestRzmp(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzp(A, Z);
    if (K.values[1] === -1) {
        K = getClosestRzmz(A, Z);
    }
    return K;
}

function getClosestRzpm(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzm(A, Z);
    if (K.values[1] === +1) {
        K = getClosestRzpz(A, Z);
    }
    return K;
}

function getClosestRzpp(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzp(A, Z);
    if (K.values[1] === +1) {
        K = getClosestRzpz(A, Z);
    }
    return K;
}

function getClosestRmzm(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzm(A, Z);
    if (K.values[0] === -1) {
        K = getClosestRmzz(A, Z);
    }
    return K;
}

function getClosestRmzp(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzp(A, Z);
    if (K.values[0] === -1) {
        K = getClosestRmzz(A, Z);
    }
    return K;
}

function getClosestRpzm(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzm(A, Z);
    if (K.values[0] === +1) {
        K = getClosestRpzz(A, Z);
    }
    return K;
}

function getClosestRpzp(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzp(A, Z);
    if (K.values[0] === +1) {
        K = getClosestRpzz(A, Z);
    }
    return K;
}

function getClosestRmmz(A: Matrix, Z: Vector): Vector {
    let K = getClosestRmzz(A, Z);
    if (K.values[1] === -1) {
        K = getClosestRzmz(A, Z);
    }
    return K;
}

function getClosestRmpz(A: Matrix, Z: Vector): Vector {
    let K = getClosestRmzz(A, Z);
    if (K.values[1] === +1) {
        K = getClosestRzpz(A, Z);
    }
    return K;
}

function getClosestRpmz(A: Matrix, Z: Vector): Vector {
    let K = getClosestRpzz(A, Z);
    if (K.values[1] === -1) {
        K = getClosestRzmz(A, Z);
    }
    return K;
}

function getClosestRppz(A: Matrix, Z: Vector): Vector {
    let K = getClosestRpzz(A, Z);
    if (K.values[1] === +1) {
        K = getClosestRzpz(A, Z);
    }
    return K;
}

// The eight vertex regions.
function getClosestRmmm(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzm(A, Z);  // K[2] = -1
    if (K.values[1] === -1) {
        K = getClosestRzmz(A, Z);  // K[1] = -1
        if (K.values[0] === -1) {
            K = getClosestRmzz(A, Z);  // K[0] = -1
        }
    }
    else if (K.values[0] === -1) {
        K = getClosestRmzz(A, Z);  // K[0] = -1
        if (K.values[1] === -1) {
            K = getClosestRzmz(A, Z);  // K[1] = -1
        }
    }
    return K;
}

function getClosestRmmp(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzp(A, Z);  // K[2] = +1
    if (K.values[1] === -1) {
        K = getClosestRzmz(A, Z);  // K[1] = -1
        if (K.values[0] === -1) {
            K = getClosestRmzz(A, Z);  // K[0] = -1
        }
    }
    else if (K.values[0] === -1) {
        K = getClosestRmzz(A, Z);  // K[0] = -1
        if (K.values[1] === -1) {
            K = getClosestRzmz(A, Z);  // K[1] = -1
        }
    }
    return K;
}

function getClosestRmpm(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzm(A, Z);  // K[2] = -1
    if (K.values[1] === +1) {
        K = getClosestRzpz(A, Z);  // K[1] = +1
        if (K.values[0] === -1) {
            K = getClosestRmzz(A, Z);  // K[0] = -1
        }
    }
    else if (K.values[0] === -1) {
        K = getClosestRmzz(A, Z);  // K[0] = -1
        if (K.values[1] === +1) {
            K = getClosestRzpz(A, Z);  // K[1] = +1
        }
    }
    return K;
}

function getClosestRmpp(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzp(A, Z);  // K[2] = +1
    if (K.values[1] === +1) {
        K = getClosestRzpz(A, Z);  // K[1] = +1
        if (K.values[0] === -1) {
            K = getClosestRmzz(A, Z);  // K[0] = -1
        }
    }
    else if (K.values[0] === -1) {
        K = getClosestRmzz(A, Z);  // K[0] = -1
        if (K.values[1] === +1) {
            K = getClosestRzpz(A, Z);  // K[1] = +1
        }
    }
    return K;
}

function getClosestRpmm(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzm(A, Z);  // K[2] = -1
    if (K.values[1] === -1) {
        K = getClosestRzmz(A, Z);  // K[1] = -1
        if (K.values[0] === +1) {
            K = getClosestRpzz(A, Z);  // K[0] = +1
        }
    }
    else if (K.values[0] === +1) {
        K = getClosestRpzz(A, Z);  // K[0] = +1
        if (K.values[1] === -1) {
            K = getClosestRzmz(A, Z);  // K[1] = -1
        }
    }
    return K;
}

function getClosestRpmp(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzp(A, Z);  // K[2] = +1
    if (K.values[1] === -1) {
        K = getClosestRzmz(A, Z);  // K[1] = -1
        if (K.values[0] === +1) {
            K = getClosestRpzz(A, Z);  // K[0] = +1
        }
    }
    else if (K.values[0] === +1) {
        K = getClosestRpzz(A, Z);  // K[0] = +1
        if (K.values[1] === -1) {
            K = getClosestRzmz(A, Z);  // K[1] = -1
        }
    }
    return K;
}

function getClosestRppm(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzm(A, Z);  // K[2] = -1
    if (K.values[1] === +1) {
        K = getClosestRzpz(A, Z);  // K[1] = +1
        if (K.values[0] === +1) {
            K = getClosestRpzz(A, Z);  // K[0] = +1
        }
    }
    else if (K.values[0] === +1) {
        K = getClosestRpzz(A, Z);  // K[0] = +1
        if (K.values[1] === +1) {
            K = getClosestRzpz(A, Z);  // K[1] = +1
        }
    }
    return K;
}

function getClosestRppp(A: Matrix, Z: Vector): Vector {
    let K = getClosestRzzp(A, Z);  // K[2] = +1
    if (K.values[1] === +1) {
        K = getClosestRzpz(A, Z);  // K[1] = +1
        if (K.values[0] === +1) {
            K = getClosestRpzz(A, Z);  // K[0] = +1
        }
    }
    else if (K.values[0] === +1) {
        K = getClosestRpzz(A, Z);  // K[0] = +1
        if (K.values[1] === +1) {
            K = getClosestRzpz(A, Z);  // K[1] = +1
        }
    }
    return K;
}

export class DistPoint3Parallelepiped3
    implements DCPQuery<Vector, Parallelepiped3,
    DistPoint3Parallelepiped3Result> {
    compute(point: Vector, ppd: Parallelepiped3):
        DistPoint3Parallelepiped3Result {
        // For a parallelepiped point X, let Y = {Dot(V0,X-C),Dot(V1,X-C),
        // Dot(V2,X-C)}. Compute the quadratic function q(Y) = (Y-Z)^T * A *
        // (Y-Z) / 2 where A = B^T * B is a symmetric matrix.
        const B = new Matrix(3, 3);
        B.setCol(0, ppd.axis[0]);
        B.setCol(1, ppd.axis[1]);
        B.setCol(2, ppd.axis[2]);
        const A = multiplyATB(B, B);

        // Transform the query point to parallelepiped coordinates,
        // Z = Inverse(B) * (P - C).
        const diff = sub(point, ppd.center);
        const Z = mulMatrix(inverse3x3(B).inverse, diff);

        // Get the minimizer for q(Y).
        const K = this.getMinimizer(A, Z);

        const closest0 = point.clone();
        const closest1 = add(ppd.center, mulMatrix(B, K));
        const delta = sub(closest0, closest1);
        const sqrDistance = dot(delta, delta);
        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            closest: [closest0, closest1]
        };
    }

    // Compute the minimizer (in parallelepiped coordinates) of the quadratic
    // function q. The domain is the cube [-1,1]^3.
    getMinimizer(A: Matrix, Z: Vector): Vector {
        const z0 = Z.values[0];
        const z1 = Z.values[1];
        const z2 = Z.values[2];
        let K: Vector;

        if (z2 < -1) {
            if (z1 < -1) {
                if (z0 < -1) {
                    K = getClosestRmmm(A, Z);
                }
                else if (z0 <= +1) {
                    K = getClosestRzmm(A, Z);
                }
                else {
                    K = getClosestRpmm(A, Z);
                }
            }
            else if (z1 <= +1) {
                if (z0 < -1) {
                    K = getClosestRmzm(A, Z);
                }
                else if (z0 <= +1) {
                    K = getClosestRzzm(A, Z);
                }
                else {
                    K = getClosestRpzm(A, Z);
                }
            }
            else {
                if (z0 < -1) {
                    K = getClosestRmpm(A, Z);
                }
                else if (z0 <= +1) {
                    K = getClosestRzpm(A, Z);
                }
                else {
                    K = getClosestRppm(A, Z);
                }
            }
        }
        else if (z2 <= +1) {
            if (z1 < -1) {
                if (z0 < -1) {
                    K = getClosestRmmz(A, Z);
                }
                else if (z0 <= +1) {
                    K = getClosestRzmz(A, Z);
                }
                else {
                    K = getClosestRpmz(A, Z);
                }
            }
            else if (z1 <= +1) {
                if (z0 < -1) {
                    K = getClosestRmzz(A, Z);
                }
                else if (z0 <= +1) {
                    // The query point is inside the parallelepiped; the
                    // minimizer is Z itself (upstream GetClosestRzzz).
                    K = Z.clone();
                }
                else {
                    K = getClosestRpzz(A, Z);
                }
            }
            else {
                if (z0 < -1) {
                    K = getClosestRmpz(A, Z);
                }
                else if (z0 <= +1) {
                    K = getClosestRzpz(A, Z);
                }
                else {
                    K = getClosestRppz(A, Z);
                }
            }
        }
        else {
            if (z1 < -1) {
                if (z0 < -1) {
                    K = getClosestRmmp(A, Z);
                }
                else if (z0 <= +1) {
                    K = getClosestRzmp(A, Z);
                }
                else {
                    K = getClosestRpmp(A, Z);
                }
            }
            else if (z1 <= +1) {
                if (z0 < -1) {
                    K = getClosestRmzp(A, Z);
                }
                else if (z0 <= +1) {
                    K = getClosestRzzp(A, Z);
                }
                else {
                    K = getClosestRpzp(A, Z);
                }
            }
            else {
                if (z0 < -1) {
                    K = getClosestRmpp(A, Z);
                }
                else if (z0 <= +1) {
                    K = getClosestRzpp(A, Z);
                }
                else {
                    K = getClosestRppp(A, Z);
                }
            }
        }

        return K;
    }
}
