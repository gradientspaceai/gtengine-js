// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox3Frustum3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection query uses the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The potential separating axes include the 3 box face normals, the 5
// distinct frustum normals (the near and far planes have the same normal),
// and cross products of normals, one from the box and one from the frustum.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only a TIQuery specialization for this pair of primitives, which becomes
// IntrOrientedBox3Frustum3TI.

import type { Frustum3 } from './Frustum3';
import type { OrientedBox } from './OrientedBox';
import type { TIQuery } from './TIQuery';
import { dot, sub } from './Vector';

// The result of IntrOrientedBox3Frustum3TI queries.
export interface IntrOrientedBox3Frustum3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrOrientedBox3Frustum3TIResult {
    return { intersect: false };
}

// Test-intersection query for a solid oriented box and a solid frustum in 3D.
export class IntrOrientedBox3Frustum3TI implements
    TIQuery<OrientedBox, Frustum3, IntrOrientedBox3Frustum3TIResult> {

    test(box: OrientedBox, frustum: Frustum3): IntrOrientedBox3Frustum3TIResult {
        const result = defaultTIResult();

        // Convenience variables.
        const axis = box.axis;
        const extent = box.extent.values;

        const diff = sub(box.center, frustum.origin);  // C-E

        const A = [0, 0, 0];      // Dot(R,A[i])
        const B = [0, 0, 0];      // Dot(U,A[i])
        const C = [0, 0, 0];      // Dot(D,A[i])
        const D = [0, 0, 0];      // (Dot(R,C-E),Dot(U,C-E),Dot(D,C-E))
        const NA = [0, 0, 0];     // dmin*Dot(R,A[i])
        const NB = [0, 0, 0];     // dmin*Dot(U,A[i])
        const NC = [0, 0, 0];     // dmin*Dot(D,A[i])
        const ND = [0, 0, 0];     // dmin*(Dot(R,C-E),Dot(U,C-E),?)
        const RC = [0, 0, 0];     // rmax*Dot(D,A[i])
        const RD = [0, 0, 0];     // rmax*(?,?,Dot(D,C-E))
        const UC = [0, 0, 0];     // umax*Dot(D,A[i])
        const UD = [0, 0, 0];     // umax*(?,?,Dot(D,C-E))
        const NApRC = [0, 0, 0];  // dmin*Dot(R,A[i]) + rmax*Dot(D,A[i])
        const NAmRC = [0, 0, 0];  // dmin*Dot(R,A[i]) - rmax*Dot(D,A[i])
        const NBpUC = [0, 0, 0];  // dmin*Dot(U,A[i]) + umax*Dot(D,A[i])
        const NBmUC = [0, 0, 0];  // dmin*Dot(U,A[i]) - umax*Dot(D,A[i])
        const RBpUA = [0, 0, 0];  // rmax*Dot(U,A[i]) + umax*Dot(R,A[i])
        const RBmUA = [0, 0, 0];  // rmax*Dot(U,A[i]) - umax*Dot(R,A[i])
        let DdD: number, radius: number, p: number;
        let fmin: number, fmax: number, tmp: number;
        let i: number, j: number;

        // M = D
        D[2] = dot(diff, frustum.dVector);
        for (i = 0; i < 3; ++i) {
            C[i] = dot(axis[i], frustum.dVector);
        }
        radius =
            extent[0] * Math.abs(C[0]) +
            extent[1] * Math.abs(C[1]) +
            extent[2] * Math.abs(C[2]);
        if (D[2] + radius < frustum.dMin || D[2] - radius > frustum.dMax) {
            result.intersect = false;
            return result;
        }

        // M = n*R - r*D
        for (i = 0; i < 3; ++i) {
            A[i] = dot(axis[i], frustum.rVector);
            RC[i] = frustum.rBound * C[i];
            NA[i] = frustum.dMin * A[i];
            NAmRC[i] = NA[i] - RC[i];
        }
        D[0] = dot(diff, frustum.rVector);
        radius =
            extent[0] * Math.abs(NAmRC[0]) +
            extent[1] * Math.abs(NAmRC[1]) +
            extent[2] * Math.abs(NAmRC[2]);
        ND[0] = frustum.dMin * D[0];
        RD[2] = frustum.rBound * D[2];
        DdD = ND[0] - RD[2];
        const MTwoRF = frustum.getMTwoRF();
        if (DdD + radius < MTwoRF || DdD > radius) {
            result.intersect = false;
            return result;
        }

        // M = -n*R - r*D
        for (i = 0; i < 3; ++i) {
            NApRC[i] = NA[i] + RC[i];
        }
        radius =
            extent[0] * Math.abs(NApRC[0]) +
            extent[1] * Math.abs(NApRC[1]) +
            extent[2] * Math.abs(NApRC[2]);
        DdD = -(ND[0] + RD[2]);
        if (DdD + radius < MTwoRF || DdD > radius) {
            result.intersect = false;
            return result;
        }

        // M = n*U - u*D
        for (i = 0; i < 3; ++i) {
            B[i] = dot(axis[i], frustum.uVector);
            UC[i] = frustum.uBound * C[i];
            NB[i] = frustum.dMin * B[i];
            NBmUC[i] = NB[i] - UC[i];
        }
        D[1] = dot(diff, frustum.uVector);
        radius =
            extent[0] * Math.abs(NBmUC[0]) +
            extent[1] * Math.abs(NBmUC[1]) +
            extent[2] * Math.abs(NBmUC[2]);
        ND[1] = frustum.dMin * D[1];
        UD[2] = frustum.uBound * D[2];
        DdD = ND[1] - UD[2];
        const MTwoUF = frustum.getMTwoUF();
        if (DdD + radius < MTwoUF || DdD > radius) {
            result.intersect = false;
            return result;
        }

        // M = -n*U - u*D
        for (i = 0; i < 3; ++i) {
            NBpUC[i] = NB[i] + UC[i];
        }
        radius =
            extent[0] * Math.abs(NBpUC[0]) +
            extent[1] * Math.abs(NBpUC[1]) +
            extent[2] * Math.abs(NBpUC[2]);
        DdD = -(ND[1] + UD[2]);
        if (DdD + radius < MTwoUF || DdD > radius) {
            result.intersect = false;
            return result;
        }

        const dRatio = frustum.getDRatio();

        // M = A[i]
        for (i = 0; i < 3; ++i) {
            p = frustum.rBound * Math.abs(A[i]) +
                frustum.uBound * Math.abs(B[i]);
            NC[i] = frustum.dMin * C[i];
            fmin = NC[i] - p;
            if (fmin < 0) {
                fmin *= dRatio;
            }
            fmax = NC[i] + p;
            if (fmax > 0) {
                fmax *= dRatio;
            }
            DdD = A[i] * D[0] + B[i] * D[1] + C[i] * D[2];
            if (DdD + extent[i] < fmin || DdD - extent[i] > fmax) {
                result.intersect = false;
                return result;
            }
        }

        // M = Cross(R,A[i])
        for (i = 0; i < 3; ++i) {
            p = frustum.uBound * Math.abs(C[i]);
            fmin = -NB[i] - p;
            if (fmin < 0) {
                fmin *= dRatio;
            }
            fmax = -NB[i] + p;
            if (fmax > 0) {
                fmax *= dRatio;
            }
            DdD = C[i] * D[1] - B[i] * D[2];
            radius =
                extent[0] * Math.abs(B[i] * C[0] - B[0] * C[i]) +
                extent[1] * Math.abs(B[i] * C[1] - B[1] * C[i]) +
                extent[2] * Math.abs(B[i] * C[2] - B[2] * C[i]);
            if (DdD + radius < fmin || DdD - radius > fmax) {
                result.intersect = false;
                return result;
            }
        }

        // M = Cross(U,A[i])
        for (i = 0; i < 3; ++i) {
            p = frustum.rBound * Math.abs(C[i]);
            fmin = NA[i] - p;
            if (fmin < 0) {
                fmin *= dRatio;
            }
            fmax = NA[i] + p;
            if (fmax > 0) {
                fmax *= dRatio;
            }
            DdD = -C[i] * D[0] + A[i] * D[2];
            radius =
                extent[0] * Math.abs(A[i] * C[0] - A[0] * C[i]) +
                extent[1] * Math.abs(A[i] * C[1] - A[1] * C[i]) +
                extent[2] * Math.abs(A[i] * C[2] - A[2] * C[i]);
            if (DdD + radius < fmin || DdD - radius > fmax) {
                result.intersect = false;
                return result;
            }
        }

        // M = Cross(n*D+r*R+u*U,A[i])
        for (i = 0; i < 3; ++i) {
            const fRB = frustum.rBound * B[i];
            const fUA = frustum.uBound * A[i];
            RBpUA[i] = fRB + fUA;
            RBmUA[i] = fRB - fUA;
        }
        for (i = 0; i < 3; ++i) {
            p = frustum.rBound * Math.abs(NBmUC[i]) +
                frustum.uBound * Math.abs(NAmRC[i]);
            tmp = -frustum.dMin * RBmUA[i];
            fmin = tmp - p;
            if (fmin < 0) {
                fmin *= dRatio;
            }
            fmax = tmp + p;
            if (fmax > 0) {
                fmax *= dRatio;
            }
            DdD = D[0] * NBmUC[i] - D[1] * NAmRC[i] - D[2] * RBmUA[i];
            radius = 0;
            for (j = 0; j < 3; ++j) {
                radius += extent[j] * Math.abs(A[j] * NBmUC[i] -
                    B[j] * NAmRC[i] - C[j] * RBmUA[i]);
            }
            if (DdD + radius < fmin || DdD - radius > fmax) {
                result.intersect = false;
                return result;
            }
        }

        // M = Cross(n*D+r*R-u*U,A[i])
        for (i = 0; i < 3; ++i) {
            p = frustum.rBound * Math.abs(NBpUC[i]) +
                frustum.uBound * Math.abs(NAmRC[i]);
            tmp = -frustum.dMin * RBpUA[i];
            fmin = tmp - p;
            if (fmin < 0) {
                fmin *= dRatio;
            }
            fmax = tmp + p;
            if (fmax > 0) {
                fmax *= dRatio;
            }
            DdD = D[0] * NBpUC[i] - D[1] * NAmRC[i] - D[2] * RBpUA[i];
            radius = 0;
            for (j = 0; j < 3; ++j) {
                radius += extent[j] * Math.abs(A[j] * NBpUC[i] -
                    B[j] * NAmRC[i] - C[j] * RBpUA[i]);
            }
            if (DdD + radius < fmin || DdD - radius > fmax) {
                result.intersect = false;
                return result;
            }
        }

        // M = Cross(n*D-r*R+u*U,A[i])
        for (i = 0; i < 3; ++i) {
            p = frustum.rBound * Math.abs(NBmUC[i]) +
                frustum.uBound * Math.abs(NApRC[i]);
            tmp = frustum.dMin * RBpUA[i];
            fmin = tmp - p;
            if (fmin < 0) {
                fmin *= dRatio;
            }
            fmax = tmp + p;
            if (fmax > 0) {
                fmax *= dRatio;
            }
            DdD = D[0] * NBmUC[i] - D[1] * NApRC[i] + D[2] * RBpUA[i];
            radius = 0;
            for (j = 0; j < 3; ++j) {
                radius += extent[j] * Math.abs(A[j] * NBmUC[i] -
                    B[j] * NApRC[i] + C[j] * RBpUA[i]);
            }
            if (DdD + radius < fmin || DdD - radius > fmax) {
                result.intersect = false;
                return result;
            }
        }

        // M = Cross(n*D-r*R-u*U,A[i])
        for (i = 0; i < 3; ++i) {
            p = frustum.rBound * Math.abs(NBpUC[i]) +
                frustum.uBound * Math.abs(NApRC[i]);
            tmp = frustum.dMin * RBmUA[i];
            fmin = tmp - p;
            if (fmin < 0) {
                fmin *= dRatio;
            }
            fmax = tmp + p;
            if (fmax > 0) {
                fmax *= dRatio;
            }
            DdD = D[0] * NBpUC[i] - D[1] * NApRC[i] + D[2] * RBmUA[i];
            radius = 0;
            for (j = 0; j < 3; ++j) {
                radius += extent[j] * Math.abs(A[j] * NBpUC[i] -
                    B[j] * NApRC[i] + C[j] * RBmUA[i]);
            }
            if (DdD + radius < fmin || DdD - radius > fmax) {
                result.intersect = false;
                return result;
            }
        }

        result.intersect = true;
        return result;
    }
}
