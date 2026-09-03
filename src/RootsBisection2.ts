// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsBisection2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Estimate a root to continuous functions F(x,y) and G(x,y) defined on a
// rectangle [xMin,xMax]x[yMin,yMax]. The requirements are that for each
// y' in [yMin,yMax], A(x) = F(x,y') satisfies A(xMin) * A(xMax) < 0, which
// guarantees A(x) has a root. Also, for each x' in [xMin,xMax],
// B(y) = G(x',y) satisfies B(yMin) * B(yMax) < 0, which guarantees B(y) has
// a root. Bisection is performed in the x-direction for A(x). Let x' be the
// root. Bisection is then performed in the y-direction for B(y). Let y' be
// the root. The function value is A(x') = F(x',y'). This effectively is a
// bisection of C(x) = F(x,h(x)) along the curve where G(x,h(x)) = 0.
//
// Port notes:
// - Upstream has two constructors selected by type traits, one for
//   floating-point Real and one for arbitrary-precision Real. As with the
//   RootsBisection1 port, only the floating-point instantiation is ported
//   (the arbitrary-precision path is deferred until the bisectors support
//   BSNumber-based numeric types), so the constructor takes just
//   xMaxIterations and yMaxIterations.
// - Upstream's operator() writes xRoot/yRoot/fAtRoot/gAtRoot to reference
//   parameters and returns the x-iteration count, so the port returns
//   { iterations, xRoot, yRoot, fAtRoot, gAtRoot }. Following the
//   RootsBisection1 precedent, the method is named find.
// - Upstream keeps the outputs as mutable members that persist between
//   calls, and the nested x-bisection can leave them unwritten (see the
//   RootsBisection1 port note about maxIterations equal to 1, and gtengine-js
//   issue #84). Because the RootsBisection1 port zero-initializes its outputs
//   in that case instead of leaving them unwritten, the members here are set
//   to zero where upstream would report a stale value from a previous call;
//   on a first call the behavior matches upstream's zero-initialized members.

import { RootsBisection1 } from './RootsBisection1.js';

export interface RootsBisection2Result {
    // The number of iterations used by the x-direction bisector, with the
    // same interpretation as RootsBisection1: 0 means F(xMin,*)*F(xMax,*) > 0
    // and it is unknown whether the rectangle contains a root, 1 means an
    // x-interval endpoint is an exact root of the x-bisection, and 2 or
    // larger is the bisection iteration count.
    iterations: number;

    // The estimated root (xRoot,yRoot) with fAtRoot = F(xRoot,yRoot) and
    // gAtRoot = G(xRoot,yRoot). See the caveats in the class comments about
    // the y-outputs when the x-bisection does not itself bisect.
    xRoot: number;
    yRoot: number;
    fAtRoot: number;
    gAtRoot: number;
}

export class RootsBisection2 {
    private mXBisector: RootsBisection1;
    private mYBisector: RootsBisection1;
    private mXRoot: number;
    private mYRoot: number;
    private mFAtRoot: number;
    private mGAtRoot: number;
    private mNoGuaranteeForRootBound: boolean;

    // Use this constructor when the numeric type is floating point (the only
    // instantiation ported; see the port notes).
    constructor(xMaxIterations: number, yMaxIterations: number) {
        this.mXBisector = new RootsBisection1(xMaxIterations);
        this.mYBisector = new RootsBisection1(yMaxIterations);
        this.mXRoot = 0;
        this.mYRoot = 0;
        this.mFAtRoot = 0;
        this.mGAtRoot = 0;
        this.mNoGuaranteeForRootBound = false;
    }

    find(F: (x: number, y: number) => number, G: (x: number, y: number) => number,
        xMin: number, xMax: number, yMin: number, yMax: number): RootsBisection2Result {
        // XFunction(x) = F(x,y), where G(x,y) = 0.
        const xFunction = (x: number): number => {
            // YFunction(y) = G(x,y)
            const yFunction = (y: number): number => G(x, y);

            // Bisect in the y-variable to find the root of YFunction(y).
            const yResult = this.mYBisector.find(yFunction, yMin, yMax);
            this.mYRoot = yResult.root;
            this.mGAtRoot = yResult.fAtRoot;

            // Upstream assigns (rather than accumulates) here, so only the
            // most recent y-bisection is reflected in the flag. The quirk is
            // preserved; see the PR "Upstream bug suspects" discussion.
            this.mNoGuaranteeForRootBound = (yResult.iterations === 0);
            return F(x, this.mYRoot);
        };

        // Bisect in the x-variable to find the root of XFunction(x).
        const xResult = this.mXBisector.find(xFunction, xMin, xMax);
        this.mXRoot = xResult.root;
        this.mFAtRoot = xResult.fAtRoot;
        this.mNoGuaranteeForRootBound =
            this.mNoGuaranteeForRootBound || (xResult.iterations === 0);

        return {
            iterations: xResult.iterations,
            xRoot: this.mXRoot,
            yRoot: this.mYRoot,
            fAtRoot: this.mFAtRoot,
            gAtRoot: this.mGAtRoot
        };
    }

    noGuaranteeForRootBound(): boolean {
        return this.mNoGuaranteeForRootBound;
    }
}
