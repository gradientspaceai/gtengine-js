# Group 39 (numerical) — C++ oracle report

Family `v39-numerical`, 46 cases: 39 compared bit for bit (`exact: true`),
0 with a tolerance, 7 declared deviations. The deep run
(`npm run oracle:deep -- 2000 v39-numerical`, 92 000 records, 5.0 s wall)
passes; outside the deviation cases it compares 1 466 537 floating-point
outputs and **every single one is bit-identical to the MSVC build**.
Nothing in this group needed a tolerance: every path is `+ - * / sqrt`, and
the one libm call on an upstream path (`std::cos` in
`RemezAlgorithm::ComputeInitialXNodes`) was measured to agree bit for bit
between the MSVC runtime and V8 on the finite set of arguments it can
receive (see "The one libm call" below).

**No port defect was found in this group.** Every main case is bit-identical
on every one of its deep-run records, including the minimizer trajectories,
the SVD and eigensolver iteration counts, the Givens sequences, the sort
orders and the Remez exchange counts.

## Coverage

Every entry's "deep run" column is measured over 2000 records per case.

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `Minimize1.h` | `Minimize1.getMinimum` | exact | 8000 reals, 100 % |
| | `Minimize1.getMinimum.defaultGuess` (the 2-argument overload) | exact | 4000 reals, 100 % |
| | `Minimize1.setEpsilonTolerance` (defaults, `SetEpsilon`, `SetTolerance`, `GetEpsilon`, `GetTolerance`) | exact | 12 000 reals, 100 % |
| | `Minimize1.getMinimum.deviation` | deviation (#298) | 1634 of 2000 records deviate |
| | `Minimize1.getMinimum.endpointGuessDeviation` | deviation (#298) | 2000 of 2000 records deviate |
| | `Minimize1.invalidInput` | exact, throw parity | 1365 of 2000 records throw on both sides; 1270 reals, 100 % |
| `MinimizeN.h` | `MinimizeN.getMinimum` | exact | 9650 reals, 100 % |
| | `MinimizeN.setEpsilon` (`SetEpsilon` clamp and `GetEpsilon`) | exact | 11 745 reals, 100 % |
| | `MinimizeN.getMinimum.deviation` | deviation (#146) | 2000 of 2000 records deviate |
| | `MinimizeN.getMinimum.reuseDeviation` (second call on the same object) | deviation (#146) | 2000 of 2000 records deviate |
| `OdeSolver.h` | covered through the four concrete solvers: `GetTDelta`, `SetTDelta` | exact | — |
| `OdeEuler.h` | `OdeEuler.update` | exact | 24 289 reals, 100 % |
| `OdeMidpoint.h` | `OdeMidpoint.update` | exact | 24 026 reals, 100 % |
| `OdeRungeKutta4.h` | `OdeRungeKutta4.update` | exact | 23 705 reals, 100 % |
| `OdeImplicitEuler.h` | `OdeImplicitEuler.update` | exact | 24 213 reals, 100 % |
| | `OdeImplicitEuler.update.singularJacobian` | exact | 11 940 reals, 100 % |
| (`Matrix.h`, the `Inverse` overload `OdeImplicitEuler` and `BlockLDLTDecomposition` resolve to) | `Matrix.inverse.isGaussianElimination` | exact | 26 000 reals, 100 %; the "resolves to Gaussian elimination" boolean is true on all 2000 records |
| `SingularValueDecomposition.h` | `SingularValueDecomposition.solve` (`Solve`, `GetSingularValues`, `GetU`, `GetV`, `GetS`, `GetUColumn`, `GetVColumn`, `GetSingularValue`, the iteration count) | exact | 99 772 reals, 100 % |
| | `.solve.rankDeficient` | exact | 100 461 reals, 100 % |
| | `.solve.tiedSingularValues` | exact | 99 327 reals, 100 % |
| | `.solve.zeroMatrix` | exact | 99 499 reals, 100 % |
| | `.solve.extremeScale` (a lattice matrix times 2^k, \|k\| <= 200) | exact | 100 569 reals, 100 % |
| | `.invalidInput` (constructor and `Solve` asserts) | exact, throw parity | 1415 of 2000 throw |
| | `.invalidIndex` (#478 accessor range checks) | exact, throw parity | 1119 of 2000 throw; 2637 reals, 100 % |
| `SymmetricEigensolver.h` | `SymmetricEigensolver.solve` (`Solve` with each `sortType`, `GetEigenvalues`, `GetEigenvectors`, `GetEigenvector`, `GetEigenvalue`, `GetEigenvectorMatrixType`) | exact | 44 960 reals, 100 % |
| | `.solve.nonConvergence` | exact | 20 176 reals, 100 % |
| | `.solve.invalidSize` (size 0, 1 or a zero budget) | exact | 2000 reals, 100 % |
| | `.solve.decoupledDeviation` | deviation (#80) | 2000 of 2000 records deviate |
| `RemezAlgorithm.h` | `RemezAlgorithm.execute.degreeOne` (no libm at all) | exact | 18 000 reals, 100 % |
| | `.execute.rational` | exact | 33 108 reals, 100 % |
| | `.execute.polynomial` | exact | 31 788 reals, 100 %, 51 of 2000 throw on both sides |
| | `.execute.singleIteration` | exact | 32 247 reals, 100 %, 37 of 2000 throw |
| | `.execute.invalidInput` | exact, throw parity | 1525 of 2000 throw; 1889 reals, 100 % |
| `CholeskyDecomposition.h` | `CholeskyDecomposition.factorAndSolve` (`Factor`, `SolveLower`, `SolveUpper`, run-time class) | exact | 26 301 reals, 100 % |
| | `.factorAndSolve.fixedSize` (the compile-time `CholeskyDecomposition<double,3>`) | exact | 21 852 reals, 100 % |
| | `.invalidSize` (`LogError` on a non-square matrix or a mismatched vector) | exact, throw parity | 1004 of 2000 throw; 4854 reals, 100 % |
| | `BlockCholeskyDecomposition.factorAndSolve` (compile-time class: `Get`, `Set`, `Factor`, `SolveLower`, `SolveUpper`) | exact | 62 110 reals, 100 % |
| | `BlockCholeskyDecomposition.invalidSize` | exact, throw parity | 1482 of 2000 throw |
| | `BlockCholeskyDecomposition.runtimeStrideDeviation` | deviation (#209) | 2000 of 2000 records deviate |
| `LDLTDecomposition.h` | `LDLTDecomposition.factorAndSolve` (`Factor`, both `Solve` overloads, run-time class) | exact | 50 242 reals, 100 % |
| | `.factorAndSolve.fixedSize` (the compile-time `LDLTDecomposition<double,4>`) | exact | 67 992 reals, 100 % |
| | `.invalidSize` | exact, throw parity | 1136 of 2000 throw; 1900 reals, 100 % |
| | `BlockLDLTDecomposition.factorAndSolve` (run-time class: all four `Convert` overloads, `Get`, `Factor`, `Solve`) | exact | 187 065 reals, 100 % |
| | `BlockLDLTDecomposition.factorAndSolve.fixedSize` (the compile-time `BlockLDLTDecomposition<double,2,2>` and its `Convert` overloads) | exact | 101 272 reals, 100 % |
| | `BlockLDLTDecomposition.getSetConvert` (`Set`/`Get` over the full matrix, vector conversions) | exact | 74 658 reals, 100 % |
| | `BlockLDLTDecomposition.invalidSize` | exact, throw parity | 1498 of 2000 throw |
| | `BlockLDLTDecomposition.convertBlockToVector.deviation` | deviation (#209) | 2000 of 2000 records deviate (upstream throws on every one) |

Generator coverage (deep-run histograms, 2000 records per case).
`Minimize1`'s objective cycles four families — a uniform quartic, a
small-lattice quartic (exact ties and flat pieces), an even quartic on an
interval symmetric about 0 (exactly equal endpoint values) and a rational
with a sharp peak — crossed with a uniform or lattice interval and with the
initial guess at the midpoint, at either endpoint or at a random interior
point; the four families are drawn 517/516/540/427 times in the main case.
`MinimizeN` runs in 2, 3 and 4 dimensions (868/614/518) on an anisotropic
quadratic with a cross term (980 records) and on a Rosenbrock valley (1020),
with `maxIterations` spread over 1-6 (777/548/179/160/176/160). The Ode
generator spreads four entry families (uniform 518, lattice 488, stiff with a
large negative diagonal 517, a zero first row 477) over dimensions 2-4
(631/682/687), linear and quadratic right-hand sides (997/1003) and 1-4
steps, and changes `tDelta` between the first and second step. The SVD
generator covers `numCols` 2-4 (669/653/678) and `numRows` 2-5
(156/410/721/713) over six matrix families and the whole `multiplier` range
0.5-24; the iteration count ranges over 0-14 with 765 non-convergent records
across the five `.solve*` cases. The symmetric generator covers sizes 2-6
(712/378/353/278/279 in the main case, 3-6 in the two decoupling cases) over
nine matrix families and all three `sortType` values (663/646/691). Remez
covers degrees 2-5 evenly plus a dedicated degree-1 case, reaches the
non-oscillatory exit on 418 of the 2000 polynomial records and 216 of the
single-iteration ones, and throws the `GetXExtreme` bracket assert on 51 and
37 of them. The Cholesky and LDLT generators cover 1 <= N <= 5 and nine
(blockSize, numBlocks) shapes with 1 <= blockSize, numBlocks <= 4, over six
matrix families drawn 296-364 times each: SPD from an integer `L0*L0^T`
(exact in integers), SPD from a uniform `L0`, positive semidefinite,
indefinite, zero, and an integer symmetric matrix whose first pivot is zero.
`Factor` succeeds on 681 of 2000 scalar Cholesky records, 642 of the
fixed-size ones, and 1000 and 998 of the two scalar LDLT cases.

## The one libm call

`RemezAlgorithm::ComputeInitialXNodes` evaluates `std::cos(j * (GTE_C_HALF_PI
/ degree))` for the odd `j` below `2*degree`. Those are the only non-`sqrt`
libm calls anywhere in this group, and for a fixed `degree` they are a fixed
finite set: the case file's generators draw the degree, not the angle. All 44
of them for `degree = 1..8` were compared between the MSVC 19.44 runtime and
V8 (Node 24) and every one is bit-identical, so the Remez cases are declared
`exact` rather than given a tolerance. Degree 1 does not touch libm at all:
its single Chebyshev cosine is at pi/2 and `ComputeInitialXNodes` overwrites
it with an exact zero; `RemezAlgorithm.execute.degreeOne` is the case that
holds that guarantee independently of the measurement.

Everything else in the family is `+ - * / sqrt fabs` on doubles. In
particular, the objective functions handed to `Minimize1`, `MinimizeN` and
`RemezAlgorithm` are arithmetic-only by construction (Horner on a quartic, a
rational with a positive denominator, a quadratic form, a Rosenbrock valley),
so the minimizer *control flow* — every V-shape test, every bracket update,
every bisection sign — is comparable, which a libm-derived objective would
not be.

## Which `Inverse` overload the C++ really calls

`OdeImplicitEuler::Update` writes `Inverse(dgMatrix)` and
`BlockLDLTDecomposition::Factor` writes `Inverse(Djj, &invertible)`, both on a
`Matrix<N,N,double>`. These are dependent calls resolved by ADL at the point
of instantiation, so the answer depends on the translation unit: with
`Matrix2x2.h`, `Matrix3x3.h` or `Matrix4x4.h` visible, partial ordering picks
their closed-form overloads, which differ from `Matrix.h`'s
Gaussian-elimination template in the last bits (v34). The port implements the
Gaussian-elimination form. The case file therefore includes none of those
headers — `Oracle.h` pulls in `Vector2/3/4.h`, `Matrix.h`, `GVector.h` and
`GMatrix.h`, none of which reaches them — and the case
`Matrix.inverse.isGaussianElimination` re-checks it on every record by
comparing `Inverse(M)` for a 3x3 and a 2x2 against an explicit
`GaussianElimination<double>()` call and emitting the result as a boolean. It
is true on all 2000 deep-run records.

## Port defects fixed

None. Every main case of this group is bit-identical to the MSVC build on
every deep-run record.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `Minimize1.getMinimum.deviation` | [#298](https://github.com/gradientspaceai/gtengine-js/issues/298) | 1634 of 2000 records. An even quartic on an interval symmetric about 0 has exactly equal endpoint values, so the vertex of the interpolating parabola lands on the middle sample to within round-off. Upstream's exact `tv == tm` test misses it, one of the asymmetric branches runs and the bracket collapses; and on the following iteration `\|denom\| <= epsilon` and upstream returns, abandoning the search. The port resolves the vertex to within 8 ulps of the bracket and refines around the midpoint instead. The selector is `Min1Differs`, a verbatim `ProbeMinimize1` in the case file run twice — once as upstream, once as the port — keeping the draws whose `(tMin, fMin)` differ bit for bit; the 366 agreeing records are draws for which the capped rejection loop (64 attempts) found no deviating configuration and kept its last candidate. |
| `Minimize1.getMinimum.endpointGuessDeviation` | [#298](https://github.com/gradientspaceai/gtengine-js/issues/298) | 2000 of 2000 records, and each one is a demonstrably wrong upstream answer. With the initial guess at an endpoint the first bracket is `{(t0,f0),(t1,f1),(t1,f1)}`, `denom` is exactly zero, and upstream returns after three evaluations reporting the better endpoint. `F(t) = c*t^2` here, whose minimum is 0 at the interior point t = 0; upstream returns, for instance, `tMin = 1, fMin = 1.2724` on `[-2,1]`. |
| `MinimizeN.getMinimum.deviation` | [#146](https://github.com/gradientspaceai/gtengine-js/issues/146) | 2000 of 2000 records. Upstream's stray `mDConjIndex = 0` writes the conjugate direction into slot 0, where the cycling shift immediately overwrites it, while slot `dimensions` keeps the first iteration's conjugate direction forever; the direction set degenerates and Powell's method stalls. The selector is a verbatim `ProbeMinimizeN` in the case file run twice, once with the stray assignment and once without, keeping the records where the reported minima differ. |
| `MinimizeN.getMinimum.reuseDeviation` | [#146](https://github.com/gradientspaceai/gtengine-js/issues/146) | 2000 of 2000 records. `mDConjIndex` is never reset by `GetMinimum`, so upstream's second call on the same object starts from whatever the first left behind. The case runs `GetMinimum` twice and compares the second result. |
| `SymmetricEigensolver.solve.decoupledDeviation` | [#80](https://github.com/gradientspaceai/gtengine-js/issues/80) | 2000 of 2000 records. When a Householder step is degenerate (`length == 0`), the reflection applied to the matrix is the identity but upstream still stores `2/Dot(v,v) == 2` as its parameter, and `GetEigenvectors`/`GetEigenvector` rebuild `I - 2*e*e^T` from it. The selector is `EigenDecouples`, a verbatim copy of upstream's `Tridiagonalize` that reports the degenerate step — an exact characterization, since the port changes only the stored value. Measured on the deep run: **upstream's Q violates `A v = lambda v` on 584 of the 1556 checked records, with worst relative residual 3.07, while the port satisfies it to 2.4e-15 on all of them.** |
| `BlockCholeskyDecomposition.runtimeStrideDeviation` | [#209](https://github.com/gradientspaceai/gtengine-js/issues/209) | 2000 of 2000 records. The run-time `BlockCholeskyDecomposition<Real,0,0>` addresses scalars inside a `BlockSize`-by-`BlockSize` block with its block-level helper `GetIndex(row,col) = col + row*NumBlocks`, so `SolveLower`, `SolveUpper`, `LowerTriangularSolver` and `SubtractiveUpdate` touch the wrong entries whenever `BlockSize != NumBlocks`; the compile-time specialization, which the port follows, uses `block(row,col)`. The case uses `NumBlocks = 2 or 3 < BlockSize = 3 or 4`, where every one of upstream's indices is still inside the block's storage (`BlockSize < NumBlocks` would read past it, which is why the generator excludes it). Measured on the deep run: **126 of the 155 records upstream factors successfully violate `L L^T = A`, with worst relative residual 4.38.** |
| `BlockLDLTDecomposition.convertBlockToVector.deviation` | [#209](https://github.com/gradientspaceai/gtengine-js/issues/209) | 2000 of 2000 records, all of them throws on the C++ side. `BlockLDLTDecomposition<T>::Convert(BlockVector, GVector&)` verifies each block vector with `current.GetSize() == NumBlocks`, but block vectors have `BlockSize` components, so upstream rejects every valid input with `BlockSize != NumBlocks`. The generator draws `NumBlocks != BlockSize` explicitly. |

The preserved upstream defects of these headers are compared bit for bit and
are reached on purpose:
[#478](https://github.com/gradientspaceai/gtengine-js/issues/478)
(`MinimizeN` stalling in a proper subspace while reporting convergence — see
the independent-reference table; `RemezAlgorithm::ComputeXExtremes` pinning
the outer nodes to `xMin`/`xMax`; `RemezAlgorithm::Execute` leaving
`GetXNodes()` one exchange ahead of `GetErrors()` and `GetCoefficients()`,
which the Remez cases emit as-is) and
[#147](https://github.com/gradientspaceai/gtengine-js/issues/147)
(`maxBracketIterations` validated and stored but never read; the cases draw
it and it changes nothing on either side).

`RemezAlgorithm::ComputePCoefficients` reads `poly[i]` past the coefficient
array when the Newton product loses its leading coefficient (the other half
of #147, fixed in the port by returning 0 there). That is upstream undefined
behavior and it is kept out of the generators: it needs `u[degree] - E *
v[degree]` to be *exactly* zero, which the drawn coefficients never produce,
and the 92 000-record deep run never hit it.

## Independent-reference checks (agreement is not correctness)

Every main case's 2000-record deep-run output was checked once against a
reference computed outside both implementations. Except where stated the
check runs on the C++ outputs, which the port reproduces bit for bit, so it
checks both sides.

| check | result |
| --- | --- |
| SVD: `U S V^T = A`, `U^T U = I`, `V^T V = I`, and the singular values in descending order (never forming `A^T A`, which would floor the measurement at sqrt(eps)) | `.solve` 1733 converged records, worst `\|USV^T - A\|/\|A\|` 9.9e-15, worst `\|U^TU - I\|` 2.2e-15, worst `\|V^TV - I\|` 2.4e-15; `.rankDeficient` 1745 records, 7.1e-15 / 2.0e-15 / 2.2e-15; `.extremeScale` 1757 records, 8.1e-15 / 2.2e-15 / 1.8e-15; `.tiedSingularValues` and `.zeroMatrix` exact (diagonal input). 0 descending-order violations over all 9235 checked records. |
| `SymmetricEigensolver`: `A v = lambda v` and `Q^T Q = I` | main case 1921 converged records, worst `\|Av - lv\|/\|A\|` 3.9e-15, worst `\|Q^TQ - I\|` 2.2e-15, 0 violations above 1e-10 — for upstream and for the port alike. |
| `SymmetricEigensolver`, the #80 deviation case, upstream versus the port | 1556 checked records: **upstream 584 violations, worst residual 3.07; the port 0 violations, worst 2.4e-15.** The port's fix corrects genuinely wrong answers, and both keep `\|Q^TQ - I\| <= 1.8e-15` (upstream's Q stays orthogonal; it is simply the wrong orthogonal matrix). |
| `CholeskyDecomposition`: `L L^T = A` and the residual `A x - b` of `SolveLower` followed by `SolveUpper` | run-time class 681 factored records, worst `\|LL^T - A\|/\|A\|` 3.5e-16, worst `\|Ax - b\|/\|b\|` 2.0e-14; compile-time class 642 records, 3.9e-16 / 8.9e-15. |
| `BlockCholeskyDecomposition` (compile-time): `L L^T = A` over the assembled matrix | 726 factored records, worst 4.3e-16. |
| `LDLTDecomposition`: `L D L^T = A` and the solve residual | run-time class 1000 factored records, worst `\|LDL^T - A\|/\|A\|` 1.7e-13, worst `\|Ax - b\|/\|b\|` 3.5e-13; compile-time class 998 records, 5.4e-13 / 1.9e-12. Every one of the 330 positive-semidefinite records is correctly rejected, as are the 331 zero and 339 zero-first-pivot ones. |
| `BlockLDLTDecomposition`: `L D L^T = A` and the solve residual | 1242 factored records. Per generator mode: SPD-integer 1.0e-15, SPD-uniform 3.0e-15, indefinite 3.7e-11, zero-first-pivot 4.2e-13 — and **positive semidefinite 1.28, with a solve residual of 9.2e15 over 57 records**. See the suspects below. |
| Ode solvers: order of convergence on `x' = -x`, `x(0) = 1`, over [0,1] with 64/128/256/512 steps, against `exp(-1)` | `OdeEuler` observed order 1.00, `OdeMidpoint` 2.00, `OdeRungeKutta4` 4.01, `OdeImplicitEuler` 1.00 — each the textbook order. |
| `Minimize1`: is `fMin` a minimum? `fMin == F(tMin)` bit for bit, and `fMin` against a 401-point grid over `[t0,t1]` | `F(tMin) == fMin` exactly on all 2000 records. The grid beats `fMin` on 201 of them; 99 of those had `maxSubdivisions <= 3`. Re-running the port with a much larger budget (24 subdivisions, 200 bisections, epsilon 1e-14) leaves only **6** beaten. So this is the drawn budget on a multimodal objective, not a defect: `Minimize1` is a bounded-effort global search. |
| `MinimizeN`: `fMin` against 2000 random probes in the box | probes beat `fMin` on 794 of 2000 records; 717 of those had `maxIterations <= 2`. With a large budget (12, 128, 60 iterations) 114 remain beaten, worst relative excess 5.5 — the documented limitation of basic Powell on a curved valley, [#478](https://github.com/gradientspaceai/gtengine-js/issues/478), preserved. |
| `RemezAlgorithm`: does the error equioscillate? | over the 1348 rational and 1022 polynomial records that ran at least three exchanges, the node errors alternate in sign on **every** record and `max\|e\| / min\|e\|` is **1.00** to two decimals on every record: the minimax characterization holds exactly. |
| `Matrix.inverse.isGaussianElimination` | the C++ side's "this `Inverse` is the Gaussian-elimination overload" boolean is true on all 2000 records. |

## Not covered

| header / entry point | reason |
| --- | --- |
| `OdeSolver.h` as a standalone case | Abstract base class with no computational entry point of its own. Its two members `SetTDelta`/`GetTDelta` are exercised through all four concrete solvers (each `.update` case reads `GetTDelta`, changes it between step 0 and step 1 and reads it again), and `Update` is pure virtual. |
| `OdeSolver`/`OdeEuler`/... with `TVector = GVector<Real>` | Upstream templates each solver on `TVector` so it can be `Vector<N,Real>` or `GVector<Real>`. The arithmetic is the same expression in both (`xIn + tDelta * F`, componentwise), and `OdeImplicitEuler` cannot use `GVector` at all because it calls `TMatrix::Identity()`, which only `Matrix<N,N,Real>` has. The cases therefore instantiate `Vector<N,double>`/`Matrix<N,N,double>`, which is what the port's single run-time-sized `Vector`/`Matrix` corresponds to. |
| `SymmetricEigensolver::GetEigenvectors` after a `Solve` that did not converge | It loops forever, on the MSVC build and in the port alike. See the suspects below; the oracle driver would not survive such a record, so every case that calls `GetEigenvectors` requires a converged `Solve` by rejection sampling, and `SymmetricEigensolver.solve.nonConvergence` covers the other branch through `GetEigenvalues`, `GetEigenvector` and `GetEigenvalue`, which are safe. |
| `SingularValueDecomposition::GetS`, `GetU`, `GetV` with a null pointer | Upstream signals "give me this output" with a pointer and asserts non-null. The port returns the array, so the null case cannot be expressed. |
| `BlockCholeskyDecomposition<Real, B, N>` with `B*N > 6` and `BlockLDLTDecomposition` with `B*N > 9` | Record size. Goldens are committed and a 12x12 block matrix is 144 doubles of input and several hundred of output per record; the shapes drawn cover every structural combination (`B < N`, `B == N`, `B > N`, `N == 1`, `B == 1`). |
| `CholeskyDecomposition<Real, N>` / `LDLTDecomposition<T, N>` / `BlockLDLTDecomposition<T, B, N>` for every compile-time `N` | One representative of each compile-time specialization is covered (`CholeskyDecomposition<double,3>`, `LDLTDecomposition<double,4>`, `BlockLDLTDecomposition<double,2,2>`) beside the run-time ones. The two specializations of each class are textually identical apart from the size source and the `LogAssert`s, and the port collapsed them into one run-time-sized class, so a second compile-time size would test the same code. |
| `BlockCholeskyDecomposition`'s run-time class on `BlockSize < NumBlocks` | Upstream's in-block `GetIndex(i,j) = j + i*NumBlocks` then exceeds the block's `BlockSize*BlockSize` storage and `GMatrix::operator[]` reads out of range — undefined behavior, kept out of the generators. `BlockSize > NumBlocks` gives the same defect with every index in range and is the deviation case. |
| `RemezAlgorithm` with `degree > 5` | The initial-node cosines were verified bit-identical up to degree 8, but the record grows as `3*degree + 5` doubles and degrees 1-5 already reach every branch (oscillatory and non-oscillatory exits, the `GetXExtreme` bracket assert, single-iteration and multi-iteration runs). |
| `Minimize1`/`MinimizeN`/`RemezAlgorithm` with objectives that call libm | Deliberate: a 1 ulp `sin`/`cos` difference between the MSVC runtime and V8 changes a bracket test and sends the search to a different local minimum, which no tolerance can repair (ORACLE.md). The objectives are arithmetic-only so that the trajectories are comparable. |

## Upstream bug suspects

Two new ones, neither already in `docs/UPSTREAM-FINDINGS.md`. In both the port
and the MSVC build agree bit for bit and both are wrong; the port preserves
the behavior and each has a pinned regression test.

**1. `SymmetricEigensolver::GetEigenvectors` loops forever after a `Solve`
that did not converge (hang).** `Solve` calls `ComputePermutation` only on
the converged path, so after exhausting the iteration budget `mPermutation`
still holds the zeros the constructor's `resize` gave it. Every accessor then
takes its "sorting was requested" branch, because the flag for "no sorting"
is `mPermutation[0] == -1`. For `GetEigenvalues` that is merely wrong (it
reports `mDiagonal[0]` for every index). For `GetEigenvectors` it is fatal:

```cpp
while ((next = mPermutation[current]) != start)
{
    ...
    current = next;   // next is always 0
}
```

At `i = 1` the loop starts with `start = 1`, `current = 1`, reads
`mPermutation[1] == 0`, sets `current = 0`, reads `mPermutation[0] == 0`,
and never terminates.

*Reproduction*: `SymmetricEigensolver<double> solver(2, 1); double A[4] =
{1,1,1,2}; solver.Solve(A, 1);` returns `0xFFFFFFFF`, `GetEigenvalues` gives
`(0.381966, 0.381966)`, and `GetEigenvectors` hangs. Measured directly: the
MSVC build was killed after 8 s with no progress, and the port
(`getEigenvectors` on the same input) was killed after 2 minutes. It is
reachable from any non-converged `Solve`: the oracle case
`SymmetricEigensolver.solve.nonConvergence` selects exactly those inputs and
finds one on the first attempt for 2000 records out of 2000, over sizes 3 to
6 and all nine matrix families, with budgets of 1 to 3 iterations.

*Port: preserved*, because the port already behaves identically and a fix
would have to invent a return value upstream does not define. A one-line fix
is available if the project wants it — set `mPermutation[0] = -1` before the
iteration loop, so a failed `Solve` leaves the accessors on their "no
sorting" path — but it changes an output on a path where upstream currently
produces none. Pinned by `test/SymmetricEigensolver.test.ts`, "upstream
(preserved): a non-converged solve leaves the permutation in a state that
makes getEigenvectors loop forever", which asserts the state that leads to
the hang without entering the loop, and by the oracle case
`SymmetricEigensolver.solve.nonConvergence`, which reproduces the
non-convergence and the degenerate `GetEigenvalues` result on 2000 records.

**2. `BlockLDLTDecomposition::Factor` reports success for a singular
diagonal block, and returns a factorization that is nowhere near `A`
(result-corrupting).** The scalar `LDLTDecomposition` tests its pivot with
the exact `Djj == 0` and therefore rejects a singular matrix reliably: over
the deep run it rejected all 330 positive-semidefinite records, all 331 zero
ones and all 339 zero-first-pivot ones. The block version instead asks
`Inverse(Djj, &invertible)`, which is `GaussianElimination`'s full-pivoting
elimination. For a block of size 3 or more the elimination has already
divided by earlier pivots, so the entries are no longer integers, the pivot
of an exactly singular block is a round-off residue instead of zero, the
`invertible` flag comes back true, and `invDjj` has entries around 1e16.

*Reproduction* (found by searching; the deep run reproduces the same class on
57 records of `BlockLDLTDecomposition.factorAndSolve`, worst relative error
1.28). `blockSize = 3`, `numBlocks = 2`, and

```
 4   0   4  -4  -4   2
 0   1  -1   1   0  -1
 4  -1   5  -5  -4   3
-4   1  -5   9   2  -1
-4   0  -4   2   7  -3
 2  -1   3  -1  -3   9
```

which is `L0*L0^T` for an integer lower-triangular `L0` with one zero on the
diagonal, hence positive semidefinite and exactly singular. The scalar
`LDLTDecomposition(6)` returns `false` on it. `BlockLDLTDecomposition(3,2)`
returns `true`, its `L D L^T` differs from `A` by 0.556 of `max|A|`, and
`Solve` reports a "solution" whose first component is 9.0e15.

*Port: preserved*, by the same judgment already recorded for
`GaussianElimination`'s denormal pivots
([#375](https://github.com/gradientspaceai/gtengine-js/issues/375)) and for
`BandedMatrix::ComputeInverse` (group 38): there is no exact condition
separating "this pivot is the round-off residue of an exact zero" from "this
pivot is legitimately small", and a tolerance would turn correct
factorizations into failures. The input also violates the class's documented
positive-definite precondition — but so does the input the scalar class
correctly rejects, and `Factor` returning `true` is what turns a precondition
violation into a silent wrong answer. Pinned by
`test/LDLTDecomposition.test.ts`, "upstream (preserved): factor reports
success for a singular diagonal block of size 3", and by the oracle case
`BlockLDLTDecomposition.factorAndSolve`, whose positive-semidefinite
generator mode reaches it.

**Not a defect, for the record.** The 201 `Minimize1` and 794 `MinimizeN`
records on which a grid or a random probe finds a lower value are the drawn
iteration budgets and the documented limitation of basic Powell, not a
failure of the implementations: raising the budget removes 97 % and 86 % of
them respectively.
