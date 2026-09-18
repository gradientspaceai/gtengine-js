# Oracle group 19 — `v19-distance`

Differential testing of the 19 headers of verify group 19 against the upstream
GTE build (MSVC 194435215, x64, `/O2 /fp:precise`, upstream commit
`d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`). See [ORACLE.md](../../ORACLE.md).

- `oracle/cpp/cases/v19-distance.cpp` — 43 cases
- `oracle/golden/v19-distance.txt` — 860 committed records (20 per case)
- `test/oracle/v19-distance.oracle.test.ts` — the TypeScript replay

**Result.** 39 ordinary cases, 780 records, 6976 floating-point outputs,
**100.00% bit-identical** to the C++ build. 4 further cases demonstrate
deliberate fixes of upstream defects. No case is skipped and no tolerance is
used anywhere in the family: every ordinary case is declared `{ exact: true }`,
which is sound because every path in these headers is `+ - * /`, `sqrt`,
`fabs`, `min`, `max` and comparisons only — no libm call.

Deep runs: `npm run oracle:deep -- 2000 v19` passes (86 000 records), and an
extra sweep at 20 000 records per case (860 000 records) also passes with no
disagreement.

Every generator mixes five modes selected from `io.index()`: uniform random
(modes 0, 1), a small integer lattice (modes 2, 3 — exact arithmetic, and
therefore exactly parallel directions, exactly coincident points, exactly
zero direction components, zero box extents and vertices exactly on a line),
and a construction aimed at a specific branch (mode 4 — a second direction
that is an exact scalar multiple of the first, a degenerate segment, a
collinear triangle, the zero direction). Every mode records the same number of
doubles, so the replay reads the inputs without knowing which mode produced
them.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `DistAlignedBoxAlignedBox.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistCircle2Circle2.h` | `compute` | exact | pass |
| `DistLine2AlignedBox2.h` | `compute`, `doQuery` | exact | pass |
| `DistLine2Circle2.h` | `compute` | exact | pass |
| `DistLine2Triangle2.h` | `compute`, `compute.deviation` | exact / deviation | pass |
| `DistLine3CanonicalBox3.h` | `compute`, `compute.deviation` | exact / deviation | pass |
| `DistLineLine.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistLineRay.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistLineSegment.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistPoint2Circle2.h` | `compute` | exact | pass |
| `DistPointCanonicalBox.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistPointLine.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistPointRay.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistPointRectangle.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistPointSegment.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistPointTriangle.h` | `compute.2d`, `compute.3d`, `useConjugateGradient.2d`, `useConjugateGradient.3d` | exact | pass |
| `DistRayRay.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistRaySegment.h` | `compute.2d`, `compute.3d`, `compute.3d.deviation` | exact / deviation | pass |
| `DistSegmentSegment.h` | `compute.2d`, `compute.3d`, `computeEndpoints.3d`, `computeRobust.2d`, `computeRobust.3d`, `computeRobustEndpoints.3d`, `computeRobust.3d.deviation` | exact / deviation | pass |

Entry points covered: every `operator()` of every header, plus
`DistPointTriangle::UseConjugateGradient`, the four-endpoint overloads
`DistSegmentSegment::operator()(P0,P1,Q0,Q1)` and
`DistSegmentSegment::ComputeRobust(P0,P1,Q0,Q1)`, `ComputeRobust` on the
two-segment overload, and `DistLine2AlignedBox2::DoQuery`. The last is
`protected` upstream (with friendship granted to
`DCPQuery<T, Line2<T>, OrientedBox2<T>>`) but the port exports it as the free
function `distLine2AlignedBox2DoQuery` for `DistLine2OrientedBox2` to call, so
it is a public entry point of the port; the C++ case reaches it through a
derived class that re-publishes it with a `using` declaration. It sets only
`parameter` and `closest` — `distance` and `sqrDistance` are left at their
`Result{}` defaults on that path and are not emitted — and it reflects its
`origin` and `direction` in place without undoing the reflection, so those
mutated in/out parameters are emitted as outputs too.

Some branch notes:

- `DistLine2AlignedBox2` and `DistLine3CanonicalBox3` dispatch on the number of
  nonzero direction components. Mode 4 passes the exactly zero direction
  (`DoQuery0D`) and the lattice modes give directions with one or two zero
  components (`DoQuery1D`, `DoQuery2D`), so all four dispatch arms are reached.
- `DistCircle2Circle2` reaches the concentric and cocircular branches through
  mode 4, and exact external/internal tangency through the integer-radius
  lattice modes.
- `DistRaySegment.compute.2d/3d` accept their inputs through a rejection loop
  (see below), so their mode choice is made from unrecorded draws: modes 2 and
  3 use lattice origins, lattice direction vectors and lattice endpoints, which
  reach `det == 0` (the parallel branch) exactly, and mode 4 makes the segment
  an exact scalar multiple of the ray direction.
- `DistPointTriangle` reaches all seven regions in both 2D and 3D through the
  uniform modes and reaches the degenerate (collinear) triangle through mode 4.

## Port defects fixed

### `src/DistLine2Triangle2.ts` — the fix for upstream issue #441 was applied too broadly

`LineIntersectsTwoEdges` computes the edge parameter. Upstream evaluates

```cpp
T s = DotPerp(D, P - V[i0]) / DotPerp(D, V[i1] - V[i0]);
```

and the port had replaced it unconditionally by the algebraically equal
quotient of the normal components already computed for the sign test,
`ncomp[i0] / (ncomp[i0] - ncomp[i1])`, because upstream's denominator can round
to exactly 0 while the two normal components are tiny and of opposite sign,
which makes every output NaN (docs/UPSTREAM-FINDINGS.md, issue
[#441](https://github.com/gradientspaceai/gtengine-js/issues/441)).

The two expressions round differently. The oracle found the port disagreeing
with the C++ build on **252 of 2000 records** (12.6%) of
`DistLine2Triangle2.compute` — every record whose vertex signs send it to
`LineIntersectsTwoEdges`, which is the ordinary line-crosses-triangle case, not
a degenerate one. The affected outputs were `parameter`, `barycentric[i0]`,
`barycentric[i1]` and both closest points; the error was a few ulps, up to
`1.2e-15` scaled. `distance` and `sqrDistance` stayed bit-identical because
both are 0 on this path.

Fix: evaluate upstream's expression and fall back to the normal components only
when upstream's denominator is exactly 0, which is the only input class on
which upstream is defective. The defect fix is preserved (the deviation case
below still disagrees with the C++ build on every one of its records) and the
generic path is now bit-identical. `test/DistLine2Triangle2.test.ts` gains a
regression test that asserts the port reproduces upstream's quotient bit for
bit on 200 constructed `++-` configurations, so the fix cannot be re-broadened
without a test failure.

Follow-up for the orchestrator: the "Suggested fix (applied in the port)"
paragraph of the `DistLine2Triangle2.h` entry in `docs/UPSTREAM-FINDINGS.md`
still describes the unconditional substitution. It should say that the port
applies it only when upstream's denominator rounds to zero. ORACLE.md does not
allow a group branch to touch that file, so the edit is not in this branch.

No other port defect was found: after this one fix, all 780 committed records
and all 860 000 deep-run records of the 39 ordinary cases agree bit for bit.

## Deliberate deviations demonstrated

Each of these cases generates only inputs on which upstream is defective; the
harness requires the port to disagree with the C++ build on at least one
record, so the fix is demonstrated against the real build rather than against a
TypeScript re-derivation.

| case | records disagreeing | upstream defect |
| --- | ---: | --- |
| `DistRaySegment.compute.3d.deviation` | 20 of 20 | The ray parameter is never clamped to `s0 >= 0` in the parallel branch and in regions 1 and 5, so upstream reports a "closest ray point" that is behind the ray origin and a distance smaller than the true one (issue [#126](https://github.com/gradientspaceai/gtengine-js/issues/126)). The generator places the segment behind the ray origin and keeps only records for which upstream's `parameter[0]` is negative. The main `compute.2d`/`compute.3d` cases keep only records for which upstream's `parameter[0] >= 0`, where the port's clamp is a no-op and the two agree bit for bit. |
| `DistSegmentSegment.computeRobust.3d.deviation` | 9 of 20 | `ComputeIntersection` replaces an endpoint t-coordinate `f/b` that falls outside `[0,1]` by `1/2`; the port clamps it to the nearest endpoint of `[0,1]` (issue [#418](https://github.com/gradientspaceai/gtengine-js/issues/418)). The generator makes `P1-P0` perpendicular to `P1-Q1`, so the line `dR/ds = 0` passes through the domain corner `(1,1)` and the ratio rounds one ulp outside the interval; a C++-side replica of `GetClampedRoot` selects the records that reach it. |
| `DistLine2Triangle2.compute.deviation` | 20 of 20 | The edge-parameter denominator rounds to exactly 0 while the vertex signs still differ, so upstream divides by zero and returns NaN (issue [#441](https://github.com/gradientspaceai/gtengine-js/issues/441)). The generator places two vertices on the line and keeps only records for which the C++ `distance` is NaN. |
| `DistLine3CanonicalBox3.compute.deviation` | 20 of 20 | For a flat (zero-extent) box the cancellation-prone accumulation `pme^2 + tmp^2 + PpE^2 + delta*parameter` can make `sqrDistance` slightly negative, so `sqrt` gives NaN; the port clamps the squared distance at zero (issue [#421](https://github.com/gradientspaceai/gtengine-js/issues/421)). The generator aims a line at a point of a flat box and keeps only records for which the C++ `sqrDistance` is negative. |

## Not covered

Nothing. All 19 headers of the group are implemented by the port and all have
cases. No header in the group is pure data.

Two narrower notes, not gaps in entry-point coverage:

- `DistLine2Circle2::DoQuery` and `DistCircle2Circle2::DoQuery` are `protected`
  and `private` upstream and are not exported by the port, so they are covered
  only through `operator()`, which reaches every branch of both.
- `DistSegmentSegment`'s four-endpoint overloads and
  `DistRaySegment`'s deviation case are exercised in 3D only; the 2D and 3D
  paths are the same `N`-dimensional template instantiated twice, and both
  dimensions are covered for the two-segment overloads.

## Upstream bug suspects

None new. The four defects exercised above are already recorded in
`docs/UPSTREAM-FINDINGS.md`. The documented unreliability of the
`det = max(a00*a11 - a01*a01, 0) > 0` parallelism test shared by
`DistLineLine`, `DistLineRay`, `DistLineSegment`, `DistRayRay`,
`DistRaySegment` and `DistSegmentSegment::operator()` (issue
[#418](https://github.com/gradientspaceai/gtengine-js/issues/418)) was
confirmed against the real build and is deliberately preserved: the lattice and
mode-4 generators feed exactly parallel direction pairs, upstream takes the
nonparallel branch on some of them, and the port takes the same branch and
produces the same bits.
