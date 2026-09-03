# Verification wave: independent review of every ported file

The first pass ported each batch with one agent and no independent review.
This wave has a *different* agent re-read every ported file against its
upstream header, add property-based tests, and fix what it finds. State lives
in GitHub: one `verify-batch` issue and one PR per group in
`plan/verify-groups.json`; `porting-status.json` records `verified` per file.

## Per-file procedure (do all of it, for every file in the group)

1. **Read both sides completely.** Upstream header at
   `D:\git\GeometricTools-upstream\GTE\Mathematics\<Name>.h` (baseline commit
   `d29e7758ae26`) and `src/<Name>.ts`. Walk the upstream file function by
   function and locate its counterpart in the port. Do not skim.
2. **Translation-hazard checklist.** For each function confirm:
   - every branch, early return and loop bound is present (off-by-one on
     `<` vs `<=`, `numElements - 1`, reversed loops);
   - integer semantics where upstream uses `int`: division (`Math.trunc` or
     `| 0`), modulo of negatives, `size_t` underflow guards, bit shifts on
     values above 2^31;
   - value vs reference: C++ copies on assignment; TS aliases. Inputs stored
     in members or results must be cloned when upstream copied them. Results
     must not alias internal scratch state that a later call overwrites;
   - `std::swap`, `std::sort` comparators (strict weak ordering, stable?),
     `std::map`/`std::set` iteration order, `std::numeric_limits` constants;
   - `T(0)`/`T(1)`/`(T)0.5` literals, `std::fabs` vs `Math.abs`,
     `std::sqrt` of tiny negatives, `atan2` argument order;
   - constructor defaults and member initialisation order; static locals;
   - output parameters that upstream fills only on some paths: the port's
     result object must have the same "unset" semantics or document a value;
   - template specialisations: was the right one (floating-point vs
     rational/exact) ported, and is the choice recorded?
   - the previous porter's "Upstream bug suspects" claims: re-derive them; a
     wrong "fix" is worse than a preserved quirk.
3. **Property-based tests.** Add a `describe('<Name> verification', ...)`
   block to `test/<Name>.test.ts` using `test/helpers/arbitraries.ts`
   (fast-check). Aim for properties that would catch a translation error,
   not restatements of the code:
   - invariants (closest points lie on their primitives; distance equals
     `|closest[0] - closest[1]|`; `sqrDistance == distance^2`; symmetric
     queries agree under argument swap; TI and FI agree on `intersect`);
   - cross-checks against an independent computation (brute-force sampling,
     a simpler formula, a different algorithm in the library, the exact
     `BSRational` path where one exists);
   - round trips (encode/decode, `Log/Exp`, forward/inverse transforms);
   - degenerate inputs (zero extents, coincident points, parallel/collinear
     configurations, empty arrays) and the documented behaviour for them;
   - upstream worked examples and any numbers quoted in upstream comments.
   Use `check(arb, pred)` (200 runs) and `seededRandom` for brute-force loops.
   A property must pass deterministically; if it needs a tolerance, justify
   it in a comment (condition number, catastrophic cancellation, ...).
4. **Fix what you find.**
   - *Port bug* (the TS disagrees with the C++): fix it, add a regression
     test that fails before and passes after, note it in the PR.
   - *Upstream bug*: follow PORTING.md "Upstream bugs" (fix if it corrupts
     results, otherwise preserve faithfully) and report under
     `## Upstream bug suspects`.
   - *API wart* that violates PORTING.md conventions: fix it if local to the
     file; otherwise describe it under `## API notes` for the orchestrator.
5. **Verdict per file**, one of:
   - `verified`: reviewed line by line, properties added, no defects;
   - `fixed`: one or more port defects fixed (list them);
   - `upstream-bug`: upstream defect found (may combine with the above).

## PR format

- Branch `verify/V<nn>`, title `Verify V<nn>: <category> (<n> files)`,
  body starts with `Closes #<issue>`.
- `## Verification findings`: a table `file | verdict | notes`.
- `## Upstream bug suspects`: as in PORTING.md (omit if none).
- `## API notes`: cross-file convention problems found (omit if none).
- Commits touch only `src/<Name>.ts` and `test/<Name>.test.ts` of files in
  the group (plus `test/helpers/*` only if you add a generic generator;
  keep those additions small and general). Never `src/index.ts`,
  `porting-status.json`, `plan/*`, or config files.
- Quality gates: `npm run typecheck` and `npm test` green.

## Orchestrator duties

- File the group's issue (label `verify-batch`) when work starts; run at most
  two verifier agents concurrently, each in its own worktree.
- After merge: `npm run update:verification -- <groupId>` marks the files
  `verified` in `porting-status.json` (commit on main), file each upstream
  suspect as an `upstream-bug` issue, and act on `## API notes`.
