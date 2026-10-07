# Critical review of PR #78: parallel evaluation through one map (#59)

Reviewed at `b3b8602` (`feature/59-parallel`). CI was still pending when this was written.

## Verdict

The engineering is careful. The pins went in first, and the PR checked that they are sensitive. The tests check that the workers really are separate processes running the source under test. The closure-frame guard has a test showing it catches the hazard it exists for. The "Not routed" section is honest. The RNG bug was found, fixed and given a test. These are things a reviewer usually has to ask for.

The weak point is the main design choice. With the default `parallel = TRUE`, **the numbers a user gets depend on whether a `future::plan()` is set.** Separately, two claims go further than the evidence: that the map works "under any plan, plant included", and that a plan never moves random numbers (true for the caller, not for the workers). Items 1–3 are worth settling before merge. The rest can be follow-ups.

## Major

### 1. A plan changes the answer by default, which undercuts the PR's own reproducibility rule

The RNG decision rests on the idea that a result should not depend on whether a plan is set. `derivative_control(parallel = TRUE)` breaks that rule for the resident Jacobian, and so for the classifier, the singularity solver, sensitivity and continuation:

- On DD99 with the Newton equilibrium at the default `equilibrium_eps`, J's error against the closed form goes from 1.3e-6 sequentially to 3.4e-4 under a plan, about 260× worse (PR comment, "Accuracy under a plan").
- On plant, the 2-resident Jacobian moves by 2.5e-4 relative (benchmark table).
- The two paths also leave the gradient closure in different warm states. Sequentially, `state` and `seed_birth_rate` end at the last stencil point. Under a plan they stay at the centre (`R/singularity.R:88-102`). So the singularity solver's next probe starts from a different place, and its trajectory and solve count can differ on an iterated equilibrium. DD99 at `1e-9` happens to give identical counts. Nothing guarantees that on plant.

A user who sets `plan(multisession)` for unrelated work will get different J, roots and borderline verdicts from a colleague who did not. The `parallel = FALSE` escape hatch only helps users who know about it.

**Suggested resolution:** start every stencil point from the centre in the sequential path too. Then:

- results are the same with or without a plan. They are already independent of the worker count: 2 and 4 workers give the same J change in the benchmark.
- `derivative_control(parallel = )` can go, along with its documentation paragraph and its test. `pip_control(parallel = )` remains a separate question.
- `points()` can always go through `regnans_map()`, which removes the duplicated merged-resident refusal logic.
- it should be no slower. The centre is `h` from every point, while the chain starts the down-step `2h` away and the next level about `1.4h` away.

The cost is a one-off re-pin of the sequential values, and that is cheap now that they are pinned. This reverses the ROADMAP decision "With no plan, the sequential path is today's, bit for bit". The PR's own accuracy analysis supports doing so: neither start is more accurate in general, so the chain buys nothing except bit-compatibility with the past.

If you keep the chain, consider defaulting `parallel = FALSE`, so that results are reproducible unless a user opts in. Whichever way you go, record the decision and its reason in ROADMAP "Decided for #59".

### 2. `future.seed = NULL` hides RNG use in workers instead of making it safe

`R/parallel.R:42`. `future.apply` documents `NULL` as "like `FALSE` but without the check whether random numbers were generated". I checked it under `plan(multicore, workers = 2)`. With `set.seed(1)` before each call, `future_lapply(1:4, function(i) runif(1), future.seed = NULL)` gave different values on each of three runs, and no warning. So for worker code the guarantee only holds if no RNG is used, and nothing enforces that.

Nothing routed in regnans draws random numbers today (I grepped `R/`). But the routed work calls harness code, and on plant it calls plant. Neither is under this repo's control. One day a harness will use a random restart or jitter, and under a plan it will then give silently irreproducible results.

**Suggested fix:** use `future.seed = TRUE` inside a save and restore of the caller's `.Random.seed`. `withr::with_preserve_seed(future_lapply(..., future.seed = TRUE))` left the caller's next `runif(1)` at 0.2655, the same as without a plan. Writing it out inline avoids adding a dependency. This gives:

- reproducible streams for each element,
- an untouched caller,
- no false alarm from the misuse check when Rcpp creates `.Random.seed`, because seeds are set explicitly.

Note that an integer `future.seed` is not an alternative. It also advances the caller's stream: I measured 0.573, the same as `TRUE`.

### 3. "Any plan, plant included" has not been tested

The PR body, ROADMAP, AGENTS.md and the roxygen docs all say the routed paths work under any plan on plant. The evidence:

- `test-plant-smoke-singularity.R` runs only `multicore`, the one plan where the pointer guard is skipped.
- It checks that `community_clear_residents(out)` holds no pointer. The full payload also carries `state` (the equilibrium solver's `demography_state`), the seed densities, and for the parameter Jacobian the caller's `parameter` closure. None of these is checked.
- No multisession worker ever loads plant and solves an SCM. The benchmark also uses `multicore` only.

**Suggested fix:** add one plant test under multisession, such as classifying the 1-resident singularity (two solves). If that is too slow for the smoke tier, change the claim to "multicore on plant; any plan on the toy harnesses". The parameter Jacobian on plant is a further limit. Any `p -> community` closure that holds a solved community is refused under a non-fork plan. The docs say so, but that path then does not support "plant included" either.

## Moderate

### 4. The TEP map and the PIP map break the PR's own rule that worker functions are top-level

- `R/pip.R:483` (TEP) maps an anonymous closure. Its frame holds `community`, `pip` and `mutual`. If `community_tep()` gets a solved plant community, the guard now refuses the call under multisession. ROADMAP lists "the TEP's pairs" under "any plan, plant included".
- `R/pip.R:236` (PIP) maps a closure whose frame holds `tf`. The linear `tf$fwd` and `tf$inv` (`R/util.R:87`) are closures over `community_trait_transform()`'s frame, which holds the whole community passed in. Plant uses the log scale (`log` and `exp` are builtins), so this does not bite today. It shows how easily a pointer can travel inside a closure's frame, and returning `identity` there fixes it.

Both closures existed before this PR, but the PR's rules now describe them. Move them to top level, or exclude them from the "plant included" claim.

### 5. The accuracy finding deserves its own issue, not a remark in a PR comment

The finite-difference error of a nested solve is about `equilibrium_eps / h`. At the defaults that is `1e-5 / 1e-3 ≈ 1e-2` relative, and it applies to the sequential path as well. This is arguably the most important thing the PR found, and it is recorded only as "it is a separate issue" with no issue number. The `derivative_control()` docs say a plan "moves the Jacobian by up to that much". Few readers will work out that "that much" can be 1% at the defaults.

File the issue (tie the stencil's equilibrium tolerance to its step, or tighten `equilibrium_eps` within stencils). Link it from ROADMAP and from the `derivative_control()` docs.

### 6. Solve counts are lost when a worker errors

In `points()` (`R/singularity.R:99`), `evaluations` is only increased after the whole map returns. If one worker's solve errors, the singularity solver's error path attaches `e$evaluations` without the solves the other workers finished. Sequentially, each finished solve is counted before an error. The tests treat solve counts as a contract, so this is a gap, though a minor one.

## Minor

- **`evaluate = attr(gradient, "points")`** (`R/derivatives.R:361`). Passing `NULL` explicitly overrides `util_jacobian()`'s default, so a gradient without the attribute fails with "attempt to apply non-function". Every current caller passes a `singularity_gradient_fn()` closure, so nothing breaks today. `attr(gradient, "points") %||% function(p) lapply(p, gradient)` would make it robust.
- **`like <- community_reset(community)`** (`R/derivatives.R:429`) sends the whole reset community to every worker, but `parameter_community()` only reads `like$trait_names`. Passing the trait names would be enough.
- **The guard serialises the whole payload in the main process** just to count pointers, and `future` then serialises it again. That is cheap next to an SCM solve. The guard also errs on the safe side in two ways: objects with harmless pointers (a `data.table`'s self-reference) are refused, and so is a `makeForkCluster()` cluster, whose workers are forks. Both are fine, but worth a line in the comment above `regnans_pointers()`.
- **`helper-parallel.R:20`** calls `pkgload::is_dev_package()` without first calling `skip_if_not_installed("pkgload")`. Each multisession use also builds a new two-process cluster and runs `load_all()` on both workers, about ten times per test run under testthat's own parallelism. Measure the wall-time cost, and consider one cluster per file.
- **Scope of "Closes #59".** #59's title covers grids and scouting. Grids are deliberately not routed (waiting on plant#656) and scouting is #60. ROADMAP records this, so closing is fine. A one-line note on #59 when it closes would save someone reopening it.
- **Speed-ups are modest** (1.4–1.8× at 2 workers for a classification). The main payoff of this PR is that #60 can be built on the map. The PR body could say so plainly, so the benchmark is not read as the case for merging.

## Done well

- Pins committed before the change, with a check that they are sensitive. Seeding every solve from the original start moved J by about 1e-6, while the oracle tolerances (1e-4) would have hidden it.
- Tests check that two workers are actually used, with separate PIDs and the dev namespace path, so a test cannot pass by falling back to `lapply()`.
- A guard test catches an anonymous `function(i) i` whose frame holds a pointer: exactly the hazard that exists.
- A real RNG bug was found and pinned with a test. `future.seed = TRUE` was moving the caller's stream, so whether a plan was set changed stochastic assembly.
- A plan with several workers but no `future.apply` installed now gives an error instead of a silent `lapply()`.
- The "Not routed, and why" section, including the upstream plant#656 finding (34× at 51 mutants).
