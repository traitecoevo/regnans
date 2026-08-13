# regnans — developer guide for Claude

R package for **assembling ecological communities with the `plant` model**: it
repeatedly invades empty/occupied trait space, solves each community to
demographic equilibrium, computes invasion-fitness landscapes, and finds
evolutionary attractors (selection gradients, singular strategies).

It is a thin evolutionary-assembly layer **on top of `plant`**. `plant` provides
the demographic model (the SCM — Size- and Patch-structured Cohort Model);
`regnans` orchestrates many `plant` runs.

## Status (June 2026)

The package was mid-way through a refactor that **ports equilibrium/demography
code out of `plant` and into this package**, while `plant` itself went through
breaking interface changes. As of now:

- `devtools::load_all()` **succeeds**.
- The **core path works and is verified**: `community_start` → `community_add` →
  `community_demography()` (default `equilibrium_iteration` solver) →
  `community_selection_gradient()`. At equilibrium resident fitness ≈ 0 and the
  selection gradient is finite — matches `vignettes/solving_attractors.Rmd`.
- **Viable-bounds path works** and is reimplemented on the community machinery:
  `community_viable_fitness_1D(community)` / `community_viable_bounds()` use an
  empty community's `fitness_function` as the fundamental-fitness function (no
  dependency on plant's removed `fundamental_fitness()`/`viable_fitness()`).
- See **Known issues** below for what is *not* yet working.

### Important: plant removed the whole fitness/equilibrium subsystem (#388)

plant NEWS (#388): "All fitness/equilibrium functionality was removed from plant;
it now lives in `regnans`." Removed from plant: `fitness_landscape()`,
`solve_max_fitness()`, `viable_fitness()`, `fundamental_fitness()`,
`equilibrium_birth_rate()`, `positive_1d()`, `max_growth_rate()`, etc. The port
is *exactly* the job of bringing these into this package. The intended approach
is to **reimplement on the community machinery** (run the SCM via the community's
`fitness_function` / `community_demography`) rather than copy plant's old code
verbatim. Done so far: equilibrium (→ `community_demography`), viable bounds
(→ `community_viable_fitness_1D`), `positive_1d`/`positive_1d_bracket` (pure
numeric helpers, in `R/community_fitness_viable.R`), and the fitness-max
functions `max_fitness()` / `max_growth_rate()` (`R/community_fitness_solve_max.R`,
operating on `community$fitness_function`). The package no longer references any
removed/undefined plant symbol (verified with `codetools::findGlobals`).

`HANDOFF.md` is an ephemeral note from the previous session (not committed);
this file supersedes it for durable guidance.

## Dev commands

```r
devtools::load_all(".")     # load package (works)
devtools::document()        # regenerate man/ + NAMESPACE from roxygen
devtools::test()            # run testthat suite (see baseline below)
```

There is a manual smoke script mirroring the core vignette path; the canonical
end-to-end exercises live in `vignettes/solving_attractors.Rmd`,
`vignettes/assembly_fitmax.Rmd`, `vignettes/assembly_stochastic.Rmd`.

## Architecture

### The `community` object (`R/community.R`)

A list with class `community`. Built by `community_start(bounds, ...)`, then
mutated through a pipeline. Key fields:

- `bounds`, `trait_names`, `traits` (matrix), `birth_rate` (per resident)
- `demography_control` — how to solve demographic equilibrium (see controls)
- `model_support` — the bridge to `plant`: `list(p = <plant Parameters>,
  plant_control = <plant control()>)`. Also caches `node_schedule_times` /
  `node_schedule_ode_times` between solves.
- `fitness_control`, `fitness_function`, `fitness_points`, `resident_fitness`,
  `selection_gradient` — populated by the fitness/gradient functions.

The typical pipeline (see `solving_attractors.Rmd`):

```r
community_start(bounds, model_support = list(p = ..., plant_control = ...)) |>
  community_add(trait_matrix(x, "lma"), birth_rate = 200) |>
  community_demography() |>            # solve to equilibrium birth rates
  community_selection_gradient()       # fitness gradient at residents
```

### plant bridge (`R/community_plant.R`)

All direct `plant` coupling lives here, in the `plant_community_*` functions.
These are no longer wired up by aliases at the bottom of the file; instead
`harness_plant()` (in `R/harness.R`) points the model-agnostic connectors at
them. See **Model harnesses** below — the plant layer is now one harness among
several, and the connectors dispatch through `community$harness`.

- `plant_default_assembly_pars(hmat, max_patch_lifetime, fixed_RA)` — base FF16
  `Parameters`. **Important:** it regenerates `node_schedule_times_default` to
  match `max_patch_lifetime` (otherwise the default schedule spans the original
  ~105 and the SCM errors with "time_max must be greater than current time").
- `plant_default_assembly_control(...)` — wrapper around `plant::control()`.
- `plant_community_parameters(community)` — builds a `plant` `Parameters` from
  the community (strategies via `plant::generate_strategy`, restores cached schedule
  times).
- `plant_community_make_demography_runner(community)` — returns a closure
  `runner(birth_rates) -> offspring_production` that runs the SCM once. It
  stashes `last_schedule_times` / `last_offspring_production` / `history` in its
  environment for the cleanup step. Resets the integration schedule to defaults
  when birth rates jump by more than
  `demography_control$equilibrium_large_birth_rate_change`.
- `plant_community_demography_runner_cleanup(community, runner, converged)` —
  reads the runner's environment and writes the equilibrium `birth_rate`,
  schedule times, and `converged`/`progress` attrs back onto the community.

### Demography / equilibrium (`R/community_demography.R`)

`community_demography(community)` dispatches on
`demography_control$equilibrium_solver_name`:

- `single_step` — one SCM step, no iteration.
- `equilibrium_iteration` — **default, working**. Fixed-point iteration of
  incoming→outgoing offspring until `equilibrium_eps` is reached.
- `equilibrium_solve_nleqslv` / `equilibrium_solve_dfsane` — root-finding via
  `util_nlsolve` (`R/util_nlsolve.R`). Verified against a known fixed point and
  against the SCM; worth reaching for when the iteration converges slowly.
- `equilibrium_hybrid` — iterate then root-find, alternating solvers, rejecting
  any solution that drove a still-viable species extinct.

After solving, `plant_community_update_fitness_function()` builds the mutant
invasion-fitness closure on the community.

### Fitness, gradients, assembly

- `R/community_fitness_landscape.R`, `community_fitness_viable.R`,
  `community_fitness_solve_max.R` — invasion-fitness landscapes (some use
  `mlr3`/Gaussian-process surrogates) and viable trait bounds.
- `R/solve_attractors.R` — `community_selection_gradient()`,
  `community_solve_singularity_1D()`. Finite-difference gradients use the
  internal `gradient_points()`/`gradient_extrapolate()` in `R/util_gradient.R`.
- `R/singularity.R` — `community_solve_singularity()` (N-D root-find on the
  selection gradient) and `community_classify_singularity()` (CSS / branching
  point / repeller / Garden of Eden, with eigen-decompositions). Their
  second-order finite differences (`util_hessian()`, `util_jacobian()`) live
  beside the gradient helpers in `R/util_gradient.R`. See **Singular
  strategies** below.
- `R/assembler.R` — `assembler_start`/`assembler_run`/`assembler_control` drive
  full assembly (births → demography → deaths) over many steps.
- `R/births*.R`, `R/deaths.R` — add/remove strategies (maximum-fitness or
  stochastic births; inviable-strategy removal).

### Control objects — keep them straight

| Object | Built by | Purpose |
|---|---|---|
| `demography_control` | `demographic_step_control()` | Equilibrium solving: `equilibrium_solver_name`, `equilibrium_eps`, `equilibrium_nsteps`, `equilibrium_large_birth_rate_change`, `equilibrium_extinct_birth_rate`, etc. Lives at `community$demography_control`. |
| `plant_control` | `plant_default_assembly_control()` / `plant::control()` | Passed straight to `run_scm()`. A plant `Control` S4 object — **cannot** hold extra fields, so all equilibrium params go in `demography_control`. Lives at `community$model_support$plant_control`. |
| `fitness_control` | (list) | How fitness landscapes are sampled (`method`, `n_evals`, …). |
| `assembler_control` | `assembler_control()` | Assembly loop: birth/death type, tolerances. |

Note: plant renamed `seed_rain` → `birth_rate`/`offspring`; control field names
here follow current plant terminology.

## plant dependency and interface

- Built against the post-#459 `plant` (installed `2.0.0.9001`). Baseline recorded
  in `.plant-interface-version`.
- Source checkouts of plant live alongside this repo, e.g.
  `../plant-dev1` (see its `NEWS.md` "Breaking changes" for old→new maps) and
  `../plant-master`.
- Key migrated calls: `scm_base_control()`→`control()`;
  `run_scm_collect(p)`→`run_scm(p, collect=TRUE)`; `build_schedule(p, ctrl)` +
  `attr(., "offspring_production")` → `run_scm(p, ctrl, refine_schedule=TRUE)`
  then `scm$parameters` / `scm$offspring_production`;
  `p$node_schedule_ode_times` (Parameters field) → `p$ode_times`.
- `run_scm(..., collect=FALSE)` returns the **SCM object** (`scm$parameters`,
  `scm$offspring_production`, `scm$net_reproduction_ratios`, `scm$run_mutant`);
  with `collect=TRUE` it returns a tidied list.
- After future plant updates, re-run the `plant-update-interface` skill (in
  `../plant-dev1/.claude/skills/`) — it reads plant's NEWS and updates
  `.plant-interface-version`.

## Test baseline

`devtools::test()` is **green: 401 pass, 0 fail, 0 skip, 0 warn**. Tests run in
parallel (`Config/testthat/parallel: true`); the `test-plant-smoke*.R` files
dominate the wall-clock as they are the only ones that run the real SCM. The
`test-harness-*.R` and `test-singularity.R` files run no SCM and are fast.

(The count has grown as the toy-harness tier has: 197 → 256 → 401. What matters
is that a change moves it up and moves nothing to FAIL.)

Note: the testthat parallel workers may fail to find `plant` on startup in some
shells; run `TESTTHAT_PARALLEL=FALSE Rscript -e 'devtools::test()'` if so.

- `test-community.R` (new) — covers the community interface: `trait_matrix`,
  `bounds`, `demographic_step_control`, `community_start/add/drop`,
  `length.community`, the `max_patch_lifetime` schedule regression, and
  integration tests for `community_demography` (empty + single resident,
  reference birth rate ≈ 0.06846) and `community_selection_gradient`.
- `test-support-fitness.R` — `positive_1d`, `bounds`/`check_bounds`/`check_point`,
  and `community_viable_fitness_1D` (viable interval ≈ [0.0405, 0.7993] for lma).
- `test-solve-attractors.R` (new, #27) — `community_solve_singularity_1D`: 1D
  attractor (≈0.1417 for lma) plus the non-bracketing `edge_ok` warning/error
  branches.
- `test-fitness-landscape.R` (new, #27) — `community_fitness_landscape` grid
  method (resident flagged, fitness ~0 at the equilibrium resident, auto-solves
  demography, rejects unknown methods). The bayesopt/surrogate method is not
  covered.
- `test-assembler.R` (new, #27) — `assembler_control` defaults/validation,
  `mutational_vcv_proportion` (diagonal log-scale vcv), the maximum-fitness and
  stochastic assembly loops (`assembler_start`/`assembler_run`), and
  `tidy_assembly` output shape.
- `helper-assembly.R` (new) — shared `assembly_model_support(max_patch_lifetime
  = 30)` used by the integration tests (previously inlined in test-community.R).
- `test-singularity.R` — `community_solve_singularity` (1-D, 2-trait, both
  solvers, trait scales, edge/validation branches) and
  `community_classify_singularity`, both against the analytic oracles tabulated
  under **Singular strategies** below.
- `test-demography-solvers.R` — all five equilibrium solvers on DD99, plus the
  genuine fixed-point tests and the `equilibrium_hybrid` extinct-species
  accept/reject branches built on `helper-harness-map.R`.
- `helper-harness-map.R` — a test-only harness whose demography runner is an
  arbitrary map `n -> map(n)`. The shipped toy harnesses return their
  equilibrium analytically, so no solver ever iterates on them; this one has a
  real, tunably-slow fixed point, which is what actually tests the root finders.
- `test-community-plots.R` — `community_plot_fitness_landscape`, forcing
  `ggplot_build()` so the aesthetics are actually evaluated.
- `test-plant-smoke-singularity.R` — the SCM anchor for the above: the
  alternative equilibrium solvers agreeing with the iteration, and the N-D
  solver plus classifier running on the real model. Deliberately
  **reference-free** (internal consistency, not pinned trait values) so it
  survives changes to the plant parameterisation; the pinned references stay in
  `test-plant-smoke.R`.

The five old test files that called plant's removed API were resolved
(drop-or-implement, not skipped): `positive_1d` was reimplemented so its test was
kept; the fundamental/viable tests were rewritten against the community machinery;
`test-x_equilibriumR.R` (superseded by `test-community.R`), `test-support-assembly.R`
(`assembly_parameters`), `test-trait_fitness.R` (dead empty loop), and the
duplicate `test-fitness-support.R` were deleted.

## Known issues / TODO

- **2D maximum-fitness births (`find_max_fitness_2d`) are unverified** — the path
  is now self-contained (uses `sys$fitness_function`) but relies on
  `fitness_slopes`/`maximize_logspace` and is commented as "not very well tested".
  Only 1D fitmax births are exercised.
- `plant_community_check_for_inviable_strategies` still has a TODO to drop its
  direct plant dependency and reuse the community fitness functions.
- `equilibrium_extinct_birth_rate` (`demographic_step_control()`, default `1e-3`)
  is an **absolute** birth rate, and it now decides who counts as extinct in
  `equilibrium_hybrid` as well as in the inviable-strategy check. Its meaning
  depends entirely on the model's units, so it needs re-tuning when the plant
  parameterisation changes the scale of equilibrium birth rates.
- The bayesopt/surrogate fitness-landscape method is still uncovered by tests.

### Resolved (see "Singular strategies" below)

- ~~`equilibrium_hybrid` solver is broken~~ — reworked for the `community`
  object. It read `eq_solution$strategies` and called `run_scm()` on the result,
  so the extinct-species re-check could not run at all (and would have errored
  if it had). It now works off `community$birth_rate` and the harness
  connectors, so it is model-agnostic, and each attempt continues from the
  previous iteration rather than restarting.
- ~~`equilibrium_solve_*` (nleqslv/dfsane) are unverified~~ — both verified,
  no code changes needed. The DD99 tests only checked dispatch (the explicit
  harness returns its equilibrium directly, so no solver ever iterates);
  `helper-harness-map.R` adds a harness whose demography runner is a genuine
  map with a known fixed point, and `test-demography-solvers.R` now shows both
  solvers reaching it in the case they exist for — where `equilibrium_iteration`
  runs out of steps first. `test-plant-smoke-singularity.R` anchors all three
  alternative solvers against the real SCM.
- ~~`community_plot_fitness_landscape` is not working~~ — fixed. It required a
  `fitness_surrogate_function` (which only the bayesopt method creates) and read
  a column literally named `x` from `fitness_points`, whose first column is
  named after the trait. It now plots grid landscapes, computes one if absent,
  follows the community's trait scale, and uses the surrogate only when present.
- ~~`fitness_control` has no default~~ — `fitness_landscape_control()` supplies
  `method = "grid"`, `n_evals`, `n_init`, and `community_start()` normalises
  whatever it is given through it. `community_fitness_landscape()` also fills in
  defaults for communities built by hand.
- ~~Two malformed `@param` roxygen tags in `R/community_plots.R`~~ —
  `devtools::document()` now runs clean.

## Singular strategies: solving and classifying (`R/singularity.R`)

`community_solve_singularity_1D()` brackets the scalar selection gradient with
`uniroot` and is strictly one-trait. Two dimension-agnostic functions sit
alongside it; both go through the harness connectors only, so they run on the
toy harnesses exactly as on the plant SCM.

- **`community_solve_singularity(community, x0, bounds, solver, ...)`** —
  multivariate root-find on `community_selection_gradient()` via `util_nlsolve`
  (`nleqslv` or `dfsane`). Searches on the community's trait scale (for `"log"`
  traits the residual is the gradient w.r.t. `log(x)`, far better conditioned).
  Discards any residents on the way in — a singular point is monomorphic — and
  returns the community *at* the root, with `attr(., "singularity")`. Candidates
  are clamped to `bounds`; landing on a bound warns (or errors with
  `edge_ok = FALSE`), mirroring the 1-D solver.
- **`community_classify_singularity(community, ...)`** — the second-order
  conditions, covering 1-D and N-D with one code path (a 1-D result is just
  1x1 matrices). Returns a `singularity_classification` object:
  - `hessian` — curvature of invasion fitness in the *mutant* direction, with
    the resident held fixed. Negative definite = ESS. Computed by
    `util_hessian()` in one vectorised call to `fitness_function`
    (`1 + 4k^2` mutant evaluations, cheap).
  - `jacobian` — derivative of the selection gradient w.r.t. the *resident*.
    Eigenvalues with negative real parts = convergence stability; a negative
    definite symmetric part = *strong* convergence stability (any mutational
    covariance). Computed by `util_jacobian()`, `2k` full equilibrium solves —
    this dominates the cost.
  - the four-way `classification`: CSS / branching point / repeller / Garden of
    Eden, plus `degenerate` when an eigenvalue is within `tol` of zero.
  - the full eigen-decompositions and, where the point is invadable,
    `branching_direction` — the leading Hessian eigenvector, i.e. the direction
    in trait space the population splits along. In a multi-trait problem that
    direction *is* the result, so it is returned, not just the verdict.

Analytic oracles used in `test-singularity.R` (all reproduced to 1e-4):

| model | Hessian | Jacobian | verdict |
|---|---|---|---|
| DD99 | `r(1/σ_C² − 1/σ_K²)` | `−r/σ_K²` | branching iff `σ_C < σ_K` |
| GK98 (symmetric 3-patch) | `−1/σ² + 2d²/(3σ⁴)` | `−1/σ²` | branching iff `d/σ > √(3/2)` |
| JJ12 | `< 0` | `< 0` | CSS at `x* = x_opt − aσ²` |
| `harness_dd99_nd` | `diag(r(1/σ_C,d² − 1/σ_K,d²))` | `diag(−r/σ_K,d²)` | branches along the narrow-kernel axis |

The two-trait DD99 case is the one that matters for multi-trait work: with
`σ_C = (0.4, 1.5)`, `σ_K = (1, 1)` the Hessian is `diag(5.25, −0.556)` and the
branching direction is `(1, 0)` — the classifier picks out *which* trait
disruptive selection acts on.

## Model harnesses: plant vs fast toy models (#33)

Running the full `plant` SCM for every fitness/equilibrium evaluation is slow and
has no closed-form answer to validate algorithms against. The pipeline is now
**model-agnostic**: a community carries a `harness` (`community$harness`) that
wires it to a backend. The pipeline only ever calls six connectors (defined in
`R/harness.R`), each forwarding to the harness's own implementation:

| connector | builds |
|---|---|
| `community_parameters()` | model parameters |
| `community_make_demography_runner()` | closure `birth_rates -> offspring` |
| `community_demography_runner_cleanup()` | writes equilibrium state back |
| `community_viable_bounds()` | viable trait region (empty community) |
| `community_check_for_inviable_strategies()` | residents to drop |
| `community_update_fitness_function()` | the invasion-fitness closure |

`harness_plant(model, version)` (the default) forwards these to the existing
`plant_community_*` code, so callers passing only `model_support` are unchanged.
It is identified by the plant physiological model (`FF16`, the current default,
or `TF24`) and the plant package `version` — a hook for supporting different
plant models/interfaces side by side. (The plant `Parameters` themselves still
come from `community_start(model_support = ...)`; `model`/`version` are recorded
metadata for now.)

`harness_explicit(fitness, equilibrium, ...)` implements all six connectors
generically from two primitives — a vectorised **invasion-fitness** function and
an **equilibrium solve** — each backed by C++ in its own file under `src/`
(`DD99.cpp`, `GK98.cpp`, `GM99.cpp`, `JJ12.cpp`; the package's first C++, Rcpp via
`LinkingTo`, `@useDynLib` in `R/zzz.R`, `src/Makevars` C++17). "Explicit" is about
the mechanism (fitness/equilibrium computed directly, not via the SCM), **not** a
claim that every quantity is closed-form. For these backends `community$birth_rate`
holds resident **abundance/density**, not a plant offspring rate. Four models
ship (named by author/year), each with test oracles:

- **`harness_dd99()`** — Dieckmann & Doebeli 1999 Gaussian competition. `x* = x0`;
  branches iff `σ_C < σ_K`, else ESS. (fitness is a per-capita rate, ~0 at resident)
- **`harness_gk98()`** — Geritz, Kisdi, Meszéna & Metz 1998 soft-selection
  (Levene) model. `x* =` capacity-weighted mean optimum (0 for the symmetric
  3-patch default); branches iff `d/σ > √(3/2) ≈ 1.2247`.
- **`harness_gm99()`** — Geritz, van der Meijden & Metz 1999 seed-size safe-site
  model (Poisson sites, size-asymmetric seedling competition; the model in
  Daniel's MATLAB at
  `OneDrive/.../Offspring-SmithFretwellReview/models/Geritz/`). Only `α·R` and
  `β·R` matter; no closed-form `x*` (numerical) and fitness is a Poisson series;
  size-asymmetric competition drives branching — Fig. 5: `αR=4.5` gives a CSS,
  `αR=7` branches at `βR=15`. Strong oracle: equilibrium invariant `W_m(m)=1`.
  **Distinct from `harness_gk98`** (1998 soft-selection).
- **`harness_jj12()`** — Johansson & Jonzén 2012 migratory-bird arrival time
  (simplified analytic form of Brännström et al. 2013 §4). CSS, no branching;
  closed form `x* = x_opt − a·σ²`. The benchmark for testing numerical
  gradient/singularity solvers vs. an exact answer.

(fitness returns a log ratio for jj12/gk98/gm99; a per-capita rate for dd99 —
both ~0 at the resident.) Tests live in
`tests/testthat/test-harness-{dd99,gk98,gm99,jj12}.R`; they assert each singular
strategy is recovered by `community_solve_singularity_1D` and that the
invasion-fitness curvature flips sign at the ESS/branching boundary.

### Mixtures, trait scale, and multi-trait models

DD99/GK98/GM99 handle **mixtures**, so the full assembler runs on them, not just
PIPs (`test-assembly.R`): DD99 packs a limiting-similarity community, GK98 (seeded
— soft selection has no empty habitat to colonise) branches to a dimorphism, GM99
assembles a seed-size polymorphism (slowest — its multi-resident fitness is a
multi-dim Poisson sum).

- **`community_start(trait_scale = "log"/"linear")`** (`community_trait_transform()`)
  controls how trait space is spaced/searched during assembly. `"log"` (default)
  suits strictly-positive traits (plant, GM99 seed size); `"linear"` suits traits
  centred at 0 (DD99, GK98). Used by the fitness-landscape grid, the
  nearest-resident distance in births/`should_move`, and `find_max_fitness_2d`.
- **`harness_dd99_nd()`** is the multi-trait DD99 (cf. Ito & Dieckmann 2007;
  product-Gaussian kernels). The explicit harness is dimension-aware
  (`explicit_resident_traits()`: vector for 1-trait models, the trait matrix for
  nD). `assembler_run()` works in nD warning-free and assembles ~27 species
  across a 2D trait plane. The nD path was hardened (part of #9):
  `assembler_append_history()` skips the 1D fitness-landscape grid for nD;
  `find_max_fitness_2d()` is a real multistart hill-climb (from each resident +
  the centre; it no longer depends on the never-populated `fitness_slopes`); and
  `maximize_scaled()` insets boundary start points so residents sitting on a
  bound don't break `nmkb`.

See also `x_misc/Revolve/doc/models.md` and the per-model write-ups +
`assembly.qmd` in `overstorey_staging/`.

## Issue & project-board conventions

Development across `plant`, `regnans`, and `overstorey` is tracked on a
shared [project board](https://github.com/orgs/traitecoevo/projects/5). New issues
are auto-added to the board with **no Status** — that's the triage queue. A maintainer
sets Status (e.g. Backlog) during triage, so you don't need to set it yourself.

When opening an issue (including whenever the user asks you to create one), always:

- **Set exactly one type label.** Only three labels exist in these repos — do not
  invent new ones:
  - `bug` — an existing feature not functioning as intended
  - `task` — a discrete task needed for a feature (the default for normal work)
  - `epic` — a new feature or capability, usually an umbrella over several tasks
- **Prefix the title with a theme tag** in square brackets so the board sorts
  cleanly. Assembly work is usually `[evol assembly]`; reuse another existing
  theme where it fits, or fall back to `[other]`:

  | Tag | Scope |
  |---|---|
  | `[evol assembly]` | Evolutionary assembly linking plant to regnans |
  | `[TF24 hydraulics]` | Hydraulics component of the TF24 strategy |
  | `[TF24 allometry]` | Flexible allometry for the TF24 model |
  | `[TF24 nsc]` | Non-structural carbohydrate storage in TF24 |
  | `[acclimation]` | Acclimation of leaf and other traits |
  | `[simplify interface]` | Consistent interface to the plant & regnans models |
  | `[Env drivers]` | Driving the model with environmental drivers |
  | `[speed]` | Performance — making the model run faster |
  | `[patch variations]` | Multiple patch setups (multi-patch, stochastic metapopulation, continuous patch) |
  | `[forecasting]` | Enabling forecasting with the plant model |
  | `[documentation]` | Documenting model capabilities (any of the three repos) |
  | `[other]` | Anything not covered above |

  A title may carry more than one tag when it genuinely spans themes.

Create issues with `gh issue create -R traitecoevo/regnans
--title "[evol assembly] …" --label task` (swap in `bug`/`epic` as appropriate).


## Plant family

`regnans` is part of the **plant family** in the [`traitecoevo`](https://github.com/traitecoevo)
org — a hub-and-spoke set of packages built around the
[`plant`](https://github.com/traitecoevo/plant) size- and trait-structured forest model.

- **Docs hub** — family user guides & theory: <https://traitecoevo.github.io/overstorey/>
- **Cross-package orientation** — how the family fits together (who depends on whom,
  source-of-truth rules, cross-repo gotchas) lives in
  [`plant-meta`](https://github.com/traitecoevo/plant-meta); start with its
  [`AGENTS.md`](https://github.com/traitecoevo/plant-meta/blob/main/AGENTS.md). Keep
  family-wide concerns there, not here.
- **Issues & board** — follow the
  [issue guide](https://github.com/traitecoevo/plant-meta/blob/main/governance/issue-guide.md);
  work is tracked on [board #5](https://github.com/orgs/traitecoevo/projects/5) (new issues
  auto-add with no Status = the triage queue). Labels: `bug` / `task` / `epic` plus `blocked`,
  `needs-info`, `cross-package`, `breaking`, `question`.
- **Commit messages** — the repo squash-merges, so a PR's title and body are copied verbatim into
  permanent history. Keep them short and durable, and put the working detail (measurements,
  alternatives rejected, what you tried first) in the first PR comment instead — see
  [`commit-messages.md`](https://github.com/traitecoevo/plant-meta/blob/main/governance/commit-messages.md).
