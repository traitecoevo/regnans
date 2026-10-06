# regnans roadmap

What this package does, how each method stands, and where it is going. This is the one place that says what is *planned*; the detail of any item lives in its issue. Worked examples are on the docs hub under [Theory → Adaptive dynamics](https://traitecoevo.github.io/overstorey/theory/adaptive-dynamics.html); the developer map is [`AGENTS.md`](AGENTS.md); the function reference is `man/`.

Status words used below: **done** (implemented and tested), **partial** (works but untested, limited, or with a known defect), **planned** (not started), **blocked** (waiting on another package — see [What regnans needs from plant](#what-regnans-needs-from-plant)).

## Scope

`regnans` is a toolkit of adaptive-dynamics methods for models whose invasion fitness and demographic equilibrium are **computed numerically** — the `plant` size-structured forest model first, but any model that can be wired to the harness interface. The ambition is a reasonably comprehensive set of the standard adaptive-dynamics analyses (invasion fitness, equilibria, selection gradients, singular strategies and their classification, evolutionary dynamics, community assembly, bifurcation structure in parameter space), built so that they work on expensive numerical models rather than only on closed-form toy models. Closed-form reference models ship alongside, but as test oracles, not as the design target.

Out of scope: the genetic basis of traits; resident attractors that are not fixed points (limit cycles, stochastic quasi-stationary states), which would need monodromy or periodic-orbit sensitivities and are deferred, as in [#43](https://github.com/traitecoevo/regnans/issues/43).

## Design principles

These constrain future code. A change that breaks one of them should say so in its PR.

### 1. Model-agnostic through the harness

A community carries a `harness` (`R/harness.R`), and the pipeline reaches the model only through its connectors: `parameters`, `make_demography_runner`, `demography_runner_cleanup`, `viable_bounds`, `check_for_inviable_strategies`, `update_fitness_function`. New methods call connectors, never `plant::` directly. Two existing violations to retire: `plant_community_check_for_inviable_strategies()` and the plant-specific default start point in `community_viable_fitness_1D()`.

### 2. Derivatives are part of the contract; algorithms are designed for exact derivatives

Models in this family are gaining automatic differentiation (plant [PR 648](https://github.com/traitecoevo/plant/pull/648), [odelia](https://github.com/traitecoevo/odelia)). Algorithms here are chosen on the assumption that exact derivatives exist — Newton rather than fixed-point iteration, implicit rather than explicit steppers, contour continuation rather than grids, gradient-based rather than derivative-free optimisation — and are not held back to accommodate models that cannot differentiate themselves. A model without derivatives is wrapped in `harness_fd()`, which supplies the same connectors by finite differences: a compatibility and verification path, not a design target.

Decided shape of the interface (implementation is Phase 0 below):

- A harness declares `h$provides`, a subset of `fitness_gradient`, `fitness_jet2`, `fitness_parameter_gradient`, `demography_jacobian_n`, `demography_jacobian_x`, `demography_jacobian_p`. `harness_fd()` fills in whatever is missing; `harness_provides(h, what)` queries it.
- `update_fitness_function()` may also set `community$fitness_derivatives`, a list of closures `gradient(y)`, `jet2(y, v)`, `parameter_gradient(y, which)` over the cached resident environment. The environment never crosses the boundary.
- The demography runner takes an optional `derivatives = c("n", "x", "p")`; the returned offspring vector then carries `jacobian_n` (n × n), `jacobian_x` (n × nk) and `jacobian_p` attributes computed in the same run, because in forward mode the tangents ride along one integration. Requesting derivatives also means "freeze the discretisation" (no schedule refinement inside a derivative).
- Connector derivatives are with respect to **raw** traits, birth rates and parameters. Trait-scale chain rules (`community_trait_transform()`, including the Hessian's `x² s'' + x s'` term on a log scale) live in regnans.
- One dispatch file, `R/derivatives.R`, is the only place that decides where a derivative comes from: `derivative_control()`, `community_fitness_gradient()`, `community_fitness_hessian()`, `community_fitness_jet2()`, `community_fitness_parameter_gradient()`, `community_demography_jacobian()`, `community_equilibrium_sensitivity()` and `community_selection_gradient_jacobian()`. Algorithms call these and never see whether the answer came from AD, analytic code or `harness_fd()`.
- The equilibrium is a fixed point `n = f(n)`, so the implicit-function sensitivity is `dn*/dθ = (I − f_n)⁻¹ f_θ`. The harness supplies `f_n`; the `I −` and the residual transforms belong to regnans.
- Terminology: the *directional second derivative* of fitness along `v` is `vᵀ H v`. A "2-jet" is the triple `(s, v·∇s, vᵀHv)` — the Taylor coefficients to second order along one line — which nested forward-mode AD returns in one pass. This file uses the plain phrase and reserves "2-jet" for the plant specification, where it names the implementation.
- The first provider is DD99 with analytic derivatives in `src/DD99.cpp`, so dispatch, chain rules and the verification contract are proven before plant supplies anything.

### 3. Every model-supplied derivative is checkable against finite differences in CI

`harness_check_derivatives(community, n_points, tol_rel, tol_abs, which = h$provides)` (exported, in `R/derivatives.R`) draws random points on the trait scale inside the bounds, evaluates each advertised derivative both ways, compares column by column, returns a data frame and errors on the worst failure. Free cross-checks come with it: the gradient at `y = x` equals the selection gradient; the third component of a jet equals `vᵀHv`; `I − f_n` is non-singular at a converged equilibrium. `tests/testthat/test-derivatives-contract.R` runs it over every shipped harness that advertises anything. A harness author in another package opts in with one `test_that()` that calls the exported function (plant will need `tol_rel` near `1e-2` for ODE tolerances). This mirrors the "AD versus FD over the same reconstruction" contract in plant's AD work.

### 4. The plant boundary stays narrow

plant exposes fitness, its derivatives with respect to the mutant trait, the resident traits, the birth rates and selected parameters, and a warm-startable resident solve that returns its state. regnans owns Newton and root-finding, the implicit-function linear algebra, continuation, branching bookkeeping, plotting, all finite-difference fallbacks and the verification harness. If regnans ever needs to know what a light profile is, the boundary has leaked.

### 5. Derivatives of a discretised trajectory are not smooth near events

An AD gradient taken where a trait moves an event in the run (a cohort crossing `hmat` between two steps) is a derivative of step placement, not of the model; plant has measured one such cell at 211 times the median of its grid. Any method that root-finds or continues on gradients carries refusal and outlier diagnostics and periodic cold finite-difference verification from the start.

### 6. Reference models have analytic answers but are exercised as numerical models

DD99, GK98, GM99, JJ12 and the multi-trait DD99 keep their closed-form singular strategies, Hessians and Jacobians as *oracles*. The pipeline must also be able to run them the numerical way: equilibrium found by the iterative or Newton solvers on the dynamics map rather than the closed-form `equilibrium()` shortcut (today no solver ever iterates on them), derivatives obtained through the same connectors the plant harness uses — the toy C++ templated on its scalar type so that AD differentiates it exactly as it does plant, with `harness_fd()` as the second route — and the analytic values used only inside `expect_equal()`. Planned as a `mode = c("analytic", "numerical")` switch on `harness_explicit()`. DD99 is Lotka–Volterra with a Gaussian kernel, so it also has closed-form coexistence boundaries for the continuation work.

## Code discipline

How the methods above get built. These are the rules a reviewer holds a PR to; they are specific to the defects this codebase has accumulated and to the cost structure of numerical models.

1. **Small, additive, bit-identical steps.** One change per PR. A refactor of numerical code (the dispatch layer, a solver rewrite) must leave every existing test value unchanged; where a value is not already pinned, pin it in a fixture *before* the refactor, in its own PR. New capability is added beside old behaviour, not threaded through it. This is the discipline plant's AD work used to land 17k lines without moving a reference number.
2. **Work for the elegant solution.** A quick-and-dirty version that works today becomes the thing every later method has to route around. Spend the extra hour on the formulation — the right unknowns, the right residual, the right place for the state to live — before writing the loop; prefer the version with fewer concepts even when it is more work to reach. If the elegant version is not yet clear, write the issue, not the code.
3. **Every number has an oracle.** A new method is first shown to recover a closed-form answer on a reference harness (`expect_equal` to `1e-4` or tighter, with the tolerance justified), then shown internally consistent on plant in a reference-free smoke test (resident fitness ≈ 0, gradient vanishes at the root, two solvers agree). "Looks plausible" is not a test. Every derivative column, whatever its source, is checked against a finite difference. The test count goes up and nothing moves to FAIL.
4. **Fail loudly; never plausible-but-wrong.** Validate at the boundary with the existing `check_*` helpers and error with a sentence that names the fix. No silent fallbacks, no seeds set inside functions, no hard-coded workarounds in place of a control field, no `NaN` that propagates where a refusal object with a reason should be returned. A method that cannot answer says why.
5. **One access point, not variants.** Seek the integrated solution: a single entry point whose behaviour is set by control variables, never a family of functions that are variations on each other. The `_1D` / `_2d` / `_nd` siblings already in the package (`community_solve_singularity_1D()` beside `community_solve_singularity()`, `find_max_fitness_1D()` / `find_max_fitness_2d()`, `community_viable_fitness_1D()`) are the pattern to retire, not to extend: one function, dimension-agnostic, with the 1-D case falling out of the general one. Likewise the model is reached only through the harness connectors; derivatives only through `R/derivatives.R`; trait-scale transforms only through `community_trait_transform()`; nonlinear solves only through `util_nlsolve()`; equilibria only through `community_demography()`. A second path to any of these is a bug, not a convenience.
6. **No dead or half-alive code.** No `browser()`, no commented-out blocks, no no-op arguments, no "option 2" sketches, no helpers nothing calls. Delete it or turn it into an issue. A dependency is declared in `DESCRIPTION` in the commit that first uses it.
7. **Controls, not hidden state.** Every tunable is a documented, validated field of a `*_control()` constructor with a default. Functions take a community and return a modified community; state that must persist between calls (schedule, warm starts) lives in named fields of the object, not in closure environments beyond the one documented runner cache and not in globals.
8. **Cost is visible.** Every equilibrium solve and every model evaluation is counted on the community (`attr(., "n_evals")` or a progress field) so a test can assert the cost table, and so a regression in cost is caught like a regression in value. Nothing re-solves an equilibrium it already has; warm starts are carried, not recomputed. Unit tests run algorithms on the toy harnesses; plant smoke tests are small (`max_patch_lifetime = 30`) and reference-free.
9. **House style.** `snake_case`; `community_<verb>_<noun>()` for pipeline steps, `harness_<model>()` for backends, `util_*` for shared numerics; pipeline-friendly signatures (community first, community back). Roxygen on every export with a correct `@title`, a runnable `@examples` on a toy harness, and `devtools::document()` clean. Comments are written for their forward value only: they say *why* when the code cannot, at the density of the surrounding file. No numbers from one run, no retractions ("this used to…"), no history, no references to the session that wrote it — git and the issue tracker hold that. Before leaving a comment, ask what a reader in two years needs from it; if the answer is nothing, delete it. Tests live in the file named for the feature and are parallel-safe.
10. **Plan before code.** Anything beyond a one-line fix has an issue naming its oracle and acceptance criteria before implementation, following the family [issue guide](https://github.com/traitecoevo/plant-meta/blob/main/governance/issue-guide.md). The PR title and body are short and durable because the repo squash-merges; measurements, rejected alternatives and what was tried first go in the first PR comment. The PR updates this file's inventory row in the same change.

## What exact derivatives change

Three tiers answer "does AD just speed things up, or does it enable new methods?"

**Tier A — exists today on finite differences; exact derivatives make it exact and cheaper.** Selection gradient (`2k` mutant evaluations become one pass), mutant Hessian (`1 + 4k²` evaluations become one second-order pass), N-D singular-strategy Newton, ESS/CSS classification, the canonical equation in one or two traits, maximum-fitness births (gradient ascent replaces Nelder–Mead).

**Tier B — only practical with exact derivatives.** The common thread is a derivative *through* the demographic equilibrium, or many traits, or a second derivative of a nested solve.

- **Newton on the resident equilibrium** with `∂f/∂n`. Fixed-point iteration converges at the rate of the spectral radius of `f_n`, which approaches one as species interact strongly — exactly the multi-species communities on which plant has been slowest. Newton converges quadratically there. With warm starts and Broyden updates between exact refreshes, an equilibrium solve inside an assembly loop or a continuation drops from 10–50 SCM runs to roughly 4–8. This is the largest single multiplier in the package.
- **Joint eco-evolutionary Newton**: solve `n = f(n, x)` and `g(x, n) = 0` together in one Newton on `(n, x)` instead of nesting an equilibrium solve inside every gradient evaluation. The augmented-system idea from the continuation work, applied to singular strategies and coalitions.
- **Implicit differentiation of the equilibrium**, `dn*/dx = (I − f_n)⁻¹ f_x`: the primitive the rest of this tier is built on.
- **Stiff integration of the canonical equation** ([#43](https://github.com/traitecoevo/regnans/issues/43)): an implicit stepper needs `∂g/∂x` every step; by finite differences that is `2k` equilibrium solves per step.
- **Continuation and bifurcation analysis**: coexistence boundaries in parameter space, singular-strategy branches along an environmental parameter, loci of branching points. Pseudo-arclength Newton on an augmented system with nested second derivatives; with finite differences the failure mode is confidently tracing the wrong branch.
- **High-dimensional trait evolution**: 10–48 TF24 traits co-evolving. Reverse mode gives all partials in about one run; finite differences cost `k` runs per gradient and derivative-free search is hopeless.
- **Parameter sensitivity of evolutionary outcomes**, `dx*/dp` for every model parameter in one adjoint sweep.
- **Active-subspace reduction** of the fitness landscape from gradient covariances, and **gradient-enhanced emulators** (each gradient adds `k` observations per run to a surrogate).

**Tier C — sampling methods: derivatives help at the margin; parallelism and caching dominate.** Grid by default, derivative-accelerated where it counts. A pairwise invasibility plot is the zero-level set of `s(y, x)` in the (resident, mutant) plane, so with `∂s/∂y` and `∂s/∂x` its boundary curves can be continued from the diagonal instead of classifying a grid, using the same machinery as the coexistence boundaries. Viable bounds become a Newton solve on `s(y) = 0`, and the N-D viable region is the continuation of that contour. Landscapes gain gradient-enhanced surrogates. Stochastic assembly and the scouting stage of continuation stay derivative-free and embarrassingly parallel.

Approximate costs in SCM-run units `R` for `k` traits and `n` residents. A mutant evaluation against a recorded resident costs `m ≈ 0.02–0.1 R`; a fixed-point equilibrium solve `E_FD ≈ 10–50 R`; a Newton equilibrium solve with exact `∂f/∂n` and a warm start `E_AD ≈ 4–8 R` for `n = 1`.

| Consumer | Finite differences today | With first derivatives (`∂f/∂n`, `∂s/∂y`, `∂·/∂p`) | Plus second derivatives (`vᵀHv`, `∂f/∂x`, mixed terms) |
|---|---|---|---|
| Resident equilibrium, one solve | `E_FD` = 10–50 R, worse as `n` grows | `E_AD` ≈ 4–8 R | `E_AD` ≈ 4–8 R |
| Selection gradient, equilibrium given | `(1 + 2k) m` | ≈ `(1 + k/2) m`, exact | ≈ `(1 + k/2) m`, exact |
| Mutant Hessian | `(1 + 4k²) · 2m` | `(1 + 4k²) · 2m`, unchanged | `k(k+1)/2` directional passes, exact |
| N-D singular strategy, ≈ 6 Newton iterations | ≈ `6(1 + 2k) E_FD`; `k = 2`: ~900 R | ≈ `6(1 + 2k) E_AD`; `k = 2`: ~180 R | `6 E_AD` + linear solve; `k = 2`: ~40 R |
| Classification (ESS + convergence stability) | Hessian + `2k E_FD`; `k = 2`: ~120 R | Hessian + `2k E_AD`; `k = 2`: ~25 R | ≈ `1 E_AD`; ~6 R |
| Canonical equation, per implicit step | `E_FD + 2k E_FD` | `E_AD + 2k E_AD` (FD across residents, AD inside) | `E_AD` |
| Parameter sensitivity, per parameter | `2 E_FD` | `1 E_AD` with a parameter tangent | `1 E_AD` with a parameter tangent |
| Continuation corrector step, residual dimension `D ≈ n(k+1) + k + 1` | ≈ `2D E_FD`; `k = 1, n = 2`: ~300 R | ≈ `(1 + D/2) E_AD`; ~25 R | ≈ `(1 + D/2) E_AD`, mixed terms exact |

Reading the table: `∂f/∂n` (Newton equilibrium plus warm start) multiplies everything; `∂s/∂y` buys exactness more than time; the exact convergence-stability Jacobian `dg/dx = ∂g/∂x|_n + ∂g/∂n (I − f_n)⁻¹ f_x` is the *hardest* quantity because it needs mixed second derivatives through the resident environment, so finite differences across residents with AD inside is the long-lived path for classification, the canonical equation and continuation, with the exact version a late upgrade.

## Methods inventory

Tier letters refer to the section above. Tests are in `tests/testthat/`; docs are overstorey pages under `theory/adaptive-dynamics/`.

### Model interface

| Method | Functions | Status | Tier | Tests / docs | Tracking |
|---|---|---|---|---|---|
| Harness contract, six connectors | `harness_plant()`, `harness_explicit()`, `R/harness.R` | done | — | `test-harness-*.R` | [#33](https://github.com/traitecoevo/regnans/issues/33) |
| Reference models with analytic oracles | `harness_dd99()`, `harness_dd99_nd()`, `harness_gk98()`, `harness_gm99()`, `harness_jj12()` | done; every reference model supplies its fitness gradient and Hessian in closed form (GM99's from the same Poisson sum as its fitness) | — | `test-harness-{dd99,gk98,gm99,jj12}.R`, `test-singularity.R`; `DD99.qmd`, `GK98.qmd`, `GM99.qmd`, `JJ12.qmd` | — |
| Reference models run as numerical models (`mode = "numerical"`, scalar-templated C++) | `harness_explicit()` | planned | — | — | [#55](https://github.com/traitecoevo/regnans/issues/55) |
| Derivative connectors, `harness_fd()`, dispatch layer | `R/derivatives.R`, `derivative_control()`, `community_fitness_gradient()`, `community_fitness_hessian()`, `community_selection_gradient_jacobian()`, `harness_provides()`, `harness_fd()` | done for the mutant direction (gradient, Hessian); resident Jacobian is finite-differenced over the gradient until plant E3/E4; demography Jacobians arrive with the Newton solver | A enabler | `test-derivatives.R` | [#50](https://github.com/traitecoevo/regnans/issues/50) |
| Derivative verification contract | `harness_check_derivatives()` | done | — | `test-derivatives-contract.R` | [#50](https://github.com/traitecoevo/regnans/issues/50) |
| Warm-started resident solve (initial guess including environment) | runner / `model_support` | blocked | B | — | [plant#650](https://github.com/traitecoevo/plant/issues/650) (E5) |

### Resident demographic equilibrium

| Method | Functions | Status | Tier | Tests / docs | Tracking |
|---|---|---|---|---|---|
| Fixed-point iteration (default), `nleqslv`, `dfsane`, hybrid, single step | `community_demography()`, `demographic_step_control()`, `util_nlsolve()` | done | — | `test-demography-solvers.R`, `test-community.R`, `test-plant-smoke-singularity.R`; `assembly_fitmax.qmd` | — |
| Newton with exact `∂f/∂n`, warm start, Broyden updates between refreshes | `equilibrium_newton` solver | planned; exact version blocked | B | — | [#56](https://github.com/traitecoevo/regnans/issues/56); [plant#650](https://github.com/traitecoevo/plant/issues/650) (E3) |
| Equilibrium sensitivity `(I − f_n)⁻¹ f_θ` | `community_equilibrium_sensitivity()` | planned; exact version blocked | B | — | [#50](https://github.com/traitecoevo/regnans/issues/50); [plant#650](https://github.com/traitecoevo/plant/issues/650) (E3, E4, E7); [odelia#39](https://github.com/traitecoevo/odelia/issues/39) |
| Joint eco-evolutionary Newton on `(n, x)` | — | planned | B | — | [#53](https://github.com/traitecoevo/regnans/issues/53) |

### Invasion analysis of a resident community

| Method | Functions | Status | Tier | Tests / docs | Tracking |
|---|---|---|---|---|---|
| Invasion fitness of mutants in the resident environment | `community$fitness_function`, `max_growth_rate()` | done | — | every `test-harness-*.R`, `test-plant-smoke.R`; `adaptive-dynamics.qmd` | — |
| Fitness landscape, grid | `community_fitness_landscape()`, `community_plot_fitness_landscape()` | done | C | `test-fitness-landscape.R`, `test-community-plots.R`; `assembly_stochastic.qmd` | — |
| Fitness landscape, Bayesian-optimisation surrogate | `fitness_landscape_control(method = "bayesopt")` | partial: untested; `DiceKriging` and `nloptr` undeclared; ignores `trait_scale` and `bounds`; sets a seed internally | C | — | [#27](https://github.com/traitecoevo/regnans/issues/27) |
| Gradient-enhanced surrogates | — | planned | B | — | [#50](https://github.com/traitecoevo/regnans/issues/50) (downstream) |
| Maximum of fitness within bounds | `max_fitness()` | done in 1-D; N-D untested | A | `test-plant-smoke.R` | [#27](https://github.com/traitecoevo/regnans/issues/27) |
| Viable trait bounds, 1-D | `community_viable_fitness_1D()`, `community_viable_bounds()` | done | C | `test-support-fitness.R`, `test-plant-smoke.R` | — |
| Viable trait region, N-D | — | planned (currently errors) | C | — | [#58](https://github.com/traitecoevo/regnans/issues/58) |
| Pairwise invasibility plots | `community_pip()` | planned (overstorey hand-rolls one per model page) | C | `DD99.qmd`, `GK98.qmd`, `GM99.qmd` | [#51](https://github.com/traitecoevo/regnans/issues/51) |
| Mutual invasibility and trait-evolution plots (two residents) | — | planned | C | — | [#51](https://github.com/traitecoevo/regnans/issues/51) |

### Selection gradients and singular strategies

| Method | Functions | Status | Tier | Tests / docs | Tracking |
|---|---|---|---|---|---|
| Selection gradient, central finite difference | `community_selection_gradient()`, `R/util_gradient.R` | done; `log_scale` argument is a no-op | A | `test-community.R`, `test-harness-jj12.R`, `test-plant-smoke.R`; `solving_attractors.qmd` | — |
| Singular strategy, 1-D bracket | `community_solve_singularity_1D()` | done | A | `test-solve-attractors.R`; `solving_attractors.qmd` | — |
| Singular strategy, N-D root-find on the trait scale | `community_solve_singularity()` | done | A | `test-singularity.R`, `test-plant-smoke-singularity.R` | [#49](https://github.com/traitecoevo/regnans/pull/49) |
| Classification: CSS, branching point, repeller, Garden of Eden; branching direction | `community_classify_singularity()`, `util_hessian()`, `util_jacobian()` | done | A | `test-singularity.R`, `test-plant-smoke-singularity.R` | [#49](https://github.com/traitecoevo/regnans/pull/49) |
| Polymorphic singular coalitions (two or more coexisting residents) | — | planned (the N-D solver discards residents) | B | — | [#53](https://github.com/traitecoevo/regnans/issues/53) |
| Parameter sensitivity `dx*/dp` and continuation of `x*` along one parameter | — | planned | B | — | [#54](https://github.com/traitecoevo/regnans/issues/54) |

### Evolutionary dynamics and community assembly

| Method | Functions | Status | Tier | Tests / docs | Tracking |
|---|---|---|---|---|---|
| Maximum-fitness assembly | `assembler_start()`, `assembler_run()`, `birth_type = "maximum"` | done in 1-D; 2-D multistart partial; live `browser()` at `R/births_maximum.R:15` | A | `test-assembler.R`, `test-assembly.R`; `assembly_fitmax.qmd`, `assembly.qmd` | [#57](https://github.com/traitecoevo/regnans/issues/57) |
| Stochastic assembly (mutation and immigration) | `birth_type = "stochastic"`, `mutational_vcv_proportion()` | done | C | `test-assembler.R`; `assembly_stochastic.qmd` | — |
| Deaths and inviable-strategy removal | `community_deaths()`, `check_for_inviable_strategies` connector | done; plant-coupled check to retire | — | `test-inviable.R` | — |
| Canonical equation, monomorphic, with branching detection and dimorphic continuation | — | planned | A then B | — | [#52](https://github.com/traitecoevo/regnans/issues/52) |
| Canonical equation, stiff / implicit, multi-trait | — | blocked | B | — | [#43](https://github.com/traitecoevo/regnans/issues/43), [odelia#35](https://github.com/traitecoevo/odelia/issues/35) |
| Tidy output and plots | `tidy_assembly()`, `plot_community()`, `plot_community_2d()` | partial | — | `test-assembler.R` | [#37](https://github.com/traitecoevo/regnans/issues/37) |

### Parameter-space structure

| Method | Functions | Status | Tier | Tests / docs | Tracking |
|---|---|---|---|---|---|
| Coexistence boundaries by continuation: scouting, ray bisection, pseudo-arclength on the augmented system, test functions and branch switching | — | planned; stage 1 (scouting + bisection) needs no derivatives, stages 2–5 blocked on plant E1–E6 | B | — | [#60](https://github.com/traitecoevo/regnans/issues/60) |
| Loci of branching points and ESS boundaries in parameter space | — | planned | B | — | [#60](https://github.com/traitecoevo/regnans/issues/60) |

### Acceleration

| Method | Functions | Status | Tier | Tests / docs | Tracking |
|---|---|---|---|---|---|
| Parallel evaluation of grids, scouting and independent branches (`future` / `mirai`) | — | planned | C | — | [#59](https://github.com/traitecoevo/regnans/issues/59) |
| Caching of plant runs via `logpile` | none in regnans | planned, in plant | C | — | [plant#651](https://github.com/traitecoevo/plant/issues/651) (E9) |
| Emulators | see landscapes above | partial | C | `scripts/gps/` survey | [#27](https://github.com/traitecoevo/regnans/issues/27) |

Caching is deliberately invisible here. plant gains `enable_logpile(path)`, after which `run_scm()`, the mutant path and the derivative endpoints look results up by the fingerprint of their inputs, so every grid, scouting pass, restart and repeated classification in regnans stops re-running identical SCMs without regnans knowing a cache exists. regnans's only obligation is that every plant call it makes is keyable: deterministic inputs, an explicit control object, a frozen schedule when caching, and no hidden state in the runner closure.

## What regnans needs from plant

plant [PR 648](https://github.com/traitecoevo/plant/pull/648) (reverse-mode AD of a stand's census by trait) supersedes [PR 553](https://github.com/traitecoevo/plant/pull/553) (the earlier prototype) and plant will add the endpoints regnans specifies, tracked in [plant#650](https://github.com/traitecoevo/plant/issues/650) (E1–E8) and [plant#651](https://github.com/traitecoevo/plant/issues/651) (E9). This section is that specification, in 648's vocabulary: run a resident with `record_trajectory = TRUE`, then ask for derivatives; results have the shape `list(value, gradient, refusal, control)`; columns are named `"<species>.<parameter>"`; a refused metric is a `NaN` row with a reason; `stand_gradient_compare()` refuses gradients taken at different `Control` values.

Notation: `n` residents with `k` traits each, `b_i` the resident birth rates, `R0_i` the net reproduction ratio of species `i`, `s(y)` the invasion fitness of a rare mutant `y` in the recorded resident environment.

Ordered by what unblocks most. E1, E2 and E3 are equal first priority.

- **E1 — Fitness as a differentiable metric.** `R0_i` and `offspring_production_i` for every resident species, for **FF16 and TF24**, with the same value/gradient/refusal/control shape as the census metrics. In 648 today the census metrics are `leaf_area`, `mass_above_ground` and `area_stem`, TF24 only, and FF16 declares none. PR 553 had `offspring_production_gradient()` for FF16, TF24 and TF24f: port the quantity, keep 648's discipline. Consumers: everything below.
- **E2 — Mutant gradient `∂s(y)/∂y` with the environment frozen.** The rare-mutant derivative: canopy held at the recorded resident, one mutant cohort integrated on the active scalar (`run_mutant()` templated). Shape `k` per mutant. Must accept a *set* of mutants against one recorded resident without re-running the resident. Consumers: the selection gradient at `y = x`, maximum-fitness births, N-D singular-strategy Newton, gradient-enhanced landscapes. Tier A: exact, and one mutant run instead of `2k`. Absent in 648 today.
- **E3 — Birth-rate axis `∂R0_i/∂b_j`**, `n × n`, as a differentiable column `"<species>.birth_rate"`. This is `f_n`, the Jacobian of the offspring map the equilibrium solvers iterate on. Consumers: the Newton equilibrium solve (the fix for slow multi-species communities), the implicit-function sensitivity `(I − f_n)⁻¹ f_x`, joint eco-evolutionary Newton. Largest multiplier in the cost table. Absent in 648; PR 553 had `birth_rate_gradient()` for FF16 with cross-species `R0` deferred.
- **E4 — Resident-feedback gradient `∂R0_i/∂x_j`**, `n × nk`: every species' trait re-shades the canopy (553's `feedback = "resident"`; 648's trajectory term applied to the fitness metric). This is `f_x`. Consumers: implicit-function sensitivity, the Jacobian of the selection gradient with respect to the resident (convergence stability), polymorphic coalitions, the canonical-equation Jacobian. Known gap to document: the frozen-grid total differed from a fully adaptive finite difference by 16–30 % on geometry traits in 553's measurements; regnans will treat the AD object as the definition and verify against a frozen-schedule finite difference.
- **E5 — Warm-startable resident solve.** `run_scm()` or a sibling that accepts an initial guess — birth rates, node schedule with ODE times and step sizes, and the environment (light-profile spline) — and returns the converged state in the same form, so that a Newton or continuation step costs about a tenth of a cold run. Whether the environment enters the unknown vector or stays an inner fixed point is plant's design decision; the continuation work recommends the former.
- **E6 — Second derivative in the mutant direction.** `∂²s/∂y²` as a `k × k` Hessian, or at minimum the directional second derivative along a supplied direction `d`: the 2-jet `(s, d·∇s, dᵀHd)` from nested forward mode in one pass. One jet suffices for one trait; `k` jets for `k` traits. Consumers: the ESS test in `community_classify_singularity()`, branching detection in the canonical equation, the `∂²r/∂y² → 0` test function in continuation. This is the `active<active<double>>` type of [#43](https://github.com/traitecoevo/regnans/issues/43). 648 is first order only.
- **E7 — Non-strategy parameter axis `∂R0_i/∂p`, `∂s/∂p`** for selected `Parameters` and `Control` entries (disturbance interval, patch lifetime, …) as differentiable columns. Consumers: parameter sensitivity `dx*/dp`, continuation in parameter space. 648 differentiates strategy parameters only. Lowest priority of the first-order items; finite differences on `p` are acceptable meanwhile.
- **E6b — Mixed second derivatives through the resident environment**, `∂(∂s/∂y)/∂x_j` and `∂(∂s/∂y)/∂b_j`. Needed only for the *exact* convergence-stability Jacobian. Lowest priority of all: finite differences across residents with E2 inside cover it. Listed so plant knows the ceiling.
- **E8 — Discipline, applying to all of the above.** Keep 648's refusal with a reason, the recorded and compared `control`, and a finite-difference check callable from R so regnans's CI can verify every column. Add an `events` hook to force a step at a known crossing (the `hmat` outlier case) so an outlying gradient can be resolved rather than discarded. Derivative runs must be on a **frozen discretisation**: schedule refinement is a discrete branch and must not sit inside any derivative, AD or FD; the API takes the schedule from the recorded resident. Mutant derivatives (E2, E6) reuse one recorded resident across many mutants.
- **E9 — `enable_logpile(path)`.** A setup call that points plant at a local `logpile` and makes `run_scm()`, the mutant path and the derivative endpoints cache-aware: fingerprint = parameters, strategies, birth rates, control, schedule, plant version and which derivatives were requested; hits return the stored result including gradient, refusal and control; misses run and store. Off by default, no behaviour change when off, and invisible to regnans, which never references `logpile`.

What stays in regnans: Newton and root-finding, the implicit-function linear algebra, continuation, branching bookkeeping, plotting, all finite-difference fallbacks, the verification harness. regnans never needs the light profile.

Related: [plant#472](https://github.com/traitecoevo/plant/issues/472) scope B, [plant#537](https://github.com/traitecoevo/plant/issues/537), [odelia#35](https://github.com/traitecoevo/odelia/issues/35), [odelia#39](https://github.com/traitecoevo/odelia/issues/39), [regnans#43](https://github.com/traitecoevo/regnans/issues/43). Open for plant to decide: whether E3 and E5 are built in plant or on odelia's steady-state work.

## Sequencing

**Phase 0 — finite differences only; lands now, gets faster later.**

1. ([#50](https://github.com/traitecoevo/regnans/issues/50)) `R/derivatives.R` dispatch and `derivative_control()`; re-point `community_selection_gradient()`, `community_classify_singularity()` and `community_solve_singularity()` (which today pays for two finite-difference Jacobians — pass `jac =` to `nleqslv`). A pure refactor; the existing tests are the regression suite.
2. `h$provides`, the `fitness_derivatives` slot, the runner's `derivatives =` convention, `harness_fd()`.
3. ([#50](https://github.com/traitecoevo/regnans/issues/50), [#55](https://github.com/traitecoevo/regnans/issues/55)) DD99 analytic derivatives as the first provider; `harness_check_derivatives()` and its contract test.
4. ([#51](https://github.com/traitecoevo/regnans/issues/51)) `community_pip()` and its plot; the two-resident mutual-invasibility plot.
5. ([#52](https://github.com/traitecoevo/regnans/issues/52), [#56](https://github.com/traitecoevo/regnans/issues/56)) The monomorphic canonical equation (explicit or `deSolve` stiff stepping with `community_selection_gradient_jacobian()` as the Jacobian), branching detection through the classifier, then dimorphic continuation. The `equilibrium_newton` solver, designed now with a Broyden/finite-difference Jacobian.
6. ([#53](https://github.com/traitecoevo/regnans/issues/53)) Polymorphic singular coalitions.
7. ([#54](https://github.com/traitecoevo/regnans/issues/54)) Parameter sensitivity and continuation of `x*` along one parameter.
8. ([#60](https://github.com/traitecoevo/regnans/issues/60)) Coexistence boundaries, stage 1 (scouting and bisection, no derivatives) and the augmented residual with a finite-difference Jacobian; gate on the measured cost of a warm-started solve.
9. Housekeeping, listed below ([#57](https://github.com/traitecoevo/regnans/issues/57), [#58](https://github.com/traitecoevo/regnans/issues/58), [#59](https://github.com/traitecoevo/regnans/issues/59)).

**Phase 1 — plant E1 and E2 land.** `harness_plant()` sets `fitness_derivatives$gradient` and `provides`; every gradient-based method becomes exact with no further change here; plant's CI runs the contract test.

**Phase 2 — plant E3, then E7.** `jacobian_n` from the runner; Newton equilibrium on by default; exact `community_equilibrium_sensitivity()`. The big multiplier.

**Phase 3 — plant E4, E6, E6b.** Directional second derivatives, `jacobian_x` and the mixed terms; exact `dg/dx`; the stiff canonical equation of [#43](https://github.com/traitecoevo/regnans/issues/43); continuation stages 2–5; high-dimensional trait evolution.

## Housekeeping

Known defects and dead ends, to clear in Phase 0:

- A live `browser()` at `R/births_maximum.R:15`, reached when a community has no `bounds` ([#57](https://github.com/traitecoevo/regnans/issues/57)).
- The Bayesian-optimisation landscape uses `DiceKriging` and `nloptr` without declaring them, ignores `trait_scale` and the `bounds` argument, and calls `set.seed(1)` internally ([#57](https://github.com/traitecoevo/regnans/issues/57)); it has no tests ([#27](https://github.com/traitecoevo/regnans/issues/27)).
- `plant_community_check_for_inviable_strategies()` reaches plant directly and hard-codes `eps_test`; `community_viable_fitness_1D()` takes its default start from `model_support$p`.
- `equilibrium_extinct_birth_rate` is an absolute threshold whose meaning depends on the model's units.
- Dead code in `R/util.R` (`maximize_logspace`, `closest_log`, `rescale`, `unrescale`, `norm2`) and commented-out surrogate plotting in `R/community_fitness_landscape.R`.
- Copy-paste `@title`s on `plot_community()` and `demographic_step_control()`.
- Dimension-specific siblings to fold into one dimension-agnostic entry point each (discipline rule 5, [#58](https://github.com/traitecoevo/regnans/issues/58)): `community_solve_singularity_1D()` into `community_solve_singularity()`, `find_max_fitness_1D()` / `find_max_fitness_2d()` into one maximiser, `community_viable_fitness_1D()` into the N-D viable region.
- Plot functions ([#37](https://github.com/traitecoevo/regnans/issues/37)) and the untested exports listed in [#27](https://github.com/traitecoevo/regnans/issues/27).

## Maintaining this file

Update the status of a row in the pull request that changes it; do not let a finished method sit as "planned". Working detail — measurements, alternatives rejected, what was tried first — goes in the issue or the first PR comment, not here. If a planned item is dropped, delete its row rather than leaving it to read as intent.
