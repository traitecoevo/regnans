# {regnans}: Modelling community assembly with the {plant} model

<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/traitecoevo/regnans/workflows/R-CMD-check/badge.svg)](https://github.com/traitecoevo/regnans/master) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/traitecoevo/regnans/branch/master/graph/badge.svg)](https://codecov.io/gh/traitecoevo/plant?branch=master) -->
<!-- badges: end -->

The [{plant}](https://github.com/traitecoevo/plant) package for R is an extensible framework for modelling size-, patch- and trait-structured demography. The [{regnans}](https://github.com/traitecoevo/regnans) package extends {plant} to model community assembly and trait evolution, by repeatedly invading trait space, solving each community to demographic equilibrium, and following selection gradients to evolutionary attractors.

> **Status:** under active development. The fitness and equilibrium machinery was
> recently moved out of {plant} into this package, so a compatible recent {plant}
> is required (see `.plant-interface-version`). Developer notes and build/test
> commands live in [`AGENTS.md`](AGENTS.md).

## What it does

- **Resident demographic equilibrium** — `community_demography()`, with fixed-point iteration, `nleqslv`/`dfsane` root-finding and a hybrid.
- **Invasion fitness and landscapes** — `community_fitness_landscape()`, `max_fitness()`, viable trait bounds (`community_viable_fitness_1D()`).
- **Selection gradients and singular strategies** — `community_selection_gradient()`, `community_solve_singularity_1D()`, the multi-trait `community_solve_singularity()`, and `community_classify_singularity()` (CSS / branching point / repeller / Garden of Eden, with the branching direction).
- **Community assembly** — `assembler_run()` with maximum-fitness or stochastic births, and inviable-strategy deaths.
- **Model harnesses** — `harness_plant()` for the `plant` SCM, plus fast reference models with analytic answers (`harness_dd99()`, `harness_gk98()`, `harness_gm99()`, `harness_jj12()`) for learning the workflow and testing the algorithms.

Worked examples are on the docs hub under [Theory → Adaptive dynamics](https://traitecoevo.github.io/overstorey/theory/adaptive-dynamics.html).

## Where it is going

The aim is a reasonably comprehensive toolkit of adaptive-dynamics methods for *numerically solved* models, designed around the exact derivatives that `plant` is gaining through automatic differentiation. Next up:

- pairwise and mutual invasibility plots;
- the canonical equation of adaptive dynamics with branching detection, then its stiff multi-trait form;
- Newton equilibrium solves and implicit sensitivities once `plant` exposes the Jacobians;
- continuation of singular strategies and coexistence boundaries through parameter space.

[`ROADMAP.md`](ROADMAP.md) has the full inventory with the status of every method, what the package needs from `plant`, and the order of work.

## The plant model family

`{regnans}` is one of three companion repositories, with work coordinated
on a shared [project board](https://github.com/orgs/traitecoevo/projects/5)
("Plant model development"):

| Repository | Role |
|---|---|
| [**plant**](https://github.com/traitecoevo/plant) | Core C++/R model: size- and trait-structured demography, the SCM solver, and the physiological strategies (FF16, TF24, …). |
| [**regnans**](https://github.com/traitecoevo/regnans) | Evolutionary community assembly on top of `plant` — invasion fitness, demographic equilibria, and selection gradients (this repo). |
| [**overstorey**](https://github.com/traitecoevo/overstorey) | The narrative documentation / field guide site (user guides, theory, worked reproductions). |

Issues from all three repositories feed into the
[project board](https://github.com/orgs/traitecoevo/projects/5), which is the
single place to see what is planned, in progress, or done across the family.

## Reporting issues

Bug reports and feature requests are welcome via the
[GitHub issue tracker](https://github.com/traitecoevo/regnans/issues). New
issues are automatically added to the
[project board](https://github.com/orgs/traitecoevo/projects/5) with status
**Backlog**.

To keep the board sortable, please:

1. **Apply one type label** — the three repositories share the same set:
   - `bug` — an existing feature not functioning as intended
   - `task` — a discrete piece of work needed for a feature
   - `epic` — a new feature or capability, usually an umbrella over several tasks
2. **Prefix the title with a theme tag** in square brackets. Assembly work is
   usually `[evol assembly]`; use another existing theme where it fits, or
   `[other]`:

   `[TF24 hydraulics]` · `[TF24 allometry]` · `[TF24 nsc]` · `[acclimation]` ·
   `[simplify interface]` · `[evol assembly]` · `[Env drivers]` · `[speed]` ·
   `[patch variations]` · `[documentation]` · `[other]`

## Installation

**Requirements**

- You must be using R 4.5.0 or newer. At this stage the package is not on CRAN. Your options for installing are described below.

**Option 1, using `remotes::install_github`**

The `{regnans}` package can be installed direct from GitHub using the [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html) package:

```r
remotes::install_github("traitecoevo/regnans", dependencies=TRUE)
```

To install a specific (older) release, decide for the version number that you want to install in https://github.com/traitecoevo/regnans/releases  e.g.

```r
remotes::install_github("traitecoevo/regnans@v1.0.0", dependencies=TRUE)
```

with `"v1.0.0"` replaced by the appropriate version number.

**Option 2, building from source**

If familiar with [git](https://git-scm.com/) you might find it easiest to build `regnans` directly from the source code. This is most useful if developing new models or strategies, or to contribute new features.

First, clone the `regnans` repository

```
git clone https://github.com/traitecoevo/regnans
```

Open an R session in the folder, then to install dependencies run

```
devtools::install_deps()
```

Then to compile the project

```
devtools::install()
```
or 

```
devtools::load_all()
```

## Plant family

`regnans` is part of the **plant family** of packages in the
[`traitecoevo`](https://github.com/traitecoevo) org, built around the
[`plant`](https://github.com/traitecoevo/plant) forest model. Docs hub:
<https://traitecoevo.github.io/overstorey/>.

**Contributing:** please skim the family
[issue guide](https://github.com/traitecoevo/plant-meta/blob/main/governance/issue-guide.md)
before filing — issues across the family are triaged on
[board #5](https://github.com/orgs/traitecoevo/projects/5), and cross-package context lives in
[`plant-meta`](https://github.com/traitecoevo/plant-meta).
