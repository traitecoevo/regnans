# regnans (development version)

## Changes to results

* Finite differences across equilibrium solves (the resident Jacobian and the parameter Jacobian, and through them classification, the singularity solver, sensitivity and continuation) solve each stencil point to the cube of the stencil's relative step, or to `equilibrium_eps` where that is tighter, so their error is set by the step rather than by `equilibrium_eps` over the step (#80). On iterated equilibria they move towards their closed forms (DD99 from 3e-4 to 5e-7, GM99 from 5e-2 to 4e-6 relative); a stencil point whose solve stops short of its tolerance warns.
* The canonical equation's RODAS stepper takes the right-hand side's Jacobian from `community_selection_gradient_jacobian()` instead of odelia's forward difference, which was too noisy for its step-size control on an iterated equilibrium: such a run now costs what it does with the model's own equilibrium (120-238 solves on the reference cases, not 427-2178).
* Of two residents in a canonical run within two finite-difference steps of each other, the rarer goes extinct, recorded as a `"merge"` event: the converse of a branching.
* `canonical_control()` defaults to `rtol = 1e-2` and `atol = 1e-4` (were 1e-3 and 1e-5): a run's end point is found by Newton, so the tighter trajectory cost RODAS 25-30% more equilibrium solves without moving it. The time at which a run reports reaching a stable point can move by up to one step.

## New features

* `community_label()` labels a community by the number of species its parameters support: the evolutionarily stable community reached by alternating evolution of the residents to their singular coalition and immigration at the maximum of invasion fitness, with a status saying whether the count is the model's (`stable`) or censored by the search (`capped`, `edge`, `unresolved`). A label needs bounds finite on the trait scale. `boundary_scout()` labels the points of a parameter design through the parallel map, its columns matched to the map's parameters by name, `boundary_design()` builds one, and `boundary_control()` holds their settings (#60).
* `boundary_bisect()` finds where the coexistence label changes along rays in parameter space, each ray bisected on a worker; `boundary_sweep()` labels a path cold and carried in both directions and flags points where the model's alternative attractors make the history decide, every route searching the cold label's box and an unresolved route leaving the flag `NA` (#60).
* `community_solve_boundary()` solves for a point on a coexistence boundary: the residents at their singular coalition and either a new fitness peak `y*` exactly marginal (`type = "invasion"`) or a resident becoming a branching point (`"branching"`), one parameter free, by Newton with a Jacobian from the derivative layer (#60).
* `community_continue_boundary()` follows a coexistence boundary through two parameters by pseudo-arclength continuation, rejecting steps whose tangent turns, monitoring test functions (a resident going extinct or ceasing to be an ESS, the marginal peak annihilating, residents merging, another peak becoming marginal, a crossing branch) and ending at the first event, located by false position, then verifying points cold against labels either side. It returns a `regnans_branch` with `print()`, `plot()` and `as.data.frame()` (cubic Hermite between points); `boundary_switch()` gives the start of the branch an event seeds, and `boundary_diagram()` collects branches (#60).
* `assembler_control()` gains `evolve` (`"none"`, `"newton"` or `"canonical"`): the residents evolve to their singular coalition between maximum-fitness births, so the assembly ends at an evolutionarily stable community; `max_residents` and `edge_tol` refuse a birth beyond a cap or on the edge of the bounds; `branch_offset` places the daughter of a resident at a branching point. An assembler counts its equilibrium solves in `$evaluations`; one saved before these fields restores with their defaults.
* `community_selection_gradient_jacobian()` returns the Jacobian of the residents' equilibrium birth rates, from the same solves, in `attr(., "birth_rate_jacobian")`.

## Dependency changes

* Migrated to plant `develop@8657b3a` (traitecoevo/plant#648), which removed `run_scm(use_ode_times =)` and `Control$save_RK45_cache`, and now replays any ODE schedule a `Parameters` object carries. regnans' equilibrium runner now integrates freely at each birth rate rather than replaying the previous iteration's steps. The resident run behind mutant fitness replays the equilibrium run exactly (times and step sizes), so resident fitness at equilibrium is 0 again. The plant interface baseline is recorded in `.plant-interface-version`.

# regnans 0.2

* Earlier changes are recorded in the git history and the pull requests on GitHub.
