# regnans (development version)

## Changes to results

* Finite differences across equilibrium solves (the resident Jacobian and the parameter Jacobian, and through them classification, the singularity solver, sensitivity and continuation) solve each stencil point to the cube of the stencil's relative step, or to `equilibrium_eps` where that is tighter, so their error is set by the step rather than by `equilibrium_eps` over the step (#80). On iterated equilibria they move towards their closed forms (DD99 from 3e-4 to 5e-7, GM99 from 5e-2 to 4e-6 relative); a stencil point whose solve stops short of its tolerance warns.
* The canonical equation's RODAS stepper takes the right-hand side's Jacobian from `community_selection_gradient_jacobian()` instead of odelia's forward difference, which was too noisy for its step-size control on an iterated equilibrium: such a run now costs what it does with the model's own equilibrium (120-238 solves on the reference cases, not 427-2178).
* Of two residents in a canonical run within two finite-difference steps of each other, the rarer goes extinct, recorded as a `"merge"` event: the converse of a branching.

## New features

* `community_selection_gradient_jacobian()` returns the Jacobian of the residents' equilibrium birth rates, from the same solves, in `attr(., "birth_rate_jacobian")`.

## Dependency changes

* Migrated to plant `develop@8657b3a` (traitecoevo/plant#648), which removed `run_scm(use_ode_times =)` and `Control$save_RK45_cache`, and now replays any ODE schedule a `Parameters` object carries. regnans' equilibrium runner now integrates freely at each birth rate rather than replaying the previous iteration's steps. The resident run behind mutant fitness replays the equilibrium run exactly (times and step sizes), so resident fitness at equilibrium is 0 again. The plant interface baseline is recorded in `.plant-interface-version`.

# regnans 0.2

* Earlier changes are recorded in the git history and the pull requests on GitHub.
