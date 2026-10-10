# regnans (development version)

## Dependency changes

* Migrated to plant `develop@8657b3a` (traitecoevo/plant#648), which removed `run_scm(use_ode_times =)` and `Control$save_RK45_cache`, and now replays any ODE schedule a `Parameters` object carries. regnans' equilibrium runner now integrates freely at each birth rate rather than replaying the previous iteration's steps. The resident run behind mutant fitness replays the equilibrium run exactly (times and step sizes), so resident fitness at equilibrium is 0 again. The plant interface baseline is recorded in `.plant-interface-version`.

# regnans 0.2

* Earlier changes are recorded in the git history and the pull requests on GitHub.
