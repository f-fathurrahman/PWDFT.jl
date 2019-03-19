Some examples for testing various `KS_solve_*` functions.

Use `run.jl` script as the driver:
```
julia run.jl atom_H/main_H_atom.jl Emin
```

`Emin` can be changed to `SCF`, `DCM`, or `TRDCM`. Not all systems supports all
`KS_solve_*` functions.
