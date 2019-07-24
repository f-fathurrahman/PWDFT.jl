Separate SCF and band calculation.

SCF is done in `run_scf.jl`.
After SCF is converged Hamiltonian is written to file `Ham.data` using
`Serialization.serialize` function available in standard library of Julia.

In `run_bands.jl`, we read Hamiltonian from `Ham.data` and replace some
variables for band calculation. In order for this to work, we need to
make sure that the same `ecutwfc` is used.

Note that we have not yet provide a clean API for band structure calculation.

