Directory structure:

- src: source code of the package

- pseudopotentials: GTH pseudopotential parameters for LDA and GGA

- examples: a collection of example scripts

- structures: a collection of several molecular and crystalline structures
  in xyz format.


Important files in src directory:

- PWDFT.jl: main package file.

- Atoms.jl: contains struct definition and functions for working with molecular and
  crystalline structures.

- PWGrid.jl: contains struct definition and functions for working plane wave basis.

- PsPot_GTH: contains struct definition and functions for working with GTH pseudopotentials.

- Electrons.jl: contains struct definition and functions for working electronic states.

- KPoints.jl: contains struct definition and functions for working with k-points

- Potentials.jl: contains struct definition and functions for working with local 
  potential terms.

- PsPotNL.jl: contains struct definition and functions for working with nonlocal pseudopotential
  terms.

- Hamiltonian.jl: contains struct definition and functions for working with Kohn-Sham
  Hamiltonian.

- op_K.jl, op_V_loc.jl, op_V_Ps_loc.jl, and op_H.jl: contain functions for applying kinetic,
  local potential, local pseudopotential, and Hamiltonian operators to wave functions.

- KS_solve_SCF.jl: contains a function to solve Kohn-Sham problem using self-consistent field
  algorithm.

- KS_solve_Emin_PCG.jl: contain a function to solve Kohn-Sham problem using direct minimization
  (preconditioned conjugate gradient) of Kohn-Sham energy.
