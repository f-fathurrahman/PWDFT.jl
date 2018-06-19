# Manual precompilation script

precompile(Atoms,(Int64,Int64,Array{Float64,2},
                  Array{Int64,1},Array{String,1},Array{String,1},
                  Array{Float64,2},Array{Float64,1}))

precompile(init_atoms_xyz, (String,))
precompile(init_atoms_xyz, (String, Bool))

precompile(init_atoms_xyz_string, (String,))
precompile(init_atoms_xyz_string, (String, Bool))

precompile(KPoints, (Int64,Array{Float64,2},Array{Float64,1},Array{Float64,2}))
precompile(KPoints, (Atoms,))

precompile(KPoints, (Atoms,Array{Int64,1},Array{Int64,1}))
precompile(KPoints, (Atoms,Array{Int64,1},Array{Int64,1}, Int64))
precompile(KPoints, (Atoms,Array{Int64,1},Array{Int64,1}, Bool))
precompile(KPoints, (Atoms,Array{Int64,1},Array{Int64,1}, Int64,Bool))

precompile(PWGrid,(Float64,Array{Float64,2}))
precompile(PWGrid,(Float64,Array{Float64,2},KPoints))

precompile(PsPot_GTH, (Nothing,))
precompile(PsPot_GTH, (String,))

precompile(PWHamiltonian, (Atoms,Float64))
precompile(PWHamiltonian, (Atoms,Array{PsPot_GTH,1},Float64))

precompile(KS_solve_Emin_PCG!, (PWHamiltonian,))
