atoms_tuple = (Natoms = 2, Nspecies = 1, positions = [0.0 2.7079775377810678; 0.0 2.7079775377810678; 0.0 2.7079775377810678], atm2species = [1, 1], atsymbs = ["Fe", "Fe"], SpeciesSymbols = ["Fe"], LatVecs = [5.4159550755621355 0.0 0.0; 0.0 5.4159550755621355 0.0; 0.0 0.0 5.4159550755621355], Zvals = [16.0], masses = [0.0]);
atoms = Atoms(atoms_tuple...);
pspfiles = ["../../pseudopotentials/GBRV_LDA/fe_lda_v1.5.uspp.F.UPF"];
ecutwfc = 20.0;
options_tuple = (dual = 5.0, Nspin_wf = 2, Nspin_dens = 2, meshk = [5, 5, 5], shiftk = [0, 0, 0], time_reversal = true, Ns = (0, 0, 0), kpoints = nothing, kpts_str = nothing, xcfunc = "VWN", use_xc_internal = false, extra_states = nothing, Nstates = 22, use_symmetry = true, use_smearing = true, smearing_kT = 0.005, starting_magn = [0.1], angle1 = nothing, angle2 = nothing, lspinorb = false, noncollinear = false);
options = HamiltonianOptions(options_tuple...);
pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
for isp in 1:atoms.Nspecies
    pspots[isp] = PsPot_UPF(pspfiles[isp]);
end
Ham = Hamiltonian(atoms, pspots, ecutwfc, options);

