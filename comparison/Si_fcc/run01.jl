using PWDFT

include("script_PWINPUT.jl")

Ham = Hamiltonian(
    atoms, pspfiles, ecutwfc,
    meshk=[meshk[1], meshk[2], meshk[3]],
    dual=dual, Nstates=nbnd,
    xcfunc=xcfunc,
    Ns_=Ns);


psiks = rand_BlochWavefunc(Ham);

electrons_scf!(Ham, psiks, NiterMax=100, use_smearing=use_smearing, kT=kT, betamix=0.1);

