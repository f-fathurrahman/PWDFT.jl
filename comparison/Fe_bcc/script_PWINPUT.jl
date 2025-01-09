# Using PAW JTH LDA not converged

atoms = Atoms(2, 1, [0.0 2.7079775377810678; 0.0 2.7079775377810678; 0.0 2.7079775377810678], [1, 1], ["Fe", "Fe"], ["Fe"], [5.4159550755621355 0.0 0.0; 0.0 5.4159550755621355 0.0; 0.0 0.0 5.4159550755621355], [16.0], [0.0])
pspfiles = ["./Fe.upf"]
ecutwfc = 20.0
ecutrho = 100.0
dual = ecutrho/ecutwfc
meshk = (3, 3, 3)
nbnd = 22
Ns = (0, 0, 0)
xcfunc = "VWN"
use_smearing = true
kT = 0.005
