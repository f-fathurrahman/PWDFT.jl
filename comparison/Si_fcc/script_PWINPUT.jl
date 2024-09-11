atoms = Atoms(2, 1, [0.0 -2.565775; 0.513155 2.565775; 0.513155 2.565775], [1, 1], ["Si", "Si"], ["Si"], [-5.13155 0.0 -5.13155; 0.0 5.13155 5.13155; 5.13155 5.13155 0.0], [0.0], [0.0])
pspfiles = ["./Si.pz-nl-kjpaw_psl.1.0.0.UPF"]
ecutwfc = 15.0
ecutrho = 60.0
dual = ecutrho/ecutwfc
meshk = (3, 3, 3)
nbnd = 4
Ns = (0, 0, 0)
xcfunc = "VWN"
use_smearing = false
kT = 0.0

Etot_pwscf =  -89.18998303*0.5 # from Ry to Ha

forces_pwscf = zeros(Float64, 3, atoms.Natoms)
forces_pwscf[:,1] .= [0.11389374, -0.14106932, -0.14106932]
forces_pwscf[:,2] .= [-0.11389374, 0.14106932, 0.14106932]
forces_pwscf[:,:] *= 0.5 # from Ry to Ha

stress_pwscf = zeros(Float64, 3, 3)
stress_pwscf[1,:] .= [0.000593151, -0.000643861, -0.000643861]
stress_pwscf[2,:] .= [-0.000643861, 0.000305365, 0.00051235]
stress_pwscf[3,:] .= [-0.000643861, 0.00051235, 0.000305365]
stress_pwscf[:,:] *= 0.5  # from Ry to Ha
