atoms = Atoms(2, 2, [2.6550652062682487 0.5310130412536498; 2.6550652062682487 0.5310130412536498; 2.6550652062682487 0.0], [1, 2], ["As", "Al"], ["As", "Al"], [5.310130412536497 5.310130412536497 0.0; 5.310130412536497 0.0 5.310130412536497; 0.0 5.310130412536497 5.310130412536497], [5.0, 3.0], [0.0, 0.0])
pspfiles = ["./As.UPF", "./Al.UPF"]
ecutwfc = 20.0
ecutrho = 100.0
dual = ecutrho/ecutwfc
meshk = (3, 3, 3)
nbnd = -1
Ns = (0, 0, 0)
xcfunc = "VWN"
use_smearing = false
kT = 0.0

Etot_pwscf = -46.44783520*0.5 # from Ry to Ha

forces_pwscf = zeros(Float64, 3, atoms.Natoms)
forces_pwscf[:,1] .= [0.11292586, 0.11292586, 0.09142182]
forces_pwscf[:,2] .= [-0.11292586, -0.11292586, -0.09142182]
forces_pwscf[:,:] *= 0.5 # from Ry to Ha

stress_pwscf = zeros(Float64, 3, 3)
stress_pwscf[1,:] .= [0.000358565, 0.00038945 , 0.000611996]
stress_pwscf[2,:] .= [0.00038945 , 0.000358565, 0.000611996]
stress_pwscf[3,:] .= [0.000611996, 0.000611996, 0.000645677]
stress_pwscf[:,:] *= 0.5  # from Ry to Ha