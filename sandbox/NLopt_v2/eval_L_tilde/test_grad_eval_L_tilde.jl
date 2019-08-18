include("eval_L_tilde.jl")

function test_grad_eval_L_tilde()
    Random.seed!(1234)
    Ham = create_Ham_atom_Pt_smearing()
    evars = rand_ElectronicVars(Ham)
    g_evars = ElectronicVars( copy(evars.psiks), copy(evars.Haux) )

    Etot = eval_L_tilde!(Ham, evars)
    println("Etot = ", Etot)

    grad_eval_L_tilde!(Ham, evars, g_evars)

    Etot = eval_L_tilde!(Ham, evars)
    println("Etot = ", Etot)
end