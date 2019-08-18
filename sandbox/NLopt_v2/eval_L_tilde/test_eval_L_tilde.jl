include("eval_L_tilde.jl")

function test_eval_L_tilde()

    Random.seed!(1234)

    Ham = create_Ham_atom_Pt_smearing()
    evars = rand_ElectronicVars(Ham)

    println("Before eval_L_tilde!")
    display(real(evars.Haux[1]))
    println()
    display(imag(evars.Haux[1]))
    println()

    Etot, Hsub = eval_L_tilde!(Ham, evars)
    @printf("Etot = %18.10f\n", Etot)

    println("After eval_L_tilde!")
    print_Haux(evars, "evars after eval_L_tilde")
    println()
    println("Hsub = ")
    display(Hsub)
    println()

end
test_eval_L_tilde()