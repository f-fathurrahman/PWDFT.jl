include("eval_L_tilde.jl")

function test_eval_L_tilde()

    Random.seed!(1234)

    Ham = create_Ham_atom_Pt_smearing()
    evars = rand_ElectronicVars(Ham)

    # prepare guess wavefunc
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64,Npoints,Nspin)
    @assert Nspin == 1
    Rhoe[:,1] = guess_rhoe( Ham )
    update!(Ham, Rhoe)
    # eigenvalues are not needed for this case
    _ = diag_LOBPCG!( Ham, evars.psiks, verbose=false, verbose_last=false, NiterMax=10 )

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
    println("real(Hsub) = ")
    display(real(Hsub[1]))
    println()
    println("imag(Hsub) = ")
    display(imag(Hsub[1]))
    println()

    Etot, Hsub = eval_L_tilde!(Ham, evars)
    @printf("Etot = %18.10f\n", Etot)

    Etot, Hsub = eval_L_tilde!(Ham, evars)
    @printf("Etot = %18.10f\n", Etot)

end
test_eval_L_tilde()