include("eval_L_tilde.jl")

function test_grad_eval_L_tilde()
    Random.seed!(1234)
    Ham = create_Ham_Al_fcc_smearing()
    evars = rand_ElectronicVars(Ham)
    g_evars = ElectronicVars( copy(evars.psiks), copy(evars.Haux) )

#=
    # prepare guess wavefunc
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Rhoe = zeros(Float64,Npoints,Nspin)
    @assert Nspin == 1
    Rhoe[:,1] = guess_rhoe( Ham )
    update!(Ham, Rhoe)
    # eigenvalues are not needed for this case
    _ = diag_LOBPCG!( Ham, evars.psiks, verbose=false, verbose_last=false, NiterMax=10 )

    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin -1)*Nkpt
        evars.Haux[i] = Hermitian( evars.psiks[i]' * ( Ham * evars.psiks[i] ) )
    end
=#

    print_Haux(evars, "evars before eval_L_tilde")
    Etot = eval_L_tilde!(Ham, evars)
    print_Haux(evars, "evars after eval_L_tilde")

    println("Etot = ", Etot)

    grad_eval_L_tilde!(Ham, evars, g_evars)
    print_Haux(g_evars, "g_evars before grad_eval_L_tilde")

    Etot = eval_L_tilde!(Ham, evars)
    println("Etot = ", Etot)
end
test_grad_eval_L_tilde()