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

    for i in 1:length(evars.Haux)
        evars.Haux[i] = Hermitian( evars.psiks[i]' * ( Ham * evars.psiks[i] ) )
    end

    Rhoe_1 = copy(Ham.rhoe)
    print_Haux(evars, "evars before eval_L_tilde!")
    
    rotate_evars!( Ham, evars )
    Etot = eval_L_tilde!(Ham, evars)
    
    print_Haux(evars, "evars after eval_L_tilde!")    
    @printf("Etot = %18.10f\n", Etot)

    Rhoe_2 = copy(Ham.rhoe)
    print_Haux(evars, "evars before eval_L_tilde!")
    
    Etot = eval_L_tilde!(Ham, evars)
    
    print_Haux(evars, "evars after eval_L_tilde!")    
    @printf("Etot = %18.10f\n", Etot)

    println("sum diff rhoe = ", sum(Rhoe_2 - Rhoe_1))

end
test_eval_L_tilde()