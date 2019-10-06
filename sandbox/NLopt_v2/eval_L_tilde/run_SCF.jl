include("eval_L_tilde.jl")

function main()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_atom_Pt()
    Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()

    KS_solve_SCF!( Ham,
        print_final_ebands=true,
        mix_method="rpulay",
        betamix=0.2,
        use_smearing=true,
        savewfc=true )

    #KS_solve_Emin_PCG!( Ham, print_final_ebands=true, i_cg_beta=4 )

end

main()
