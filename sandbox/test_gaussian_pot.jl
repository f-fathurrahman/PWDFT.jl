function init_V_gauss_G(a::Atoms, A, alpha)
end

function test_main()
    # Initialize positions
    atoms = init_xyz("pos.xyz")
    # set cell
    atoms.cell = 16.0*diagm(ones(3))
    center!(atoms)
    #
    #pw_calc = PWDFT( atoms, ecutwfc, )

    # init potential
    V_ps_loc = init_V_gauss_G( atoms, A, alpha )
end