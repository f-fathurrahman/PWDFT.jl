using PWDFT

function test_LDA_VWN()
    Npoints = 5
    Rhoe = zeros( Float64, Npoints )
    Rhoe = [0.1, 0.2, 0.3, 0.4, 0.5]

    #
    """
    ccall( (:calc_Vxc_VWN, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, Vxc )
    #
    ccall( (:calc_epsxc_VWN, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, epsxc )
    """
    epsxc = calc_epsxc_VWN( Rhoe )
    Vxc = calc_Vxc_VWN( Rhoe )

    for ip = 1:Npoints
        @printf("%3d %18.10f %18.10f %18.10f\n", ip, Rhoe[ip], epsxc[ip], Vxc[ip])
    end
end


function test_GGA_PBE()
    ecutwfc_Ry = 30.0
    LatVecs = 16.0*eye(3)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )

    srand(1234)
    Ngwx = pw.gvecw.Ngwx
    dVol = pw.Ω/prod(pw.Ns)

    Nstates = 4
    Focc = 2.0*ones(Nstates)

    psi = ortho_gram_schmidt( rand(Complex128,Ngwx,Nstates) )

    Rhoe = calc_rhoe( pw, Focc, psi )
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( pw, Rhoe )
    @printf("sum Vxc = %18.10f\n", sum(Vxc))

    epsxc = calc_epsxc_PBE( pw, Rhoe )
    E_xc = dot( Rhoe, epsxc ) * dVol
    @printf("sum E_xc = %18.10f\n", E_xc)
end

#test_LDA_VWN()

test_GGA_PBE()
