using LinearAlgebra
using Random
using Printf
using Libxc
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("Libxc_GGA.jl")

include("XC_x_slater.jl")
include("XC_c_pw.jl")
include("XC_c_vwn.jl")

include("XC_x_slater_spin.jl")
include("XC_c_pw_spin.jl")
include("XC_c_vwn_spin.jl")

include("XC_x_pbe.jl")
include("XC_c_pbe.jl")

include("XC_c_pbe_spin.jl")

function test_LDA_VWN()

    @printf("---------------\n")
    @printf("Testing LDA VWN\n")
    @printf("---------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    pw = Ham.pw
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    Rhoe = guess_rhoe_atomic(Ham)
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_VWN( Ham.xc_calc, Rhoe )

    epsxc = calc_epsxc_VWN( Ham.xc_calc, Rhoe )
    E_xc = dot( Rhoe, epsxc ) * dVol

    etxc = 0.0
    
    Vxc_v2 = zeros(Npoints)
    for ip in 1:Npoints
        ss_x, vx = XC_x_slater( Rhoe[ip] )
        ss_c, vc = XC_c_vwn( Rhoe[ip] )
        etxc = etxc + (ss_x + ss_c)*Rhoe[ip]
        Vxc_v2[ip] = vx + vc
    end


    @printf("E_xc v1 = %18.10f\n", E_xc)
    @printf("E_xc v2 = %18.10f\n", etxc*dVol)
    
    @printf("sum Vxc v1 = %18.10f\n", sum(Vxc))
    @printf("sum Vxc v2 = %18.10f\n", sum(Vxc_v2))

end


function test_LDA_VWN_spinpol()

    @printf("-----------------------\n")
    @printf("Testing LDA VWN spinpol\n")
    @printf("-----------------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Nspin=2 )

    pw = Ham.pw
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    Rhoe = guess_rhoe_atomic(Ham, starting_magnetization=[0.1])
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_VWN( Ham.xc_calc, Rhoe )

    epsxc = calc_epsxc_VWN( Ham.xc_calc, Rhoe )
    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    E_xc = dot( Rhoe_total, epsxc ) * dVol

    etxc = 0.0
    Vxc_v2 = zeros(Npoints,2)
    
    for ip in 1:Npoints    
        ρ = Rhoe[ip,1] + Rhoe[ip,2]
        ζ = (Rhoe[ip,1] - Rhoe[ip,2])/ρ

        ss_x, vxup, vxdn = XC_x_slater_spin( ρ, ζ )
        ss_c, vcup, vcdn = XC_c_vwn_spin( ρ, ζ )
        
        etxc = etxc + (ss_x + ss_c)*ρ

        Vxc_v2[ip,1] = vxup + vcup
        Vxc_v2[ip,2] = vxdn + vcdn
    end

    @printf("E_xc v1 = %18.10f\n", E_xc)
    @printf("E_xc v2 = %18.10f\n", etxc*dVol)

    for ispin = 1:2
        @printf("sum Vxc v1 = %18.10f\n", sum(Vxc[:,ispin]))
        @printf("sum Vxc v2 = %18.10f\n", sum(Vxc_v2[:,ispin]))        
    end

end



function test_GGA_PBE()

    @printf("---------------\n")
    @printf("Testing GGA PBE\n")
    @printf("---------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth", "Si-q4.gth")]
    ecutwfc = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[1,1,1], xcfunc="PBE" )

    pw = Ham.pw
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/prod(pw.Ns)

    Rhoe = guess_rhoe_atomic(Ham)
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( Ham.xc_calc, pw, Rhoe )

    epsxc = calc_epsxc_PBE( Ham.xc_calc, pw, Rhoe )
    
    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe[:,1] )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = dot( gRhoe[:,ip], gRhoe[:,ip] )
    end

    #
    Vxc_v2 = zeros(Npoints) # should depend also to Nspin
    # h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
    h = zeros(3,Npoints)
    #
    dh = zeros(Npoints)

    etxc = 0.0
    etgxc = 0.0
    for ip in 1:Npoints

        ss_x, vx = XC_x_slater( Rhoe[ip,1] )
        gss_x, v1x, v2x = XC_x_pbe( Rhoe[ip,1], gRhoe2[ip] )
        
        ss_c, vc = XC_c_pw( Rhoe[ip,1] )
        gss_c, v1c, v2c = XC_c_pbe( Rhoe[ip,1], gRhoe2[ip] )

        Vxc_v2[ip] = vx + vc + v1x + v1c
        h[1:3,ip] = (v2x + v2c)*gRhoe[1:3,ip]

        etxc = etxc + ( ss_x + ss_c )*Rhoe[ip,1]
        etgxc = etgxc + ( gss_x + gss_c )

    end

    println("sum abs gRhoe 1 = ", sum(abs.(gRhoe[1,:])))
    println("sum abs gRhoe 2 = ", sum(abs.(gRhoe[2,:])))
    println("sum abs gRhoe 3 = ", sum(abs.(gRhoe[3,:])))

    println("sum gRhoe2 = ", sum(gRhoe2))

    dh[:] = op_nabla_dot(pw, h)

    println("sum abs h  = ", sum(abs.(h)))
    println("sum abs dh = ", sum(abs.(dh)))
    for ip in 1:Npoints
        Vxc_v2[ip] = Vxc_v2[ip] - dh[ip]
    end

    @printf("E_xc v1 = %18.10f\n", dot( Rhoe, epsxc ) * dVol)
    @printf("E_xc v2 = %18.10f\n", (etxc + etgxc)*dVol )

    @printf("sum Vxc v1 = %18.10f\n", sum(Vxc))
    @printf("sum Vxc v2 = %18.10f\n", sum(Vxc_v2))


end



function test_GGA_PBE_spinpol()
    
    @printf("-----------------------\n")
    @printf("Testing GGA PBE spinpol\n")
    @printf("-----------------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[1,1,1], xcfunc="PBE", Nspin=2 )

    pw = Ham.pw
    dVol = pw.CellVolume/prod(pw.Ns)

    Rhoe = guess_rhoe_atomic(Ham, starting_magnetization=[0.1])

    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( Ham.xc_calc, pw, Rhoe )

    Npoints = prod(Ham.pw.Ns)
    #f = open("TEMP_Vxc_spinpol.dat", "w")
    #for ip = 1:Npoints
    #    @printf(f, "%d %18.10f %18.10f\n", ip, Vxc[ip,1], Vxc[ip,2])
    #end
    #close(f)

    epsxc = calc_epsxc_PBE( Ham.xc_calc, pw, Rhoe )
    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    E_xc = dot( Rhoe_total, epsxc ) * dVol


    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] )
    gRhoe2_up = zeros( Float64, Npoints )

    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )
    gRhoe2_dn = zeros( Float64, Npoints )

    gRhoe = op_nabla( pw, Rhoe_total )
    gRhoe2 = zeros( Float64, Npoints )

    for ip = 1:Npoints
        gRhoe2_up[ip] = dot( gRhoe_up[:,ip], gRhoe_up[:,ip] )
        gRhoe2_dn[ip] = dot( gRhoe_dn[:,ip], gRhoe_dn[:,ip] )
        gRhoe2[ip] = dot( gRhoe[:,ip],  gRhoe[:,ip] )
    end

    #
    Vxc_v2 = zeros(Npoints,2) # should depend also to Nspin
    # h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
    h_up = zeros(3,Npoints)
    h_dn = zeros(3,Npoints)
    #
    dh_up = zeros(Npoints)
    dh_dn = zeros(Npoints)

    etxc = 0.0
    etgxc = 0.0
    for ip in 1:Npoints

        ρ_up = Rhoe[ip,1]
        ρ_dn = Rhoe[ip,2]
        ρ = ρ_up + ρ_dn
        ζ = (ρ_up - ρ_dn)/ρ

        ss_x, vxup, vxdn = XC_x_slater_spin( ρ, ζ )
        ss_c, vcup, vcdn = XC_c_pw_spin( ρ, ζ )

        gss_xup, v1xup, v2xup = XC_x_pbe( 2*ρ_up, 4*gRhoe2_up[ip] )
        gss_xdn, v1xdn, v2xdn = XC_x_pbe( 2*ρ_dn, 4*gRhoe2_dn[ip] )

        gss_x = 0.5 * (gss_xup + gss_xdn)
        v2xup = 2.0 * v2xup
        v2xdn = 2.0 * v2xdn

        gss_c, v1cup, v1cdn, v2c = XC_c_pbe_spin( ρ, ζ, gRhoe2[ip] )
        v2cup = v2c
        v2cdn = v2c
        v2cud = v2c

        Vxc_v2[ip,1] = vxup + vcup + v1xup + v1cup
        Vxc_v2[ip,2] = vxdn + vcdn + v1xdn + v1cdn

        for i in 1:3
           grup = gRhoe_up[i,ip]
           grdn = gRhoe_dn[i,ip]
           h_up[i,ip] = ( v2xup + v2cup ) * grup + v2cud * grdn
           h_dn[i,ip] = ( v2xdn + v2cdn ) * grdn + v2cud * grup
        end

        etxc = etxc + ( ss_x + ss_c )*ρ # use total Rhoe
        etgxc = etgxc + ( gss_x + gss_c )

    end

    dh_up[:] = op_nabla_dot(pw, h_up)
    dh_dn[:] = op_nabla_dot(pw, h_dn)

    println("sum abs h  up = ", sum(abs.(h_up)))
    println("sum abs dh up = ", sum(abs.(dh_up)))

    println("sum abs h  dn = ", sum(abs.(h_dn)))
    println("sum abs dh dn = ", sum(abs.(dh_dn)))

    for ip in 1:Npoints
        Vxc_v2[ip,1] = Vxc_v2[ip,1] - dh_up[ip]
        Vxc_v2[ip,2] = Vxc_v2[ip,2] - dh_dn[ip]        
    end

    @printf("E_xc    = %18.10f\n", E_xc)
    @printf("E_xc v2 = %18.10f\n", (etxc + etgxc)*dVol )

    @printf("avg Vxc v1 up = %18.10f\n", sum(Vxc[:,1])/Npoints)
    @printf("avg Vxc v2 up = %18.10f\n", sum(Vxc_v2[:,1])/Npoints)

    @printf("avg Vxc v1 dn = %18.10f\n", sum(Vxc[:,2])/Npoints)
    @printf("avg Vxc v2 dn = %18.10f\n", sum(Vxc_v2[:,2])/Npoints)

end


#test_LDA_VWN()
#test_LDA_VWN_spinpol()
#@time test_LDA_VWN_spinpol()
#@time test_LDA_VWN_spinpol()
#test_GGA_PBE()
test_GGA_PBE_spinpol()

#test_spinpol(xc="VWN", Nspin=1)
#test_spinpol(xc="VWN", Nspin=2)
#test_spinpol(xc="PBE", Nspin=1)
#test_spinpol(xc="PBE", Nspin=2)
