using LinearAlgebra
using Random
using Printf
using Libxc
using PWDFT

import PWDFT: guess_rhoe_atomic

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

    Vxc = calc_Vxc_VWN( Rhoe )

    epsxc = calc_epsxc_VWN( Rhoe )
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

    Vxc = calc_Vxc_VWN( Rhoe )

    epsxc = calc_epsxc_VWN( Rhoe )
    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    E_xc = dot( Rhoe_total, epsxc ) * dVol

    etxc = 0.0
    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3) 
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

    Vxc = calc_Vxc_PBE( pw, Rhoe )

    epsxc = calc_epsxc_PBE( pw, Rhoe )

    etxc = 0.0
    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    
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

    etgxc = 0.0
    ex_only = 0.0
    ec_only = 0.0
    SMALL = 1e-10
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

    Rhoe = guess_rhoe_atomic(Ham, starting_magnetization=[0.0])

    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( pw, Rhoe )

    Npoints = prod(Ham.pw.Ns)
    f = open("TEMP_Vxc_spinpol.dat", "w")
    for ip = 1:Npoints
        @printf(f, "%d %18.10f %18.10f\n", ip, Vxc[ip,1], Vxc[ip,2])
    end
    close(f)

    epsxc = calc_epsxc_PBE( pw, Rhoe )
    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    E_xc = dot( Rhoe_total, epsxc ) * dVol
    @printf("E_xc = %18.10f\n", E_xc)

    psiks = rand_BlochWavefunc(Ham)
    update!( Ham, Rhoe )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    Nstates = Ham.electrons.Nstates

    for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin-1)*Nkpt
            psiks[ikspin] = ones(Ngw[ik], Nstates)/Ngw[ik]
            psiks[ikspin] = op_V_loc(Ham, psiks[ikspin])
            ss = sum(psiks[ikspin])
            @printf("sum psiks[ikspin] = [%18.10f,%18.10f]\n", real(ss), imag(ss))
        end
    end

end


#test_LDA_VWN()
#@time test_LDA_VWN_spinpol()
#@time test_LDA_VWN_spinpol()
test_GGA_PBE()
#test_GGA_PBE_spinpol()

#test_spinpol(xc="VWN", Nspin=1)
#test_spinpol(xc="VWN", Nspin=2)
#test_spinpol(xc="PBE", Nspin=1)
#test_spinpol(xc="PBE", Nspin=2)
