using Printf
using LinearAlgebra
using PWDFT

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function init_Ham_Si_fcc( a::Float64, meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(a))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    #return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk, Ns_=(32,32,32) )
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
end

function init_V_loc_from_radial( pw::PWGrid, pspot::PsPot_GTH )

    #r0 = 0.44642857142857e-3
    #r_ratio = 1.0247
    #cclog = log(r_ratio)
    #Nradial = 495

    r0 = 1e-5
    r_ratio = 1.0125
    cclog = log(r_ratio)
    Nradial = 2100

    
    r = zeros(Nradial)
    
    Vloc_radial = zeros(Nradial)
    
    r[1] = r0
    
    Vloc_radial[1] = PWDFT.eval_Vloc_R( pspot, r[1])
    
    @printf("%5d %18.10f %18.10f\n", 1, r[1], Vloc_radial[1])

    for i = 2:Nradial
        r[i] = r_ratio*r[i-1]
        Vloc_radial[i] = PWDFT.eval_Vloc_R( pspot, r[i])
        @printf("%5d %18.10f %18.10f\n", i, r[i], Vloc_radial[i])
    end

    CellVolume = pw.CellVolume
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2

    psi = zeros(Nradial)
    vps = zeros(Ng)

    for i = 1:Nradial
      psi[i] = Vloc_radial[i]*r[i]^3
    end 
    integ = sum(psi)*cclog
    vps[1] = integ*4.0*pi/CellVolume

    for ig = 2:Ng
        Gm = sqrt(G2[ig])
        for i = 1:Nradial
            arg = r[i]*Gm
            psi[i] = Vloc_radial[i]*sin(arg)/arg * r[i]^3
        end
        integ = sum( psi ) * cclog
        vps[ig] = integ*4.0*pi/CellVolume
    end

    vps_2 = eval_Vloc_G( pspot, G2 )/CellVolume

    println("")
    for ig = 1:100
        @printf("%18.10f %18.10f\n", vps[ig], vps_2[ig])
    end


end


function main()
    LATCONST = 10.2631
    Ham = init_Ham_Si_fcc( LATCONST, [3,3,3] )

    init_V_loc_from_radial( Ham.pw, Ham.pspots[1] )
end

main()