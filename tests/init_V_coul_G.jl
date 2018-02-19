using PWDFT
using PWDFT.PW

function init_V_coulomb_G( pw::PWGrid, strf::Array{Complex128,2}, Znucls::Array{Float64} )

    Nsp1 = size(strf)[2]
    Nsp2 = size(Znucls)[1]
    #
    # This is the restriction for this function
    #
    if Nsp1 != Nsp2
        @printf("ERROR: Nsp1 /= Nsp2: %d %d\n", Nsp1, Nsp2)
        exit()
    end
    Nspecies = Nsp1

    Npoints = prod(pw.Ns)
    Ω = pw.Ω

    Vg = zeros(Complex128, Npoints)
    G2 = pw.gvec.G2
    V  = zeros(Float64, Npoints)
    #
    for isp = 1:Nspecies
        prefactor = -4*pi*Znucls[isp]/Ω
        #
        # V[ig=0] = 0.0
        #
        for ig = 2:Npoints
            Vg[ig] = prefactor/G2[ig]*strf[ig,isp]
        end
        V[:] = V[:] + real( G_to_R(pw.Ns, Vg) ) * Npoints
    end
    return V
end


# Calculate structure factor
function calc_strfact( atpos::Array{Float64,2}, Nspecies::Int,
    atm2species::Array{Int,1}, G::Array{Float64,2} )
    #
    Ng = size(G)[2]
    Na = size(atpos)[2]
    Sf = zeros(Complex128,Ng,Nspecies)
    for ia = 1:Na
        isp = atm2species[ia]
        for ig = 1:Ng
            GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
            Sf[ig,isp] = Sf[ig,isp] + cos(GX) - im*sin(GX)
        end
    end
    return Sf
end

function calc_strfact( atoms::Atoms, pw::PWGrid )
    return calc_strfact( atoms.positions, atoms.Nspecies, atoms.atm2species, pw.gvec.G )
end



function test_main()
    const LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0*0.5
    pw = PWGrid( ecutwfc_Ry, LatVecs )
    println(pw)
    #
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)
    #
    strf = calc_strfact( atoms, pw )
    #
    V = init_V_coulomb_G( pw, strf, [1.0] )
    println("sum V = ", sum(V))
end

test_main()
