function calc_forces_Ps_loc( Ham )
    F_Ps_loc = zeros(Float64, 3, Ham.atoms.Natoms)
    calc_forces_Ps_loc!( Ham, F_Ps_loc )
    return F_Ps_loc
end

function calc_forces_Ps_loc!( Ham::Hamiltonian, F_Ps_loc::Array{Float64,2} )
    calc_forces_Ps_loc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe, F_Ps_loc )
    return
end

function calc_forces_Ps_loc!(
    atoms::Atoms, pw::PWGrid, pspots::Array{PsPot_GTH,1},
    Rhoe::Array{Float64,2},
    F_Ps_loc::Array{Float64,2}
)
    
    Ng = pw.gvec.Ng
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    G = pw.gvec.G
    G2 = pw.gvec.G2

    idx_g2r = pw.gvec.idx_g2r
    
    Npoints = prod(pw.Ns)
    
    RhoeG = zeros(ComplexF64, Npoints)
    Nspin = size(Rhoe,2)
    # XXX Use sum or reduce
    for ispin in 1:Nspin, ip in 1:Npoints
        RhoeG[ip] = RhoeG[ip] + Rhoe[ip,ispin]
    end
    R_to_G!(pw, RhoeG)
    @views RhoeG[:] .= RhoeG[:]/Npoints # this normalization is required
    #
    # Precalculate Vg for each species
    # Using G-shells
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Vgl = zeros(Float64, Ngl)
    idx_g2shells = pw.gvec.idx_g2shells
    Vgl = zeros(Float64, Ngl, Nspecies)
    for isp in 1:Nspecies
        psp = pspots[isp]
        for igl in 1:Ngl
            Vgl[igl,isp] = eval_Vloc_G( psp, G2_shells[igl] )
        end
    end
    # We don't need factor of CellVolume.
    # It will be cancelled by multiplication in the expression for F_Ps_loc
    
    X = atoms.positions
    #
    for ia in 1:Natoms
        #
        isp = atm2species[ia]
        psp = pspots[isp]
        #
        # G=0 contribution is zero, G[:,ig=0] = [0, 0, 0]
        for ig in 2:Ng
            GX = G[1,ig]*X[1,ia] + G[2,ig]*X[2,ia] + G[3,ig]*X[3,ia]
            Sf = cos(GX) - im*sin(GX)
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            @views F_Ps_loc[:,ia] .+= real(im*G[:,ig]*Vgl[igl,isp]*conj(RhoeG[ip])*Sf)
        end
    end
    #
    return
end