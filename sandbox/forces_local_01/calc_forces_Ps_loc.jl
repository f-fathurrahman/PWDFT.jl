function calc_forces_Ps_loc(Ham::Hamiltonian)
    
    Ng = Ham.pw.gvec.Ng
    Natoms = Ham.atoms.Natoms

    G = Ham.pw.gvec.G
    G2 = Ham.pw.gvec.G2

    V_Ps_loc = Ham.potentials.Ps_loc
    idx_g2r = Ham.pw.gvec.idx_g2r
    
    Npoints = prod(Ham.pw.Ns)
    Rhoe = zeros(Float64,Npoints)
    
    for ispin = 1:Ham.electrons.Nspin
        for ip = 1
            Rhoe[ip] = Rhoe[ip] + Ham.rhoe[ip,ispin]
        end
    end
    RhoeG = R_to_G(Ham.pw, Rhoe)*Npoints

    Ω = Ham.pw.CellVolume
    
    F_Ps_loc = zeros(ComplexF64,3,Natoms)

    atpos = Ham.atoms.positions

    #Vg = R_to_G(Ham.pw, V_Ps_loc)
    
    for ia = 1:Natoms
        
        isp = Ham.atoms.atm2species[ia]
        psp = Ham.pspots[isp]
        
        # G=0 contribution is zero, G[:,ig=0] = [0, 0, 0]
        for ig = 2:Ng
            
            GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
            Sf = cos(GX) - im*sin(GX)

            Vg_ig = eval_Vloc_G( psp, G2[ig] )

            ip = idx_g2r[ig]
            F_Ps_loc[:,ia] = F_Ps_loc[:,ia] +
                             im*G[:,ig]*Vg_ig*conj(RhoeG[ip])*Sf*Natoms
            #F_Ps_loc[:,ia] = F_Ps_loc[:,ia] +
            #G[:,ig]*Vg_ig*(sin(GX)*real(RhoeG[ip]) + cos(GX)*imag(RhoeG[ip]))
        end
        println(F_Ps_loc[:,ia])
    end

    return real(F_Ps_loc)
end

function calc_forces_Ps_loc_finite_diff(Ham::Hamiltonian)
    Δ = 0.001

    pos_orig = copy(Ham.atoms.positions)
    Natoms = Ham.atoms.Natoms

    F_Ps_loc = zeros(3,Natoms)

    for ia = 1:Natoms
    for idir = 1:3
        
        # set to original positions
        Ham.atoms.positions[:,:] = pos_orig[:,:]
        
        Ham.atoms.positions[idir,ia] = pos_orig[idir,ia] + 0.5*Δ
        KS_solve_Emin_PCG!(Ham)
        Eplus = Ham.energies.Ps_loc
        
        Ham.atoms.positions[idir,ia] = pos_orig[idir,ia] - 0.5*Δ
        KS_solve_Emin_PCG!(Ham)
        Eminus = Ham.energies.Ps_loc

        F_Ps_loc[idir,ia] = -(Eplus - Eminus)/Δ
    end
    end

    Ham.atoms.positions[:,:] = pos_orig[:,:]

    return F_Ps_loc
end


