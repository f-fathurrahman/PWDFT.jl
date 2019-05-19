using Printf
using PWDFT

function test_LDA_VWN_small()
    Npoints = 5
    Rhoe = Array{Float64}(undef,Npoints)
    Rhoe[:] = [0.1, 0.2, 0.3, 0.4, 0.5]

    epsxc = calc_epsxc_VWN( Rhoe )
    Vxc = calc_Vxc_VWN( Rhoe )

    for ip = 1:Npoints
        @printf("%3d %18.10f %18.10f %18.10f\n", ip, Rhoe[ip], epsxc[ip], Vxc[ip])
    end
end


function test_LDA_VWN_spinpol_small()
    
    Npoints = 5
    Nspin = 2
    
    Rhoe = Array{Float64}(undef,Npoints,Nspin)
    
    Rhoe[:,1] = [0.1, 0.2, 0.3, 0.4, 0.5]
    Rhoe[:,2] = [0.1, 0.2, 0.3, 0.4, 0.5]

    epsxc = calc_epsxc_VWN( Rhoe )
    Vxc = calc_Vxc_VWN( Rhoe )

    for ip = 1:Npoints
        @printf("%3d | %18.10f | %18.10f | [%18.10f %18.10f]\n",
                ip, Rhoe[ip], epsxc[ip], Vxc[ip,1], Vxc[ip,2])
    end

end

function test_spinpol( ; xc="VWN", Nspin=1 )

    @assert( xc=="VWN" || xc=="PBE" )
    
    @assert Nspin <= 2
    
    @printf("\n")
    @printf("---------------------\n")
    @printf("Testing LibXC: %s\n", xc)
    if Nspin == 1
        @printf("Non spin polarized\n")
    else
        @printf("Spin polarized\n")
    end
    @printf("---------------------\n")

    atoms = init_atoms_xyz_string(
        """
        1

        H   0.0   0.0   0.0
        """
    )
    atoms.LatVecs = gen_lattice_fcc(5.0)

    meshk = [2,2,2]
    shiftk = [0,0,0]
    kpoints = KPoints( atoms, meshk, shiftk )

    ecutwfc_Ry = 30.0
    LatVecs = gen_lattice_sc(16.0)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs, kpoints=kpoints )

    Random.seed!(1234)
    Ngw = pw.gvecw.Ngw
    dVol = pw.CellVolume/prod(pw.Ns)

    Nkpt = kpoints.Nkpt
    Nkspin = Nkpt*Nspin

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    Nstates = 4

    if Nspin == 2
        Focc = 1.0*ones(Nstates,Nkspin)
        #idx_up = 1:Nkpt
        #idx_dn = Nkpt+1:2*Nkpt
        #Focc[Nstates,idx_dn] = 0.0
    else
        Focc = 2.0*ones(Nstates,Nkspin)
        #Focc[Nstates,:] = 1.0
    end

    for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            psiks[ikspin] = ortho_gram_schmidt( rand(ComplexF64,Ngw[ik],Nstates) )
        end
    end

    Npoints = prod(pw.Ns)
    Rhoe = zeros(Npoints,Nspin)
    
    for ispin = 1:Nspin
        idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
        Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
    end

    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    if Nspin == 2
        diffRhoe = Rhoe[:,1] - Rhoe[:,2]
        @printf("Integrated diff rhoe = %18.10f\n", sum(diffRhoe)*dVol)
    end

    if xc == "VWN"
        Vxc = calc_Vxc_VWN( Rhoe )
    else
        Vxc = calc_Vxc_PBE( pw, Rhoe )
    end

    for ispin = 1:Nspin
        @printf("ispin = %d, sum Vxc = %18.10f\n", ispin, sum(Vxc[:,ispin]))
    end

    if xc == "VWN"
        epsxc = calc_epsxc_VWN( Rhoe )
    else
        epsxc = calc_epsxc_PBE( pw, Rhoe )
    end

    if Nspin == 2
        Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    else
        Rhoe_total = Rhoe[:,1]
    end

    E_xc = dot( Rhoe_total, epsxc ) * dVol
    @printf("sum E_xc = %18.10f\n", E_xc)
end