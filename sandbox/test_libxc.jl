using LinearAlgebra
using Random
using Printf
using PWDFT

function test_LDA_VWN()
    Npoints = 5
    Rhoe = Array{Float64}(undef,Npoints)
    Rhoe[:] = [0.1, 0.2, 0.3, 0.4, 0.5]

    epsxc = calc_epsxc_VWN( Rhoe )
    Vxc = calc_Vxc_VWN( Rhoe )

    for ip = 1:Npoints
        @printf("%3d %18.10f %18.10f %18.10f\n", ip, Rhoe[ip], epsxc[ip], Vxc[ip])
    end
end


function test_LDA_VWN_spinpol()
    
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


function test_GGA_PBE()

    @printf("---------------\n")
    @printf("Testing GGA PBE\n")
    @printf("---------------\n")

    ecutwfc_Ry = 30.0
    LatVecs = gen_lattice_sc(16.0)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )

    srand(1234)
    Ngwx = pw.gvecw.Ngwx
    dVol = pw.CellVolume/prod(pw.Ns)
    
    Nkpt = 1
    psik = Array{Array{ComplexF64,2},1}(undef,Nkpt)

    Nstates = 4
    Focc = 2.0*ones(Nstates,Nkpt)
    psik[1] = ortho_gram_schmidt( rand(ComplexF64,Ngwx,Nstates) )

    Rhoe = calc_rhoe( pw, Focc, psik )
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( pw, Rhoe )
    @printf("sum Vxc = %18.10f\n", sum(Vxc))

    epsxc = calc_epsxc_PBE( pw, Rhoe )
    E_xc = dot( Rhoe, epsxc ) * dVol
    @printf("sum E_xc = %18.10f\n", E_xc)
end



function test_GGA_PBE_spinpol()
    
    @printf("-----------------------\n")
    @printf("Testing GGA PBE spinpol\n")
    @printf("-----------------------\n")

    ecutwfc_Ry = 30.0
    LatVecs = gen_lattice_sc(16.0)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )

    srand(1234)
    Ngwx = pw.gvecw.Ngwx
    dVol = pw.CellVolume/prod(pw.Ns)
    
    Nkpt = 2
    psik = Array{Array{ComplexF64,2},1}(undef,Nkpt)

    Nstates = 4
    Focc = 1.0*ones(Nstates,Nkpt)

    psik[1] = ortho_gram_schmidt( rand(ComplexF64,Ngwx,Nstates) )
    psik[2] = ortho_gram_schmidt( rand(ComplexF64,Ngwx,Nstates) )

    Npoints = prod(pw.Ns)
    Rhoe = zeros(Npoints,2)
    
    Rhoe[:,1] = calc_rhoe( pw, Focc[:,1:1], psik[1:1] )
    Rhoe[:,2] = calc_rhoe( pw, Focc[:,2:2], psik[2:2] )

    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( pw, Rhoe )
    @printf("sum Vxc = %18.10f\n", sum(Vxc))

    epsxc = calc_epsxc_PBE( pw, Rhoe )
    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    E_xc = dot( Rhoe_total, epsxc ) * dVol
    @printf("sum E_xc = %18.10f\n", E_xc)
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

    srand(1234)
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

test_LDA_VWN()
test_LDA_VWN_spinpol()
test_GGA_PBE()
test_GGA_PBE_spinpol()

test_spinpol(xc="VWN", Nspin=1)
test_spinpol(xc="VWN", Nspin=2)

test_spinpol(xc="PBE", Nspin=1)
test_spinpol(xc="PBE", Nspin=2)
