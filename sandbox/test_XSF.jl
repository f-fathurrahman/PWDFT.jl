using Printf
using LinearAlgebra
using PWDFT

function test_main()
    LatVecs = gen_lattice_sc(16.0)
    ecutwfc_Ry = 30.0
    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )
    println( pw )

    Natoms = 1
    atpos = zeros( 3, Natoms )
    atpos[:,1] = [ 8.0, 2.0, 7.0 ]

    ik = 1
    Npoints = prod( pw.Ns )
    idx_gw2g = pw.gvecw.idx_gw2g
    idx_gw2r = pw.gvecw.idx_gw2r
    gwave = pw.gvec.G[:,idx_gw2g[ik]]
    Ngwx = pw.gvecw.Ngwx
    println("Ngwx = ", Ngwx)

    l = 1
    ia = 1
    psi = Array{ComplexF64}(undef,Ngwx)
    ctmp =Array{ComplexF64}(undef,Npoints)
    for m = -l:l
        ctmp .= 0.0 + im*0.0 # need to zero out ctmp first
        for ig = 1:Ngwx
            g = gwave[:,ig]
            Gm = norm(g)
            GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
            Sf = cos(GX) - im*sin(GX)
            psi[ig] = (1.0*im)^l * Ylm_real(l,m,g) * exp( -0.5*Gm^2 ) * Sf
        end
        ctmp[idx_gw2r[ik]] = psi
        psiR_real = real( G_to_R( pw, ctmp ) ) * Npoints
        filxsf = "TEMP_psi_l_" * string(l) * "_m_" * string(m) * ".xsf"
        write_xsf( filxsf, LatVecs/ANG2BOHR, atpos/ANG2BOHR, molecule=false )
        write_xsf_data3d_crystal(  filxsf, pw.Ns, LatVecs/ANG2BOHR, psiR_real )
        @printf("File %s is written\n", filxsf)
    end

end

test_main()


