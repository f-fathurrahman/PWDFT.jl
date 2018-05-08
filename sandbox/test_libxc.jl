using PWDFT

function test_LDA_VWN()
    Npoints = 5
    Rhoe = zeros( Float64, Npoints )
    Rhoe = [0.1, 0.2, 0.3, 0.4, 0.5]

    epsxc = calc_epsxc_VWN( Rhoe )
    Vxc = calc_Vxc_VWN( Rhoe )

    for ip = 1:Npoints
        @printf("%3d %18.10f %18.10f %18.10f\n", ip, Rhoe[ip], epsxc[ip], Vxc[ip])
    end
end


function test_LDA_VWN_spinpol()
    Npoints = 5
    Rhoe = zeros( Float64, Npoints, 2 )
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
    LatVecs = 16.0*eye(3)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )

    srand(1234)
    Ngwx = pw.gvecw.Ngwx
    dVol = pw.Ω/prod(pw.Ns)
    
    Nkpt = 1
    psik = Array{Array{Complex128,2},1}(Nkpt)

    Nstates = 4
    Focc = 2.0*ones(Nstates,Nkpt)
    psik[1] = ortho_gram_schmidt( rand(Complex128,Ngwx,Nstates) )

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
    LatVecs = 16.0*eye(3)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )

    srand(1234)
    Ngwx = pw.gvecw.Ngwx
    dVol = pw.Ω/prod(pw.Ns)
    
    Nkpt = 2
    psik = Array{Array{Complex128,2},1}(Nkpt)

    Nstates = 4
    Focc = 1.0*ones(Nstates,Nkpt)

    psik[1] = ortho_gram_schmidt( rand(Complex128,Ngwx,Nstates) )
    psik[2] = ortho_gram_schmidt( rand(Complex128,Ngwx,Nstates) )

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

#test_LDA_VWN()
#test_LDA_VWN_spinpol()
test_GGA_PBE()
test_GGA_PBE_spinpol()
