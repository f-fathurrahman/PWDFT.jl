using LinearAlgebra
using Random
using Printf
using PWDFT

function test_GGA_PBE()

    @printf("---------------\n")
    @printf("Testing GGA PBE\n")
    @printf("---------------\n")

    ecutwfc_Ry = 30.0
    LatVecs = gen_lattice_sc(16.0)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )

    Random.seed!(1234)
    Ngwx = pw.gvecw.Ngwx
    dVol = pw.CellVolume/prod(pw.Ns)
    
    Nkpt = 1
    psik = Array{Array{ComplexF64,2},1}(undef,Nkpt)

    Nstates = 4
    Focc = 2.0*ones(Nstates,Nkpt)
    psik[1] = ortho_sqrt( rand(ComplexF64,Ngwx,Nstates) )

    Rhoe = calc_rhoe( pw, Focc, psik )
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( pw, Rhoe )
    @printf("sum Vxc = %18.10f\n", sum(Vxc))

    epsxc = calc_epsxc_PBE( pw, Rhoe )
    E_xc = dot( Rhoe, epsxc ) * dVol
    @printf("sum E_xc = %18.10f\n", E_xc)
end

function norm1_gauss(dr, sigma)
    dr2 = dot(dr,dr)
    c1 = 2*sigma^2
    cc1 = sqrt(2*pi*sigma^2)^3
    return exp(-dr2/c1)/cc1
end

function diff_norm1_gauss(dr, sigma, i::Int64)
    dr2 = dot(dr,dr)
    c1 = 2*sigma^2
    cc1 = sqrt(2*pi*sigma^2)^3
    return 2*dr[i]*exp(-dr2/c1)/(cc1*c1)
end

function test_op_nabla()
    ecutwfc_Ry = 1.0
    LatVecs = gen_lattice_sc(16.0)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )

    Npoints = prod(pw.Ns)
    f = zeros(Npoints)
    f_x = zeros(Npoints)
    f_y = zeros(Npoints)
    f_z = zeros(Npoints)
    center = [8.0, 8.0, 8.0]
    #
    sigma = 0.5
    for ip = 1:Npoints
        dr = pw.r[:,ip] - center[:]
        f[ip] = norm1_gauss(dr, sigma)
        f_x[ip] = diff_norm1_gauss(dr, sigma, 1)
        f_y[ip] = diff_norm1_gauss(dr, sigma, 2)
        f_z[ip] = diff_norm1_gauss(dr, sigma, 3)
    end

    dVol = pw.CellVolume/Npoints

    integ_f = sum(f)*dVol
    println("integ_f = ", integ_f)

    grad_f = op_nabla(pw, f)
    println("size grad_f = ", size(grad_f))

    println("integ grad_f = ", sum(grad_f)*dVol)
    println("integ f_x = ", sum(f_x)*dVol)

    println("diff x = ", sum(grad_f[1,:] - f_x[:]) )
    println("diff y = ", sum(grad_f[2,:] - f_y[:]) )
    println("diff z = ", sum(grad_f[3,:] - f_z[:]) )
end

function test_op_nabla_dot()
    ecutwfc_Ry = 40.0
    LatVecs = gen_lattice_sc(16.0)
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )

    Npoints = prod(pw.Ns)
    ff = zeros(3,Npoints)
    center = [8.0, 8.0, 8.0]
    #
    dff_x = zeros(Npoints)
    dff_y = zeros(Npoints)
    dff_z = zeros(Npoints)
    # P*i + Q*j + R*k
    #
    sigma_x = 0.5
    sigma_y = 0.65
    sigma_z = 0.75
    for ip = 1:Npoints
        dr = pw.r[:,ip] - center[:]
        ff[1,ip] = norm1_gauss(dr, sigma_x) # func of x, y, and z
        ff[2,ip] = norm1_gauss(dr, sigma_y) # func of x, y, and z
        ff[3,ip] = norm1_gauss(dr, sigma_z) # func of x, y, and z
        #
        dff_x[ip] = diff_norm1_gauss(dr, sigma_x, 1)
        dff_y[ip] = diff_norm1_gauss(dr, sigma_y, 2)
        dff_z[ip] = diff_norm1_gauss(dr, sigma_z, 3)
    end
    div_ff_analytic = dff_x + dff_y + dff_z

    dVol = pw.CellVolume/Npoints

    integ_f_x = sum(ff[1,:])*dVol
    println("integ_f_x = ", integ_f_x)

    integ_f_y = sum(ff[2,:])*dVol
    println("integ_f_y = ", integ_f_y)

    integ_f_z = sum(ff[3,:])*dVol
    println("integ_f_z = ", integ_f_z)

    div_ff = op_nabla_dot(pw, ff)
    println("size div_ff = ", size(div_ff))

    println("diff = ", sum(div_ff - div_ff_analytic))
end


test_op_nabla()

#test_op_nabla_dot()

#@time test_GGA_PBE()
#@time test_GGA_PBE()