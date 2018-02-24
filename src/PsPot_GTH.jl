struct PsPot_GTH
    atsymb::String
    zval::Int64
    rlocal::Float64
    rc::Array{Float64,1}  # indexed by l, originally [0:3]
    c::Array{Float64,1}   # coefficient in local pseudopotential
    h::Array{Float64,3}   # l,1:3,1:3
    lmax::Int64           # l = 0, 1, 2, 3 (s, p, d, f)
    Nproj_l::Array{Int64,1}  # originally 0:3
    rcut_NL::Array{Int64,1}  # originally 0:3, needed for real space evaluation
end


function PsPot_GTH( filename )
    c = zeros(Float64,4)
    rlocal = 0.0
    rc = zeros(Float64,4)
    h = zeros(Float64,4,3,3)
    Nproj_l = zeros(Int64,4)
    rcut_NL = zeros(Float64,4)

    f = open(filename)

    line = readline(f)
    l = split(line)
    atsymb = l[1]

    line = readline(f)
    l = split(line)
    n_s = parse( Int64, l[1] )
    n_p = parse( Int64, l[2] )
    n_d = parse( Int64, l[3] )
    n_f = parse( Int64, l[4] )

    zval = n_s + n_p + n_d + n_f

    line = readline(f)
    l = split(line)
    rlocal = parse( Float64, l[1] )
    n_c_local = parse( Int64, l[2] )

    line = readline(f)
    l = split(line)
    for i = 1:n_c_local
        c[i] = parse( Float64, l[i] )
    end

    line = readline(f)
    l = split(line)
    lmax = parse( Int64, l[1] ) - 1


    # l = 1 -> s
    # l = 2 -> p
    # l = 3 -> d
    # l = 4 -> f

    # lmax is however use the usual physics convention

    for i = 1:lmax+1
        line = readline(f)
        l = split(line)
        rc[i] = parse( Float64, l[1] )
        Nproj_l[i] = parse( Float64, l[2] )
        for ii = 1:Nproj_l[i]
            line = readline(f)
            l = split(line)
            iread = 0
            for iii = ii:Nproj_l[i]
                iread = iread + 1
                h[i,ii,iii] = parse( Float64, l[iread] )
            end
        end
    end

    close(f)

    for k = 1:4
        for i = 1:3
            for j = i+1:3
                h[k,j,i] = h[k,i,j]
            end
        end
    end

    return PsPot_GTH(atsymb, zval, rlocal, rc, c, h, lmax, Nproj_l, rcut_NL)

end


# Display information about GTH pseudopotential
import Base.println
function println( psp::PsPot_GTH )

    const ANGMOM = ["s", "p", "d", "f"]

    Nproj_l = psp.Nproj_l
    c = psp.c
    rlocal = psp.rlocal
    rc = psp.rc
    h = psp.h

    @printf("\nLocal pseudopotential info:\n\n")
    @printf("rloc: %f, c: %f, %f, %f, %f\n", rlocal, c[1], c[2], c[3], c[4])
    @printf("\n")
    @printf("Nonlocal pseudopotential info:\n\n")
    for i=1:psp.lmax+1
        @printf("Angular momentum: %s, rc = %f\n", ANGMOM[i], rc[i])
        @printf("h = \n")
        for pi = 1:Nproj_l[i]
            for pj = 1:Nproj_l[i]
                @printf("%18.10f ", h[i,pi,pj])
            end
            @printf("\n")
        end
        @printf("\n")
    end

end


"""
Evaluate GTH local pseudopotential in R-space
"""
function eval_Vloc_R( psp::PsPot_GTH, r::Array{Float64,2} )

    Npoints = size(r)[1]
    Vloc = zeros(Npoints)

    term1 = psp.c[1]
    for ip = 1:Npoints
        rrloc = r[ip]/psp.rloc
        for i = 2:4
            term1 = term1 + psp.c[i]*rrloc^(2.0*(i-1))
        end
        Vloc[ip] = -psp.zval/r[ip] * erf( rrloc/sqrt(2.0) ) + exp(-0.5*rrloc^2)*term1
    end
    return Vloc
end


"""
Evaluate GTH local pseudopotential in G-space
"""
function eval_Vloc_G( psp::PsPot_GTH, G2::Float64, Ω::Float64 )

    rloc = psp.rlocal
    zval = psp.zval
    c1 = psp.c[1]
    c2 = psp.c[2]
    c3 = psp.c[3]
    c4 = psp.c[4]

    pre1 = -4*pi*zval/Ω
    pre2 = sqrt(8*pi^3)*rloc^3/Ω
    Gr = sqrt(G2)*rloc
    expGr2 = exp(-0.5*Gr^2)

    const SMALL = eps()

    if sqrt(G2) > SMALL
        Vg = pre1/G2*expGr2 + pre2*expGr2 * (c1 + c2*(3-Gr^2) +
             c3*(15 - 10*Gr^2 + Gr^4) + c4*(105 - 105*Gr^2 + 21*Gr^4 - Gr^6) )
    else
        # println("Small G2 = ", G2)
        # limiting value
        Vg = 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15*c3 + 105*c4)
    end

    return Vg
end



"""
Evaluate GTH local pseudopotential in G-space
"""
function eval_Vloc_G( psp::PsPot_GTH, G2::Array{Float64,1}, Ω::Float64 )

    Ng = size(G2)[1]
    Vg = zeros(Ng)

    rloc = psp.rlocal
    zval = psp.zval
    c1 = psp.c[1]
    c2 = psp.c[2]
    c3 = psp.c[3]
    c4 = psp.c[4]

    pre1 = -4*pi*zval/Ω
    pre2 = sqrt(8*pi^3)*rloc^3/Ω
    #
    for ig=2:Ng
        Gr = sqrt(G2[ig])*rloc
        expGr2 = exp(-0.5*Gr^2)
        Vg[ig] = pre1/G2[ig]*expGr2 + pre2*expGr2 * (c1 + c2*(3-Gr^2) +
                    c3*(15 - 10*Gr^2 + Gr^4) + c4*(105 - 105*Gr^2 + 21*Gr^4 - Gr^6) )
    end
    # limiting value, with minus sign ?
    Vg[1] = 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15*c3 + 105*c4)

    return Vg
end


"""
Evaluate GTH projector function in R-space.
"""
function eval_proj_R( psp::PsPot_GTH, l, i, r::Float64 )
    x = sqrt( gamma( l + (4*i-1)/2.0 ) )
    if l==0 & i==1
        rr = 1.0
    else
        rr = r^(l + 2*(i-1))
    end
    fprj = sqrt(2) * rr * exp(-r^2/(2*psp.rc[l+1])) /
           ( psp.rc[l+1]^(l + (4*i-1)/2) * x )
    return fprj
end



"""
Evaluate GTH projector function in G-space.
"""
function eval_proj_G( psp::PsPot_GTH, l, iproj, Gm, Ω )

    # NOTICE that Gm is the magnitudes of G-vectors
    Ng = size(Gm)[1]

    Vprj = zeros(Ng)

    rrl = psp.rc[l+1]

    # s-channel
    if l == 0
        if iproj==1
            for ig = 1:Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = exp( -0.5*Gr2 )
            end
        elseif iproj==2
            for ig = 1:Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = 2.0/sqrt(15.0) * exp( -0.5*Gr2 ) * ( 3.0 - Gr2 )
            end
        elseif iproj==3
            for ig = 1:Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = (4.0/3.0)/sqrt(105.0) * exp( -0.5*Gr2 ) * (15.0 - 10.*Gr2 + Gr2^2)
            end
        end  # if iproj

    # p-channel
    elseif l == 1
        if iproj == 1
            for ig = 1:Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = (1.0/sqrt(3.0)) * exp(-0.5*Gr2) * Gm[ig]
            end
        elseif iproj == 2
            for ig = 1:Ng
                Gr2 = (Gm[ig]*rrl)^2
                Vprj[ig] = (2./sqrt(105.)) * exp(-0.5*Gr2) * Gm[ig]*(5. - Gr2)
            end
        elseif iproj == 3
            for ig = 1:Ng
                Gr2 = ( Gm[ig]*rrl)^2
                Vprj[ig] = (4./3.)/sqrt(1155.) * exp(-0.5*Gr2) * Gm[ig] * (35. - 14.*Gr2 + Gr2^2)
            end
        end # if iproj

    # d-channel
    elseif l == 2
        if iproj == 1
            for ig = 1:Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = (1.0/sqrt(15.0)) * exp(-0.5*Gr2) * Gm[ig]^2
            end
        elseif iproj == 2
            for ig = 1:Ng
                Gr2 = (Gm[ig]*rrl)^2
                Vprj[ig] = (2./3.)/sqrt(105.) * exp(-0.5*Gr2) * Gm[ig]^2 * (7.-Gr2)
            end
        end # if iproj

    # f-channel
    elseif l == 3
        # XXX only one projector
        for ig = 1:Ng
            Gr2 = ( Gm[ig]*rrl )^2
            Vprj[ig] = Gm[ig]^3 * exp(-0.5*Gr2)
        end

    end  # if l

    pre =  4 * pi^(5./4.) * sqrt( 2.^(l+1) * rrl^(2*l+3) / Ω )
    Vprj[:] = pre * Vprj[:]
    return Vprj
end
