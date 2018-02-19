using PWDFT
using PWDFT.Utilities

using SpecialFunctions: erf

mutable struct PsPot_HGH
    itype::Int
    atsymb::String
    zval::Float64
    lloc::Int
    lmax::Int
    rloc::Float64
    rc::Array{Float64}   # indexed (l+1), l=0,1,2,3
    c::Array{Float64}    # indexed 1,2,3,4
    h::Array{Float64,3}  # originally indexed [0:3,1:3,1:3]
    k::Array{Float64,3}  # indexed [0:3,1:3,1:3]
    nprj::Array{Int}
    snprj::Int           # snprj = sum(nprj)
    lll::Array{Int}
    ipr::Array{Int}
end


# Constructor
function PsPot_HGH( itype::Int, atsymb::String, filename::String; verbose=false )

    if verbose
        @printf("\nFile = %s\n", filename)
    end
    file = open( filename )

    comment = readline(file)

    lines = split( readline(file) )
    zval = parse( Float64, lines[1] )
    zion = parse( Float64, lines[2] )

    lines  = split( readline(file) )
    pspcod = parse( Int, lines[1] )
    pspxc  = parse( Int, lines[2] )
    lmax = parse( Int, lines[3] )
    lloc = parse( Int, lines[4] )
    mmax   = parse( Int, lines[5] )
    r2well = parse( Int, lines[6] )

    rc = zeros(Float64,4)
    c = zeros(Float64,4)
    h = zeros(Float64,4,3,3)
    k = zeros(Float64,4,3,3)

    lines = split( readline(file) )
    rloc = parse( Float64, lines[1] )
    c[1] = parse( Float64, lines[2] )
    c[2] = parse( Float64, lines[3] )
    c[3] = parse( Float64, lines[4] )
    c[4] = parse( Float64, lines[5] )

    const ANGMOM = ["s", "p", "d", "f"]

    l = 0  # s
    lines = split( readline(file) )
    #
    rc[1] = parse( lines[1] )
    #
    h[1,1,1] = parse( lines[2] )
    h[1,2,2] = parse( lines[3] )
    h[1,3,3] = parse( lines[4] )

    for l = 1:3
        lines = split( readline(file) )
        #
        rc[l+1] = parse( lines[1] )
        #
        h[l+1,1,1] = parse( lines[2] )
        h[l+1,2,2] = parse( lines[3] )
        h[l+1,3,3] = parse( lines[4] )
        #
        lines = split( readline(file) )
        k[l+1,1,1] = parse( lines[1] )
        k[l+1,2,2] = parse( lines[2] )
        k[l+1,3,3] = parse( lines[3] )
        #
    end

    # additional info
    lines = split(readline(file))
    nprj = zeros(Int,4)
    nprj[1] = parse( Int, lines[1] )
    nprj[2] = parse( Int, lines[2] )
    nprj[3] = parse( Int, lines[3] )
    nprj[4] = parse( Int, lines[4] )

    snprj  = sum(nprj)
    lll = zeros(Int,snprj)
    ipr = zeros(Int,snprj)

    iv = 0
    for il = 0:lmax
        for ii = 1:nprj[il+1]
            iv = iv + 1
            lll[iv] = il
            ipr[iv] = ii
        end
    end

    close(file)

    const M_HALF = 0.5
    const M_THREE = 3.0
    const M_FIVE = 5.0
    const M_ONE = 1.0

    # from Octopus code, see also the appendix of HGH paper

    h[0+1, 1, 2] = -M_HALF    * sqrt(M_THREE/M_FIVE) * h[0+1, 2, 2]
    h[0+1, 1, 3] =  M_HALF    * sqrt(M_FIVE/21.0)    * h[0+1, 3, 3]
    h[0+1, 2, 3] = -M_HALF    * sqrt(100.0/63.0)     * h[0+1, 3, 3]
    h[1+1, 1, 2] = -M_HALF    * sqrt(M_FIVE/7.0)     * h[1+1, 2, 2]
    h[1+1, 1, 3] =  M_ONE/6.0 * sqrt(35.0/11.0)      * h[1+1, 3, 3]
    h[1+1, 2, 3] = -M_ONE/6.0 * ( 14.0 / sqrt(11.0)) * h[1+1, 3, 3]
    h[2+1, 1, 2] = -M_HALF    * sqrt(7.0/9.0)        * h[2+1, 2, 2]
    h[2+1, 1, 3] =  M_HALF    * sqrt(63.0/143.0)     * h[2+1, 3, 3]
    h[2+1, 2, 3] = -M_HALF    * (18.0/sqrt(143.0))   * h[2+1, 3, 3]

    k[0+1, 1, 2] = -M_HALF    * sqrt(M_THREE/M_FIVE) * k[0+1, 2, 2]
    k[0+1, 1, 3] =  M_HALF    * sqrt(M_FIVE/21.0)    * k[0+1, 3, 3]
    k[0+1, 2, 3] = -M_HALF    * sqrt(100.0/63.0)     * k[0+1, 3, 3]
    k[1+1, 1, 2] = -M_HALF    * sqrt(M_FIVE/7.0)     * k[1+1, 2, 2]
    k[1+1, 1, 3] =  M_ONE/6.0 * sqrt(35.0/11.0)      * k[1+1, 3, 3]
    k[1+1, 2, 3] = -M_ONE/6.0 * (14.0 / sqrt(11.0))  * k[1+1, 3, 3]
    k[2+1, 1, 2] = -M_HALF    * sqrt(7.0/9.0)        * k[2+1, 2, 2]
    k[2+1, 1, 3] =  M_HALF    * sqrt(63.0/143.0)     * k[2+1, 3, 3]
    k[2+1, 2, 3] = -M_HALF    * (18.0/sqrt(143.0))   * k[2+1, 3, 3]

    # Parameters are symmetric.
    # be careful of index !!!
    for ik = 1:4
        for i = 1:3
            for j = i+1:3
                h[ik, j, i] = h[ik, i, j]
                k[ik, j, i] = k[ik, i, j]
            end
        end
    end

    psp = PsPot_HGH( itype, atsymb, zval, lloc, lmax, rloc,
                     rc, c, h, k, nprj, snprj, lll, ipr )

    return psp
end

# Display information about HGH pseudopotential
function info_PsPot_HGH( psp::PsPot_HGH )

    const ANGMOM = ["s", "p", "d", "f"]

    @printf("\nLocal pseudopotential info:\n\n")
    @printf("rloc: %f, c: %f, %f, %f, %f\n", psp.rloc, psp.c[1], psp.c[2], psp.c[3], psp.c[4])
    @printf("\n")
    @printf("Nonlocal pseudopotential info:\n\n")
    for i=1:4
        @printf("Angular momentum: %s, rc = %f\n", ANGMOM[i], psp.rc[i])
        @printf("h = \n")
        PrintMatrix( reshape(psp.h[i,:,:],(3,3) ) )
        @printf("\n")
    end

end


"""
Evaluate HGH local pseudopotential in R-space
"""
function HGH_eval_Vloc_R( psp, r::Array{Float64,2} )

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
Evaluate HGH local pseudopotential in G-space
"""
function HGH_eval_Vloc_G( psp, G2, Ω )

    Ng = size(G2)[1]
    Vg = zeros(Ng)

    rloc = psp.rloc
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
Evaluate HGH projector function in R-space.
"""
function HGH_eval_proj_R( psp, l, i, r::Float64 )
    x = sqrt( gamma( l + (4*i-1)/2.0 ) )
    if l==0 & i==1
        rr = 1.0
    else
        rr = r^(l + 2*(i-1))
    end
    fprj = sqrt(2) * rr * exp(-r^2/(2*psp.rc[l+1])) /
           ( p.rc[l+1]^(l + (4*i-1)/2) * x )
    return fprj
end



"""
Evaluate HGH projector function in G-space.
"""
function HGH_eval_proj_G( psp, l, iproj, Gm, Ω )

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
