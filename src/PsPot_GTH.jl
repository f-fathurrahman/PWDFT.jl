using SpecialFunctions: erf, gamma

struct PsPot_GTH
    pspfile::String
    atsymb::String
    zval::Int64
    rlocal::Float64
    rc::Array{Float64,1}  # indexed by l, originally [0:3]
    c::Array{Float64,1}   # coefficient in local pseudopotential
    h::Array{Float64,3}   # l,1:3,1:3
    lmax::Int64           # l = 0, 1, 2, 3 (s, p, d, f)
    Nproj_l::Array{Int64,1}  # originally 0:3
    rcut_NL::Array{Float64,1}  # originally 0:3, needed for real space evaluation
end

# Dummy PsPot_GTH
function PsPot_GTH()
    pspfile = "nothing"
    atsymb = "X"
    zval = 1.0  # to avoid error when generating electron info
    rlocal = 0.0
    rc = zeros(Float64,4)
    c = zeros(Float64,4) 
    h = zeros(Float64,4,3,3)
    lmax = -1
    Nproj_l = zeros(Int64,4)
    rcut_NL = zeros(Float64,4)
    return PsPot_GTH(pspfile, atsymb, zval, rlocal, rc, c, h, lmax, Nproj_l, rcut_NL)
end


"""
Initialize PsPot_GTH with parameters given in file given by path
`filename`.
"""
function PsPot_GTH( filename::String )
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

    return PsPot_GTH(filename, atsymb, zval, rlocal, rc, c, h, lmax, Nproj_l, rcut_NL)

end

"""
Evaluate GTH local pseudopotential in R-space
"""
function eval_Vloc_R( psp::PsPot_GTH, r::Array{Float64,1} )

    Npoints = size(r)[1]
    Vloc = zeros(Npoints)

    term1 = psp.c[1]
    for ip = 1:Npoints
        rrloc = r[ip]/psp.rlocal
        for i = 2:4
            term1 = term1 + psp.c[i]*rrloc^(2.0*(i-1))
        end
        Vloc[ip] = -psp.zval/r[ip] * erf( rrloc/sqrt(2.0) ) + exp(-0.5*rrloc^2)*term1
    end
    return Vloc
end


function eval_Vloc_R( psp::PsPot_GTH, r::Float64 )
    term1 = psp.c[1]
    rrloc = r/psp.rlocal
    for i = 2:4
        term1 = term1 + psp.c[i]*rrloc^(2.0*(i-1))
    end
    Vloc = -psp.zval/r * erf( rrloc/sqrt(2.0) ) + exp(-0.5*rrloc^2)*term1
    return Vloc
end


"""
Evaluate GTH local pseudopotential in G-space
"""
function eval_Vloc_G( psp::PsPot_GTH, G2::Array{Float64,1} )
    Ng = size(G2)[1]
    Vg = Array{Float64}(undef,Ng)
    for ig = 1:Ng
        Vg[ig] = eval_Vloc_G( psp, G2[ig] )
    end
    return Vg
end


function eval_Vloc_G( psp::PsPot_GTH, G2::Float64 )
    rloc = psp.rlocal
    zval = psp.zval
    c1 = psp.c[1]
    c2 = psp.c[2]
    c3 = psp.c[3]
    c4 = psp.c[4]

    pre1 = -4*pi*zval
    pre2 = sqrt(8*pi^3)*rloc^3
    Gr = sqrt(G2)*rloc
    expGr2 = exp(-0.5*Gr^2)

    SMALL = 1.0e-8
    #SMALL = eps()

    if sqrt(G2) > SMALL
        Vg = pre1/G2*expGr2 + pre2*expGr2 * (c1 + c2*(3.0 - Gr^2) +
             c3*(15.0 - 10.0*Gr^2 + Gr^4) + c4*(105.0 - 105.0*Gr^2 + 21.0*Gr^4 - Gr^6) )
    else
        Vg = 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15.0*c3 + 105.0*c4)  # to match QE
        #Vg = 0.0  # E_pspcore needs to be added later
    end

    return Vg
end


function eval_Vloc_G_short( psp::PsPot_GTH, G2::Float64 )
    rloc = psp.rlocal
    zval = psp.zval
    c1 = psp.c[1]
    c2 = psp.c[2]
    c3 = psp.c[3]
    c4 = psp.c[4]

    pre2 = sqrt(8*pi^3)*rloc^3
    Gr = sqrt(G2)*rloc
    expGr2 = exp(-0.5*Gr^2)

    if sqrt(G2) > eps()
        Vg = pre2*expGr2 * (c1 + c2*(3.0 - Gr^2) +
             c3*(15.0 - 10.0*Gr^2 + Gr^4) + c4*(105.0 - 105.0*Gr^2 + 21.0*Gr^4 - Gr^6) )
    else
        Vg = 0.0
    end

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
    fprj = sqrt(2) * rr * exp( -r^2 / (2*psp.rc[l+1]^2) ) /
           ( psp.rc[l+1]^(l + (4*i-1)/2) * x )
    return fprj
end



"""
Evaluate GTH projector function in G-space.
"""
function eval_proj_G( psp::PsPot_GTH, l::Int64, iproj::Int64, Gm::Array{Float64,1}, CellVolume::Float64 )

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
                Vprj[ig] = (4.0/3.0)/sqrt(105.0) * exp( -0.5*Gr2 ) * (15.0 - 10.0*Gr2 + Gr2^2)
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
                Vprj[ig] = (2.0/sqrt(105.0)) * exp(-0.5*Gr2) * Gm[ig]*(5.0 - Gr2)
            end
        elseif iproj == 3
            for ig = 1:Ng
                Gr2 = ( Gm[ig]*rrl)^2
                Vprj[ig] = (4.0/3.0)/sqrt(1155.0) * exp(-0.5*Gr2) * Gm[ig] * (35.0 - 14.0*Gr2 + Gr2^2)
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
                Vprj[ig] = (2.0/3.0)/sqrt(105.0) * exp(-0.5*Gr2) * Gm[ig]^2 * (7.0-Gr2)
            end
        end # if iproj

    # f-channel
    elseif l == 3
        # XXX only one projector
        for ig = 1:Ng
            Gr2 = ( Gm[ig]*rrl )^2
            Vprj[ig] = Gm[ig]^3 * exp(-0.5*Gr2)/sqrt(105.0)
        end

    end  # if l

    pre =  4 * pi^(5.0/4.0) * sqrt( 2.0^(l+1) * rrl^(2*l+3) / CellVolume )
    Vprj[:] = pre * Vprj[:]
    return Vprj
end



function eval_proj_G( psp::PsPot_GTH, l::Int64, iproj::Int64, Gm::Float64, CellVolume::Float64 )

    Vprj = 0.0

    rrl = psp.rc[l+1]

    # s-channel
    if l == 0
        if iproj==1
            Gr2 = ( Gm*rrl )^2
            Vprj = exp( -0.5*Gr2 )
        elseif iproj==2
            Gr2 = ( Gm*rrl )^2
            Vprj = 2.0/sqrt(15.0) * exp( -0.5*Gr2 ) * ( 3.0 - Gr2 )
        elseif iproj==3
            Gr2 = ( Gm*rrl )^2
            Vprj = (4.0/3.0)/sqrt(105.0) * exp( -0.5*Gr2 ) * (15.0 - 10.0*Gr2 + Gr2^2)
        end  # if iproj

    # p-channel
    elseif l == 1
        if iproj == 1
            Gr2 = ( Gm*rrl )^2
            Vprj = (1.0/sqrt(3.0)) * exp(-0.5*Gr2) * Gm
        elseif iproj == 2
            Gr2 = ( Gm*rrl )^2
            Vprj = (2.0/sqrt(105.0)) * exp(-0.5*Gr2) * Gm * (5.0 - Gr2)
        elseif iproj == 3
            Gr2 = (  Gm*rrl )^2
            Vprj = (4.0/3.0)/sqrt(1155.0) * exp(-0.5*Gr2) * Gm * (35.0 - 14.0*Gr2 + Gr2^2)
        end # if iproj

    # d-channel
    elseif l == 2
        if iproj == 1
            Gr2 = ( Gm*rrl )^2
            Vprj = (1.0/sqrt(15.0)) * exp(-0.5*Gr2) * Gm^2
        elseif iproj == 2
            Gr2 = ( Gm*rrl )^2
            Vprj = (2.0/3.0)/sqrt(105.0) * exp(-0.5*Gr2) * Gm^2 * (7.0 - Gr2)
        end # if iproj

    # f-channel
    elseif l == 3
        # XXX only one projector
        Gr2 = ( Gm*rrl )^2
        Vprj = Gm^3 * exp(-0.5*Gr2)/sqrt(105.0)

    end  # if l
    
    pre =  4.0 * pi^(5.0/4.0) * sqrt( 2.0^(l+1) * rrl^(2*l+3) / CellVolume )
    Vprj = pre * Vprj
    return Vprj
end

include("PsPot_GTH_io.jl")
