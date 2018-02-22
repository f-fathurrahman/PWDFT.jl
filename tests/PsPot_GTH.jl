using PWDFT: print_matrix


struct PsPot_GTH
    atsymb::String
    zval::Int64
    rlocal::Float64
    rc::Array{Float64,1}  # indexed by l, originally [0:3]
    c::Array{Float64,1}   # coefficient in local pseudopotential
    h::Array{Float64,3}   # l,1:3,1:3
    lmax::Int64
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
    println("zval = ", zval)

    line = readline(f)
    l = split(line)
    rlocal = parse( Float64, l[1] )
    n_c_local = parse( Int64, l[2] )

    println("rlocal = ", rlocal)
    println("n_c_local = ", n_c_local)

    line = readline(f)
    l = split(line)
    for i = 1:n_c_local
        c[i] = parse( Float64, l[i] )
    end
    println("c = ", c[1:n_c_local])

    line = readline(f)
    l = split(line)
    lmax = parse( Int64, l[1] ) - 1

    println("lmax = ", lmax)

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
        print_matrix(h[k,:,:])
    end

    return PsPot_GTH(atsymb, zval, rlocal, rc, c, h, lmax, Nproj_l, rcut_NL)

end


# Display information about HGH pseudopotential
import Base.println
function println( psp::PsPot_GTH )

    const ANGMOM = ["s", "p", "d", "f"]

    @printf("\nLocal pseudopotential info:\n\n")
    @printf("rloc: %f, c: %f, %f, %f, %f\n", psp.rlocal, psp.c[1], psp.c[2], psp.c[3], psp.c[4])
    @printf("\n")
    @printf("Nonlocal pseudopotential info:\n\n")
    for i=1:4
        @printf("Angular momentum: %s, rc = %f\n", ANGMOM[i], psp.rc[i])
        @printf("h = \n")
        print_matrix( reshape(psp.h[i,:,:],(3,3) ) )
        @printf("\n")
    end

end
