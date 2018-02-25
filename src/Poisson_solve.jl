#
# Given electron density in real space, return Hartree potential in reciprocal
# space
#
function Poisson_solve( pw::PWGrid, rhoR )
    #
    G2 = pw.gvec.G2
    Ns = pw.Ns
    Npoints = prod(Ns)
    #
    ctmp = 4.0*pi*R_to_G( pw, rhoR )
    #
    ctmp[1] = 0.0
    for ip = 2:Npoints
        ctmp[ip] = ctmp[ip]/G2[ip]
    end
    return ctmp
end
