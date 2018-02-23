#
# Given electron density in real space, return Hartree potential in reciprocal
# space
#
function Poisson_solve( PW::PWGrid, rhoR )
    #
    G2 = PW.gvec.G2
    Ns = PW.Ns
    Npoints = prod(Ns)
    #
    ctmp = 4.0*pi*c_R_to_G( Ns, rhoR )
    #
    ctmp[1] = 0.0
    for ip = 2:Npoints
        ctmp[ip] = ctmp[ip]/G2[ip]
    end
    return ctmp
end
