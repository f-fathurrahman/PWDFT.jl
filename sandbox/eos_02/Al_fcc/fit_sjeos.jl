using Printf
import Polynomials
using DelimitedFiles: readdlm
using PWDFT: Ry2eV, ANG2BOHR

function main()
    
    dat = readdlm("TEMP_EOS_data.dat")
    
    # calculate primitive fcc volumes and convert it to angstrom^3
    volumes = 0.25*dat[:,1].^3 / (ANG2BOHR^3)
    
    # energies in eV
    energies = dat[:,2]*2*Ry2eV
    energies_abinit = dat[:,3]*2*Ry2eV
    energies_pwscf = dat[:,4]*2*Ry2eV

    V0, E0, B_GPa = do_fit_sjeos( volumes, energies )
    @printf("Bulk modulus (PWDFT)  = %18.10f GPa\n", B_GPa)

    V0, E0, B_GPa = do_fit_sjeos( volumes, energies_abinit )
    @printf("Bulk modulus (ABINIT) = %18.10f GPa\n", B_GPa)

    V0, E0, B_GPa = do_fit_sjeos( volumes, energies_pwscf )
    @printf("Bulk modulus (PWSCF)  = %18.10f GPa\n", B_GPa)

end

function do_fit_sjeos( volumes, energies )

    # fit to 3rd order polynomial in V^(-1/3)
    pfit0 = Polynomials.polyfit( volumes.^(-1/3), energies, 3 )
    pfit1 = Polynomials.polyder(pfit0)
    pfit2 = Polynomials.polyder(pfit1)

    rs = Polynomials.roots(pfit1)
    t = 0.0
    for r in rs
        if (r > 0.0) && (pfit2(r) > 0.0)
            t = r
            break
        end
    end

    V0 = t^(-3)
    E0 = pfit0(t)
    B = t^5 * pfit2(t)/9.0

    kJ = 6.241509125883258e+21

    B_GPa = B/kJ*1.0e24

    return V0, E0, B_GPa

end

main()