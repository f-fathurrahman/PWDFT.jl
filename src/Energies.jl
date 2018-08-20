"""
The type for collecting all energy terms: kinetic, local and
nonlocal pseudopotential, Hartree, XC, nuclei-nuclei terms
and electronic entropy term.
"""
mutable struct Energies
    Kinetic::Float64
    Ps_loc::Float64
    Ps_nloc::Float64
    Hartree::Float64
    XC::Float64
    NN::Float64
    mTS::Float64
end

"""
Creates an instance of `Energies` with value of zeros for all fields.
"""
function Energies()
    return Energies(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

import Base: println, sum

function println( energies::Energies; use_smearing=false )
    
    @printf("Kinetic    energy: %18.10f\n", energies.Kinetic )
    @printf("Ps_loc     energy: %18.10f\n", energies.Ps_loc )
    @printf("Ps_nloc    energy: %18.10f\n", energies.Ps_nloc )
    @printf("Hartree    energy: %18.10f\n", energies.Hartree )
    @printf("XC         energy: %18.10f\n", energies.XC )
    if use_smearing
        @printf("-TS              : %18.10f\n", energies.mTS)
    end
    @printf("-------------------------------------\n")
    E_elec = energies.Kinetic + energies.Ps_loc + energies.Ps_nloc +
             energies.Hartree + energies.XC + energies.mTS
    @printf("Electronic energy: %18.10f\n", E_elec)
    @printf("NN         energy: %18.10f\n", energies.NN )
    @printf("-------------------------------------\n")
    E_total = E_elec + energies.NN
    if use_smearing
        @printf("Total free energy: %18.10f\n", E_total)
        @printf("\n")
        @printf("Total energy (extrapolated to T=0): %18.10f\n", E_total - 0.5*energies.mTS)
    else
        @printf("Total      energy: %18.10f\n", E_total )
    end
end

"""
Get total energy by summing all of its components.
"""
function sum( energies::Energies )
    # irrespective of whether we are using smearing or not
    return energies.Kinetic + energies.Ps_loc + energies.Ps_nloc +
           energies.Hartree + energies.XC + energies.NN + energies.mTS
end

