import Base: show
function show( io::IO, energies::Energies; use_smearing=false, is_paw=false )

    @printf(io, "Kinetic    energy: %18.10f\n", energies.Kinetic )
    @printf(io, "Ps_loc     energy: %18.10f\n", energies.Ps_loc )
    @printf(io, "Ps_nloc    energy: %18.10f\n", energies.Ps_nloc )
    @printf(io, "Hartree    energy: %18.10f\n", energies.Hartree )
    @printf(io, "XC         energy: %18.10f\n", energies.XC )

    is_paw && 
    @printf(io, "EHxc_paw   energy: %18.10f\n", energies.EHxc_paw)

    if ( abs(energies.mTS) > eps() ) || use_smearing
        @printf(io, "-TS              : %18.10f\n", energies.mTS)
    end

    @printf(io, "-------------------------------------\n")
    
    E_elec = energies.Kinetic + energies.Ps_loc + energies.Ps_nloc +
             energies.Hartree + energies.XC + energies.mTS
    
    is_paw && ( E_elec += energies.EHxc_paw )

    @printf(io, "Electronic energy: %18.10f\n", E_elec)
    @printf(io, "NN         energy: %18.10f\n", energies.NN )
    @printf(io, "-------------------------------------\n")
    
    E_total = E_elec + energies.NN
    
    if use_smearing
        @printf(io, "Total free energy: %18.10f\n", E_total)
        @printf(io, "\n")
        @printf(io, "Total energy (extrapolated to T=0): %18.10f\n", E_total - 0.5*energies.mTS)
    else
        @printf(io, "Total      energy: %18.10f\n", E_total )
    end

    # core energy for PAW case is not yet accounted for

end
show( energies::Energies; use_smearing=false, is_paw=false ) = show( stdout,
    energies, use_smearing=use_smearing, is_paw=is_paw )

import Base: println
println( energies::Energies; use_smearing=false, is_paw=false ) = show(
    energies, use_smearing=use_smearing, is_paw=is_paw )