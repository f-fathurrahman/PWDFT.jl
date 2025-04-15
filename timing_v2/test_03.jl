#=
NOTES:
Most difference in time for initializing Hamiltonian is here.
USPP and PAW gave significant longer times than NCPP.
=#

function test_03( Ham )

    if eltype(Ham.pspots) == PsPot_UPF
        is_gga = false # XXX FIX THIS !!!!
        Nspin = 1
        res = @be PsPotNL_UPF( Ham.atoms, Ham.pw, Ham.pspots, is_gga=is_gga, Nspin=Nspin )
        display(res)
    elseif eltype(Ham.pspots) == PsPot_GTH
        res = @be PsPotNL( Ham.atoms, Ham.pw, Ham.pspots, check_norm=false )
        display(res)
    else
        error("Not supporting mixed pseudopotential types: GTH and UPF")
    end

    return
end
