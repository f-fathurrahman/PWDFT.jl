#
# Various utilities to read PWDFT.jl log
#

function read_pwdft_Natoms( filename::String )
    Natoms = 0
    f = open(filename, "r")
    while !eof(f)
        l = readline(f)
        if occursin("Natoms", l)
            Natoms = parse( Int64, split(l, "=")[2] )
        end
    end
    close(f)
    return Natoms
end

function read_pwdft_energies( filename::String )
    E_kin      = 0.0
    E_hartree  = 0.0
    E_xc       = 0.0
    E_ewald    = 0.0
    E_pspCore  = 0.0
    E_Ps_loc   = 0.0
    E_Ps_nloc  = 0.0
    mTS        = 0.0
    E_total    = 0.0

    f = open(filename, "r")
    while !eof(f)
        l = readline(f)
        if occursin("Kinetic    energy", l)
            E_kin = parse( Float64, split(l, ":")[2] )
        end
        if occursin("Hartree    energy", l)
            E_hartree = parse( Float64, split(l, ":")[2] )
        end
        if occursin("XC         energy", l)
            E_xc = parse( Float64, split(l, ":")[2] )
        end        
        if occursin("NN         energy", l)
            E_ewald = parse( Float64, split(l, ":")[2] )
        end
        if occursin("PspCore    energy", l)
            E_pspCore = parse( Float64, split(l, ":")[2] )
        end
        if occursin("Ps_loc     energy", l)
            E_Ps_loc = parse( Float64, split(l, ":")[2] )
        end
        if occursin("Ps_nloc    energy", l)
            E_Ps_nloc = parse( Float64, split(l, ":")[2] )
        end
        if occursin("-TS ", l)
            mTS = parse( Float64, split(l, ":")[2] )
        end

        # only of these will be parsed
        if occursin("Total free energy:", l)
            E_total = parse( Float64, split(l, ":")[2] )
        end
        if occursin("Total      energy:", l)
            E_total = parse( Float64, split(l, ":")[2] )
        end
    end
    close(f)

    energies = Energies( E_kin, E_Ps_loc, E_Ps_nloc, 
                         E_hartree, E_xc, E_ewald, E_pspCore, mTS )
    return energies

end