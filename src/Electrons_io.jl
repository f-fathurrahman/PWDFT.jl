import Base: print

"""
Display some information about an instance of `Electrons`.
Occupation numbers are by default will only displayed for
several states only, except when `all_states=true`.
All occupation numbers will be displayed for all kpoints.
"""
function print( io::IO, electrons::Electrons; header=true, all_states=false )

    Nspin = electrons.Nspin
    Focc = electrons.Focc
    Nelectrons = electrons.Nelectrons
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nkspin = size(Focc)[2]
    Nkpt = Nkspin/Nspin

    if header
        @printf(io, "\n")
        @printf(io, "                                    ---------\n")
        @printf(io, "                                    Electrons\n")
        @printf(io, "                                    ---------\n")
        @printf(io, "\n")
    end
    @printf(io, "Nspin         = %8d\n", Nspin)
    @printf(io, "Nkpt          = %8d\n", Nkpt)
    @printf(io, "Nelectrons    =  %18.10f\n", Nelectrons)
    @printf(io, "Nstates       = %8d\n", Nstates)
    @printf(io, "Nstates_occ   = %8d\n", Nstates_occ)
    @printf(io, "Nstates_empty = %8d\n\n", Nstates - Nstates_occ)

    if Nspin == 1
        @printf(io, "Occupation numbers: (spin-paired)\n\n")
    else
        @printf(io, "Occupation numbers: (spin-polarized)\n\n")
    end

    if Nstates < 8
        all_states = true
    end

    Nk_per_line = 8

    if all_states
        for ist = 1:Nstates
            @printf(io, "state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf(io, "%8.5f ", Focc[ist,iks])
                if (iks % Nk_per_line) == 0
                    @printf(io, "\n")
                    @printf(io, "              ")
                end
            end
            #
            @printf(io, "\n")
        end
    else
        for ist = 1:4
            @printf(io, "state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf(io, "%8.5f ", Focc[ist,iks])
                if (iks % Nk_per_line) == 0
                    @printf(io, "\n")
                    @printf(io, "              ")
                end
            end
            @printf(io, "\n")
        end
        @printf(io, ".....\n\n")
        #
        for ist = Nstates-3:Nstates
            @printf(io, "state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf(io, "%8.5f ", Focc[ist,iks])
                if (iks % Nk_per_line) == 0
                    @printf(io, "\n")
                    @printf(io, "              ")
                end
            end
            @printf(io, "\n")
        end
    end
end
print( electrons::Electrons; header=true, all_states=false ) = print( stdout, electrons, header=header, all_states=all_states )


function print_ebands( io::IO, electrons::Electrons, kpoints::KPoints; unit="hartree" )
    if unit == "eV"
        ebands = electrons.ebands*2*Ry2eV
    else
        ebands = electrons.ebands
    end

    Nspin = electrons.Nspin
    Focc = copy(electrons.Focc)
    Nkspin = size(Focc)[2]
    Nkpt = kpoints.Nkpt
    wk = kpoints.wk
        
    Nstates = electrons.Nstates

    for ik = 1:Nkpt
        @printf(io, "ik = %5d, k = (%18.10f,%18.10f,%18.10f)\n\n",
                ik, kpoints.k[1,ik], kpoints.k[2,ik], kpoints.k[3,ik])
        if electrons.Nspin == 2
            for ist = 1:Nstates
                @printf(io, "%8d %13.10f %18.10f -- %13.10f %18.10f\n", ist,
                        Focc[ist,ik], ebands[ist,ik], Focc[ist,ik+Nkpt], ebands[ist,ik+Nkpt])
            end
        else
            for ist = 1:Nstates
                @printf(io, "%8d %13.10f %18.10f\n", ist,
                        Focc[ist,ik], ebands[ist,ik])
            end
        end
        @printf(io, "\n")
    end

    for ik = 1:Nkpt
        Focc[:,ik] = wk[ik]*Focc[:,ik]
    end
    if Nspin == 2
        for ik = 1:Nkpt
            Focc[:,ik+Nkpt] = wk[ik]*Focc[:,ik+Nkpt]
        end
    end
    @printf(io, "sum(weighted Focc) = %18.10f\n", sum(Focc))
end
print_ebands( electrons::Electrons, kpoints::KPoints; unit="hartree" ) = print_ebands( stdout, electrons, kpoints, unit=unit )
