import Base: println

"""
Display some information about an instance of `Electrons`.
Occupation numbers are by default will only displayed for
several states only, except when `all_states=true`.
All occupation numbers will be displayed for all kpoints.
"""
function println( electrons::Electrons; header=true, all_states=false )

    Nspin = electrons.Nspin
    Focc = electrons.Focc
    Nelectrons = electrons.Nelectrons
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nkspin = size(Focc)[2]
    Nkpt = Nkspin/Nspin

    if header
        @printf("\n")
        @printf("                                    ---------\n")
        @printf("                                    Electrons\n")
        @printf("                                    ---------\n")
        @printf("\n")
    end
    @printf("Nspin         = %8d\n", Nspin)
    @printf("Nkpt          = %8d\n", Nkpt)
    @printf("Nelectrons    =  %18.10f\n", Nelectrons)
    @printf("Nstates       = %8d\n", Nstates)
    @printf("Nstates_occ   = %8d\n", Nstates_occ)
    @printf("Nstates_empty = %8d\n\n", Nstates - Nstates_occ)

    if Nspin == 1
        @printf("Occupation numbers: (spin-paired)\n\n")
    else
        @printf("Occupation numbers: (spin-polarized)\n\n")
    end

    if Nstates < 8
        all_states = true
    end

    Nk_per_line = 8

    if all_states
        for ist = 1:Nstates
            @printf("state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf("%8.5f ", Focc[ist,iks])
                if (iks % Nk_per_line) == 0
                    @printf("\n")
                    @printf("              ")
                end
            end
            #
            @printf("\n")
        end
    else
        for ist = 1:4
            @printf("state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf("%8.5f ", Focc[ist,iks])
                if (iks % Nk_per_line) == 0
                    @printf("\n")
                    @printf("              ")
                end
            end
            @printf("\n")
        end
        @printf(".....\n\n")
        #
        for ist = Nstates-3:Nstates
            @printf("state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf("%8.5f ", Focc[ist,iks])
                if (iks % Nk_per_line) == 0
                    @printf("\n")
                    @printf("              ")
                end
            end
            @printf("\n")
        end
    end
end


function print_ebands( electrons::Electrons, kpoints::KPoints; unit="hartree" )
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
        @printf("ik = %5d, k = (%18.10f,%18.10f,%18.10f)\n\n",
                ik, kpoints.k[1,ik], kpoints.k[2,ik], kpoints.k[3,ik])
        if electrons.Nspin == 2
            for ist = 1:Nstates
                @printf("%8d %13.10f %18.10f -- %13.10f %18.10f\n", ist,
                        Focc[ist,ik], ebands[ist,ik], Focc[ist,ik+Nkpt], ebands[ist,ik+Nkpt])
            end
        else
            for ist = 1:Nstates
                @printf("%8d %13.10f %18.10f\n", ist,
                        Focc[ist,ik], ebands[ist,ik])
            end
        end
        @printf("\n")
    end

    for ik = 1:Nkpt
        Focc[:,ik] = wk[ik]*Focc[:,ik]
    end
    if Nspin == 2
        for ik = 1:Nkpt
            Focc[:,ik+Nkpt] = wk[ik]*Focc[:,ik+Nkpt]
        end
    end
    @printf("sum(weighted Focc) = %18.10f\n", sum(Focc))
end
