mutable struct Electrons
    Nelectrons::Float64
    Nstates::Int64
    Nstates_occ::Int64
    Focc::Array{Float64,2}
    ebands::Array{Float64,2}
    Nspin::Int64
end

# dummy Electrons
function Electrons()
    Nelectrons = 1
    Nstates = 1
    Nstates_occ = 1
    Focc = zeros(Nstates,1) # Nkpt=1
    ebands = zeros(Nstates,1) # use Nkpt=1
    Nspin = 1
    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )
end


function Electrons( atoms::Atoms, Pspots::Array{PsPot_GTH,1};
                    Nspin=1, Nkpt=1,
                    Nstates=nothing, Nstates_empty=0 )
    
    @assert( Nspin <= 2 )

    Nelectrons = get_Nelectrons(atoms,Pspots)

    is_odd = round(Int64,Nelectrons)%2 == 1

    # If Nstates is not specified and Nstates_empty == 0, we calculate
    # Nstates manually from Nelectrons
    if (Nstates == nothing)
        Nstates = round( Int64, Nelectrons/2 )
        if Nstates*2 < Nelectrons
            Nstates = Nstates + 1
        end
        if Nstates_empty > 0
            Nstates = Nstates + Nstates_empty
        end
    end

    Focc = zeros(Float64,Nstates,Nkpt*Nspin)
    ebands = zeros(Float64,Nstates,Nkpt*Nspin)
    
    Nstates_occ = Nstates - Nstates_empty
    
    if Nspin == 1
        for ist = 1:Nstates_occ-1
            Focc[ist,:] .= 2.0
        end
        if is_odd
            Focc[Nstates_occ,:] .= 1.0
        else
            Focc[Nstates_occ,:] .= 2.0
        end
    else
        for ist = 1:Nstates_occ-1
            Focc[ist,:] .= 1.0
        end
        idx_up = 1:Nkpt
        if is_odd
            # assign only to the spin up
            Focc[Nstates_occ,idx_up] .= 1.0
        else
            Focc[Nstates_occ,:] .= 1.0
        end
    end

    sFocc = sum(Focc)/Nkpt
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        @printf("ERROR: diff sum(Focc) and Nelectrons is not small\n")
        @printf("sum Focc = %f, Nelectrons = %f\n", sFocc, Nelectrons)
        exit()
    end

    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )
end


function Electrons( atoms::Atoms, Zvals::Array{Float64,1};
                    Nspin=1, Nkpt=1,
                    Nstates=nothing, Nstates_empty=0 )

    @assert( Nspin <= 2 )

    Nelectrons = 0.0
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    for ia = 1:Natoms
        isp = atm2species[ia]
        Nelectrons = Nelectrons + Zvals[isp]
    end

    is_odd = round(Int64,Nelectrons)%2 == 1

    # If Nstates is not specified and Nstates_empty == 0, we calculate
    # Nstates manually from Nelectrons
    if (Nstates == nothing)
        Nstates = round( Int64, Nelectrons/2 )
        if Nstates*2 < Nelectrons
            Nstates = Nstates + 1
        end
        if Nstates_empty > 0
            Nstates = Nstates + Nstates_empty
        end
    end

    Focc = zeros(Float64,Nstates,Nkpt*Nspin)
    ebands = zeros(Float64,Nstates,Nkpt*Nspin)

    Nstates_occ = Nstates - Nstates_empty
    
    if Nspin == 1
        for ist = 1:Nstates_occ-1
            Focc[ist,:] .= 2.0
        end
        if is_odd
            Focc[Nstates_occ,:] .= 1.0
        else
            Focc[Nstates_occ,:] .= 2.0
        end
    else
        for ist = 1:Nstates_occ-1
            Focc[ist,:] .= 1.0
        end
        idx_up = 1:Nkpt
        if is_odd
            # assign only to the spin up
            Focc[Nstates_occ,idx_up] .= 1.0
        else
            Focc[Nstates_occ,:] .= 1.0
        end
    end

    sFocc = sum(Focc)/Nkpt
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        @printf("ERROR diff sum(Focc) and Nelectrons is not small\n")
        @printf("sum Focc = %f, Nelectrons = %f\n", sFocc, Nelectrons)
        exit()
    end

    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )


end


function get_Nelectrons( atoms::Atoms, Pspots::Array{PsPot_GTH,1} )
    Nelectrons = 0.0
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    for ia = 1:Natoms
        isp = atm2species[ia]
        Nelectrons = Nelectrons + Pspots[isp].zval
    end
    return Nelectrons
end


function get_Zvals( PsPots::Array{PsPot_GTH,1} )
    Nspecies = size(PsPots)[1]
    Zvals = zeros(Float64, Nspecies)
    for isp = 1:Nspecies
        Zvals[isp] = PsPots[isp].zval
    end
    return Zvals
end


import Base.println
function println( electrons::Electrons, all_states=false )

    Nspin = electrons.Nspin
    Focc = electrons.Focc
    Nelectrons = electrons.Nelectrons
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nkspin = size(Focc)[2]

    @printf("\n")
    @printf("                                     ---------\n")
    @printf("                                     Electrons\n")
    @printf("                                     ---------\n")
    @printf("\n")
    @printf("Nspin         = %8d\n", Nspin)
    @printf("Nkpt          = %8d\n", Nkspin/Nspin)
    @printf("Nelectrons    =  %18.10f\n", Nelectrons)
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
    if all_states
        for ist = 1:Nstates
            @printf("state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf("%8.5f ", Focc[ist,iks])
            end
            @printf("\n")
        end
    else
        for ist = 1:4
            @printf("state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf("%8.5f ", Focc[ist,iks])
            end
            @printf("\n")
        end
        @printf(".....\n")
        #
        for ist = Nstates-3:Nstates
            @printf("state #%4d = ", ist)
            for iks = 1:Nkspin
                @printf("%8.5f ", Focc[ist,iks])
            end
            @printf("\n")
        end
    end
end

