mutable struct ElectronsInfo
    Nelectrons::Float64
    Nstates::Int64
    Nstates_occ::Int64
    Focc::Array{Float64,1}
end

# dummy ElectronsInfo
function ElectronsInfo()
    Nelectrons = 1
    Nstates = 1
    Nstates_occ = 1
    Focc = zeros(Nstates)
    return ElectronsInfo( Nelectrons, Nstates, Nstates_occ, Focc )
end


function ElectronsInfo( atoms::Atoms, Pspots::Array{PsPot_GTH,1};
                        Nstates=nothing, Nstates_empty=0 )

    Nelectrons = get_Nelectrons(atoms,Pspots)

    is_odd = round(Int64,Nelectrons)%2 == 1

    if Nstates == nothing
        Nstates = round( Int64, Nelectrons/2 )
        if is_odd
            Nstates = Nstates + 1
        end
    end

    Focc = zeros( Float64, Nstates )
    Nstates_occ = Nstates - Nstates_empty
    for ist = 1:Nstates_occ-1
        Focc[ist] = 2.0
    end

    if is_odd
        Focc[Nstates_occ] = 1.0
    else
        Focc[Nstates_occ] = 2.0
    end

    # Check if the generated Focc is consistent
    if abs( sum(Focc) - Nelectrons ) > eps()
        @printf("ERROR diff sum(Focc) and Nelectrons is not small\n")
        @printf("sum Focc = %f, Nelectrons = %f\n", sum(Focc), Nelectrons)
        exit()
    end

    return ElectronsInfo( Nelectrons, Nstates, Nstates_occ, Focc )
end


function ElectronsInfo( atoms::Atoms, Zvals::Array{Float64,1};
                        Nstates=nothing, Nstates_empty=0 )

    Nelectrons = 0.0
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    for ia = 1:Natoms
        isp = atm2species[ia]
        Nelectrons = Nelectrons + Zvals[isp]
    end

    is_odd = round(Int64,Nelectrons)%2 == 1

    if Nstates == nothing
        Nstates = round( Int64, Nelectrons/2 )
        if is_odd
            Nstates = Nstates + 1
        end
    end

    Focc = zeros( Float64, Nstates )
    Nstates_occ = Nstates - Nstates_empty
    for ist = 1:Nstates_occ-1
        Focc[ist] = 2.0
    end

    if is_odd
        Focc[Nstates_occ] = 1.0
    else
        Focc[Nstates_occ] = 2.0
    end

    # Check if the generated Focc is consistent
    if abs( sum(Focc) - Nelectrons ) > eps()
        @printf("ERROR diff sum(Focc) and Nelectrons is not small\n")
        @printf("sum Focc = %f, Nelectrons = %f\n", sum(Focc), Nelectrons)
        exit()
    end

    return ElectronsInfo( Nelectrons, Nstates, Nstates_occ, Focc )


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
function println( electrons::ElectronsInfo )
    @printf("\n")
    @printf("                       -------------\n")
    @printf("                       ElectronsInfo\n")
    @printf("                       -------------\n")
    @printf("\n")
    @printf("Nelectrons    =  %18.10f\n", electrons.Nelectrons)
    @printf("Nstates_occ   = %8d\n", electrons.Nstates_occ)
    @printf("Nstates_empty = %8d\n\n", electrons.Nstates - electrons.Nstates_occ)
    @printf("Occupation numbers:\n\n")
    for ist = 1:electrons.Nstates
        @printf("state #%8d = %8.5f\n", ist, electrons.Focc[ist])
    end
end

