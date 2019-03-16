const ZATOMS = Dict(
    "H"  => 1,
    "He" => 2,
    "Li" => 3,
    "Be" => 4,
    "B"  => 5,
    "C"  => 6,
    "N"  => 7,
    "O"  => 8,
    "F"  => 9,
    "Ne" => 10,
    "Na" => 11,
    "Mg" => 12,
    "Al" => 13,
    "Si" => 14,
    "P"  => 15,
    "S"  => 16,
    "Cl" => 17,
    "Ar" => 18,
    "K"  => 19,
    "Ca" => 20,
    "Sc" => 21,
    "Ti" => 22,
    "V"  => 23,
    "Cr" => 24,
    "Mn" => 25,
    "Fe" => 26,
    "Co" => 27,
    "Ni" => 28,
    "Cu" => 29,
    "Zn" => 30,
    "Ga" => 31,
    "Ge" => 32,
    "As" => 33,
    "Se" => 34,
    "Br" => 35,
    "Kr" => 36,
    "Rb" => 37,
    "Sr" => 38,
    "Y"  => 39,
    "Zr" => 40,
    "Nb" => 41,
    "Mo" => 42,
    "Tc" => 43,
    "Ru" => 44,
    "Rh" => 45,
    "Pd" => 46,
    "Ag" => 47,
    "Cd" => 48,
    "In" => 49,
    "Sn" => 50,
    "Sb" => 51,
    "Te" => 52,
    "I"  => 53,
    "Xe" => 54,
    "Cs" => 55,
    "Ba" => 56,
    "La" => 57,
    "Ce" => 58,
    "Pr" => 59,
    "Nd" => 60,
    "Pm" => 61,
    "Sm" => 62,
    "Eu" => 63,
    "Gd" => 64,
    "Tb" => 65,
    "Dy" => 66,
    "Ho" => 67,
    "Er" => 68,
    "Tm" => 69,
    "Yb" => 70,
    "Lu" => 71,
    "Hf" => 72,
    "Ta" => 73,
    "W"  => 74,
    "Re" => 75,
    "Os" => 76,
    "Ir" => 77,
    "Pt" => 78,
    "Au" => 79,
    "Hg" => 80,
    "Tl" => 81,
    "Pb" => 82,
    "Bi" => 83,
    "Po" => 84,
    "At" => 85,
    "Rn" => 86,)


"""
Returns an array of atomic numbers (with size `Nspecies`) given an instance
of `Atoms`.
"""
function get_Zatoms( atoms::Atoms )
    Nspecies = atoms.Nspecies
    SpeciesSymbols = atoms.SpeciesSymbols
    Zatoms = zeros(Float64,Nspecies)
    for isp = 1:Nspecies
        Zatoms[isp] = ZATOMS[SpeciesSymbols[isp]]
    end
    return Zatoms
end


##
## Helper functions
##

"""
Given all atomic symbols present in a system, `atsymbs`,
this function returns number of species present by searching for
unique symbols. This function use a naive algorithm.
Similar functionality can be obtained by using the built-in `unique` function.
```
Nspecies = length(unique(str))
```
"""
function get_Nspecies( atsymbs::Array{String,1} )
    # Determine number of species
    Nspecies = 0
    Natoms = size(atsymbs)[1]
    for ia = 1:Natoms
        k2 = 0
        for k1 = 1:ia-1
            if atsymbs[k1] == atsymbs[ia]
                k2 = 1
            end
        end
        # find different
        if k2 == 0
            Nspecies = Nspecies + 1
        end
    end
    return Nspecies
end


"""
Given all atomic symbols present in a system, `atsymbs`,
this function returns unique atomic species symbol.
This function use a naive algorithm.
Similar functionality can be obtained by using the built-in `unique` function.
```
SpeciesSymbols = unique(str)
```
"""
function get_SpeciesSymbols( Nspecies, atsymbs )
    
    SpeciesSymbols = Array{String}(undef,Nspecies)
    Natoms = size(atsymbs)[1]

    idx1 = 0
    for ia = 1:Natoms
        k2 = 0
        for k1 = 1:ia-1
            if atsymbs[k1] == atsymbs[ia]
                k2 = 1
            end
        end
        # Found different species
        if k2==0
            idx1 = idx1 + 1
            SpeciesSymbols[idx1] = atsymbs[ia]
        end
    end

    return SpeciesSymbols

end

"""
Get `atm2species` which described mapping between `atsymbs` and `SpeciesSymbols`.
"""
function get_atm2species( atsymbs, SpeciesSymbols )

    Natoms = size(atsymbs)[1]
    Nspecies = size(SpeciesSymbols)[1]

    atm2species = Array{Int64}(undef,Natoms)

    for ia = 1:Natoms
        for isp = 1:Nspecies
            if atsymbs[ia] == SpeciesSymbols[isp]
                atm2species[ia] = isp
            end
        end
    end

    return atm2species
end