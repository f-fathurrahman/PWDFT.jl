using PWDFT

# To avoid using Glob (ugh ...)
const ATOM_LIST = [
    "Ag", "Al", "Ar", "As", "Au",
    "Ba", "Be", "Bi", "B", "Br",
    "Ca", "Cd", "C", "Cl", "Co", "Cr",  "Cs", "Cu",
    "Fe", "F", "Ga", "Ge", "He", "Hf", "Hg", "H",
    "I", "In", "Ir", "K", "Kr",
    "Li", "Lu", "Mg", "Mn", "Mo",
    "Na", "Nb", "Ne", "Ni", "N", "O", "Os",
    "Pb", "Pd", "P", "Po", "Pt",
    "Rb", "Re", "Rh", "Rn", "Ru",
    "Sb", "Sc", "Se", "Si", "S", "Sn", "Sr",
    "Ta", "Tc", "Te", "Ti", "Tl",
    "V", "W", "Xe", "Y", "Zn", "Zr"
]

# FIXME: Some data don't have primitive counterpart (error when processing with cif2cell)
const ATOM_LIST_PRIMITIVE = [
    "Ag", "Al", "Ar", "Au",
    "Ba", "Be", "Br",
    "Ca", "Cd", "C", "Cl", "Co", "Cr",  "Cs", "Cu",
    "Fe", "F", "Ga", "Ge", "He", "Hf", "Hg", "H",
    "I", "In", "Ir", "K", "Kr",
    "Lu", "Mg", "Mn", "Mo",
    "Nb", "Ne", "Ni", "N", "O", "Os",
    "Pb", "Pd", "P", "Po", "Pt",
    "Rb", "Re", "Rh", "Rn", "Ru",
    "Sc", "Se", "Si", "Sn", "Sr",
    "Ta", "Tc", "Te", "Ti", "Tl",
    "V", "W", "Xe", "Y", "Zn", "Zr"
]



include("read_pwscf_input.jl")

function main()

    for a in ATOM_LIST_PRIMITIVE
        atoms, meshk = read_pwscf_input("PWINPUT_primitive/"*a*".in")
        println(atoms)
    end

    for a in ATOM_LIST
        atoms, meshk = read_pwscf_input("PWINPUT/"*a*".in")
        println(atoms)
    end

end

main()