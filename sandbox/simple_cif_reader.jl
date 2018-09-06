using Printf
using PWDFT

function simple_cif_reader( filename::String; verbose=false )
    
    f = open(filename, "r")
    lines = readlines(f)
    close(f)
    
    phase_name = "X"
    a = 0.0
    b = 0.0
    c = 0.0
    α = 0.0
    β = 0.0
    γ = 0.0
    for l in lines
        spl_str = split(l)
        if length(spl_str) >=2
            if spl_str[1] == "_pd_phase_name"
                phase_name = replace( replace(spl_str[2], "'" => ""), " " => "" )
            elseif spl_str[1] == "_cell_length_a"
                a = parse(Float64, spl_str[2])
            elseif spl_str[1] == "_cell_length_b"
                b = parse(Float64, spl_str[2])
            elseif spl_str[1] == "_cell_length_c"
                c = parse(Float64, spl_str[2])
            elseif spl_str[1] == "_cell_angle_alpha"
                α = parse(Float64, spl_str[2])
            elseif spl_str[1] == "_cell_angle_beta"
                β = parse(Float64, spl_str[2])
            elseif spl_str[1] == "_cell_angle_gamma"
                γ = parse(Float64, spl_str[2])
            end
        end
    end

    # assumption: this is the exact number of columns for a line that
    # contains atomic site information   
    Natoms = 0
    for l in lines
        spl_str = split(l)
        if length(spl_str) == 8
            Natoms = Natoms + 1
        end
    end

    atsymbs = Array{String}(undef,Natoms)
    atpos = zeros(Float64,3,Natoms)

    ia = 0
    for l in lines
        spl_str = split(l)
        if length(spl_str) == 8
            ia = ia + 1
            atsymbs[ia] = spl_str[8]
            atpos[1,ia] = parse(Float64, spl_str[3])
            atpos[2,ia] = parse(Float64, spl_str[4])
            atpos[3,ia] = parse(Float64, spl_str[5])
        end
    end
    
    LatVecs = ANG2BOHR*gen_lattice_triclinic(a, b, c, α, β, γ)
    atpos = LatVecs*atpos

    if verbose
        println("phase_name = ", phase_name)
        println("a = ", a)
        println("b = ", b)
        println("c = ", c)
        println("α = ", α)
        println("β = ", β)
        println("γ = ", γ)
        for ia = 1:Natoms
            @printf("%8s %18.10f %18.10f %18.10f\n", atsymbs[ia], atpos[1,ia], atpos[2,ia], atpos[3,ia])
        end
    end

    #
    # !!! THIS IS A SPECIAL CASE !!!
    #
    Nspecies = 1
    atm2species = ones(Int64,Natoms)
    SpeciesSymbols = [phase_name]
    Zvals = [0.0]

    return Atoms( Natoms, Nspecies, atpos, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )

end

const CIFS_HOME = "../_compare/DeltaCodes/DeltaCodesDFT-master/CIFs/"

function test_main()
    atoms = simple_cif_reader(CIFS_HOME*"H.cif", verbose=true)
    write_xsf( "TEMP_H.xsf", atoms )
end

test_main()
