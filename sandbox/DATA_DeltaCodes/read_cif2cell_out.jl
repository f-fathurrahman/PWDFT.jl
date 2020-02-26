using PWDFT

function main()
    f = open("CIF2CELL/Cl.OUT", "r")
    
    a = 0.0
    b = 0.0
    c = 0.0
    α = 0.0
    β = 0.0
    γ = 0.0

    line = readline(f)
    line = replace(line, "," => "")
    line = replace(line, "[" => "")
    line = replace(line, "]" => "")

    Natoms = length(split(line)) 

    atsymbs = ["" for i in 1:Natoms]
    atpos = zeros(Float64, 3, Natoms)

    while !eof(f)
        line = readline(f)

        # This should be encountered first
        if occursin("Lattice parameters", line)
            #
            line = readline(f) # skip one line
            line = readline(f)
            ll = split(line, " ", keepempty=false)
            a = parse(Float64, ll[1])
            b = parse(Float64, ll[2])
            c = parse(Float64, ll[3])
            #
            line = readline(f)
            line = readline(f)
            ll = split(line, " ", keepempty=false)
            α = parse(Float64, ll[1])
            β = parse(Float64, ll[2])
            γ = parse(Float64, ll[3])
        end

        if occursin("Representative sites", line)
            line = readline(f) # skip one line
            for ia in 1:Natoms
                line = readline(f)
                ll = split(line, " ", keepempty=false)
                atsymbs[ia] = ll[1]
                for i in 1:3
                    atpos[i,ia] = parse(Float64, ll[i+1])
                end
            end
        end
    end

    close(f)

    println("a = ", a)
    println("b = ", b)
    println("c = ", c)

    println("α = ", α)
    println("β = ", β)
    println("γ = ", γ)

    LatVecs = gen_lattice_triclinic(a, b, c, α, β, γ)
    
    # purge very small components
    if LatVecs[1,2] < 10*eps() LatVecs[1,2] = 0.0 end
    if LatVecs[1,3] < 10*eps() LatVecs[1,3] = 0.0 end
    if LatVecs[2,3] < 10*eps() LatVecs[2,3] = 0.0 end
    
    LatVecs[:,:] = ANG2BOHR*LatVecs[:,:]
    atpos[:,:] = LatVecs[:,:]*atpos[:,:]

    SpeciesSymbols = unique(atsymbs)
    Nspecies = length(SpeciesSymbols)

    # Mapping of atoms to species index
    atm2species = PWDFT.get_atm2species( atsymbs, SpeciesSymbols )

    Zvals = zeros(Nspecies)

    atoms = Atoms(Natoms, Nspecies, atpos, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )
    println(atoms)

    write_xsf( "TEMP_atoms.xsf", atoms )

end

main()