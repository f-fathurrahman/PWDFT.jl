# Overload println for Atoms
import Base: println
function println( a::Atoms; header=true )
    if header
        @printf("\n")
        @printf("                                     -----\n")
        @printf("                                     Atoms\n")
        @printf("                                     -----\n")
        @printf("\n")
    end
    @printf("Natoms   = %5d\n", a.Natoms)
    @printf("Nspecies = %5d\n", a.Nspecies)
    @printf("\n")
    for isp = 1:a.Nspecies
        @printf("Species %3d: %s\n", isp, a.SpeciesSymbols[isp])
    end
    
    @printf("\n")
    @printf("Cartesian coordinates in bohr:\n")
    @printf("\n")
    for ia = 1:a.Natoms
        @printf("%5s %18.10f %18.10f %18.10f\n", a.atsymbs[ia],
                a.positions[1,ia], a.positions[2,ia], a.positions[3,ia])
    end
    
    @printf("\n")
    @printf("Direct lattice vectors: (v1, v2, and v3 are given by column)\n")
    @printf("\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", a.LatVecs[i,1], a.LatVecs[i,2], a.LatVecs[i,3])
    end

    frac_pos = inv(a.LatVecs)*a.positions
    @printf("\n")
    @printf("Fractional coordinates:\n")
    @printf("\n")
    for ia = 1:a.Natoms
        @printf("%5s %18.10f %18.10f %18.10f\n", a.atsymbs[ia],
                frac_pos[1,ia], frac_pos[2,ia], frac_pos[3,ia])
    end

end