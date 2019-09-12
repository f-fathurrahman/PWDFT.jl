# Overload println for Atoms
import Base: show
function show( io::IO, a::Atoms; header=true )
    if header
        @printf(io, "\n")
        @printf(io, "                                     -----\n")
        @printf(io, "                                     Atoms\n")
        @printf(io, "                                     -----\n")
        @printf(io, "\n")
    end
    @printf(io, "Natoms   = %5d\n", a.Natoms)
    @printf(io, "Nspecies = %5d\n", a.Nspecies)
    @printf(io, "\n")
    for isp = 1:a.Nspecies
        @printf(io, "Species %3d: %s\n", isp, a.SpeciesSymbols[isp])
    end
    
    @printf(io, "\n")
    @printf(io, "Cartesian coordinates in bohr:\n")
    @printf(io, "\n")
    for ia = 1:a.Natoms
        @printf(io, "%5s %18.10f %18.10f %18.10f\n", a.atsymbs[ia],
                a.positions[1,ia], a.positions[2,ia], a.positions[3,ia])
    end
    
    @printf(io, "\n")
    @printf(io, "Direct lattice vectors: (v1, v2, and v3 are given by column)\n")
    @printf(io, "\n")
    for i = 1:3
        @printf(io, "%18.10f %18.10f %18.10f\n", a.LatVecs[i,1], a.LatVecs[i,2], a.LatVecs[i,3])
    end

    frac_pos = inv(a.LatVecs)*a.positions
    @printf(io, "\n")
    @printf(io, "Fractional coordinates:\n")
    @printf(io, "\n")
    for ia = 1:a.Natoms
        @printf(io, "%5s %18.10f %18.10f %18.10f\n", a.atsymbs[ia],
                frac_pos[1,ia], frac_pos[2,ia], frac_pos[3,ia])
    end

end
show( a::Atoms; header=true ) = show( stdout, a, header=header )