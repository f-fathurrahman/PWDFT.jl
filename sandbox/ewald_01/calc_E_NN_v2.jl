using SpecialFunctions: erfc

function calc_E_NN_v2( atoms::Atoms )
    return calc_E_NN_v2( atoms.LatVecs, atoms, atoms.Zvals )
end

function calc_E_NN_v2( atoms::Atoms, Zvals::Array{Float64,1} )
    return calc_E_NN_v2( atoms.LatVecs, atoms, Zvals )
end

function calc_E_NN_v2( pw::PWGrid, atoms::Atoms, Zvals::Array{Float64,1} )
    return calc_E_NN_v2( pw.LatVecs, atoms, Zvals )
end

"""
Calculates repulsive interaction energy of ions in periodic unit cell
given the following input

- `LatVecs`: lattice vectors

- `atoms`: an instance of `Atoms`

- `Zvals`: array of charges of ions

This is the actual method that does the calculation.
The simplified method `calc_E_NN(atoms::Atoms)`, the input
parameter `LatVecs` is taken from `atoms.LatVecs` and
the charges are calculated based on `atoms.Zvals`. Value of
`atoms.Zvals` is usually calculated inside `Hamiltonian` constructor.
So, in the usual case, one can simply call this method and
save the result in an `Energies` instance and the proceed to
calculate electronic energy components.

The code inspired from Prof. Natalie Holzwarth:
http://users.wfu.edu/natalie/s18phy712/computerprograms/ewaldsum.f90

The original code is written in Rydberg unit.

Revised to improve readability and using η to match the definition
used in Eq. (F.5) in Prof. Richard Martins' book.
"""
function calc_E_NN_v2( LatVecs::Array{Float64,2}, atoms::Atoms, Zvals::Array{Float64,1} )

    t1 = LatVecs[:,1]
    t2 = LatVecs[:,2]
    t3 = LatVecs[:,3]
  
    Ω = abs(det(LatVecs))

    RecVecs = 2*pi*inv(LatVecs')
    g1 = RecVecs[:,1]
    g2 = RecVecs[:,2]
    g3 = RecVecs[:,3]

    t1m = sqrt(dot(t1,t1))
    t2m = sqrt(dot(t2,t2))
    t3m = sqrt(dot(t3,t3))

    g1m = sqrt(dot(g1,g1))
    g2m = sqrt(dot(g2,g2))
    g3m = sqrt(dot(g3,g3))

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # Atomic positions
    tau = atoms.positions

    # Parameters
    gcut = 2.0
    ebsl = 1e-8

    glast2 = gcut*gcut
    gexp = -log(ebsl)    
    η = sqrt(glast2/gexp)/2

    x = 0.0
    totalcharge = 0.0
    for ia = 1:Natoms
        isp = atm2species[ia]
        x = x + Zvals[isp]^2
        totalcharge = totalcharge + Zvals[isp]
    end

    ewald = -2*η*x/sqrt(pi) - pi*(totalcharge^2)/(Ω*η^2)

    tmax = sqrt(0.5*gexp)/η

    mmm1 = round(Int64, tmax/t1m + 1.5)
    mmm2 = round(Int64, tmax/t2m + 1.5)
    mmm3 = round(Int64, tmax/t3m + 1.5)

    dtau = zeros(Float64,3)
    G = zeros(Float64,3)
    T = zeros(Float64,3)

    for ia = 1:Natoms
    for ja = 1:Natoms
        dtau[1] = tau[1,ia] - tau[1,ja]
        dtau[2] = tau[2,ia] - tau[2,ja]
        dtau[3] = tau[3,ia] - tau[3,ja]
        isp = atm2species[ia]
        jsp = atm2species[ja]
        ZiZj = Zvals[isp]*Zvals[jsp]
        for i = -mmm1:mmm1
        for j = -mmm2:mmm2
        for k = -mmm3:mmm3
            if (ia != ja) || ( (abs(i) + abs(j) + abs(k)) != 0 )
                T[1] = i*t1[1] + j*t2[1] + k*t3[1]
                T[2] = i*t1[2] + j*t2[2] + k*t3[2]
                T[3] = i*t1[3] + j*t2[3] + k*t3[3]
                rmag2 = sqrt( (dtau[1] - T[1])^2 +
                              (dtau[2] - T[2])^2 +
                              (dtau[3] - T[3])^2 )
                ewald = ewald + ZiZj*erfc(rmag2*η)/rmag2
            end
        end
        end
        end
    end
    end

    mmm1 = round(Int64, gcut/g1m + 1.5)
    mmm2 = round(Int64, gcut/g2m + 1.5)
    mmm3 = round(Int64, gcut/g3m + 1.5)
      
    for i = -mmm1:mmm1
    for j = -mmm2:mmm2
    for k = -mmm3:mmm3
        if ( abs(i) + abs(j) + abs(k) ) != 0
            G[1] = i*g1[1] + j*g2[1] + k*g3[1]
            G[2] = i*g1[2] + j*g2[2] + k*g3[2]
            G[3] = i*g1[3] + j*g2[3] + k*g3[3]        
            G2 = G[1]^2 + G[2]^2 + G[3]^2
            x = 4*pi/Ω * exp(-0.25*G2/η^2)/G2
            for ia = 1:Natoms
            for ja = 1:Natoms
                isp = atm2species[ia]
                jsp = atm2species[ja]
                ZiZj = Zvals[isp]*Zvals[jsp]
                dtau[1] = tau[1,ia] - tau[1,ja]
                dtau[2] = tau[2,ia] - tau[2,ja]
                dtau[3] = tau[3,ia] - tau[3,ja]
                Gtau = G[1]*dtau[1] + G[2]*dtau[2] + G[3]*dtau[3]
                ewald = ewald + x*ZiZj*cos(Gtau)
            end # ja
            end # ia
        end # if
    end
    end
    end

    return ewald*0.5 # Convert to Hartree
end
