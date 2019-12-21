using SpecialFunctions: erfc

#=
function calc_E_NN_mod( atoms::Atoms )
    return calc_E_NN_mod( atoms.LatVecs, atoms, atoms.Zvals )
end

function calc_E_NN_mod( atoms::Atoms, Zvals::Array{Float64,1} )
    return calc_E_NN_mod( atoms.LatVecs, atoms, Zvals )
end

function calc_E_NN_mod( pw::PWGrid, atoms::Atoms, Zvals::Array{Float64,1} )
    return calc_E_NN_mod( pw.LatVecs, atoms, Zvals )
end
=#

function calc_E_NN_mod( pw::PWGrid, atoms::Atoms )
    return calc_E_NN_mod( pw, pw.LatVecs, atoms, atoms.Zvals )
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
"""
function calc_E_NN_mod( pw::PWGrid,
    LatVecs::Array{Float64,2}, atoms::Atoms, Zvals::Array{Float64,1}
)

    t1 = LatVecs[:,1]
    t2 = LatVecs[:,2]
    t3 = LatVecs[:,3]
  
    t1m = sqrt(dot(t1,t1))
    t2m = sqrt(dot(t2,t2))
    t3m = sqrt(dot(t3,t3))

    CellVolume = pw.CellVolume
    tau = atoms.positions

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # determine eta
    TOL  = 1e-8
    Gcut = 2*pw.ecutwfc/(2*pi)
    gexp = -log(TOL)
    eta  = 0.5*Gcut^2/gexp

    tpi = 2.0*pi
    con = CellVolume/(4.0*pi)
    con2 = (4.0*pi)/CellVolume

    x = 0.0
    totalcharge = 0.0
    for ia = 1:Natoms
        isp = atm2species[ia]
        x = x + Zvals[isp]^2
        totalcharge = totalcharge + Zvals[isp]
    end

    #ewald = -cccc*x - 4.0*pi*(totalcharge^2)/(CellVolume*eta)
    # x = \sum_{I=1}^{Natoms} Z_{I}^{2}
    ewald = -eta/sqrt(pi)*x - totalcharge^2*pi/(CellVolume*2*eta^2)

    tmax = sqrt(2.0*gexp/eta^2)

    mmm1 = round(Int64, tmax/t1m + 1.5)
    mmm2 = round(Int64, tmax/t2m + 1.5)
    mmm3 = round(Int64, tmax/t3m + 1.5)

    v = zeros(Float64,3)
    w = zeros(Float64,3)

    for ia = 1:Natoms
    for ja = 1:Natoms
        #
        v[:] = tau[:,ia] - tau[:,ja]
        #
        isp = atm2species[ia]
        jsp = atm2species[ja]
        prd = Zvals[isp]*Zvals[jsp]*0.5
        for i = -mmm1:mmm1
        for j = -mmm2:mmm2
        for k = -mmm3:mmm3
            if (ia != ja) || ( (abs(i) + abs(j) + abs(k)) != 0 )
                w[:] = v[:] + i*t1 + j*t2 + k*t3
                rmag = sqrt(dot(w,w))
                arg = rmag*eta
                ewald = ewald + prd*erfc(arg)/rmag
            end
        end
        end
        end
    end
    end

    return ewald
end
