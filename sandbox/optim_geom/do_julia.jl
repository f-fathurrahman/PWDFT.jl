using LinearAlgebra

function main()
    Natoms = 2
    H = diagm( 0 => ones(Natoms*3)*70 )
    H[2,2] = 22;
    display(H); println()

    forces = zeros(3,Natoms)
    forces[:,1] = [3.0, 2.0, 3.0]
    forces[:,2] =  [3.1, 2.1, 4.0] 
    display(forces'); println()

    f = vec(forces)
    omega, V = eigen(H)

    display(omega); println()

    dr = V' * (V*f ./ abs.(omega))
    #dr = np.dot( V, np.dot(f, V) / np.fabs(omega) ).reshape((-1, 3))
    display(dr); println()

    dr = reshape(dr, (3,Natoms))
    steplengths = sqrt.(sum( dr.^2, dims=1 ))
    display(steplengths); println();
end

main()