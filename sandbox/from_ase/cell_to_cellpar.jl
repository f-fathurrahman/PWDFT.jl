function cell_to_cellpar( cell; radians=false )

    @assert length(cell) == 9

    lengths = zeros(3)
    for i in 1:3
        lengths[i] = norm(cell[:,i])
    end

    # α: between 3rd and 2nd vector
    ll = lengths[2] * lengths[3]
    if ll > 1e-16
        x = dot(cell[:,2], cell[:,3]) / ll
        α = 180.0 / pi * acos(x)
    else
        α = 90.0
    end

    # β
    ll = lengths[1] * lengths[3]
    if ll > 1e-16
        x = dot(cell[:,1], cell[:,3]) / ll
        β = 180.0 / pi * acos(x)
    else
        β = 90.0
    end

    # γ
    ll = lengths[1] * lengths[2]
    if ll > 1e-16
        x = dot(cell[:,1], cell[:,2]) / ll
        γ = 180.0 / pi * acos(x)
    else
        γ = 90.0
    end


#=
    angles = zeros(3)
    for i in 0:2
        j = i - 1 + 1
        k = i - 2 + 1

        if j < 1
            j = j + 3
        end

        if k < 1
            k = k + 3
        end

        println()
        println("i = ", i+1)
        println("j = ", j)
        println("k = ", k)

        ll = lengths[j] * lengths[k]
        if ll > 1e-16
            x = dot(cell[:,j], cell[:,k]) / ll
            angle1 = 180.0 / pi * acos(x)
        else
            angle1 = 90.0
        end
        angles[i+1] = angle1
    end
=#

    println(lengths)
    println(α, " ", β, " ", γ)

end

#=
def cell_to_cellpar(cell, radians=False):
    """Returns the cell parameters [a, b, c, alpha, beta, gamma].

    Angles are in degrees unless radian=True is used.
    """
    lengths = [np.linalg.norm(v) for v in cell]
    angles = []
    for i in range(3):
        j = i - 1
        k = i - 2
        ll = lengths[j] * lengths[k]
        if ll > 1e-16:
            x = np.dot(cell[j], cell[k]) / ll
            angle = 180.0 / pi * arccos(x)
        else:
            angle = 90.0
        angles.append(angle)
    if radians:
        angles = [angle * pi / 180 for angle in angles]
    return np.array(lengths + angles)
=#
