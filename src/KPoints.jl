"""
The type for describing Bloch wave vector of electronic states.
"""
mutable struct KPoints
    Nkpt::Int64
    mesh::Tuple{Int64,Int64,Int64}
    k::Array{Float64,2}
    wk::Array{Float64,1}
    RecVecs::Array{Float64,2} # copy of reciprocal vectors
end

# for compatibility purpose, `mesh` is defaulting to (0,0,0)
# This should be useful for band structure calculation and
# manual specification of kpoints
#=function KPoints( atoms::Atoms, k::Array{Float64,2},
                  wk::Array{Float64,1}, RecVecs::Array{Float64,2} )
    Nkpt = size(k)[1]
    return KPoints( Nkpt, (0,0,0), k, wk, RecVecs )
end
=#

function KPoints( Nkpt::Int64, k::Array{Float64,2},
                  wk::Array{Float64,1}, RecVecs::Array{Float64,2} )
    return KPoints( Nkpt, (0,0,0), k, wk, RecVecs )
end


"""
Creates a 'dummy' instance of `KPoints` with only one kpoint.
"""
function KPoints( atoms::Atoms )
    Nkpt = 1
    k = zeros(3,1)
    wk = [1.0]
    RecVecs = 2*pi*inv(atoms.LatVecs')
    return KPoints( Nkpt, (0,0,0), k, wk, RecVecs )
end


function KPoints( atoms::Atoms, mesh::Tuple{Int64,Int64,Int64}, is_shift::Array{Int64,1};
                  time_reversal=1 )
    return KPoints( atoms, (mesh[1], mesh[2], mesh[3]), is_shift, time_reversal=time_reversal )
end


"""
Generate an instance of `KPoints` with uniform grid-points reduced by symmetry
(using routines from `spglib`) given the following inputs:

- `atoms`: an instance of `Atoms`

- `mesh`: an array of three integers specifying sampling points in the k-grid
  for each directions.

- `is_shift`: an array of three integers (0 or 1) specifying whether the grid
  is shifted or not for each directions.
  **WARNING**: only `is_shift=[0,0,0]` is tested.

- `time_reversal`: an integer specifying whether time-reversal symmetry should
  be used for reducing k-points.

- `verbose`: if `true`, additional information will be printed out.
 """
function KPoints( atoms::Atoms, mesh::Array{Int64,1}, is_shift::Array{Int64,1};
                  time_reversal=1 )

    Nkpt, kgrid, mapping =
    spg_get_ir_reciprocal_mesh( atoms, mesh, is_shift, is_time_reversal=time_reversal )
    
    # search for unique mapping
    umap = unique(mapping)

    # Total number of grid points (unreduced)
    NkptTotal = prod(mesh)

    list_ir_k = []
    for ikk = 1:Nkpt
        for ik = 1:NkptTotal
            if umap[ikk] == mapping[ik]
                append!( list_ir_k, [kgrid[:,ik]] )
                break
            end
        end
    end

    kred = zeros(Float64,3,Nkpt)
    for ik = 1:Nkpt
        kred[1,ik] = list_ir_k[ik][1] / mesh[1]
        kred[2,ik] = list_ir_k[ik][2] / mesh[2]
        kred[3,ik] = list_ir_k[ik][3] / mesh[3]
    end
    
    # count for occurence of each unique mapping
    kcount = zeros(Int64,Nkpt)
    for ik = 1:Nkpt
        kcount[ik] = count( i -> ( i == umap[ik] ), mapping )
    end

    # calculate the weights
    wk = kcount[:]/sum(kcount)

    # need to calculate this here because a `PWGrid` instance is
    # not given in the inputs.
    RecVecs = 2*pi*inv(Matrix(atoms.LatVecs'))

    # convert to cartesian unit
    kred = RecVecs*kred

    return KPoints( Nkpt, (mesh[1], mesh[2], mesh[3]), kred, wk, RecVecs )

end

# unit in crystal coord
function parse_kpts_string(kpts_string::String)

    lines = split(kpts_string, "\n", keepempty=false)
    Nkpt = parse(Int64, lines[1])

    kpts = zeros(3,Nkpt)
    wk = zeros(Nkpt)
    for ik = 1:Nkpt
        ll = split(lines[ik+1], " ", keepempty=false)
        k1 = parse(Float64, ll[1])
        k2 = parse(Float64, ll[2])
        k3 = parse(Float64, ll[3])
        w  = parse(Float64, ll[4])
        kpts[:,ik] = [k1, k2, k3]
        wk[ik] = w
    end

    # normalize sum(wk) to 1
    ss = sum(wk)
    wk = wk/ss

    return Nkpt, kpts, wk

end


function kpoints_from_string( atoms::Atoms, kpts_string::String )
    Nkpt, kred, wk = parse_kpts_string(kpts_string)
    # kpts need to be converted to Cartesian form
    RecVecs = 2*pi*inv(Matrix(atoms.LatVecs'))
    kpt = RecVecs*kred
    return KPoints(Nkpt, (0,0,0), kpt, wk, RecVecs)
end


# Temporary work around for generating kpoints for band structure calculations
function kpoints_from_file( atoms::Atoms, filename::String )
    file = open(filename)
    str = readline(file)
    Nkpt = parse( Int, str )
    kred = zeros( Float64, 3,Nkpt )
    for ik = 1:Nkpt
        str = split(readline(file))
        kred[1,ik] = parse( Float64, str[1] )
        kred[2,ik] = parse( Float64, str[2] )
        kred[3,ik] = parse( Float64, str[3] )
    end
    close(file)
    # kpts need to be converted to Cartesian form
    RecVecs = 2*pi*inv(Matrix(atoms.LatVecs'))
    kpt = RecVecs*kred
    #
    wk = ones(Nkpt) # not used for non-scf calculations
    #
    return KPoints(Nkpt, (0,0,0), kpt, wk, RecVecs)
end

function kpath_from_file( atoms::Atoms, filename::String )
    file = open(filename)
    str = readline(file)
    Nkpt = parse( Int64, str )
    kred = zeros( Float64, 3,Nkpt )
    for ik = 1:Nkpt
        str = split(readline(file))
        kred[1,ik] = parse( Float64, str[1] )
        kred[2,ik] = parse( Float64, str[2] )
        kred[3,ik] = parse( Float64, str[3] )
    end

    # Kpath
    Nkpt_spec = parse(Int64, readline(file))
    println("Nkpt_spec = ", Nkpt_spec)
    kpt_spec_red = zeros(Float64, 3,Nkpt_spec)
    kpt_spec_labels = Array{String}(undef,Nkpt_spec)
    for ik = 1:Nkpt_spec
        str = split(readline(file))
        kpt_spec_red[1,ik] = parse(Float64, str[1])
        kpt_spec_red[2,ik] = parse(Float64, str[2])
        kpt_spec_red[3,ik] = parse(Float64, str[3])
        kpt_spec_labels[ik] = str[4]
    end

    close(file)
    # kpts need to be converted to Cartesian form
    RecVecs = 2*pi*inv(atoms.LatVecs')
    kpt = RecVecs*kred
    kpt_spec = RecVecs*kpt_spec_red
    #
    wk = ones(Nkpt) # not used for non-scf calculations
    #
    return KPoints(Nkpt, (0,0,0), kpt, wk, RecVecs), kpt_spec, kpt_spec_labels
end

function gen_MonkhorstPack( mesh::Array{Int64,1} )
    ik = 0
    kpts = zeros(Float64,3,prod(mesh))
    for k = 1:mesh[3]
    for j = 1:mesh[2]
    for i = 1:mesh[1]
        ik = ik + 1
        kpts[1,ik] = (2*i - mesh[1] - 1)/(2*mesh[1])
        kpts[2,ik] = (2*j - mesh[2] - 1)/(2*mesh[2])
        kpts[3,ik] = (2*k - mesh[3] - 1)/(2*mesh[3])
    end
    end
    end
    return kpts
end


# adapted from ASE
function get_special_kpoints(lattice::String)
    if lattice == "cubic"
        return Dict("G" => [0.0, 0.0, 0.0],
                    "M" => [1/2, 1/2, 0.0],
                    "R" => [1/2, 1/2, 1/2],
                    "X" => [0.0, 1/2, 0.0],
                    "G1" => [0.0, 0.0, 0.0])
    elseif lattice == "fcc"
        return Dict("G" => [0.0, 0.0, 0.0],
                    "K" => [3/8, 3/8, 3/4],
                    "L" => [1/2, 1/2, 1/2],
                    "U" => [5/8, 1/4, 5/8],
                    "W" => [1/2, 1/4, 3/4],
                    "X" => [1/2, 0.0, 1/2],
                    "G1" => [1.0, 1.0, 1.0])
    elseif lattice == "bcc"
        return Dict("G" => [0.0,  0.0, 0.0],
                    "H" => [1/2, -1/2, 1/2],
                    "P" => [1/4,  1/4, 1/4],
                    "N" => [0.0,  0.0, 1/2],
                    "G1" => [0.0, 0.0, 0.0])
    elseif lattice == "tetragonal"
        return Dict("G" => [0.0, 0.0, 0.0],
                    "A" => [1/2, 1/2, 1/2],
                    "M" => [1/2, 1/2, 0.0],
                    "R" => [0.0, 1/2, 1/2],
                    "X" => [0.0, 1/2, 0.0],
                    "Z" => [0.0, 0.0, 1/2],
                    "G1" => [0.0, 0.0, 0.0])
    elseif lattice == "orthorhombic"
        return Dict("G" => [0.0, 0.0, 0.0],
                    "R" => [1/2, 1/2, 1/2],
                    "S" => [1/2, 1/2, 0.0],
                    "T" => [0.0, 1/2, 1/2],
                    "U" => [1/2, 0.0, 1/2],
                    "X" => [1/2, 0.0, 0.0],
                    "Y" => [0.0, 1/2, 0.0],
                    "Z" => [0.0, 0.0, 1/2],
                    "G1" => [0.0, 0.0, 0.0])
    elseif lattice == "hexagonal"
        return Dict("G" => [0.0, 0.0, 0.0],
                    "A" => [0.0, 0.0, 1/2],
                    "H" => [1/3, 1/3, 1/2],
                    "K" => [1/3, 1/3, 0.0],
                    "L" => [1/2, 0.0, 1/2],
                    "M" => [1/2, 0.0, 0.0],
                    "G1" => [0.0, 0.0, 0.0])
    else
        error(@sprintf("Unknown lattice = %s\n", lattice))
    end
end


include("KPoints_io.jl")