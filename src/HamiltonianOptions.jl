# Originally the fields are keyword arguments to `Hamiltonian`
mutable struct HamiltonianOptions
    dual::Float64
    Nspin::Int64
    meshk::Vector{Int64}
    shiftk::Vector{Int64}
    time_reversal::Bool
    Ns::Tuple{Int64,Int64,Int64}
    kpoints::Union{KPoints,Nothing}
    kpts_str::String
    xcfunc::String
    use_xc_internal::Bool
    extra_states::Int64
    Nstates::Int64
    use_symmetry::Bool
end

function HamiltonianOptions()
    dual = 4.0
    Nspin = 1
    meshk = [1,1,1]
    shiftk = [0,0,0]
    time_reversal = true
    Ns = (0,0,0)
    kpoints = nothing
    kpts_str = ""
    xcfunc = "VWN"
    use_xc_internal = false
    extra_states = -1
    Nstates = -1
    use_symmetry = true
    return HamiltonianOptions(
        dual, Nspin, meshk, shiftk, time_reversal, Ns,
        kpoints, kpts_str, xcfunc, use_xc_internal,
        extra_states, Nstates, use_symmetry
    )
end