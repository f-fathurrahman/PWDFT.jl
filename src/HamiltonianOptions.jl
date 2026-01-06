# Originally the fields are keyword arguments to `Hamiltonian`
mutable struct HamiltonianOptions
    dual::Float64
    Nspin_channel::Int64 # for wavefunction
    Nspin_comp::Int64 # for density, potentials
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
    use_smearing::Bool
    starting_magn::Union{Vector{Float64},Nothing}
    angle1::Union{Vector{Float64},Nothing}
    angle2::Union{Vector{Float64},Nothing}
    lspinorb::Bool # spin-orbit coupling
    noncollinear::Bool # noncollinear magn
end

function HamiltonianOptions()
    dual = 4.0
    Nspin_channel = 1
    Nspin_comp = 1
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
    use_smearing = false
    starting_magn = nothing
    angle1 = nothing
    angle2 = nothing
    lspinorb = false
    noncollinear = false
    return HamiltonianOptions(
        dual, Nspin_channel, Nspin_comp, meshk, shiftk, time_reversal, Ns,
        kpoints, kpts_str, xcfunc, use_xc_internal,
        extra_states, Nstates, use_symmetry, use_smearing,
        starting_magn, angle1, angle2,
        lspinorb, noncollinear
    )
end

function get_domag(options::HamiltonianOptions)
    domag = false # this is only have meaning in case of lspinorb (with noncollinear)
    if options.noncollinear
        if !isnothing(options.starting_magn) && !isapprox(sum(options.starting_magn), 0.0)
            domag = true
        end
    end
    return domag
end