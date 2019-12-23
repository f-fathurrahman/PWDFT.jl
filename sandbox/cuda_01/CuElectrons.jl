mutable struct CuElectrons
    Nelectrons::Float64
    Nstates::Int64
    Nstates_occ::Int64
    Focc_gpu::CuArray{Float64,2}
    Focc::Array{Float64,2}    
    ebands_gpu::CuArray{Float64,2}
    ebands::Array{Float64,2}    
    Nspin::Int64
end

function CuElectrons()
    Nelectrons = 1
    Nstates = 1
    Nstates_occ = 1
    
    Focc_gpu = CuArrays.zeros(Float64,Nstates,1) # Nkpt=1
    Focc = zeros(Nstates,1)

    ebands_gpu = CuArrays.zeros(Float64,Nstates,1) # use Nkpt=1
    ebands = zeros(Nstates,1)

    Nspin = 1
    return CuElectrons( Nelectrons, Nstates, Nstates_occ, Focc_gpu, Focc, ebands_gpu, ebands, Nspin )
end

function CuElectrons( ele::Electrons )
    return CuElectrons(
        ele.Nelectrons,
        ele.Nstates,
        ele.Nstates_occ,
        CuArray(ele.Focc),
        ele.Focc,
        CuArray(ele.ebands),
        ele.ebands,
        ele.Nspin )
end
