using BenchmarkTools

include("PWDFT_cuda.jl")

function main_CPU(; ecutwfc=15.0)
    pw = PWGrid(ecutwfc, gen_lattice_sc(20.0))

    Npoints = prod(pw.Ns)
    Rhoe = rand(Float64,Npoints)

    RhoeG = R_to_G( pw, Rhoe )

    @printf("%8d %8d %8d CPU: ", pw.gvec.Ng, pw.Ns[1], Npoints)
    
    #@btime R_to_G( $pw, $Rhoe )
    @btime R_to_G( $pw.Ns, $Rhoe )
end

function main_GPU(; ecutwfc=15.0)
    pw = CuPWGrid(ecutwfc, gen_lattice_sc(20.0))

    Npoints = prod(pw.Ns)
    Rhoe = CuArray(rand(Float64,Npoints))

    RhoeG = R_to_G( pw, Rhoe )
    
    @printf("%8d %8d %8d GPU: ", pw.gvec.Ng, pw.Ns[1], Npoints)

    @btime begin
        #CuArrays.@sync R_to_G( $pw, $Rhoe )
        CuArrays.@sync R_to_G( $pw.Ns, $Rhoe )
    end
end

for ecut in 30.0:5.0:60.0
    println("ecut = ", ecut)
    main_CPU(ecutwfc=ecut)
    main_GPU(ecutwfc=ecut)
end