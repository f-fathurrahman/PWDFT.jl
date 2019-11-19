using BenchmarkTools

include("XC_c_pw.jl")
include("XC_c_pbe.jl")
include("XC_c_pbe_v2.jl")

function main()
    Rhoe = 1.2
    gRhoe2 = 0.1

    @btime begin
        for i in 1:5000
            XC_c_pbe($Rhoe, $gRhoe2)
        end
    end

    @btime begin
        for i in 1:5000
            XC_c_pbe_v2($Rhoe, $gRhoe2)
        end
    end

end

main()