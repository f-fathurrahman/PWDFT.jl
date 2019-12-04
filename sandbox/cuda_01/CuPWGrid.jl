using CuArrays
using LinearAlgebra

using PWDFT

function CuPWGrid( ecutwfc::Float64, LatVecs::Array{Float64,2}; kpoints=nothing, Ns_=(0,0,0) )

    ecutrho = 4.0*ecutwfc
    #
    RecVecs = 2*pi*inv(Matrix(LatVecs'))
    #RecVecs = 2*pi*invTrans_m3x3(LatVecs)

    CellVolume = abs(det(LatVecs))
    #
    LatVecsLen = Array{Float64}(undef,3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])

    Ns1 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    if any(Ns_ .== 0)
        Ns1 = good_fft_order(Ns1)
        Ns2 = good_fft_order(Ns2)
        Ns3 = good_fft_order(Ns3)
        Ns = (Ns1,Ns2,Ns3)
    else
        Ns = Ns_[:]
    end

    return

end



function main()
    CuPWGrid( 15.0, gen_lattice_fcc(10.0) )

    println("Pass here")
end

main()