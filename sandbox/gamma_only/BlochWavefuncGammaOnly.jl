mutable struct BlochWavefuncGammaOnly
    data::Array{ComplexF64,2}
end

import LinearAlgebra: dot

function dot(v1::BlochWavefuncGammaOnly, v2::BlochWavefuncGammaOnly)    
    res = 2*dot(v1.data, v2.data)
    for ist in 1:size(v1.data,2)
        res = res - v1.data[1,ist] * v2.data[1,ist]
    end
    return res
end