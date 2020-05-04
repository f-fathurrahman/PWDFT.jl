mutable struct BlochWavefuncGammaOnly
    data::Array{ComplexF64,2}
end

import LinearAlgebra: dot
function dot(v1::BlochWavefuncGammaOnly, v2::BlochWavefuncGammaOnly)
    res = 2*dot(v1.data, v2.data) - v1[1]*v2[1]
    return BlochWavefuncGammaOnly(res)
end