mutable struct BlochWavefuncGamma
    data::Array{ComplexF64,2}
end

import LinearAlgebra: dot
function dot( v1::BlochWavefuncGamma, v2::BlochWavefuncGamma )
    return 2*dot(v1.data, v2.data)
end




#
# Old definition, the DC component is not zero, but a pure real number
#
#function dot(v1::BlochWavefuncGamma, v2::BlochWavefuncGamma)
#    res = 2*dot(v1.data, v2.data)
#    for ist in 1:size(v1.data,2)
#        res = res - v1.data[1,ist] * v2.data[1,ist]
#    end
#    return res
#end