__precompile__()

module PWDFT

#
# Essentials: FFTW
#
export c_G_to_R
export c_R_to_G
export R_to_G
export G_to_R
include("FFT/wrappers_fft.jl")


#
# Submodules begin here
#

module Utilities
    export PrintMatrix
    include("Utilities/PrintMatrix.jl")
end

module PW
    export PWGrid
    export GVectors
    export GVectorsW
    include("PW/PWGrid_v03.jl")
end

module PsPot
    export PsPot_HGH
    export HGH_eval_Vloc_G
    export HGH_eval_proj_G
    export info_PsPot_HGH
    include("PsPot/PsPot_HGH.jl")
end

#
# End of submodules
#

end
