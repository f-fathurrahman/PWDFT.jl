struct MinimizeParams
    N_α_adjust_max::Int64
    αt_start::Float64
    αt_min::Float64
    αt_reduceFactor::Float64
    αt_increaseFactor::Float64
    updateTestStepSize::Bool
end

function MinimizeParams()
    N_α_adjust_max = 3
    αt_start = 1.0
    αt_min = 1e-10
    αt_reduceFactor = 0.1 #0.1
    αt_increaseFactor = 3.0
    updateTestStepSize = true
    return MinimizeParams( N_α_adjust_max, αt_start, αt_min, αt_reduceFactor,
        αt_increaseFactor, updateTestStepSize )
end
