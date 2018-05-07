# F = E - TS
# S = -k \sum( .... )
function calc_entropy( Focc::Array{Float64,2}, wk::Array{Float64,1},
                       kT::Float64; is_spinpol=false )
    #
    const SMALL = 1.e-10
    Nstates = size(Focc)[1]
    Nkpt = size(Focc)[2]
    e = 0.0
    if is_spinpol
        for ik = 1:Nkpt
            for ist = 1:Nstates
                if Focc[ist] > SMALL
                    e = e + Focc[ist,ik]*log(Focc[ist,ik])*wk[ik]
                else
                    e = e + (1.0 - Focc[ist,ik])*log(1.0 - Focc[ist,ik])*wk[ik]
                end
            end
        end
    else
        for ik = 1:Nkpt
            for ist = 1:Nstates
                # spin-degenerate case
                if Focc[ist] > SMALL
                    e = e + 0.5*Focc[ist,ik]*log(0.5*Focc[ist,ik])*wk[ik]
                else
                    e = e + (1.0 - 0.5*Focc[ist,ik])*log(1.0 - 0.5*Focc[ist,ik])*wk[ik]
                end
            end
            e = 2*e
        end
    end
    # double negative
    return kT*e
end


# F = E - TS
# S = -k \sum( .... )
function calc_entropy( Focc::Array{Float64,1}, kT::Float64; is_spinpol=false )
    const SMALL = 1.e-10
    Nstates = length(Focc)
    e = 0.0
    if is_spinpol
        for ist = 1:Nstates
            if Focc[ist] > SMALL
                e = e + Focc[ist]*log(Focc[ist])
            else
                e = e + (1.0 - Focc[ist])*log(1.0 - Focc[ist])
            end
        end
    else
        for ist = 1:Nstates
            # spin-degenerate case
            if Focc[ist] > SMALL
                e = e + 0.5*Focc[ist]*log(0.5*Focc[ist])
            else
                e = e + (1.0 - 0.5*Focc[ist])*log(1.0 - 0.5*Focc[ist])
            end
        end
        e = 2*e
    end
    # double negative
    return kT*e
end

