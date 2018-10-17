# F = E - TS
# S = -k \sum( .... )
function calc_entropy( Focc::Array{Float64,2}, wk::Array{Float64,1},
                       kT::Float64; Nspin=1 )
    #
    SMALL = eps()
    
    Nstates = size(Focc)[1]
    Nkspin = size(Focc)[2]
    Nkpt = round(Int64,Nkspin/Nspin)

    ent = 0.0
    if Nspin == 2
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            for ist = 1:Nstates
                if Focc[ist,ik] > SMALL
                    ent = ent + Focc[ist,ik]*log(Focc[ist,ik])*wk[ik]
                else
                    ent = ent + (1.0 - Focc[ist,ik])*log(1.0 - Focc[ist,ik])*wk[ik]
                end
            end
        end
        end
    else
        for ik = 1:Nkpt
            for ist = 1:Nstates
                # spin-degenerate case
                if Focc[ist,ik] > SMALL
                    ent = ent + 0.5*Focc[ist,ik]*log(0.5*Focc[ist,ik])*wk[ik]
                else
                    ent = ent + (1.0 - 0.5*Focc[ist,ik])*log(1.0 - 0.5*Focc[ist,ik])*wk[ik]
                end
            end
        end
        ent = 2*ent
    end

    # double negative
    return kT*ent
end


function calc_entropy(
    Focc::Array{Float64,2},
    wk::Array{Float64,1},
    kT::Float64,
    evals::Array{Float64,2},
    E_fermi::Float64;
    Nspin=1
)

    Nstates = size(Focc)[1]
    Nkspin = size(Focc)[2]
    Nkpt = round(Int64,Nkspin/Nspin)

    ent = 0.0
    if Nspin == 2
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            for ist = 1:Nstates
                if evals[ist,ikspin] <= E_fermi
                    ent = ent + Focc[ist,ikspin]*log(Focc[ist,ikspin])*wk[ik]
                else
                    ent = ent + (1.0 - Focc[ist,ikspin])*log(1.0 - Focc[ist,ikspin])*wk[ik]
                end
            end
        end
        end
    else
        for ik = 1:Nkpt
            for ist = 1:Nstates
                # spin-degenerate case
                if evals[ist,ik] <= E_fermi
                    ent = ent + 0.5*Focc[ist,ik]*log(0.5*Focc[ist,ik])*wk[ik]
                else
                    ent = ent + (1.0 - 0.5*Focc[ist,ik])*log(1.0 - 0.5*Focc[ist,ik])*wk[ik]
                end
            end
        end
        ent = 2*ent
    end

    # double negative
    return kT*ent
end


