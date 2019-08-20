# modify evars.psiks and evars.Haux
# modify Ham.electrons.ebands
function rotate_evars!( Ham, evars; skip_ortho=false )
    psiks = evars.psiks
    Haux = evars.Haux
    if !skip_ortho
        for i in length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end
    U_Haux = copy(Haux)
    for i in 1:length(U_Haux)
        Ham.electrons.ebands[:,i], U_Haux[i] = eigen( Haux[i] )
        Haux[i] = diagm( 0 => Ham.electrons.ebands[:,i] ) # rotate Haux
        psiks[i] = psiks[i]*U_Haux[i] # rotate psiks
    end
    return
end

# modify evars.psiks and evars.Haux
# modify Ham.electrons.ebands
function rotate_evars!( Ham, evars, d; skip_ortho=false )
    psiks = evars.psiks
    Haux = evars.Haux
    if !skip_ortho
        for i in length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end
    U_Haux = copy(Haux)
    d_U_Haux = copy(Haux)
    for i in 1:length(U_Haux)
        Ham.electrons.ebands[:,i], U_Haux[i] = eigen( Haux[i] )
        Haux[i] = diagm( 0 => Ham.electrons.ebands[:,i] ) # rotate Haux
        psiks[i] = psiks[i]*U_Haux[i] # rotate psiks
        # also rotate search direction
        λ, d_U_Haux[i] = eigen( d.Haux[i] )
        d.Haux[i] = diagm( 0 => λ )
        #ortho_sqrt!(d.psiks[i])  # need this?
        d.psiks[i] = d.psiks[i]*d_U_Haux[i]
    end
    return
end


function subspace_rotation!( Ham, psiks; rotate_psiks=true )

    Nkspin = length(psiks)
    Nstates = Ham.electrons.Nstates
    Hsub = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        Hsub[i] = zeros(ComplexF64,Nstates,Nstates)
    end

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    for ispin = 1:Nspin, ik = 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin - 1)*Nkpt
        Hr = Hermitian(psiks[i]' * op_H(Ham, psiks[i]))
        Ham.electrons.ebands[:,i], evecs = eigen(Hr)
        if rotate_psiks
            psiks[i] = psiks[i]*evecs # also rotate
        end
        Hsub[i] = psiks[i]' * ( op_H(Ham, psiks[i]) )
    end

    return Hsub
end