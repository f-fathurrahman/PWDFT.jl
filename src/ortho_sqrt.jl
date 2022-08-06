function ortho_sqrt( psi::Array{ComplexF64,2} )
    Udagger = inv(sqrt(psi'*psi))
    return psi*Udagger
end

function ortho_sqrt!( psi::Array{ComplexF64,2} )
    Udagger = inv(sqrt(psi'*psi))
    psi[:,:] = psi*Udagger
    return
end