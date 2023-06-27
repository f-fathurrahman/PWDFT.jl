function ortho_sqrt( psi::Array{ComplexF64,2} )
    Udagger = inv(sqrt(psi'*psi))
    return psi*Udagger
end

function ortho_sqrt!( psi::Array{ComplexF64,2} )
    Udagger = inv(sqrt(psi'*psi))
    psi[:,:] = psi*Udagger
    return
end

"""
    ortho_sqrt(Ham, psi)

Orthonormalize psi::Matrix{Complex} such that

psi' * S * psi = I

where S is overlap operator.

If overlap operator is not needed
(e.g. in case of norm-conserving pseudopotentials), you can pass
psi only.

Multiplication of overlap operator and psi is done by

op_S(Ham, psi)

The `ik` and `ispin` must be set properly before calling this function.

```
for ispin in 1:Nspin, ik in 1:Nkpt
    Ham.ik = ik
    Ham.ispin = ispin
    ikspin = ik + (ispin - 1)*Nkpt
    psiks[ispin] = ortho_sqrt(Ham, psiks[ikspin])
end
```
"""
function ortho_sqrt( Ham::Hamiltonian, psi::Matrix{ComplexF64} )
    if Ham.need_overlap
        O = psi' * op_S(Ham, psi)
    else
        O = psi' * psi
    end
    Udagger = inv(sqrt(O))
    return psi*Udagger
end




"""
    ortho_sqrt!(Ham, psi)

In-place version of `ortho_sqrt`.
"""
function ortho_sqrt!( Ham::Hamiltonian, psi::Matrix{ComplexF64} )
    if Ham.need_overlap
        O = psi' * op_S(Ham, psi)
    else
        O = psi' * psi
    end
    Udagger = inv(sqrt(O))
    psi[:,:] = psi*Udagger
    return
end