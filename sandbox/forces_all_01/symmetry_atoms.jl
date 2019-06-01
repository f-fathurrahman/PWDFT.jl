function init_irt( atoms::Atoms, sym_info::SymmetryInfo )
    
    Natoms = atoms.Natoms
    tau = atoms.positions
    RecVecs = inv(Matrix(atoms.LatVecs'))  # !! We don't include 2*pi factor here
    atm2species = atoms.atm2species
    
    xau = zeros(3,Natoms)
    rau = zeros(3,Natoms)
    
    for ia = 1:Natoms
        xau[:,ia] = RecVecs[1,:]*tau[1,ia] + RecVecs[2,:]*tau[2,ia] + RecVecs[3,:]*tau[3,ia]
        println(xau[:,ia])
    end

    # checking for supercell is skipped (probably not needed here)

    Nsyms = sym_info.Nsyms
    s = sym_info.s
    ft = sym_info.ft

    irt = zeros(Int64,Nsyms,Natoms)

    Nequal = 0
    for isym = 1:Nsyms
        for ia = 1:Natoms
            rau[:,ia] = s[1,:,isym] * xau[1,ia] + 
                        s[2,:,isym] * xau[2,ia] + 
                        s[3,:,isym] * xau[3,ia]
        end
    end

    for isym = 1:Nsyms
        for ia = 1:Natoms
            for ib = 1:Natoms
                if atm2species[ia] == atm2species[ib]
                    is_equal = eqvect( rau[:,ia], xau[:,ib], ft , 1e-5 )
                    #println("is_equal = ", is_equal)
                    if is_equal
                        Nequal = Nequal + 1
                        #println("is_equal = ", is_equal)
                        irt[isym,ia] = ib
                        break
                    end
                end
            end
        end
    end



    #for isym = 1:Nsyms
    #    for ia = 1:Natoms
    #        @printf("%d ", irt[isym,ia])
    #    end
    #    @printf("\n")
    #end
    #println("Nequal = ", Nequal)

    return irt

end

function eqvect(x, y, f, accep)
  
  # This function test if the difference x-y-f is an integer.
  # x, y = 3d vectors in crystal axis, f = fractional translation
  # adapted from QE-6.4

  res = ( abs(x[1]-y[1]-f[1] - round(Int64, x[1]-y[1]-f[1])) < accep ) &&
        ( abs(x[2]-y[2]-f[2] - round(Int64, x[2]-y[2]-f[2])) < accep ) &&
        ( abs(x[3]-y[3]-f[3] - round(Int64, x[3]-y[3]-f[3])) < accep )
  
  return res
end
