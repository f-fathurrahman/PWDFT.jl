function gen_kpath( atoms::Atoms, path_str::String, lattice::String; Δk = 0.02 )
    
    dict_spec_kpts = get_special_kpoints(lattice)

    Nkpt_spec = length(path_str)
    kpt_spec = zeros(3,Nkpt_spec)
    kpt_spec_labels = []
    for ik = 1:Nkpt_spec
        label = string(path_str[ik])
        kvec = dict_spec_kpts[label]
        append!( kpt_spec_labels,label )
        kpt_spec[:,ik] = kvec
    end

    # distance between two adjacent special points in path_str
    d = zeros(Nkpt_spec-1)
    for ik = 1:Nkpt_spec-1
        d[ik] = norm( kpt_spec[:,ik+1] - kpt_spec[:,ik] )
    end

    # estimate number of kpoints between two adjacent special points in kpath
    # based on Δk
    Nk = zeros(Int64,Nkpt_spec-1)
    for ik = 1:Nkpt_spec-1
        Nk[ik] = round(Int64, d[ik]/Δk)
    end

    ipk = 0
    Nkpt_on_path = sum(Nk .- 1) + 1 # total no. of kpt on the path
    kpt = zeros(3,Nkpt_on_path)
    for ik = 1:Nkpt_spec-1
        dvec = kpt_spec[:,ik+1] - kpt_spec[:,ik]
        dk = dvec./(Nk[ik]-1)
        for iik = 1:Nk[ik]-1
            kvec = kpt_spec[:,ik] + (iik-1).*dk
            ipk = ipk + 1
            kpt[:,ipk] = kvec
        end
    end
    
    # The last point
    kpt[:,Nkpt_on_path] = kpt_spec[:,Nkpt_spec]

    # convert to Cartesian form
    RecVecs = 2*pi*inv(atoms.LatVecs')
    kpt_cart = RecVecs*kpt
    kpt_spec_cart = RecVecs*kpt_spec
    #
    wk = ones(Nkpt_on_path) # not used for non-scf calculations
    #
    return KPoints(Nkpt_on_path, kpt_cart, wk, RecVecs), kpt_spec_cart, kpt_spec_labels

end
