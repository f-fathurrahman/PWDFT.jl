function dump_bandstructure( ebands, kpath, kpt_spec, kpt_spec_labels;
                             filename="TEMP_bands.dat")
    Nkpt = size(kpath, 2)
    Nstates = size(ebands, 1)

    Nspin = round( Int64, size(ebands, 2)/Nkpt )
    Nkspin = Nkpt*Nspin

    Xcoords = Array{Float64}(undef,Nkpt)
    Xcoords[1] = 0.0
    
    ikpt_spec = 2 # start searching for 2nd kpt_spec

    SMALL = 1e-10

    Nkpt_spec = size(kpt_spec, 2)
    Xticks = zeros(Float64, Nkpt_spec)
    Xticks[1] = 0.0

    for ik in 2:Nkpt
        dk = kpath[:,ik] - kpath[:,ik-1]
        Xcoords[ik] = Xcoords[ik-1] + norm( dk )
        dnorm_spec = norm(kpath[:,ik]-kpt_spec[:,ikpt_spec])
        #@printf("%d %18.10e\n", ik, dnorm_spec)
        if dnorm_spec < SMALL
            Xticks[ikpt_spec] = Xcoords[ik]
            @printf("Found ikpt_spec: %d, xtick=%f\n", ikpt_spec, Xticks[ikpt_spec])
            ikpt_spec = ikpt_spec + 1
        end
    end
    println("ikpt_spec last = ", ikpt_spec)
    if ikpt_spec < Nkpt_spec
        @warn "Not finding enough ticks for special k-points"
    end

    f = open(filename,"w")

    # write Xtick and labels as comment
    @printf(f, "#%d\n", Nkpt_spec)
    for ik = 1:Nkpt_spec
        @printf(f, "#%18.10f %s\n", Xticks[ik], kpt_spec_labels[ik])
    end
    
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin -1)*Nkpt
        @printf(f, "%18.10f ", Xcoords[ik])
        for ist in 1:Nstates
            @printf(f, "%18.10f ", ebands[ist,ikspin])
        end
        @printf(f, "\n")
    end
    close(f)
    return
end
