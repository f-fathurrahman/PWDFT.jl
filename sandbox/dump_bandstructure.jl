function dump_bandstructure( evals, kpath, kpt_spec, kpt_spec_labels;
                             filename="TEMP_bands.dat")
    Nkpt = size(kpath)[2]
    Nstates = size(evals)[1]

    Xcoords = Array{Float64}(undef,Nkpt)
    Xcoords[1] = 0.0
    
    ikpt_spec = 2 # start searching for 2nd kpt_spec

    SMALL = 1e-10

    Nkpt_spec = size(kpt_spec)[2]
    Xticks = Array{Float64}(undef,Nkpt_spec)
    Xticks[1] = 0.0

    for ik = 2:Nkpt
        dk = kpath[:,ik] - kpath[:,ik-1]
        Xcoords[ik] = Xcoords[ik-1] + norm( dk )
        dnorm_spec = norm(kpath[:,ik]-kpt_spec[:,ikpt_spec])
        #@printf("%d %18.10e\n", ik, dnorm_spec)
        if dnorm_spec < SMALL
            #@printf("Found ikpt_spec: %d\n", ikpt_spec)
            Xticks[ikpt_spec] = Xcoords[ik]
            ikpt_spec = ikpt_spec + 1
        end
    end

    f = open(filename,"w")

    # write Xtick and labels as comment
    @printf(f, "#%d\n", Nkpt_spec)
    for ik = 1:Nkpt_spec
        @printf(f, "#%18.10f %s\n", Xticks[ik], kpt_spec_labels[ik])
    end
    
    for ik = 1:Nkpt
        @printf(f, "%18.10f ", Xcoords[ik])
        for ist = 1:Nstates
            @printf(f, "%18.10f ", evals[ist,ik])
        end
        @printf(f, "\n")
    end
    close(f)
end
