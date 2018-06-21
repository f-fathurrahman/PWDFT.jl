function dump_bandstructure( evals, kpath; filename="TEMP_bands.dat")
    Nkpt = size(kpath)[2]
    Nstates = size(evals)[1]

    Xcoords = Array{Float64}(undef,Nkpt)
    Xcoords[1] = 0.0
    for ik = 2:Nkpt
        dk = kpath[:,ik] - kpath[:,ik-1]
        Xcoords[ik] = Xcoords[ik-1] + norm( dk )
    end

    f = open(filename,"w")
    
    for ik = 1:Nkpt
        @printf(f, "%18.10f ", Xcoords[ik])
        for ist = 1:Nstates
            @printf(f, "%18.10f ", evals[ist,ik])
        end
        @printf(f, "\n")
    end
    close(f)
end
