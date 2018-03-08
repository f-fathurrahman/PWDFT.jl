using PWDFT


function mm_to_nn(mm::Int64,S::Int64)
    if mm > S/2
        return mm - S
    else
        return mm
    end
end

function init_grid_G_new( Ns, RecVecs, ecutrho )

    Ng = prod(Ns)

    G  = Array{Float64}(3,Ng)
    G2 = Array{Float64}(Ng)

    ig = 0
    Ngs = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ig = ig + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G[1,ig] = RecVecs[1,1]*gi + RecVecs[2,1]*gj + RecVecs[3,1]*gk
        G[2,ig] = RecVecs[1,2]*gi + RecVecs[2,2]*gj + RecVecs[3,2]*gk
        G[3,ig] = RecVecs[1,3]*gi + RecVecs[2,3]*gj + RecVecs[3,3]*gk
        G2[ig] = G[1,ig]^2 + G[2,ig]^2 + G[3,ig]^2
        if 0.5*G2[ig] < ecutrho
            Ngs = Ngs + 1
        end
    end
    end
    end

    return Ngs
    #idx_gw2r = findn( 0.5*G2 .< ecutwfc )
end


function test_main()
    LatVecs = 16.0*diagm(ones(3))
    ecutwfc_Ry = 40.0
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )
    println(pw)

    Ns = pw.Ns
    RecVecs = pw.RecVecs
    ecutrho = pw.ecutrho

    Ngs = init_grid_G_new( Ns, RecVecs, ecutrho )
    println("Ngs = ", Ngs)
end

test_main()
