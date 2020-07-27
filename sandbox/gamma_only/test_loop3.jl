using Printf

function main()
    #Ns = (3,4,5)
    Ns = (5,8,12)
    #Ns = (5,4,3)
    ip = 0
    for k in 0:Ns[3]-1
        for j in 0:Ns[2]-1
            for i in 0:Ns[1]-1
                ip = ip + 1
                ii = i + 1
                jj = j + 1
                kk = k + 1 
                ip1 = i + 1 + Ns[1]*( j + Ns[2]*k )
                @printf("[%4d,%4d,%4d] -- %4d -- %4d\n", i, j, k, ip, ip1)
            end
            println()
        end
    end

    println("Ns = ", Ns)
    println("Npoints = ", prod(Ns))
end

main()