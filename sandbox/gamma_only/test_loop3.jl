using Printf

function main()
    Ns = (3,4,5)
    ip = 0
    for k in 0:Ns[3], j in 0:Ns[2], i in 0:Ns[1]
        ip = ip + 1
        ip1 = i + 1 + j*Ns[2] + k*Ns[2]*Ns[3]
        @printf("%4d %4d\n", ip, ip1)
    end
end

main()