using CuArrays
using PWDFT
using PWDFT_cuda

function main()

    pw = PWGrid(15.0, gen_lattice_fcc(10.0))
    
    cu_gvec = CuGVectors( pw.gvec )
    cu_gvecw = CuGVectorsW( pw.gvecw )

    println("Pass here")
end

@time main()
@time main()