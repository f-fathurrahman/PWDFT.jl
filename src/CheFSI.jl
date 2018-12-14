function diag_CheFSI!( Ham::Hamiltonian, psiks::BlochWavefunc, cheby_degree::Int64)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        ub, lb = get_ub_lb_lanczos( Ham, Nstates*2 )
        psiks[ikspin] = chebyfilt( Ham, psiks[ikspin], cheby_degree, lb, ub)
        psiks[ikspin] = ortho_sqrt( psiks[ikspin] )
    end
    end
end


function chebyfilt( Ham::Hamiltonian, X, degree, lb, ub)
    Ngw_ik  = size(X)[1]
    Nstates = size(X)[2]
    #
    ee = (ub - lb)/2
    c = (ub + lb)/2
    sigma = ee/(lb-ub)
    sigma1 = sigma
    #
    Y = zeros(ComplexF64,Ngw_ik,Nstates)
    Y1 = zeros(ComplexF64,Ngw_ik,Nstates)
    #
    Y = op_H(Ham, X) - X*c
    Y = Y*sigma1/ee
    #
    for i = 2:degree
        sigma2 = 1/(2/sigma1 - sigma)
        Y1 = ( op_H(Ham,Y) - Y*c)*2 * sigma2/ee - X*(sigma*sigma2)
        X = Y
        Y = Y1
        sigma = sigma2
    end
    return Y
end


function get_ub_lb_lanczos( Ham::Hamiltonian, nlancz::Int64 )
    #
    pw = Ham.pw
    ik = Ham.ik
    #
    Ngw_ik = pw.gvecw.Ngw[ik]
    V = zeros(ComplexF64,Ngw_ik,nlancz)
    HV = zeros(ComplexF64,Ngw_ik,nlancz)
    T = zeros(Float64,nlancz,nlancz)
    f = zeros(ComplexF64,Ngw_ik)
    s = zeros(ComplexF64,nlancz)
    h = zeros(ComplexF64,nlancz)
    #
    V[:,1] = randn(Ngw_ik) + im*randn(Ngw_ik)
    beta = norm(V[:,1])
    V[:,1] = V[:,1] ./ beta
    #
    HV[:,1] = op_H( Ham, V[:,1] )
    h[1] = real( V[:,1]' * HV[:,1] )
    #
    T[1,1] = h[1]
    # One-step of reorthogonalization
    f[:] = HV[:,1] - V[:,1]*h[1]
    s[1] = V[:,1]' * f[:]
    h[1] = h[1] + s[1]
    f[:] = f[:] - V[:,1]*s[1]
    # MAIN LOOP
    for j = 2:nlancz
        #@printf("iter lanczos = %d\n", j)
        beta = norm(f)
        T[j,j-1] = beta
        V[:,j] = f[:]/beta
        HV[:,j] = op_H( Ham, V[:,j] )
        #
        for jj = 1:j
            h[jj] = V[:,jj]' * HV[:,j]
        end
        f[:] = HV[:,j] - V[:,1:j]*h[1:j]

        # One-step of reorthogonalization
        #s = V' * f
        #h = h + s
        for jj=1:j
            s[jj] = V[:,jj]' * f[:]
            h[jj] = h[jj] + s[jj]
        end
        f[:] = f[:] - V[:,1:j]*s[1:j]
        #
        T[1:j,j] = real(h[1:j])
    end
    #
    evalsT = eigvals(T)
    #lb = evalsT[Nstates+2]
    #ub = evalsT[2*Nstates]
    lb = evalsT[Int64(nlancz/2)]
    ub = norm_matrix_induced(T) + norm(f)
    #
    return lb, ub
end


function norm_matrix_induced(A::Array{Float64,2})
    N = size(A)[1]
    # FIXME no check for matrix form

    # unit-norm vector
    d = 1/sqrt(N)
    v1 = ones(N)*d
    #
    v = A*v1
    return norm(v)
end

