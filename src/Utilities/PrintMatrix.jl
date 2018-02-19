function PrintMatrix( A )
  Nrows = size(A)[1]
  Ncols = size(A)[2]
  for ir = 1:Nrows
    for ic = 1:Ncols
      @printf("%15.8f ", A[ir,ic])
    end
    @printf("\n")
  end
end
