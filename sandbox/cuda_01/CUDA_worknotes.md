# Important packages

Several important packages:
- `CUDAnative`
- `CuArrays`

`CUDANative` is a low level library which allows us to write CUDA kernel in Julia.

`CuArrays` is high level library. Some operations are automatically
parallelized. `CuArrays` also wraps several CUDA libraries such as
CUBLAS, CUFFT, CUSOLVER, etc.

# Initialization


The function `cu` can be used to copy array from CPU to GPU, but I am having
several difficulties because the default type would be Int32 or Float32 instead
of Int64 or Float64.

The constructor CuArray(A) can be used instead.


# Using FFT

Initializing an array with random elements and copy it to the GPU:
```julia
A = rand(ComplexF64, 64, 64, 64)
d_A = cu(A)
```

We can calculate FFT of `d_A` directly using `fft` function from FFTW.
```julia
d_A_G = fft(d_A)
```

The function \jlinline{collect} can be used copy the result back to CPU.
```
A_G = collect(d_A_G)
```

# About initializing CuArrays

```
d_x = CuArrays.zeros(ComplexF64, 10, 10)
```

