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
several difficulties because the resulting type is not what I expected:
```
julia> v = zeros(Float64,10);

julia> d_v = cu(v);

julia> typeof(v), typeof(d_v)
(Array{Float64,1}, CuArray{Float32,1,Nothing})
```

However:
```
julia> v = rand(ComplexF64,10);

julia> d_v = cu(v);

julia> typeof(v), typeof(d_v)
(Array{Complex{Float64},1}, CuArray{Complex{Float64},1,Nothing})
```


I think it is advisable to use `CuArray` constructor directly.
```
julia> d_v = CuArray(v);

julia> typeof(v), typeof(d_v)
(Array{Float64,1}, CuArray{Float64,1,Nothing})
```

Using `cu` and `CuArray` involve copy from CPU to GPU. We can use
`CuArrays.fill` to initialize device array directly: (does not involve copy?)
```julia
d_v = CuArrays.fill( zero(Float64), 5) # or
d_v = CuArrays.fill( 0.0, 5 )

d_v = CuArrays.fill( zero(ComplexF64), (5,4,2) ) # or
d_v = CuArrays.fill( 0.0 + im*0.0, (5,4,2) )
```

The function collect can be used to copy the result back to CPU.
```
julia> d_v = d_v .+ 1.1;

julia> v = collect(d_v);
```


More functions:
```julia
d_x = CuArrays.zeros(ComplexF64, 10, 10)
d_c = CuArrays.rand(ComplexF64, 10, 10)
```


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


Use `CUDAnative.function_name` to replace various math functions.