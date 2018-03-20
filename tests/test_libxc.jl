const LIBXC_SO_PATH = "/home/efefer/WORKS/my_github_repos/PWDFT.jl/src/extlibs/libxc_interface.so"

# XXX: Need to set environment variables LD_LIBRARY_PATH to directory containing
# libxc.so.4

function test_main()
    Npoints = 5
    Rhoe = zeros( Float64, Npoints )
    Rhoe = [0.1, 0.2, 0.3, 0.4, 0.5]

    Vxc = zeros( Float64, Npoints )
    epsxc = zeros( Float64, Npoints )
    #
    ccall( (:calc_Vxc_VWN, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, Vxc )
    #
    ccall( (:calc_epsxc_VWN, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, epsxc )

    for ip = 1:Npoints
        @printf("%3d %18.10f %18.10f %18.10f\n", ip, Rhoe[ip], epsxc[ip], Vxc[ip])
    end
end

test_main()
