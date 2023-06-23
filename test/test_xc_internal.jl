ρ = [0.1, 0.2, 0.3, 0.4, 0.5]
@testset "LDA_VWN xc internal" begin
    
    Npoints = size(ρ,1)
    epsxc = zeros(Float64, Npoints)
    Vxc = zeros(Float64, Npoints)
    for i in 1:Npoints
        epsxc[i], Vxc[i] = PWDFT.XC_x_slater(ρ[i]) .+ PWDFT.XC_c_vwn(ρ[i])
    end

    xc_calc = LibxcXCCalculator(x_id=1, c_id=7)
    @test calc_epsxc_VWN(xc_calc, ρ) ≈ epsxc atol=1e-5
    @test calc_Vxc_VWN(xc_calc, ρ) ≈ Vxc atol=1e-5
end

