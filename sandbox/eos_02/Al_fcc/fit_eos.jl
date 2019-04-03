using Printf
using DelimitedFiles: readdlm
using PWDFT: Ry2eV, ANG2BOHR
using LsqFit

function main()
    dat = readdlm("TEMP_EOS_data.dat")
    # calculate primitive fcc volumes and convert it to angstrom^3
    volumes = 0.25*dat[:,1].^3 / (ANG2BOHR^3)
    # energies in eV
    energies = dat[:,2]*2*Ry2eV
    energies_abinit = dat[:,3]*2*Ry2eV
    energies_pwscf = dat[:,4]*2*Ry2eV

    params1 = do_fit_parabolic(volumes, energies)
    params2 = do_fit_sjeos(volumes, energies)
end

function do_fit_sjeos(volumes, energies)
    # 
    @. model(x, p) = p[1] + p[2]*x^(-1/3) + p[3]*x^(-2/3) + p[4]*x^(-1)

    p0 = [0.5, 0.5, 0.5, 0.5]

    fit = curve_fit(model, volumes, energies, p0, show_trace=true)
    
    println("Residuals:")
    println(fit.resid)
    
    dump(fit)

    p = fit.param
    Ndata = length(volumes)
    for i in 1:Ndata
        V = volumes[i]
        E = energies[i]
        E_m = model(volumes[i],p)
        diffE = abs(E - E_m)  # this is already calculated in residuals
        @printf("%18.10f %18.10f %18.10f %18.10e\n", V, E, E_m, diffE)
    end
    return p
end


function do_fit_parabolic(volumes, energies)
    # parabolic model
    @. model(x, p) = p[1] + p[2]*x + p[3]*x^2

    p0 = [1.0, 1.0, 1.0]

    fit = curve_fit(model, volumes, energies, p0, show_trace=true)
    
    println("Residuals:")
    println(fit.resid)
    
    dump(fit)

    p = fit.param

    Ndata = length(volumes)
    for i in 1:Ndata
        V = volumes[i]
        E = energies[i]
        E_m = model(volumes[i],p)
        diffE = abs(E - E_m)  # this is already calculated in residuals
        @printf("%18.10f %18.10f %18.10f %18.10e\n", V, E, E_m, diffE)
    end
    return p
end

main()