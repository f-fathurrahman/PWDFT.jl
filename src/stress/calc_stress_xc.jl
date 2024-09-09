function calc_stress_xc!( pw, potentials, Rhoe, Exc, stress_xc )

    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    Vxc = potentials.XC
    Evtxc = sum(Vxc[:,1] .* Rhoe)*dVol
    # FIXME: This is only for Nspin=1

    println("Exc = ", Exc)
    println("Evtxc = ", Evtxc)

    fill!(stress_xc, 0.0)
    for l in 1:3
        stress_xc[l,l] = -(Exc - Evtxc)/pw.CellVolume
    end

    return
end