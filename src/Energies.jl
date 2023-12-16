"""
The type for collecting all energy terms: kinetic, local and
nonlocal pseudopotential, Hartree, XC, nuclei-nuclei terms
and electronic entropy term.
"""
mutable struct Energies
    Kinetic::Float64
    Ps_loc::Float64
    Ps_nloc::Float64
    Hartree::Float64
    XC::Float64
    NN::Float64
    mTS::Float64
    EHxc_paw::Float64
    # EHxc_paw is a copy of EHxc_paw in PAWVariables
end

"""
Creates an instance of `Energies` with value of zeros for all fields.
"""
function Energies()
    return Energies(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

import Base: sum

"""
Get total energy by summing all of its components.
"""
function sum( energies::Energies )
    # irrespective of whether we are using smearing or not
    return energies.Kinetic + energies.Ps_loc + energies.Ps_nloc +
           energies.Hartree + energies.XC + energies.NN + energies.mTS + energies.EHxc_paw
end


import Base: copy
function copy( e1::Energies )
    return Energies(
        e1.Kinetic,
        e1.Ps_loc ,
        e1.Ps_nloc,
        e1.Hartree,
        e1.XC     , 
        e1.NN     , 
        e1.mTS    , 
        e1.EHxc_paw
    )
end


import Base: +
function +( e1::Energies, e2::Energies )
    return Energies(
        e1.Kinetic + e2.Kinetic,
        e1.Ps_loc  + e2.Ps_loc,
        e1.Ps_nloc + e2.Ps_nloc,
        e1.Hartree + e2.Hartree,
        e1.XC      + e2.XC,
        e1.NN      + e2.NN,
        e1.mTS     + e2.mTS,
        e1.EHxc_paw + e2.EHxc_paw
    )
end

import Base: -
function -( e1::Energies, e2::Energies )
    return Energies(
        e1.Kinetic - e2.Kinetic,
        e1.Ps_loc  - e2.Ps_loc,
        e1.Ps_nloc - e2.Ps_nloc,
        e1.Hartree - e2.Hartree,
        e1.XC      - e2.XC,
        e1.NN      - e2.NN,
        e1.mTS     - e2.mTS,
        e1.EHxc_paw - e2.EHxc_paw
    )
end

include("Energies_io.jl")
