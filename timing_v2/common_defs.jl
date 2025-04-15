using Revise
using Chairmarks
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

const DIR_PSP_LDA_GTH = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_PSP_PBE_GTH = joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth")

const DIR_PSP_LDA_ONCV = joinpath(DIR_PWDFT, "pseudopotentials", "ONCV_v0.4.1_LDA")
const DIR_PSP_PBE_ONCV = joinpath(DIR_PWDFT, "pseudopotentials", "ONCV_v0.4.1_PBE")

const DIR_PSP_LDA_GBRV = joinpath(DIR_PWDFT, "pseudopotentials", "GBRV_LDA")
const DIR_PSP_PBE_GBRV = joinpath(DIR_PWDFT, "pseudopotentials", "GBRV_PBE")

const DIR_PSP_LDA_PAW_JTH = joinpath(DIR_PWDFT, "pseudopotentials", "PAW_JTH_LDA")
const DIR_PSP_PBE_PAW_JTH = joinpath(DIR_PWDFT, "pseudopotentials", "PAW_JTH_PBE")

const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_TIMING = joinpath(DIR_PWDFT, "timing")

includet("prepare_Ham_various.jl")

