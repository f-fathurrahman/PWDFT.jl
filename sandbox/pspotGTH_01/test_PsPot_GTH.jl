using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function test_read()
    psp = PsPot_GTH(joinpath(DIR_PSP,"Pt-q10.gth"))
    println(psp)
    psp = PsPot_GTH(joinpath(DIR_PSP,"Pt-q18.gth"))
    println(psp)
    psp = PsPot_GTH(joinpath(DIR_PSP,"Li-q3.gth"))
    println(psp)
    psp = PsPot_GTH(joinpath(DIR_PSP,"C-q4.gth"))
    println(psp)
end

function test_write_psp10()
    psp = PsPot_GTH(joinpath(DIR_PSP,"Pt-q10.gth"))
    write_psp10(psp)
end

test_read()
test_write_psp10()
