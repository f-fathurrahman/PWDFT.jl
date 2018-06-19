using Printf
using PWDFT

const ALL_PS_PADE_GTH = """
Ag-q11.gth  Bi-q5.gth   Cs-q1.gth   Ge-q4.gth   Kr-q8.gth   Nb-q5.gth   P-q5.gth    Sc-q11.gth  Te-q6.gth   Zn-q12.gth
Ag-q19.gth  B-q3.gth    Cs-q9.gth   He-q2.gth   La-q11.gth  Nd-q14.gth  Pr-q13.gth  Sc-q3.gth   Ti-q12.gth  Zn-q20.gth
Ag-q1.gth   Br-q7.gth   Cu-q11.gth  Hf-q12.gth  Li-q1.gth   Ne-q8.gth   Pt-q10.gth  Se-q6.gth   Ti-q4.gth   Zn-q2.gth
Al-q3.gth   Ca-q10.gth  Cu-q19.gth  Hg-q12.gth  Li-q3.gth   Ni-q10.gth  Pt-q18.gth  Si-q4.gth   Tl-q13.gth  Zr-q12.gth
Ar-q8.gth   Ca-q2.gth   Cu-q1.gth   Hg-q2.gth   Lu-q25.gth  Ni-q18.gth  Rb-q1.gth   Sm-q16.gth  Tl-q3.gth   Zr-q4.gth
As-q5.gth   Cd-q12.gth  Dy-q20.gth  Ho-q21.gth  Mg-q10.gth  N-q5.gth    Rb-q9.gth   Sn-q4.gth   Tm-q23.gth  
At-q7.gth   Cd-q2.gth   Er-q22.gth  H-q1.gth    Mg-q2.gth   O-q6.gth    Re-q15.gth  S-q6.gth    V-q13.gth   
Au-q11.gth  Ce-q12.gth  Eu-q17.gth  In-q13.gth  Mn-q15.gth  Os-q16.gth  Re-q7.gth   Sr-q10.gth  V-q5.gth    
Au-q19.gth  Cl-q7.gth   Fe-q16.gth  In-q3.gth   Mn-q7.gth   Os-q8.gth   Rh-q17.gth  Sr-q2.gth   W-q14.gth   
Au-q1.gth   Co-q17.gth  Fe-q8.gth   I-q7.gth    Mo-q14.gth  Pb-q4.gth   Rh-q9.gth   Ta-q13.gth  W-q6.gth    
Ba-q10.gth  Co-q9.gth   F-q7.gth    Ir-q17.gth  Mo-q6.gth   Pd-q10.gth  Rn-q8.gth   Ta-q5.gth   Xe-q8.gth   
Ba-q2.gth   C-q4.gth    Ga-q13.gth  Ir-q9.gth   Na-q1.gth   Pd-q18.gth  Ru-q16.gth  Tb-q19.gth  Yb-q24.gth  
Be-q2.gth   Cr-q14.gth  Ga-q3.gth   K-q1.gth    Na-q9.gth   Pm-q15.gth  Ru-q8.gth   Tc-q15.gth  Y-q11.gth   
Be-q4.gth   Cr-q6.gth   Gd-q18.gth  K-q9.gth    Nb-q13.gth  Po-q6.gth   Sb-q5.gth   Tc-q7.gth   Y-q3.gth
"""

const ALL_PS_PBE_GTH = """
Ag-q11.gth  Bi-q5.gth   Cs-q9.gth   Ge-q4.gth   K-q9.gth    Nd-q14.gth  Pm-q15.gth  Rn-q8.gth   Ta-q13.gth  W-q6.gth
Ag-q19.gth  B-q3.gth    Cu-q11.gth  He-q2.gth   Kr-q8.gth   Ne-q8.gth   Po-q6.gth   Ru-q16.gth  Ta-q5.gth   Xe-q8.gth
Al-q3.gth   Br-q7.gth   Cu-q19.gth  Hf-q12.gth  La-q11.gth  Ni-q18.gth  P-q5.gth    Ru-q8.gth   Tb-q19.gth  Yb-q24.gth
Ar-q8.gth   Ca-q10.gth  Dy-q20.gth  Hg-q12.gth  Li-q3.gth   N-q5.gth    Pr-q13.gth  Sb-q5.gth   Tc-q15.gth  Y-q11.gth
As-q5.gth   Cd-q12.gth  Er-q22.gth  Ho-q21.gth  Lu-q25.gth  O-q6.gth    Pt-q10.gth  Sc-q11.gth  Te-q6.gth   Zn-q12.gth
At-q7.gth   Ce-q12.gth  Eu-q17.gth  H-q1.gth    Mg-q10.gth  Os-q16.gth  Pt-q18.gth  Se-q6.gth   Ti-q12.gth  Zn-q20.gth
Au-q11.gth  Ce-q30.gth  Fe-q16.gth  In-q13.gth  Mg-q2.gth   Os-q8.gth   Rb-q9.gth   Si-q4.gth   Tl-q13.gth  Zr-q12.gth
Au-q19.gth  Cl-q7.gth   F-q7.gth    In-q3.gth   Mn-q15.gth  Pb-q14.gth  Re-q15.gth  Sm-q16.gth  Tl-q3.gth
Ba-q10.gth  Co-q17.gth  Ga-q13.gth  I-q7.gth    Mo-q14.gth  Pb-q4.gth   Re-q7.gth   Sn-q4.gth   Tm-q23.gth
Be-q4.gth   C-q4.gth    Ga-q3.gth   Ir-q17.gth  Na-q9.gth   Pd-q10.gth  Rh-q17.gth  S-q6.gth    V-q13.gth
Bi-q15.gth  Cr-q14.gth  Gd-q18.gth  Ir-q9.gth   Nb-q13.gth  Pd-q18.gth  Rh-q9.gth   Sr-q10.gth  W-q14.gth
"""

function test_read_all()
    
    ps_name = split(ALL_PS_PADE_GTH)
    for p in ps_name
        filename = "../pseudopotentials/pade_gth/"*p
        psp = PsPot_GTH(filename)
        println(psp)
    end

    ps_name = split(ALL_PS_PBE_GTH)
    for p in ps_name
        filename = "../pseudopotentials/pbe_gth/"*p
        psp = PsPot_GTH(filename)
        println(psp)
    end

end

function test_get_symbol()
    psp_dict = Dict()
    ps_name = split(ALL_PS_PADE_GTH)
    println("Dict(")
    for p in ps_name
        atsymb = split(p,"-")[1]
        chg = parse(Int64, split(split(p,"-")[2],".")[1][2:end] )
        #println("atsymb = ", atsymb)
        #merge!(psp_dict, Dict(atsymb => p))
        #println("chg = ", chg)
        println("\"",atsymb,"\" => ", "\"", p,"\",")
    end
    #println(psp_dict["H"])
    println(")")
end

test_get_symbol()
