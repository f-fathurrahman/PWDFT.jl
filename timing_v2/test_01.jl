function time_PWGrid( Ham )
    # arguments
    ecutwfc = Ham.pw.ecutwfc
    LatVecs = Ham.atoms.LatVecs
    dual = Ham.pw.ecutrho/ecutwfc
    kpoints = Ham.pw.gvecw.kpoints
    Ns = Ham.pw.Ns

    # Timing PWGrid
    res = @be PWGrid( ecutwfc, LatVecs,
        kpoints=kpoints,
        Ns_=Ns,
        dual=dual
    )
    display(res)

    return
end

#=

Ham_gth = create_Ham_atom_Pt_gth()

Ham_oncv = create_Ham_atom_Pt_oncv()

Ham_gbrv = create_Ham_atom_Pt_gbrv()

=#
