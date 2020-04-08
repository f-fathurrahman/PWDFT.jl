struct PsPot_UPF
    pspfile::String
    atsymb::String
    zval::Float64
    Nr::Int64
    r::Array{Float64,1}
    rab::Array{Float64,1}
    V_local::Array{Float64,1}
    Nproj::Int64
    proj_l::Array{Int64,1}
    rcut_l::Array{Float64,1}
    proj_func::Array{Float64,2}
    Dij::Array{Float64,2}
end

function PsPot_UPF( upf_file::String )
    
    xdoc = LightXML.parse_file(upf_file)

    # get the root element
    xroot = LightXML.root(xdoc)  # an instance of XMLElement
    
    pp_header = LightXML.get_elements_by_tagname(xroot, "PP_HEADER")

    atsymb = LightXML.attributes_dict(pp_header[1])["element"]
    zval = Int64( parse( Float64, LightXML.attributes_dict(pp_header[1])["z_valence"] ) )
    lmax = parse( Int64, LightXML.attributes_dict(pp_header[1])["l_max"] )
    Nr = parse(Int64,LightXML.attributes_dict(pp_header[1])["mesh_size"])

    pp_mesh = LightXML.get_elements_by_tagname(xroot, "PP_MESH")

    #Nr = parse(Int64,LightXML.attributes_dict(pp_mesh[1])["mesh"])
    
    pp_r = LightXML.get_elements_by_tagname(pp_mesh[1], "PP_R")
    pp_r_str = LightXML.content(pp_r[1])
    pp_r_str = replace(pp_r_str, "\n" => " ")
    spl_str = split(pp_r_str, keepempty=false)

    @assert(length(spl_str) == Nr)
    
    r = zeros(Float64, Nr)
    for i = 1:Nr
        r[i] = parse(Float64, spl_str[i])
    end

    pp_rab = LightXML.get_elements_by_tagname(pp_mesh[1], "PP_RAB")
    pp_rab_str = LightXML.content(pp_rab[1])
    pp_rab_str = replace(pp_rab_str, "\n" => " ")
    spl_str = split(pp_rab_str, keepempty=false)

    @assert(length(spl_str) == Nr)

    rab = zeros(Float64, Nr)
    for i = 1:Nr
        rab[i] = parse(Float64, spl_str[i])
    end


    pp_local = LightXML.get_elements_by_tagname(xroot, "PP_LOCAL")
    pp_local_str = LightXML.content(pp_local[1])
    pp_local_str = replace(pp_local_str, "\n" => " ")
    spl_str = split(pp_local_str, keepempty=false)

    @assert(length(spl_str) == Nr)

    V_local = zeros(Float64, Nr)
    for i = 1:Nr
        V_local[i] = parse(Float64, spl_str[i])*0.5 # convert to Hartree
    end

    Nproj = parse(Int64,LightXML.attributes_dict(pp_header[1])["number_of_proj"])

    pp_nonlocal = LightXML.get_elements_by_tagname(xroot, "PP_NONLOCAL")

    proj_func = zeros(Float64,Nr,Nproj)
    proj_l = zeros(Int64,Nproj)
    rcut_l = zeros(Float64,Nproj)

    for iprj = 1:Nproj
        pp_beta = LightXML.get_elements_by_tagname(pp_nonlocal[1], "PP_BETA."*string(iprj))
        #
        proj_l[iprj] = parse( Int64, LightXML.attributes_dict(pp_beta[1])["angular_momentum"] )
        rcut_l[iprj] = parse( Float64, LightXML.attributes_dict(pp_beta[1])["cutoff_radius"] )
        #
        pp_beta_str = LightXML.content(pp_beta[1])
        pp_beta_str = replace(pp_beta_str, "\n" => " ")
        spl_str = split(pp_beta_str, keepempty=false)
        for i = 1:Nr
            proj_func[i,iprj] = parse(Float64,spl_str[i])*0.5 # Convert to Hartree
        end
    end

    Dij = zeros(Nproj,Nproj)
    Dij_temp = zeros(Nproj*Nproj)

    pp_dij = LightXML.get_elements_by_tagname(pp_nonlocal[1], "PP_DIJ")
    pp_dij_str = LightXML.content(pp_dij[1])
    pp_dij_str = replace(pp_dij_str, "\n" => " ")
    spl_str = split(pp_dij_str, keepempty=false)
    for i = 1:Nproj*Nproj
        Dij_temp[i] = parse(Float64,spl_str[i])
    end
    Dij = reshape(Dij_temp,(Nproj,Nproj))*2  # convert to Hartree

    return PsPot_UPF(upf_file, atsymb, zval, Nr, r, rab, V_local, Nproj, proj_l, rcut_l, proj_func, Dij)
end


import Base: println
function println(psp::PsPot_UPF)
    @printf("Nr = %d\n", psp.Nr)
    for i = 1:psp.Nproj
        @printf("l=%d rcut=%f\n", psp.proj_l[i], psp.rcut_l[i])
    end
    display(psp.Dij)
    println()
end