using DelimitedFiles
import LightXML

function traverse_upf(upf_file)
    xdoc = LightXML.parse_file(upf_file)
    # get the root element
    xroot = LightXML.root(xdoc)  # an instance of XMLElement
    # print its name
    println(LightXML.name(xroot))
    for c in LightXML.child_nodes(xroot)  # c is an instance of XMLNode
        println(LightXML.nodetype(c))
        if LightXML.is_elementnode(c)
            e = LightXML.XMLElement(c)  # this makes an XMLElement instance
            println(LightXML.name(e))
        end
    end
end


function parse_upf(upf_file)
    
    xdoc = LightXML.parse_file(upf_file)

    # get the root element
    xroot = LightXML.root(xdoc)  # an instance of XMLElement
    
    pp_mesh = LightXML.get_elements_by_tagname(xroot, "PP_MESH")

    pp_header = LightXML.get_elements_by_tagname(xroot, "PP_HEADER")

    Nr = parse(Int64,LightXML.attributes_dict(pp_mesh[1])["mesh"])
    println("Nr = ", Nr)
    
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
        V_local[i] = parse(Float64, spl_str[i])
    end

    Nproj = parse(Int64,LightXML.attributes_dict(pp_header[1])["number_of_proj"])
    println("Nproj = ", Nproj)

    pp_nonlocal = LightXML.get_elements_by_tagname(xroot, "PP_NONLOCAL")

    proj_func = zeros(Float64,Nr,Nproj)
    for iprj = 1:Nproj
        pp_beta = LightXML.get_elements_by_tagname(pp_nonlocal[1], "PP_BETA."*string(iprj))
        pp_beta_str = LightXML.content(pp_beta[1])
        pp_beta_str = replace(pp_beta_str, "\n" => " ")
        spl_str = split(pp_beta_str, keepempty=false)
        for i = 1:Nr
            proj_func[i,iprj] = parse(Float64,spl_str[i])
        end
        #plot_func_radial(r, proj_func[:,iprj], "TEMP_PP_BETA."*string(iprj)*".pdf")
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
    Dij = reshape(Dij_temp,(Nproj,Nproj))*2

    display(Dij)
    println()

    println("Pass here")

    return
end

using PGFPlotsX

function plot_func_radial(r, fr, plotfile::String)
    tableData = Table( ["x" => r, "y" => fr] )
    fig =
    @pgf Axis(
        {
            xlabel = "r",
            ylabel = "fr",
            width = "10cm",
            xmin = 0.0,
            xmax = 2.0,
            ymax = 3.0,
        },
        PlotInc({
                mark = "none",
            },
            tableData),
    )
    pgfsave(plotfile, fig)
    return
end


function parse_gpaw_data(gpaw_data)
    xdoc = LightXML.parse_file(gpaw_data)
    # get the root element
    xroot = LightXML.root(xdoc)  # an instance of XMLElement
    radial_grid = LightXML.get_elements_by_tagname(xroot, "radial_grid");
    #pp_r_str = content(get_elements_by_tagname(pp_mesh[1], "PP_R")[1]);
    #r = readdlm(IOBuffer(pp_r_str))
end

#@code_native traverse_upf("Ni.pbe-nd-rrkjus.UPF")
#traverse_upf("Ni.pbe-nd-rrkjus.UPF")
#@time traverse_upf("Ni.PBE")
@time parse_upf("Ni.pbesol-n-kjpaw_psl.0.1.UPF")

@time parse_upf("Pt.pz-hgh.UPF")

#traverse_upf("Si.pz-hgh.UPF")

#parse_upf("Si.pz-hgh.UPF")
#parse_upf("Ni.pbe-nd-rrkjus.UPF")

