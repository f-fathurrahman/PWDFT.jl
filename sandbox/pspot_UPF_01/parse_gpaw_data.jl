function parse_gpaw_data(gpaw_data)
    xdoc = LightXML.parse_file(gpaw_data)
    # get the root element
    xroot = LightXML.root(xdoc)  # an instance of XMLElement
    radial_grid = LightXML.get_elements_by_tagname(xroot, "radial_grid");
    #pp_r_str = content(get_elements_by_tagname(pp_mesh[1], "PP_R")[1]);
    #r = readdlm(IOBuffer(pp_r_str))
end