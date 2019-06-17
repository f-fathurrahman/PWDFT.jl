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