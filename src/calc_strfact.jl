# Calculate structure factor
function calc_strfact( atpos::Array{Float64,2}, Nspecies::Int,
    atm2species::Array{Int,1}, G::Array{Float64,2} )
    #
    Ng = size(G)[2]
    Na = size(atpos)[2]
    Sf = zeros(ComplexF64,Ng,Nspecies)
    for ia = 1:Na
        isp = atm2species[ia]
        for ig = 1:Ng
            GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
            Sf[ig,isp] = Sf[ig,isp] + cos(GX) - im*sin(GX)
        end
    end
    return Sf
end

function calc_strfact( atoms::Atoms, pw::PWGrid )
    return calc_strfact( atoms.positions, atoms.Nspecies, atoms.atm2species, pw.gvec.G )
end
