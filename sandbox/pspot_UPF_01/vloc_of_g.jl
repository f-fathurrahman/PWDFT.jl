# This routine computes the Fourier transform of the local
# part of an atomic pseudopotential, given in numerical form.
# A term erf(r)/r is subtracted in real space (thus making the
# function short-ranged) and added again in G space (for G<>0)
# The G=0 term contains \int (V_loc(r)+ Ze^2/r) 4pi r^2 dr.
# This is the "alpha" in the so-called "alpha Z" term of the energy.
#
# Adapted from the file vloc_of_g.f90
#
# Copyright (C) 2001-2007 Quantum ESPRESSO group
# 

function init_Vloc_G(
	Nr::Int64,
    r::Array{Float64,1},
	rab::Array{Float64,1},
	Vloc_at::Array{Float64},
	Zval::Float64,
	Ngl::Int64, gl::Array{Float64,1},
	CellVolume::Float64
)
 	
  	Vloc_G = zeros(Float64, Ngl)
  	aux = zeros(Float64, Nr)
  	aux1 = zeros(Float64, Nr)

  	if gl[1] < 1e-8
    	# first the G=0 term
        for ir in 1:Nr
           aux[ir] = r[ir] * ( r[ir] * Vloc_at[ir] + Zval )
        end
    	Vloc_G[1] = integ_simpson(Nr, aux, rab)
    	igl0 = 2
  	else
    	igl0 = 1
  	end
    
    println("igl0 = ", igl0)

  	# here the G != 0 terms, we first compute the part of the integrand 
  	# function independent of |G| in real space
  	for ir in 1:Nr
  	   aux1[ir] = r[ir] * Vloc_at[ir] + Zval * erf( r[ir] )
  	end
 
  	for igl in igl0:Ngl
     	Gx = sqrt( gl[igl] )
     	for ir in 1:Nr
     	  	aux[ir] = aux1[ir] * sin(Gx*r[ir])/Gx
     	end
     	Vloc_G[igl] = integ_simpson( Nr, aux, rab )
	end

  	Vloc_G[:] = Vloc_G[:] * 4*pi #/ CellVolume
  	return Vloc_G
end
