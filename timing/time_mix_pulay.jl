function time_mix_pulay()

  Npoints = 500

	# create mock data

  Rhoe = reshape( mix_pulay!(
                reshape(Rhoe,(Npoints*Nspin)),
                reshape(Rhoe_new,(Npoints*Nspin)), betamix, XX, FF, iter, MIXDIM, x_old, f_old
                ), (Npoints,Nspin) )

end
