using DelimitedFiles
using Printf
using LaTeXStrings
using PGFPlotsX
using PWDFT: Ry2eV

include("../common/read_special_kpts.jl")

function test_plot_bands(; Nspin=1)

    filebands = "TEMP_bands.dat"

    symb_kpts_spec, x_kpts_spec = read_special_kpts(filebands)

    ebands = readdlm(filebands, comments=true)

    Nbands = size(ebands,2) - 1
    Nkspin = size(ebands,1)
    Nkpt = round(Int64, Nkspin/Nspin)

    E_f = 0.2492559573
    
    plt_inc = Array{Plot,1}(undef,Nbands*Nspin)
    
    for ispin in 1:Nspin, iband in 1:Nbands
        if ispin == 2
            idx_k = Nkpt+1:2*Nkpt
            line_style = "dashed"
        else
            idx_k = 1:Nkpt
            line_style = "solid"
        end
        kcart = ebands[idx_k,1]
        ene = (ebands[idx_k,iband+1] .- E_f)*2*Ry2eV
        iband_spin = iband + (ispin - 1)*Nbands
        
        if ispin == 2
            plt_inc[iband_spin] = @pgf PlotInc( {
                "solid",
                mark = "*",
                color = "red",
                mark_options = "fill=red",
                mark_size = "1pt",
            }, Coordinates(kcart, ene) )
        else
            plt_inc[iband_spin] = @pgf PlotInc( {
                "solid",
                mark = "*",
                color = "blue",
                mark_options = "fill=blue",
                mark_size = "1pt",
            }, Coordinates(kcart, ene) )
        end
    end

    kmin = 0.0
    kmax = maximum(ebands[:,1])

    fig = @pgf Axis( {
             title = "Ni band structure",
             ylabel = "Energy (eV)",
             height = "10cm",
             width = "8cm",
             xtick = x_kpts_spec,
             xticklabels = symb_kpts_spec,
             "xmajorgrids",
             "ymajorgrids",
             xmin=kmin,
             xmax=kmax,
             ymin = -11.5,
             ymax = 5.1,
          },
          plt_inc...)
    pgfsave("TEMP_bands.tex", fig)
    pgfsave("TEMP_bands.pdf", fig)

end

test_plot_bands(Nspin=2)

