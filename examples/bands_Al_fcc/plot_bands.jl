using DelimitedFiles
using LaTeXStrings
using PGFPlotsX
using PWDFT: Ry2eV

include("../common/read_special_kpts.jl")

function test_plot_bands()

    filebands = "TEMP_bands.dat"

    symb_kpts_spec, x_kpts_spec = read_special_kpts(filebands)

    ebands = readdlm(filebands, comments=true)
    Nbands = size(ebands)[2] - 1
    kcart = ebands[:,1]

    E_f = 0.3598017974

    plt_inc = Array{Plot,1}(undef,Nbands)
    for iband = 1:Nbands
        plt_inc[iband] = @pgf PlotInc( {
            "solid",
            mark = "*",
            color = "black",
            mark_options = "fill=cyan",
            mark_size = "1pt",
        }, Coordinates(kcart, (ebands[:,iband+1] .- E_f)*2*Ry2eV) )
    end

    kmin = 0.0
    kmax = maximum(kcart)

    fig = @pgf Axis( {
             title = "Al band structure",
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

test_plot_bands()

