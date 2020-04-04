using DelimitedFiles
import LightXML
using Printf

include("PsPot_UPF.jl")

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


function test_main()
    psp = PsPot_UPF("Pt.pz-hgh.UPF")
    println(psp)
    
    psp = PsPot_UPF("Si.pz-hgh.UPF")
    println(psp)

    psp = PsPot_UPF("Si_ONCV_PBE-1.0.upf")
    println(psp)

end

test_main()
