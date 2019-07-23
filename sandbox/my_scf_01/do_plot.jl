using DelimitedFiles
using PGFPlotsX
using LaTeXStrings

function do_plot(sys_name::String, plot_title)

    betamix = 0.1:0.1:0.9

    prefix = sys_name*"_scf_history_betamix_"
    dat = readdlm(prefix*string(betamix[1]) * ".dat")

    plt_inc = Array{Plot}(undef,1)
    leg_inc = Array{LegendEntry}(undef,1)

    plt_inc[1] = @pgf PlotInc( {
            "solid",
        },
        Coordinates( dat[:,1], log10.(dat[:,3]) ) )
    leg_inc[1] = LegendEntry(L"\beta_{\mathrm{mix}}="*string(betamix[1]))

    for i in 3:2:9
        dat = readdlm(prefix*string(betamix[i])*".dat", comments=true)
        p = @pgf PlotInc( {
            "solid",
            },
            Coordinates( dat[:,1], log10.(dat[:,3]) ) )
        l = LegendEntry(L"\beta_{\mathrm{mix}}="*string(betamix[i]))
        push!(plt_inc, p)
        push!(leg_inc, l)
    end

    fig = @pgf Axis({
        title = plot_title,
        ylabel = L"\log(\Delta E)",
        xlabel = "SCF iteration",
        xmin = 0.0,
        xmax = 150,
        ymin = -8,
        ymax = 2,
        height = "6cm",
        width = "8cm",
        "xmajorgrids",
        "ymajorgrids",
        },
        plt_inc[1], leg_inc[1],
        plt_inc[2], leg_inc[2],
        plt_inc[3], leg_inc[3],
        plt_inc[4], leg_inc[4],
        plt_inc[5], leg_inc[5] )

    pgfsave("CONV_"*sys_name*".pdf", fig)
end

do_plot("H2", L"Convergence of $\mathrm{H}_{2}$")
do_plot("NH3", L"Convergence of $\mathrm{NH}_{3}$")
do_plot("Si_fcc", "Convergence of Si fcc")
do_plot("GaAs", "Convergence of GaAs")