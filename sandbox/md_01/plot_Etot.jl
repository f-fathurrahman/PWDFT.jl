using DelimitedFiles

import PyPlot
const plt = PyPlot

function main(prefix::String; dt=1.0)
    data = readdlm(prefix)
    Ndata = size(data,1)
    @views Etot = data[2:Ndata,2]
    Etot0 = Etot[1]
    plt.clf()
    plt.plot(data[2:Ndata,1], Etot .- Etot0, marker="o")
    plt.grid()
    plt.title(prefix)
    plt.savefig("IMG_"*prefix*".pdf")
end

#main("ETOT_H2O.dat")
main("ETOT_CO2_v3.dat")
#main("ETOT_H2O_v4.dat")
