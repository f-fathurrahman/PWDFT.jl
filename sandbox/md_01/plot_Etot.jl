using DelimitedFiles

import PyPlot
const plt = PyPlot

function main(prefix::String; dt=1.0)
    Etot = readdlm(prefix)
    Etot0 = Etot[1]

    Ndata = length(Etot)
    plt.clf()
    plt.plot(range(0,stop=Ndata-1)*dt, Etot .- Etot0, marker="o")
    plt.grid()
    plt.savefig("IMG_"*prefix*".pdf")
end

#main("LOG_ETOT1")
#main("LOG_ETOT2")
#main("LOG_ETOT3", dt=0.5)
main("LOG_ETOT4")