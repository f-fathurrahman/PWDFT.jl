# x = (ϵ - μ)/(2.0*width)

function smear_fermi(μ, ϵ, width)
   x = (ϵ - μ)/(2.0*width)
   return 0.5*(1.0 - tanh(x))
end

function smear_gauss(μ, ϵ, width)
    x = (ϵ - μ)/(2.0*width)
    return 0.5*erfc(x)
end

function smear_cold()
    x = (ϵ - μ)/(2.0*width)
    return 0.5*erfc(x+sqrt(0.5)) + exp( -(x + sqrt(0.5))^2 )/sqrt(2*pi)
end

function smear_fermi_entropy(μ, ϵ, width)
    x = (ϵ - μ)/(2.0*width)
    f = 0.5*( 1. - tanh(x) )
    S = 0.0
    if f > 1e-300
        S = S - f*log(f);
    if (1-f) >1e-300
        S = S - (1-f)*log(1-f)
    end
end
# case SmearingGauss: return exp(-x*x) / sqrt(M_PI);
# case SmearingCold: return exp(-std::pow(x+sqrt(0.5),2)) * (1.+x*sqrt(2.)) / sqrt(M_PI);



function find_mu(const std::vector<diagMatrix>& eps, double nElectrons, double& Bz) const
{   
    const bool& verbose = e->cntrl.shouldPrintMuSearch; // orig

    if(verbose) logPrintf("\nBisecting to find mu(nElectrons=%.15le)\n", nElectrons);
    
    //Find a range which is known to bracket the result:
    const double absTol = 1e-10, relTol = 1e-14;
    double nTol = std::max(absTol, relTol*fabs(nElectrons));
    double muMin=-0.1, muMax=+0.0;
    
    while(nElectronsCalc(muMin,eps,Bz)>=nElectrons+nTol) {
        muMin-=(muMax-muMin);
    }

    while(nElectronsCalc(muMax,eps,Bz)<=nElectrons-nTol) {
        muMax += (muMax-muMin);
    }
    
    //Bisect:
    double muTol = std::max(absTol*smearingWidth, relTol*std::max(fabs(muMin),fabs(muMax)));
    while(muMax-muMin>=muTol)
    {
        double mu = 0.5*(muMin + muMax);
        double N = nElectronsCalc(mu, eps, Bz);
        if(verbose) logPrintf("MUBISECT: mu = [ %.15le %.15le %.15le ]  N = %le\n", muMin, mu, muMax, N);
        if(N>nElectrons) muMax = mu;
        else muMin = mu;
    }
    return 0.5*(muMin + muMax);
}