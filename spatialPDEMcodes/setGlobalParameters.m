function globalParameters = setGlobalParameters()
    N = 200;

    globalParameters.mu = 0.000001;
    globalParameters.A0 = 1;
    globalParameters.C0 = 0.5;
    globalParameters.r = 0.9;
    globalParameters.reducedKill = 0.1;
    
    globalParameters.yield = 1;
    globalParameters.K = 0.25;
    globalParameters.V = 1;
    
    globalParameters.killing = 15;
    globalParameters.S0 = 0.01;
    
    globalParameters.antibioticPos = 7/10;
    
    globalParameters.diffusionRate = 0.000001;
    globalParameters.ABdiffusionRate = 0.0001;

    globalParameters.N = N;
    
    globalParameters.X = (0:(N-1))/(N-1);

    globalParameters.abBinding = 0.005;
    globalParameters.abDecay = 0.06;

end