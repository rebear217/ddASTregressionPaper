function v = f(u,globalParameters)

    N = globalParameters.N;
    mu = globalParameters.mu;
    V = globalParameters.V;
    K = globalParameters.K;
    r = globalParameters.r;
    ab = globalParameters.abBinding;
    abD = globalParameters.abDecay;
    rK = globalParameters.reducedKill;
    killing = globalParameters.killing;
    
    h = 1/(N-1);
    
    diffusionRate = globalParameters.diffusionRate;
    ABdiffusionRate = globalParameters.ABdiffusionRate;
    yield = globalParameters.yield;

    S = u(1:N);
    R = u(N+1:2*N);
    A = u(2*N+1:3*N);
    C = u(3*N+1:4*N);
    
    uptakeRate = V*C./(C+K);
    gRate = yield*uptakeRate;
    
    dC = ABdiffusionRate*D(C) - uptakeRate.*(S+R);
    dA = ABdiffusionRate*D(A) - A*abD - (S+R).*A*ab;
    dR = r*R.*gRate - rK*killing*A.*R + mu*S + diffusionRate*D(R);
    dS = S.*gRate - killing*A.*S + diffusionRate*D(S) - mu*S;
    
    v = [dS ; dR ; dA; dC];
    
end
