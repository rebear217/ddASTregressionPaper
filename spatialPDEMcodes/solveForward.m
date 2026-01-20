function outState = solveForward(inState,T,parameters)

    N = parameters.N;
    options = odeset('NonNegative',ones(4*N,1));

    [~,soln] = ode113(@(t,x)f(x,parameters),[0 T],inState,options);
    outState = soln(end,:);
    outState = outState(:);

end