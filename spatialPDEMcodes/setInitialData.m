function [iC,globalParameters] = setInitialData(globalParameters,option)

    if nargin < 2
        option = 1;
    end

    N = globalParameters.N;
    antibioticPos = globalParameters.antibioticPos;

    C = ones(N,1)*globalParameters.C0;
    A = zeros(N,1);
    R = zeros(N,1);
    S = ones(N,1);
    
    M = floor(antibioticPos*N);

    A(M:end) = globalParameters.A0;
    
    switch option
        case 1
            S = R;
            S(1:10) = globalParameters.S0/10;
            R(1:10) = S(1:10)/1000;
            globalParameters.abDecay = 0.01;
        case 2
            S = globalParameters.S0*S;
        case 3
            S = globalParameters.S0*S;
            R = S;
        otherwise
            S = globalParameters.S0*S;
            A = 0*A;
            R = 0*R;
            S(M:end) = 0;
    end
    
    iC = [S ; R ; A ; C];
end