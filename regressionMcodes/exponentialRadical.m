function [er,positivityPattern] = exponentialRadical(p,r)
    
    positivityPattern = [0,1,0];
    er = p(1) + exp(p(3) + sqrt(1+abs(p(2)).*r.^2 )) .* ( -1 + sqrt(1+abs(p(2))*r.^2) );

    % different formulation: p(3) > 0 for the following but not for the
    % above:
    % er = @(p,r)p(1) + p(3)*exp(sqrt(1+p(2).*r.^2 )) .* ( -1 + sqrt(1+p(2)*r.^2) );

end