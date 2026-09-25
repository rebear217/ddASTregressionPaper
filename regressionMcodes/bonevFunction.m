function [b,positivityPattern] = bonevFunction(p,r)
    positivityPattern = [0,0];
    b = p(1).*exp(-r.^2*p(2));
end