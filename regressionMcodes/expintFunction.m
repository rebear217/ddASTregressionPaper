function [fexpint,positivityPattern] = expintFunction(b,r)

    positivityPattern = [1,1,0];
    fexpint = abs(b(1)) + b(3) ./ expint(abs(b(2))*r.^2);

    %I = @(x)exp(-x)./x;
    %myGamma = @(a)integral(I,a,inf);
    %fexpint = @(b,r)(abs(b(1)) + b(3) ./ expint(abs(b(2))*r.^2));
    %Esuper = @(r)(0.5*exp(-r).*log(1+2./r));
    %fsuper = @(b,r)(abs(b(1)) + b(3) ./ Esuper(abs(b(2))*r.^2));
    %Esub = @(r)(exp(-r).*log(1+1./r));
    %fsub = @(b,r)(abs(b(1)) + b(3) ./ Esub(abs(b(2))*r.^2));    
end