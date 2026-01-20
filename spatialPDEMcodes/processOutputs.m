function [sl,sl2] = processOutputs(inputs,modelflag)

figure(1)
set(1,'pos',[1897         475         579         463])
N = 100;

if nargin < 2
    modelflag = 0;
end

switch modelflag
    case 1
        p0 =[0.8 3 0.005 2.1 0.0006];
        fitLab = 'Log';
    case 2
        p0 = [8.4 0.86 0.13 0.23 -0.9];
        fitLab = 'ExpInt';
    case 3
        p0 = [1 0.66 17.2 0.23 -0.9];
        fitLab = 'Radical Exp';
    case 4
        p0 = [1 0.66 17.2 -4.7];        
        fitLab = 'Bonev';
end

for j = 1:length(inputs)
    input = inputs{j};

    A0s = input.A0s;
    popDensity = input.popDensity;
    T = input.T;

    sl = semilogx(A0s,popDensity,'-','LineWidth',2,'DisplayName',['T: ',num2str(T),'h'],'HandleVisibility','off');
    hold on
    if modelflag > 0
        %popSize = sqrt(popDensity(1:N));

        popSize = popDensity(1:N);
        fitDR = fitnlm(popSize,A0s(1:N),@(par,Pop)DR(par,Pop,modelflag),p0)
        p = fitDR.Coefficients.Estimate;
        ps = 0.1:1e-4:1.02*max(popSize);
        sl2 = plot(fitDR.feval(ps),ps,'-','DisplayName',[fitLab,' fit'],'LineWidth',3);
    end
end

axis tight
xlabel('antibiotic dose (\mug/mL)')
ylabel('total bacterial density (OD)')
legend('boxoff')    
legend('location','northeast')

end

function Ac = DR(p,popSize,flag)
    %Ac = p(1) - Z(p(2:end),popSize,flag);

    Ac = Z(p(3:end),p(2)*abs(p(1)-popSize).^(1/2),flag);
end

function z = Z(p,r,flag)

    if flag == 1
        logF = p(3) + p(1)*exp((-1 + (1+p(2)*r.^2).^(1/2))/2 + ...
            log(((-1 + (1+p(2)*r.^2).^(1/2))))/2 + ...
            (p(2)/2)*r.^2./(-1 + (1+p(2)*r.^2).^(1/2)));
        z = logF;
        %b0F = [0.03 0.4 1.75];
    elseif flag == 2
        fexpint = p(3) + p(1) ./ expint(abs(p(2))*r.^2);
        z = fexpint;
        %b0ei = [0.5 0.02 1];
    elseif flag == 3
        fexpint = p(3) + p(1) .* exp(sqrt(1+abs(p(2))*r.^2)).*(-1 + sqrt(1+abs(p(2))*r.^2));
        z = fexpint;        
    else
        bonevF = p(2)*exp(-r.^2*(p(1)));    
        z = bonevF;
        %b0b = [1.3 -0.03];
    end

end