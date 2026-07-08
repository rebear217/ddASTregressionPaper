function [sl,sl2] = processOutputsCell(inputs,modelflag,plotflag)

    figure(1)
    set(1,'pos',[1897         475         579         463])
    N = 165;
    
    COL1 = [0,1,1]*0.75;
    COL2 = [1,1,0]*0.75;
    
    if nargin < 2
        modelflag = 0;
    end
    if nargin < 3
        plotflag = 1;
    end
    
    switch modelflag
        case 1
            p0 =[0.8 3 0.005 2.1 0.0006];
            fitLab = 'Log';
        case 2
            p0 = [8.4 0.86 0.13 0.23 -0.9];
            fitLab = 'ExpInt';
        case 3
            p0 = [2 0.3 1.5 0.23 0.9];
            fitLab = 'Radical Exp';
        case 4
            p0 = ones(size([1 0.66 17.2 -4.7]));        
            fitLab = 'Bonev';
    end
    
    weight = @(y)(abs(y)).^(-2);
    M = length(inputs);
    for j = 1:M
        s = (j-1)/(M-1);
        COL = s*COL1 + (1-s)*COL2;
    
        input = inputs{j};
    
        A0s = input.A0s;
        popDensity = input.popDensity;
        T = input.T;
    
        if plotflag
            sl = semilogx(A0s,popDensity,'-','LineWidth',3,'DisplayName',...
                [num2str(T),'h dose response'],'color',COL);
            hold on
        else
            sl = [];
        end
        hold on
        popSize = popDensity(1:N);
        ps = 0.000:1e-4:1.3*max(popSize);
        if modelflag > 0
            if modelflag < 4
                fitDR = fitnlm(popSize,A0s(1:N),@(par,Pop)DR(par,Pop,modelflag),p0,'Weights',weight)
                sl2 = plot(fitDR.feval(ps),ps,'-','DisplayName',...
                    [fitLab,' fit (adj R^2\approx ',num2str(fitDR.Rsquared.Adjusted,3),')'],'LineWidth',3);
                %p = fitDR.Coefficients.Estimate;
            else
                fitDR = fitlm(log(A0s(1:N)),popSize,'Weights',weight)
                sl2 = plot(A0s,fitDR.feval(log(A0s)),'-','DisplayName',...
                    [fitLab,' fit (adj R^2\approx ',num2str(fitDR.Rsquared.Adjusted,3),')'],'LineWidth',3);            
            end
        end
    end
    
    axis tight
    xlabel('antibiotic dose (\mug/mL)')
    ylabel('total bacterial density (integrated S+R)')
    legend('boxoff')    
    legend('location','northeast')
    xlim([A0s(1) A0s(end)])

end

function Ac = DR(p,popSize,flag)
    Ac = Z(p(3:end),p(2)*abs(p(1)-popSize).^(1/2),flag);
end

function z = Z(p,r,flag)

    if flag == 1
        z = p(3) + p(1)*exp((-1 + (1+p(2)*r.^2).^(1/2))/2 + ...
            log(((-1 + (1+p(2)*r.^2).^(1/2))))/2 + ...
            (p(2)/2)*r.^2./(-1 + (1+p(2)*r.^2).^(1/2)));
        %b0F = [0.03 0.4 1.75];
    elseif flag == 2
        z = p(3) + p(1) ./ expint(abs(p(2))*r.^2);
        %b0ei = [0.5 0.02 1];
    elseif flag == 3
        z = p(1) + p(3)*exp(sqrt(1+abs(p(2)).*r.^2 )) .* ( -1 + sqrt(1+abs(p(2))*r.^2) );
    else
        z = p(2)*exp(-r.^2*(p(1)));    
        %b0b = [1.3 -0.03];
    end

end