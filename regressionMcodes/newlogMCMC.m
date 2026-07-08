clear all
close all

%%  data from literature
[conc,zoi,~] = defineMicrococcusNisinData();

weights = @(yhat) 1./(abs(yhat).^2);
b0 = [0.03 0.4 1.75];

logF = @(p,r)p(1) + p(3)*exp(sqrt(1+p(2).*r.^2 )) .* ( -1 + sqrt(1+p(2)*r.^2) );

fit = fitnlm(zoi(1:end-2),conc(1:end-2),logF,b0,'Weights',weights)
out = MCMCanalysis(conc(1:end-2),zoi(1:end-2),logF,fit.Coefficients.Estimate,weights);

figure(5)
semilogy(zoi,conc,'.k','markersize',30,'HandleVisibility','off');

figure(6)
MCMCshadedPlot(out,conc,zoi,'r');

%%

for j = [1 2 3 4 5 6]
    figure(j)
    if j == 4
        %subplot(2,2,1);
        hold on
        plot([0.625,1.25],[0,0],'-k','LineWidth',4)
        text(0.8,0.1,'W','Color','k','FontSize',25)
    end
    if j == 5
        ylim([0.4 90])
        xlim([0 11])        
        title('Predictive envelopes of the radical exponential model')
    end
    if j == 6
        hold on
        plot([0.625,1.25],[0,0],'-k','LineWidth',4,'DisplayName','W interval')
        text(0.75,0.5,'W','Color','k','FontSize',25)
        xlim([0.1 140])
        MIC = fit.Coefficients.Estimate(1);
        text(10,3,['MIC\approx',num2str(MIC,3),'\mug/mL']); 
    end

    exportgraphics(gcf,['./figures/newlogMCMC',num2str(j),'200k.pdf']);
end
