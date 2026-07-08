clear all
close all

%%  data from literature
[conc,zoi,~] = defineMicrococcusNisinData();

I = @(x)exp(-x)./x;

myGamma = @(a)integral(I,a,inf);
weights = @(yhat) 1./(abs(yhat).^2);
fexpint = @(b,r)(b(1) + b(3) ./ expint(b(2)*r.^2));

fit = fitnlm(zoi(1:end-2),conc(1:end-2),fexpint,[0.5 0.02 1],'Weights',weights)
out = MCMCanalysis(conc(1:end-2),zoi(1:end-2),fexpint,fit.Coefficients.Estimate,weights);

figure(5)
semilogy(zoi,conc,'.k','markersize',30,'HandleVisibility','off');

figure(6)
MCMCshadedPlot(out,conc,zoi,'b');

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
        title('Predictive envelopes of the exponential integral model')
    end
    if j == 6
        hold on
        plot([0.625,1.25],[0,0],'-k','LineWidth',4,'DisplayName','W interval')
        text(0.75,0.5,'W','Color','k','FontSize',25)
        xlim([0.1 140])        
        MIC = fit.Coefficients.Estimate(1);
        text(10,3,['MIC\approx',num2str(MIC,3),'\mug/mL']);         
    end

    exportgraphics(gcf,['./figures/ExpIntMCMC',num2str(j),'200k.pdf']);
end
