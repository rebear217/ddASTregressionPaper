clear all
close all

%%  data from literature

%Micrococcus luteus (ATCC 10240) and nisin
%data taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC4576963/pdf/fsn30003-0394.pdf
%Figure 7 (see legend of A-J)

conc = [125 75 50 25 18.75 12.5 6.25 2.5 1.25 0.625];
zoi = [10.7 10.5 10.3 9.5 9 8 6.5 4 0.0001 0.000001];

weights = @(yhat) 1./(1 + abs(yhat).^2);
b0 = [0.03 0.4 1.75];

logF = @(p,r)p(1) + p(3)*exp(sqrt(1+p(2).*r.^2 )) .* ( sqrt(1+p(2)*r.^2) );

fit = fitnlm(zoi(3:end),conc(3:end),logF,b0,'Weights',weights)

MCMCanalysis(conc(3:end),zoi(3:end),logF,fit.Coefficients.Estimate,weights);
figure(5)
semilogy(zoi,conc,'.k','markersize',30,'HandleVisibility','off');

for j = [1 2 3 4 5]
    figure(j)
    if j == 4
        %subplot(2,2,1);
        hold on
        plot([0.625,1.25],[0,0],'-k','LineWidth',4)
        text(0.8,0.14,'W','Color','k','FontSize',25)
    end
    if j == 5
        ylim([0.4 90])
        xlim([0 11])        
        title('Predictive envelopes of the radical exponential model')
    end

    exportgraphics(gcf,['./figures/newlogMCMC',num2str(j),'200k.pdf']);
end
