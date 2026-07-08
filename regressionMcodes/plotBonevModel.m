%clear all
close all
clc

JACK = [0.2 0.5 0.2];
bonevF = @(p,r)p(1).*exp(-r.^2*p(2));

%% get data from literature

[conc,zoi,R] = defineMicrococcusNisinData();
disp('Figure 3A')

figure(1)
semilogx(conc,zoi,'.k','markersize',38,'DisplayName','ZoI data')
hold on

b0 = [1 0];
weights = @(yhat) 1./(abs(yhat).^2);

%weights = @(yhat) ones(size(yhat));
fit2 = fitnlm(zoi(1:end-2),conc(1:end-2),bonevF,b0,'Weights',weights)

M = 1;
for j = M:9
    Z = setdiff(zoi(M:end-2),zoi(j));
    C = setdiff(conc(M:end-2),conc(j));
    fit0 = fitnlm(Z,C,bonevF,fit2.Coefficients.Estimate,'Weights',weights);   
    if j == M
        plot(fit0.feval(R),R,'-','DisplayName','Jacknife frequentist fits','linewidth',2,'color',JACK)
    else
        plot(fit0.feval(R),R,'-','linewidth',2,'color',JACK,'HandleVisibility','off')
    end
end

MIC = abs(fit2.Coefficients.Estimate(1));

adjR2 = fit2.Rsquared.Adjusted;
plot(fit2.feval(R),R,'-k','DisplayName',['Bonev frequentist fit (adj R^2\approx',num2str(adjR2,3),')'],'linewidth',3)

semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
text(20,3,['MIC\approx',num2str(MIC,3),'\mug/mL']);

%semilogx(conc(1:2),zoi(1:2),'ok','markersize',14,'DisplayName','excised data','linewidth',1)
set(gca,'Ytick',0:1:max(R))

W1 = conc(end);
W2 = conc(end-1);
plot([W1,W2],[0,0],'-k','linewidth',6,'DisplayName','MIC ground truth (W)');
text(0.8,0.5,'W','FontSize',22);

MMSEMmicroNisinBonevF = MCMCupdatePointEstimates(conc(1:end-2),zoi(1:end-2),bonevF,fit2.Coefficients.Estimate,weights);

ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','northwest')
axis tight
xlim([0 150])
ylim([0 max(R)])

%% get data from literature

[conc,zoi,R] = defineSarcinaCloxData();
disp('Figure 2A')

fit = fitnlm(zoi,conc,bonevF,b0,'Weights',weights)

figure(2)
semilogx(conc,zoi,'.k','markersize',38,'DisplayName','ZoI data')
hold on

set(gca,'Ytick',0:5:max(R))

for j = 1:5
    Z = setdiff(zoi,zoi(j));
    C = setdiff(conc,conc(j));
    fitJ = fitnlm(Z,C,bonevF,fit.Coefficients.Estimate,'Weights',weights);   
    if j == 1
        plot(fitJ.feval(R),R,'-','DisplayName','Jacknife frequentist fits','linewidth',2,'color',JACK)
    else
        plot(fitJ.feval(R),R,'-','linewidth',2,'color',JACK,'HandleVisibility','off')
    end
end

adjR2 = fit.Rsquared.Adjusted;
plot(fit.feval(R),R,'-k','DisplayName',['Bonev frequentist fit (adj R^2\approx',num2str(adjR2,3),')'],'linewidth',3)

semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
MIC = abs(fit.Coefficients.Estimate(1));
text(20,5,['MIC\approx',num2str(MIC,3),'\mug/mL']);

MMSEMsarcinaCloxBonevF = MCMCupdatePointEstimates(conc,zoi,bonevF,fit.Coefficients.Estimate,weights);

ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','northwest')
axis tight
xlim([0.9 100])
ylim([0 max(R)])

%%

plotON = 1;
if plotON
    figure(1)
    exportgraphics(gcf,'./figures/Micrococcus_luteus_bonev.pdf')
    figure(2)
    exportgraphics(gcf,'./figures/Bacillus_subtilis_bonev.pdf')
end
