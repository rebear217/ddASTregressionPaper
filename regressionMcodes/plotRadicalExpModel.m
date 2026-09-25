%clear all
close all
clc

JACK = [0.2 0.5 0.2];
expRad = @exponentialRadical;

%%  data from literature

[conc,zoi,R] = defineMicrococcusNisinData();
b0 = [1.75 0.4 0.03];

figure(1)
disp('Figure 3C')

semilogx(conc,zoi,'.k','markersize',34,'DisplayName','ZoI data')
hold on

weights = @(yhat) 1./(abs(yhat).^2);
%weights = @(yhat) ones(size(yhat));
fit2 = fitnlm(zoi(1:end-2),conc(1:end-2),expRad,b0,'Weights',weights(conc(1:end-2)))

M = 1;
for j = M:9
    %Jackknife regressions:
    Z = myRemoveDatum(zoi(M:end-2),zoi(j));
    C = myRemoveDatum(conc(M:end-2),conc(j));

    fit0 = fitnlm(Z,C,expRad,fit2.Coefficients.Estimate,'Weights',weights(C));   
    if j == M
        plot(fit0.feval(R),R,'-','DisplayName','Jackknife frequentist fits','linewidth',2,'color',JACK)
    else
        plot(fit0.feval(R),R,'-','linewidth',2,'color',JACK,'HandleVisibility','off')
    end
end

MIC = abs(fit2.Coefficients.Estimate(1));

adjR2 = fit2.Rsquared.Adjusted;
plot(fit2.feval(R),R,'-r','DisplayName',['Radical exp frequentist fit (adj R^2\approx',num2str(adjR2,3),')'],'linewidth',3)

semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
text(20,3,['MIC\approx',num2str(MIC,3),'\mug/mL']);

%semilogx(conc(1:2),zoi(1:2),'ok','markersize',14,'DisplayName','excised data','linewidth',1)
set(gca,'Ytick',0:1:max(R))

W1 = conc(end);
W2 = conc(end-1);
plot([W1,W2],[0,0],'-k','linewidth',6,'DisplayName','MIC ground truth (W)');
text(0.8,0.5,'W','FontSize',22);

MMSEMmicroNisinLogF = MCMCupdatePointEstimates(conc(1:end-2),zoi(1:end-2),...
    expRad,fit2.Coefficients.Estimate,weights(conc(1:end-2)));

ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','northwest')
axis tight
xlim([0 150])
%ylim([0 max(R)])
ylim([0 13])

%% data from literature

[conc,zoi,R] = defineSarcinaCloxData();

figure(2)
disp('Figure 2C')

semilogx(conc,zoi,'.k','markersize',38,'DisplayName','ZoI data')
hold on

fit = fitnlm(zoi,conc,expRad,[1.85 0.04 0.09],'Weights',weights(conc))

set(gca,'Ytick',0:5:max(R))

for j = 1:5
    %Jackknife regressions:
    Z = myRemoveDatum(zoi,zoi(j));
    C = myRemoveDatum(conc,conc(j));

    fitJ = fitnlm(Z,C,expRad,fit.Coefficients.Estimate,'Weights',weights(C));   
    if j == 1
        plot(fitJ.feval(R),R,'-','DisplayName','Jackknife frequentist fits','linewidth',2,'color',JACK)
    else
        plot(fitJ.feval(R),R,'-','linewidth',2,'color',JACK,'HandleVisibility','off')
    end
end

adjR2 = fit.Rsquared.Adjusted;
plot(fit.feval(R),R,'-r','DisplayName',['Radical Exp frequentist fit (adj R^2\approx',num2str(adjR2,3),')'],'linewidth',3)

semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
MIC = abs(fit.Coefficients.Estimate(1));
text(20,5,['MIC\approx',num2str(MIC,3),'\mug/mL']);

MMSEMsarcinaCloxLogF = MCMCupdatePointEstimates(conc,zoi,expRad,...
    fit.Coefficients.Estimate,weights(conc));

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
    exportgraphics(gcf,'./figures/Micrococcus_luteus_radExp.pdf')
    figure(2)
    exportgraphics(gcf,'./figures/Bacillus_subtilis_radExp.pdf')
end