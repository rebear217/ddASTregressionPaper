%clear all
close all
clc

JACKcolour = [0.2 0.5 0.2];
fexpint = @expintFunction;

%%  data from literature

%Micrococcus luteus (ATCC 10240) and nisin
%data taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC4576963/pdf/fsn30003-0394.pdf

[conc,zoi,R] = defineMicrococcusNisinData();

pGuess = [1 0.02 0.5];
weights = @(yhat) 1./(abs(yhat).^2);
%weights = @(yhat) ones(size(yhat));
%weights = @(yhat) 1./(1 + abs(yhat).^3);

figure(1)
disp('Figure 3B')

semilogx(conc,zoi,'.k','markersize',34,'DisplayName','ZoI data')
hold on

fitAll = fitnlm(zoi(1:end-2),conc(1:end-2),fexpint,pGuess,'Weights',weights(conc(1:end-2)))
%fit1 = fitnlm(zoi,conc,fexpint,[0.5 0.02 1],'Weights',weights)
%fit2 = fitnlm(zoi,conc,fexpint,[0.5 0.02 1],'Weights',weights)

Nz = length(zoi(1:end-2));
for j = 1:Nz
    %Jackknife regressions:
    Z = myRemoveDatum(zoi(1:end-2),zoi(j));
    C = myRemoveDatum(conc(1:end-2),conc(j));

    fit0 = fitnlm(Z,C,fexpint,fitAll.Coefficients.Estimate,'Weights',weights(C));   
    if j == 1
        plot(fit0.feval(R),R,'-','DisplayName','Jackknife frequentist fits','linewidth',2,'color',JACKcolour)
    else
        plot(fit0.feval(R),R,'-','linewidth',2,'color',JACKcolour,'HandleVisibility','off')
    end
end

MIC = abs(fitAll.Coefficients.Estimate(1));
adjR2 = fitAll.Rsquared.Adjusted;

%plot(fitAll.feval(R),R,'-b','DisplayName','fit: all data','linewidth',1)
plot(fitAll.feval(R),R,'-b','DisplayName',['ExpInt frequentist fit (adj R^2\approx',num2str(adjR2,3),')'],'linewidth',3)
semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
text(20,3,['MIC\approx',num2str(MIC,3),'\mug/mL']);

%semilogx(conc(1:2),zoi(1:2),'ok','markersize',14,'DisplayName','excised data','linewidth',1)

W1 = conc(end);
W2 = conc(end-1);
plot([W1,W2],[0,0],'-k','linewidth',6,'DisplayName','MIC ground truth (W)');
text(0.8,0.5,'W','FontSize',22);

MMSEMmicroNisinExpInt = MCMCupdatePointEstimates(conc(1:end-2),zoi(1:end-2),fexpint,...
    fitAll.Coefficients.Estimate,weights(conc(1:end-2)));

ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','northwest')
xlim([0 150])
ylim([0 max(R)])
set(gca,'Ytick',0:1:max(R))

%% data from literature

[conc,zoi,R] = defineSarcinaCloxData();

fit = fitnlm(zoi,conc,fexpint,pGuess,'Weights',weights(conc))

figure(2)
disp('Figure 2B')

semilogx(conc,zoi,'.k','markersize',38,'DisplayName','ZoI data')
hold on
set(gca,'Ytick',0:5:max(R))

for j = 1:length(zoi)
    %Jackknife regressions:
    Z = myRemoveDatum(zoi,zoi(j));
    C = myRemoveDatum(conc,conc(j));
    fitJ = fitnlm(Z,C,fexpint,fit.Coefficients.Estimate,'Weights',weights(C));   
    if j == 1
        plot(fitJ.feval(R),R,'-','DisplayName','Jackknife frequentist fits','linewidth',2,'color',JACKcolour)
    else
        plot(fitJ.feval(R),R,'-','linewidth',2,'color',JACKcolour,'HandleVisibility','off')
    end
end

adjR2 = fit.Rsquared.Adjusted;
plot(fit.feval(R),R,'-b','DisplayName',['ExpInt frequentist fit (adj R^2\approx',num2str(adjR2,3),')'],'linewidth',3)

semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
MIC = abs(fit.Coefficients.Estimate(1));
text(20,5,['MIC\approx',num2str(MIC,3),'\mug/mL']);

MMSEMsarcinaCloxExpInt = MCMCupdatePointEstimates(conc,zoi,fexpint,...
    fit.Coefficients.Estimate,weights(conc));

ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','northwest')
axis tight
xlim([0.5 100])
ylim([0 max(R)])

%%

plotON = 1;
if plotON
    figure(1)
    exportgraphics(gcf,'./figures/Micrococcus_luteus_expint.pdf')
    figure(2)
    exportgraphics(gcf,'./figures/Bacillus_subtilis_expint.pdf')
end

