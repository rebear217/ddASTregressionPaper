%clear all
close all
clc

JACK = [0.2 0.5 0.2];

%%  data from literature

%Micrococcus luteus (ATCC 10240) and nisin
%data taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC4576963/pdf/fsn30003-0394.pdf

[conc,zoi,R] = defineMicrococcusNisinData();

I = @(x)exp(-x)./x;

myGamma = @(a)integral(I,a,inf);
%c0 = 0.5; c1 = 0.02;

Esub = @(r)(exp(-r).*log(1+1./r));
Esuper = @(r)(0.5*exp(-r).*log(1+2./r));

fexpint = @(b,r)(abs(b(1)) + b(3) ./ expint(abs(b(2))*r.^2));
fsuper = @(b,r)(abs(b(1)) + b(3) ./ Esuper(abs(b(2))*r.^2));
fsub = @(b,r)(abs(b(1)) + b(3) ./ Esub(abs(b(2))*r.^2));
pGuess = [1 0.02 0.5];

%fexpint = @(b,r)abs(b(3)) + b(1).*(log(abs(r./b(2)))).^(1/2)

weights = @(yhat) 1./(abs(yhat).^2);
%weights = @(yhat) ones(size(yhat));
%weights = @(yhat) 1./(1 + abs(yhat).^3);

figure(1)
disp('Figure 3B')

semilogx(conc,zoi,'.k','markersize',34,'DisplayName','ZoI data')
hold on

fitAll = fitnlm(zoi(1:end-2),conc(1:end-2),fexpint,pGuess,'Weights',weights)
%fit1 = fitnlm(zoi,conc,fexpint,[0.5 0.02 1],'Weights',weights)
%fit2 = fitnlm(zoi,conc,fexpint,[0.5 0.02 1],'Weights',weights)
M = 1;
for j = M:9
    Z = setdiff(zoi(M:end-2),zoi(j));
    C = setdiff(conc(M:end-2),conc(j));
    fit0 = fitnlm(Z,C,fexpint,fitAll.Coefficients.Estimate,'Weights',weights);   
    if j == M
        plot(fit0.feval(R),R,'-','DisplayName','Jacknife frequentist fits','linewidth',2,'color',JACK)
    else
        plot(fit0.feval(R),R,'-','linewidth',2,'color',JACK,'HandleVisibility','off')
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

MMSEMmicroNisinExpInt = MCMCupdatePointEstimates(conc(1:end-2),zoi(1:end-2),fexpint,fitAll.Coefficients.Estimate,weights);

ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','northwest')
xlim([0 150])
ylim([0 max(R)])
set(gca,'Ytick',0:1:max(R))

%% data from literature

[conc,zoi,R] = defineSarcinaCloxData();

fit = fitnlm(zoi,conc,fexpint,pGuess,'Weights',weights)

figure(2)
disp('Figure 2B')

semilogx(conc,zoi,'.k','markersize',38,'DisplayName','ZoI data')
hold on
set(gca,'Ytick',0:5:max(R))

for j = 1:5
    Z = setdiff(zoi,zoi(j));
    C = setdiff(conc,conc(j));
    fitJ = fitnlm(Z,C,fexpint,fit.Coefficients.Estimate,'Weights',weights);   
    if j == 1
        plot(fitJ.feval(R),R,'-','DisplayName','Jacknife frequentist fits','linewidth',2,'color',JACK)
    else
        plot(fitJ.feval(R),R,'-','linewidth',2,'color',JACK,'HandleVisibility','off')
    end
end

adjR2 = fit.Rsquared.Adjusted;
plot(fit.feval(R),R,'-b','DisplayName',['ExpInt frequentist fit (adj R^2\approx',num2str(adjR2,3),')'],'linewidth',3)

semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
MIC = abs(fit.Coefficients.Estimate(1));
text(20,5,['MIC\approx',num2str(MIC,3),'\mug/mL']);

MMSEMsarcinaCloxExpInt = MCMCupdatePointEstimates(conc,zoi,fexpint,fit.Coefficients.Estimate,weights);

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

