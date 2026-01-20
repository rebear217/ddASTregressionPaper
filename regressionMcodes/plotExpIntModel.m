clear all
close all
clc

%%  data from literature

%Micrococcus luteus (ATCC 10240) and nisin
%data taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC4576963/pdf/fsn30003-0394.pdf
%Figure 7 (see legend of A-J)

conc = [125 75 50 25 18.75 12.5 6.25 2.5 1.25 0.625];
zoi = [10.7 10.5 10.3 9.5 9 8 6.5 4 0.0001 0.000001];

I = @(x)exp(-x)./x;
R = 0:0.01:12;

myGamma = @(a)integral(I,a,inf);
c0 = 0.5;
c1 = 0.02;

Esub = @(r)(exp(-r).*log(1+1./r));
Esuper = @(r)(0.5*exp(-r).*log(1+2./r));

fexpint = @(b,r)(abs(b(3)) + b(1) ./ expint(abs(b(2))*r.^2));
fsuper = @(b,r)(abs(b(3)) + b(1) ./ Esuper(abs(b(2))*r.^2));
fsub = @(b,r)(abs(b(3)) + b(1) ./ Esub(abs(b(2))*r.^2));

%fexpint = @(b,r)abs(b(3)) + b(1).*(log(abs(r./b(2)))).^(1/2)

figure(1)

semilogx(conc,zoi,'.k','markersize',34,'DisplayName','ZoI data')
hold on

weights = @(yhat) 1./(1 + abs(yhat).^2);
%weights = @(yhat) ones(size(yhat));

fitAll = fitnlm(zoi(1:end-2),conc(1:end-2),fexpint,[0.5 0.02 1],'Weights',weights)
fit1 = fitnlm(zoi(2:end-2),conc(2:end-2),fexpint,[0.5 0.02 1],'Weights',weights)
fit2 = fitnlm(zoi(3:end-2),conc(3:end-2),fexpint,[0.5 0.02 1],'Weights',weights)

M = 1;

for j = M:9
    Z = setdiff(zoi(M:end-1),zoi(j));
    C = setdiff(conc(M:end-1),conc(j));
    fit0 = fitnlm(Z,C,fexpint,fitAll.Coefficients.Estimate,'Weights',weights);   
    if j == M
        plot(fit0.feval(R),R,'-','DisplayName','Jacknife fits','linewidth',2,'color',[1 1 1]*0.7)
    else
        plot(fit0.feval(R),R,'-','linewidth',2,'color',[1 1 1]*0.7,'HandleVisibility','off')
    end
end

MIC = abs(fit2.Coefficients.Estimate(3));
adjR2 = fit2.Rsquared.Adjusted;

%plot(fitAll.feval(R),R,'-b','DisplayName','fit: all data','linewidth',1)
plot(fit2.feval(R),R,'-b','DisplayName',['ExpInt fit (adj R^2\approx',num2str(adjR2),')'],'linewidth',3)
semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
text(MIC,0.4,num2str(MIC));

%semilogx(conc(1:2),zoi(1:2),'ok','markersize',14,'DisplayName','excised data','linewidth',1)

ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','southeast')
xlim([0 150])
ylim([0 max(R)])
set(gca,'Ytick',0:1:max(R))

%% data from literature

%data taken from
%https://pmc.ncbi.nlm.nih.gov/articles/PMC546645/pdf/applmicro00228-0032.pdf
%PAP is unclear (not stated in paper): Bacillus subtilis or Sarcina lutea?
%(Table 2)

conc = [3 6 12 24 48];
zoi = [11.68 16.62 20.96 24.86 28.21];
R = 0:0.01:32;

fit = fitnlm(zoi,conc,fexpint,[0.5 0.02 1],'Weights',weights)

figure(2)

semilogx(conc,zoi,'.k','markersize',38,'DisplayName','ZoI data')
hold on
set(gca,'Ytick',0:5:max(R))

for j = 1:5
    Z = setdiff(zoi,zoi(j));
    C = setdiff(conc,conc(j));
    fitJ = fitnlm(Z,C,fexpint,fit.Coefficients.Estimate,'Weights',weights);   
    if j == 1
        plot(fitJ.feval(R),R,'-','DisplayName','Jacknife fits','linewidth',2,'color',[1 1 1]*0.7)
    else
        plot(fitJ.feval(R),R,'-','linewidth',2,'color',[1 1 1]*0.7,'HandleVisibility','off')
    end
end

adjR2 = fit.Rsquared.Adjusted;
plot(fit.feval(R),R,'-b','DisplayName',['ExpInt fit (adj R^2\approx',num2str(adjR2),')'],'linewidth',3)

semilogx(conc,zoi,'.k','markersize',34,'HandleVisibility','off')
MIC = abs(fit.Coefficients.Estimate(3));
text(MIC,0.8,num2str(MIC));

ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','southeast')
axis tight
xlim([0.5 100])
ylim([0 max(R)])

%%

figure(1)
exportgraphics(gcf,'./figures/Micrococcus_luteus_expint.pdf')
figure(2)
exportgraphics(gcf,'./figures/Bacillus_subtilis_expint.pdf')
