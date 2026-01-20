clear all
close all
clc

gr = [1,1,1]/2;

%%  data from literature

correctLo = @(m,s) m - 1.96*abs(s);
correctHi = @(m,s) m + 1.96*abs(s);

%Micrococcus luteus (ATCC 10240) and nisin
%data taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC4576963/pdf/fsn30003-0394.pdf
%Figure 7 (see legend of A-J)

conc = [125 75 50 25 18.75 12.5 6.25 2.5 1.25 0.625];
zoi = [10.7 10.5 10.3 9.5 9 8 6.5 4 0.0001 0.000001];

R = 0:0.01:12;

%old version:
%logF = @(p,r)p(3) + p(1)*exp((-1 + (1+p(2)*r.^2).^(1/2))/2 + ...
%        log(((-1 + (1+p(2)*r.^2).^(1/2))))/2 + ...
%        (p(2)/2)*r.^2./(-1 + (1+p(2)*r.^2).^(1/2)));
%new version:
logF = @(p,r)p(3) + exp(p(1) + sqrt(1+p(2).*r.^2 )) .* ( sqrt(1+p(2)*r.^2) );

b0F = [0.03 0.4 1.75];

fexpint = @(b,r)(abs(b(3)) + b(1) ./ expint(abs(b(2))*r.^2));
b0ei = [0.5 0.02 1];

bonevF = @(p,r)p(1).*exp(-r.^2*p(2));
b0b = [1.3 -0.03];

weights = @(yhat) 1./(1 + abs(yhat).^2);

S = 3;
fitF = fitnlm(zoi(S:end-2),conc(S:end-2),logF,b0F,'Weights',weights)
fitEI = fitnlm(zoi(S:end-2),conc(S:end-2),fexpint,b0ei,'Weights',weights)
fitB = fitnlm(zoi(S:end-2),conc(S:end-2),bonevF,b0b,'Weights',weights)
fitLM = fitlm(zoi(S:end-2),log(conc(S:end-2)),'Weights',weights(log(conc(S:end-2))))

figure(1)

semilogx(conc(:),zoi(:),'.k','markersize',36,'DisplayName','ZoI data')
hold on
plot(fitB.feval(R),R,'-','linewidth',2,'color','k','DisplayName',['Bonev (adj R^2\approx',num2str(fitB.Rsquared.Adjusted),')'])
plot(fitEI.feval(R),R,'-','linewidth',2,'color','b','DisplayName',['ExpInt (adj R^2\approx',num2str(fitEI.Rsquared.Adjusted),')'])
plot(fitF.feval(R),R,'-','linewidth',2,'color','r','DisplayName',['Radical Exp (adj R^2\approx',num2str(fitF.Rsquared.Adjusted),')'])
plot(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',gr,'DisplayName',['Log-linear (adj R^2\approx',num2str(fitLM.Rsquared.Adjusted),')'])

L1 = correctLo(fitB.Coefficients.Estimate(1),fitB.Coefficients.SE(1));
L2 = correctHi(fitB.Coefficients.Estimate(1),fitB.Coefficients.SE(1));
plot([L1,L2],[7.25,7.25],'-','linewidth',6,'color','k','DisplayName','Bonev MIC 95% CI')

L1 = correctLo(fitEI.Coefficients.Estimate(3),fitEI.Coefficients.SE(3));
L2 = correctHi(fitEI.Coefficients.Estimate(3),fitEI.Coefficients.SE(3));
plot([max(1e-5,L1),L2],[7,7],'-','linewidth',6,'color','b','DisplayName','ExpInt MIC 95% CI')

L1 = correctLo(fitF.Coefficients.Estimate(3),fitF.Coefficients.SE(3));
L2 = correctHi(fitF.Coefficients.Estimate(3),fitF.Coefficients.SE(3));
plot([L1,L2],[6.75,6.75],'-','linewidth',6,'color','r','DisplayName','Radical Exp MIC 95% CI')

L1 = correctLo(fitLM.Coefficients.Estimate(1),fitLM.Coefficients.SE(1));
L2 = correctHi(fitLM.Coefficients.Estimate(1),fitLM.Coefficients.SE(1));
plot(exp([L1,L2]),[6.5,6.5],'-','linewidth',6,'color',gr,'DisplayName','Log-linear MIC 95% CI')

text(exp((L1+L2)/2),0.5,num2str(exp(fitLM.Coefficients.Estimate(1))))

set(gca,'Ytick',0:1:max(R))
ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','southeast')
axis tight
xL = xlim;
xlim([0.25 xL(2)])
exportgraphics(gcf,'./figures/3ModelCompare1.pdf')

%%

conc = [3 6 12 24 48];
zoi = [11.68 16.62 20.96 24.86 28.21];
R = 0:0.01:32;

fitF = fitnlm(zoi,conc,logF,[0.09 0.04 1.85],'Weights',weights)
fitEI = fitnlm(zoi,conc,fexpint,[0.5 0.02 1],'Weights',weights)
fitB = fitnlm(zoi,conc,bonevF,[1 0],'Weights',weights)
fitLM = fitlm(zoi,log(conc),'Weights',weights(log(conc)))

figure(2)
semilogx(conc(:),zoi(:),'.k','markersize',36,'DisplayName','ZoI data')
hold on
plot(fitB.feval(R),R,'-','linewidth',2,'color','k','DisplayName',['Bonev (adj R^2\approx',num2str(fitB.Rsquared.Adjusted),')'])
plot(fitEI.feval(R),R,'-','linewidth',2,'color','b','DisplayName',['ExpInt (adj R^2\approx',num2str(fitEI.Rsquared.Adjusted),')'])
plot(fitF.feval(R),R,'-','linewidth',2,'color','r','DisplayName',['Radical Exp (adj R^2\approx',num2str(fitF.Rsquared.Adjusted),')'])
plot(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',gr,'DisplayName',['Log-linear (adj R^2\approx',num2str(fitLM.Rsquared.Adjusted),')'])

L1 = correctLo(fitB.Coefficients.Estimate(1),fitB.Coefficients.SE(1));
L2 = correctHi(fitB.Coefficients.Estimate(1),fitB.Coefficients.SE(1));
plot([L1,L2],[15.5,15.5],'-','linewidth',6,'color','k','DisplayName','Bonev MIC 95% CI')

L1 = correctLo(fitEI.Coefficients.Estimate(3),fitEI.Coefficients.SE(3));
L2 = correctHi(fitEI.Coefficients.Estimate(3),fitEI.Coefficients.SE(3));
plot([L1,L2],[15,15],'-','linewidth',6,'color','b','DisplayName','ExpInt MIC 95% CI')

L1 = correctLo(fitF.Coefficients.Estimate(3),fitF.Coefficients.SE(3));
L2 = correctHi(fitF.Coefficients.Estimate(3),fitF.Coefficients.SE(3));
plot([L1,L2],[14.5,14.5],'-','linewidth',6,'color','r','DisplayName','Radical Exp MIC 95% CI')

L1 = correctLo(fitLM.Coefficients.Estimate(1),fitLM.Coefficients.SE(1));
L2 = correctHi(fitLM.Coefficients.Estimate(1),fitLM.Coefficients.SE(1));
plot(exp([L1,L2]),[14,14],'-','linewidth',6,'color',gr,'DisplayName','Log-linear MIC 95% CI')

text(exp((L1+L2)/2),1,num2str(exp(fitLM.Coefficients.Estimate(1))))

%set(gca,'Ytick',0:1:max(R))
ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
legend('Location','southeast')
axis tight

exportgraphics(gcf,'./figures/3ModelCompare2.pdf')
