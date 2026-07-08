%clear all
close all
clc

gr = [1,1,1]/2;

correctLo = @(m,s) m - 1.96*abs(s);
correctHi = @(m,s) m + 1.96*abs(s);

%%  data from literature
[conc,zoi,R] = defineMicrococcusNisinData();

% radical old version (an exp & log radical function):
% logF = @(p,r)p(3) + p(1)*exp((-1 + (1+p(2)*r.^2).^(1/2))/2 + ...
%        log(((-1 + (1+p(2)*r.^2).^(1/2))))/2 + ...
%        (p(2)/2)*r.^2./(-1 + (1+p(2)*r.^2).^(1/2)));

% radical new version:
logF = @(p,r)p(1) + exp(p(3) + sqrt(1+abs(p(2)).*r.^2 )) .* ( -1 + sqrt(1+abs(p(2))*r.^2) );
b0F = [1.8 0.024 -2];

% expint function
fexpint = @(b,r)(b(1) + b(3) ./ expint(abs(b(2))*r.^2));
b0ei = [1 0.02 2];

% bonev function
bonevF = @(p,r)p(1).*exp(-r.^2*p(2));
b0b = [2 -0.03];

weights = @(yhat) 1./(abs(yhat).^2);

S = 1;
fitF = fitnlm(zoi(S:end-2),conc(S:end-2),logF,b0F,'Weights',weights)
fitEI = fitnlm(zoi(S:end-2),conc(S:end-2),fexpint,b0ei,'Weights',weights)
fitB = fitnlm(zoi(S:end-2),conc(S:end-2),bonevF,b0b,'Weights',weights)
fitLM = fitlm(zoi(S:end-2),log(conc(S:end-2)),'Weights',weights(conc(S:end-2)))

figure(1)
set(1,'pos',[543   721   696   461])

semilogx(conc(:),zoi(:),'.k','markersize',36,'DisplayName','ZoI data')
hold on

YS = 4;
L1 = correctLo(fitB.Coefficients.Estimate(1),fitB.Coefficients.SE(1));
L2 = correctHi(fitB.Coefficients.Estimate(1),fitB.Coefficients.SE(1));
plot([L1,L2],YS + [7.25,7.25],'-','linewidth',6,'color','k','DisplayName','Bonev MIC 95% CI')

L1 = correctLo(fitEI.Coefficients.Estimate(1),fitEI.Coefficients.SE(1));
L2 = correctHi(fitEI.Coefficients.Estimate(1),fitEI.Coefficients.SE(1));
plot([max(1e-5,L1),L2],YS + [7,7],'-','linewidth',6,'color','b','DisplayName','ExpInt MIC 95% CI')

L1 = correctLo(fitF.Coefficients.Estimate(1),fitF.Coefficients.SE(1));
L2 = correctHi(fitF.Coefficients.Estimate(1),fitF.Coefficients.SE(1));
plot([L1,L2],YS + [6.75,6.75],'-','linewidth',6,'color','r','DisplayName','Radical Exp MIC 95% CI')

L1 = correctLo(fitLM.Coefficients.Estimate(1),fitLM.Coefficients.SE(1));
L2 = correctHi(fitLM.Coefficients.Estimate(1),fitLM.Coefficients.SE(1));
plot(exp([L1,L2]),YS + [6.5,6.5],'-','linewidth',6,'color',gr,'DisplayName','Log-linear MIC 95% CI')

text(-0.3 + exp((L1+L2)/2),-0.5,['LLR-MIC\approx',num2str(exp(fitLM.Coefficients.Estimate(1)),2)]);

set(gca,'Ytick',0:1:max(R))
ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')

semilogx(MMSEMmicroNisinLogF.Xdose,MMSEMmicroNisinLogF.Yresponse,'--','LineWidth',1,...
    'color','r','DisplayName','Radical Exp MMSEM')
semilogx(MMSEMmicroNisinBonevF.Xdose,MMSEMmicroNisinBonevF.Yresponse,'--','LineWidth',1,...
    'color','k','DisplayName','Bonev MMSEM')
semilogx(MMSEMmicroNisinExpInt.Xdose,MMSEMmicroNisinExpInt.Yresponse,'--','LineWidth',1,...
    'color','b','DisplayName','ExpInt MMSEM')

%plot(fitB.feval(R),R,'-','linewidth',2,'color','k','DisplayName',['Bonev (adj R^2\approx',num2str(fitB.Rsquared.Adjusted,3),')'])
%plot(fitEI.feval(R),R,'-','linewidth',2,'color','b','DisplayName',['ExpInt (adj R^2\approx',num2str(fitEI.Rsquared.Adjusted,3),')'])
%plot(fitF.feval(R),R,'-','linewidth',2,'color','r','DisplayName',['Radical Exp (adj R^2\approx',num2str(fitF.Rsquared.Adjusted,3),')'])
%plot(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',gr,'DisplayName',['Log-linear (adj R^2\approx',num2str(fitLM.Rsquared.Adjusted,3),')'])

plot(fitB.feval(R),R,'-','linewidth',2,'color','k','DisplayName','Bonev frequentist fit')
plot(fitEI.feval(R),R,'-','linewidth',2,'color','b','DisplayName','ExpInt frequentist fit')
plot(fitF.feval(R),R,'-','linewidth',2,'color','r','DisplayName','Radical Exp frequentist fit')
plot(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',gr,'DisplayName','Log-linear')

axis tight
xL = xlim;
xlim([0.185 101])
%legend('Location','northeastoutside')
legend('Location','southeast')
legend('boxoff')

ylim([0 12])

exportgraphics(gcf,'./figures/3ModelCompare1.pdf')

%%

[conc,zoi,R] = defineSarcinaCloxData();

fitF = fitnlm(zoi,conc,logF,[0.09 0.04 1.85],'Weights',weights)
fitEI = fitnlm(zoi,conc,fexpint,[0.5 0.02 1],'Weights',weights)
fitB = fitnlm(zoi,conc,bonevF,[1 0],'Weights',weights)
fitLM = fitlm(zoi,log(conc),'Weights',weights((conc)))

figure(2)
set(2,'pos',[543   761   643   421])

semilogx(conc(:),zoi(:),'.k','markersize',36,'DisplayName','ZoI data')
hold on
%plot(fitB.feval(R),R,'-','linewidth',2,'color','k','DisplayName',['Bonev (adj R^2\approx',num2str(fitB.Rsquared.Adjusted,3),')'])
%plot(fitEI.feval(R),R,'-','linewidth',2,'color','b','DisplayName',['ExpInt (adj R^2\approx',num2str(fitEI.Rsquared.Adjusted,3),')'])
%plot(fitF.feval(R),R,'-','linewidth',2,'color','r','DisplayName',['Radical Exp (adj R^2\approx',num2str(fitF.Rsquared.Adjusted,3),')'])
%plot(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',gr,'DisplayName',['Log-linear (adj R^2\approx',num2str(fitLM.Rsquared.Adjusted,3),')'])

YS = 12;
L1 = correctLo(fitB.Coefficients.Estimate(1),fitB.Coefficients.SE(1));
L2 = correctHi(fitB.Coefficients.Estimate(1),fitB.Coefficients.SE(1));
plot([L1,L2],YS+[15.5,15.5],'-','linewidth',6,'color','k','DisplayName','Bonev MIC 95% CI')

L1 = correctLo(fitEI.Coefficients.Estimate(1),fitEI.Coefficients.SE(1));
L2 = correctHi(fitEI.Coefficients.Estimate(1),fitEI.Coefficients.SE(1));
plot([L1,L2],YS+[15,15],'-','linewidth',6,'color','b','DisplayName','ExpInt MIC 95% CI')

L1 = correctLo(fitF.Coefficients.Estimate(1),fitF.Coefficients.SE(1));
L2 = correctHi(fitF.Coefficients.Estimate(1),fitF.Coefficients.SE(1));
plot([L1,L2],YS+[14.5,14.5],'-','linewidth',6,'color','r','DisplayName','Radical Exp MIC 95% CI')

L1 = correctLo(fitLM.Coefficients.Estimate(1),fitLM.Coefficients.SE(1));
L2 = correctHi(fitLM.Coefficients.Estimate(1),fitLM.Coefficients.SE(1));
plot(exp([L1,L2]),YS+[14,14],'-','linewidth',6,'color',gr,'DisplayName','Log-linear MIC 95% CI')

text(-0.3+exp((L1+L2)/2),-0.95,['LLR-MIC\approx',num2str(exp(fitLM.Coefficients.Estimate(1)),2)]);

semilogx(MMSEMsarcinaCloxLogF.Xdose,MMSEMsarcinaCloxLogF.Yresponse,'--','LineWidth',1,...
    'color','r','DisplayName','Radical Exp MMSEM')
semilogx(MMSEMsarcinaCloxBonevF.Xdose,MMSEMsarcinaCloxBonevF.Yresponse,'--','LineWidth',1,...
    'color','k','DisplayName','Bonev MMSEM')
semilogx(MMSEMsarcinaCloxExpInt.Xdose,MMSEMsarcinaCloxExpInt.Yresponse,'--','LineWidth',1,...
    'color','b','DisplayName','ExpInt MMSEM')

plot(fitB.feval(R),R,'-','linewidth',2,'color','k','DisplayName','Bonev frequentist fit')
plot(fitEI.feval(R),R,'-','linewidth',2,'color','b','DisplayName','ExpInt frequentist fit')
plot(fitF.feval(R),R,'-','linewidth',2,'color','r','DisplayName','Radical Exp frequentist fit')
plot(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',gr,'DisplayName','Log-linear')


%set(gca,'Ytick',0:1:max(R))
ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')

axis tight
xL = xlim;
legend('Location','southeast')
legend('boxoff')

xlim([0.3 220])
ylim([0 30])

exportgraphics(gcf,'./figures/3ModelCompare2.pdf')
