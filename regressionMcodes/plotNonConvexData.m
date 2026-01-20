clear all
close all
clc

gr = [1,1,1]/2;

%%

clc

%dataset 1 from https://europepmc.org/article/PMC/4576963, Table 1
%dataset 3 from https://pmc.ncbi.nlm.nih.gov/articles/PMC5761473, Table 3
dataLabels = {'ZoI data','ZoI data','ZoI data'};

concs = {[0.625 , 0.312 , 0.156, 0.08, 0.04],[0.625 , 0.312 , 0.156],[2.56, 3.2, 4.0, 5.0 6.25]};
zois = {[3.44 , 2.4 , 1.65, 1.2, 0.35],[3.44 , 2.6 , 1.85],[19.6, 21.8, 22.9, 24.5, 26.7]};
stes = {[1 , 1 ,1, 1, 1]/10,[0.08 , 0.12 , 0.19],[1,1,1,1,1]/10};
xlims = {[0.01 1.1],[0.01 1.1],[0.1 10]};

b0bons = {[0.093048      -0.1587],[0.088749     -0.16473],[0.88613   -0.0028004],[0.74632    -0.003168]};
b0expints = {[0.21149     0.062875  -2.3977e-12],[0.12871     0.087392  -3.2225e-13],[1.1478    0.0015803   7.4595e-11],[0.27383    0.0029442     -0.98827]};
b0logs = {[-1.3212     0.061806     -0.70121],[-3.9648       0.4753    -0.020914],[-2.3013     0.011425      0.18092],[-4.8894     0.034636        1.287]};

Rs = {0:0.001:4,0:0.001:4,0:0.001:30};
ylabels = {'zone of inhibition radius (mm)','zone of inhibition radius (mm)','zone of inhibition diameter (mm)'};
xlabels = {'erythromycin concentration (\mug/mL)','nisin concentration (\mug/mL)','Levofloxacin concentration (\mug/mL)'};
N = 3;

%logF = @(p,r)abs(p(3)) + p(1)*exp((-1 + (1+p(2)*r.^2).^(1/2))/2 + ...
%        log(((-1 + (1+p(2)*r.^2).^(1/2))))/2 + ...
%        (p(2)/2)*r.^2./(-1 + (1+p(2)*r.^2).^(1/2)));
logF = @(p,r)p(3) + exp(p(1) + sqrt(1+p(2).*r.^2 )) .* ( sqrt(1+p(2)*r.^2) );
%logF = @(p,r)p(3) + exp(p(1) + p(2).*r) .* p(2).*r;

expintF = @(b,r)(abs(b(3)) + b(1) ./ expint(abs(b(2))*r.^2));
bonevF = @(p,r)p(1).*exp(-r.^2*p(2));

weights = @(yhat) 1./(1 + abs(yhat).^2);
%weights = @(yhat) ones(size(yhat));

opt = statset('fitnlm');
opt.MaxIter = 20000;
opt.TolFun = 1e-10;
opt.TolX = 1e-10;

correctLo = @(m,s) m - 1.96*s;
correctHi = @(m,s) m + 1.96*s;

for D = 1:4
    
    d = D;
    if D == 4
        d = 3;
    end

    conc = concs{d};
    zoi = zois{d};
    ste = stes{d};

    if D >= 4
        conc = conc(1:end-1);
        zoi = zoi(1:end-1);
        ste = ste(1:end-1);
    end

    thisxlim = xlims{d};
    b0log = b0logs{D};
    b0expint = b0expints{D};
    b0bon = b0bons{D};
    R = Rs{d};
    
    if N > 1
        conc = repmat(conc,N,1);
        ste = ones(N,1)*ste;
        ste = ste .* randn(size(ste));
        zoi = repmat(zoi,N,1) + ste;
    end

    zoi = zoi(:);
    conc = conc(:);

    figure(D)

    fitBon = fitnlm(zoi,conc,bonevF,b0bon,'Weights',weights,'Options',opt)
    fitBon.Coefficients.Estimate'

    fitEI = fitnlm(zoi,conc,expintF,b0expint,'Weights',weights,'Options',opt)
    fitEI.Coefficients.Estimate'

    fitlog = fitnlm(zoi,conc,logF,b0log,'Weights',weights,'Options',opt)
    fitlog.Coefficients.Estimate'

    fitLM = fitlm(zoi,log(conc),'Weights',weights(log(conc)));
    
    semilogx(fitEI.feval(R),R,'-','linewidth',2,'color','b','DisplayName',['ExpInt (adj R^2\approx',num2str(fitEI.Rsquared.Adjusted),')'])
    hold on
    semilogx(fitlog.feval(R),R,'-','linewidth',2,'color','r','DisplayName',['Radical Exp (adj R^2\approx',num2str(fitlog.Rsquared.Adjusted),')'])
    semilogx(fitBon.feval(R),R,'-','linewidth',2,'color','k','DisplayName',['Bonev (adj R^2\approx',num2str(fitBon.Rsquared.Adjusted),')'])
    semilogx(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',gr,'DisplayName',['Log-linear (adj R^2\approx',num2str(fitLM.Rsquared.Adjusted),')'])

    %{
    L1 = fitlog.Coefficients.Estimate(3) - 2*fitlog.Coefficients.SE(3);
    L2 = fitlog.Coefficients.Estimate(3) + 2*fitlog.Coefficients.SE(3);
    plot([L1,L2],[0,0],'-','linewidth',4,'color','r')
    %}

    L1 = fitBon.Coefficients.Estimate(1) - 1.96*fitBon.Coefficients.SE(1);
    L2 = fitBon.Coefficients.Estimate(1) + 1.96*fitBon.Coefficients.SE(1);
    %plot([L1,L2],[0,0],'-','linewidth',6,'color','k','DisplayName','Bonev MIC 95% CI')

    L1 = correctLo(fitLM.Coefficients.Estimate(1),fitlog.Coefficients.SE(1));
    L2 = correctHi(fitLM.Coefficients.Estimate(1),fitlog.Coefficients.SE(1));
    %plot(exp([L1,L2]),[-0.2,-0.2],'-','linewidth',6,'color',gr,'DisplayName','Log-linear MIC 95% CI')

    semilogx(conc,zoi,'.k','markersize',20,'DisplayName',dataLabels{d},'color',[1 1 1]*0.7)
    semilogx(concs{d},zois{d},'.k','markersize',44,'DisplayName',[dataLabels{d},' replicate means'])

    if D == 4
        c = concs{d};
        z = zois{d};

        semilogx(c(end),z(end),'or','markersize',18,'DisplayName','excised datapoint')
    end

    if D == 4
        %d == -1 never happens
        L1 = correctLo(fitEI.Coefficients.Estimate(3),fitEI.Coefficients.SE(3));
        L2 = correctHi(fitEI.Coefficients.Estimate(3),fitEI.Coefficients.SE(3));
        %plot([L1,L2],[0.0,0.0],'-','linewidth',6,'color','b','DisplayName','ExpInt MIC 95% CI')

        L1 = correctLo(fitlog.Coefficients.Estimate(3),fitlog.Coefficients.SE(3));
        L2 = correctHi(fitlog.Coefficients.Estimate(3),fitlog.Coefficients.SE(3));
        %plot([L1,L2],[0.6,0.6],'-','linewidth',6,'color','r','DisplayName','Radical Exp MIC 95% CI')
    end    

    legend('location','northwest')
    axis tight
    xlim(thisxlim)
    ylim([0 R(end)])

    xlabel(xlabels{d})
    ylabel(ylabels{d})
    
    figure(D)
    exportgraphics(gcf,['./figures/ThreeModelCompareConvex',num2str(D),'.pdf'])

end
