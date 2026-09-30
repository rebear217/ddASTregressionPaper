clear all
close all
clc

gr = [1,1,1]/2;
JACKcolour = [0.2 0.5 0.2];

%%

clc

% dataset: https://pmc.ncbi.nlm.nih.gov/articles/PMC5761473, Table 3
conc = [2.56, 3.2, 4.0, 5.0 6.25];
zoi = [19.6, 21.8, 22.9, 24.5, 26.7];

xL = [0.05 10];
R = 0:0.001:28;

radExpF = @exponentialRadical;
expintF = @expintFunction;
bonevF = @bonevFunction;
Fmodels = {radExpF,expintF,bonevF};

% Due to difficulties getting fitlm convergence, here we constrain to positive
% MICvalues using the smooth nonlinear term MIC = p(1)^2:
% remember to take square roots later if the MIC parameter value is
% needed...:
sqP = @(p)[p(1)^2,p(2:end)];
radExpFsquared = @(p,r)radExpF(sqP(p),r);
expintFsquared = @(p,r)expintF(sqP(p),r);
bonevFsquared = @(p,r)bonevF(sqP(p),r);
FmodelsSquared = {radExpFsquared,expintFsquared,bonevFsquared};

p0radExp = [0.46585    0.0070733      -1.0188];
p0expint = [0.46585    0.0015803   1.1478];
p0bonev = [1   -0.0028004];
Iguesses = {p0radExp,p0expint,p0bonev};

weights = @(yhat) 1./(abs(yhat).^2);

opt = statset('fitnlm');
opt.MaxIter = 100000;
opt.TolFun = 1e-10;
opt.TolX = 1e-10;

correctLo = @(m,s) m - 1.96*s;
correctHi = @(m,s) m + 1.96*s;

cols = {'r','b','k'};
labels = {'Radical Exp','ExpInt','Bonev'};
numbers = {'1st','2nd','3rd','4th','5th'};

for ml = 1:3

    figure(ml)
    semilogx(conc,zoi,'.k','markersize',44,'DisplayName','ZoI data')
    hold on

    if ml == 3
        fitLM = fitlm(zoi,log(conc),'Weights',weights(conc));
        semilogx(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',[1,1,1]*0.7,...
            'DisplayName',['Log-linear (adj R^2\approx',num2str(fitLM.Rsquared.Adjusted,2),')'])
        MICllr = exp(fitLM.feval(0));
        text(0.065,8.5,['LLR-MIC\approx',num2str(MICllr,3),'\mug/mL']);
    end

    for D = 1:6
        
        if D == 6
            concExcise = conc;
            zoiExcise = zoi;
            lw = 3;
            lab = labels{ml};
            col = cols{ml};
        else
            %concExcise = setdiff(conc,conc(D));
            %zoiExcise = setdiff(zoi,zoi(D));
            concExcise = myRemoveDatum(conc,conc(D));
            zoiExcise = myRemoveDatum(zoi,zoi(D));
            
            lw = 2;
            lab = ['Jackknife ',numbers{D}];
            col = JACKcolour;
        end
    
        fitModel = fitnlm(zoiExcise,concExcise,FmodelsSquared{ml},Iguesses{ml},'Weights',...
            weights(concExcise),'Options',opt)
        fitModel.Coefficients.Estimate'
        
        semilogx(fitModel.feval(R),R,'-','linewidth',2,'color',col,'DisplayName',...
            [lab,' (adj R^2\approx',num2str(fitModel.Rsquared.Adjusted,2),')'],...
            'linewidth',lw)
        hold on

        if D == 6
            MIC = fitModel.feval(0);
            Coefficients = fitModel.Coefficients.Estimate;

            % use MCMC WITHOUT the MIC-transformed formulation:
            Coefficients = sqP(Coefficients')';

            % MCMC lower bound on the MIC:
            MIClb = 0;
            % MCMC upper bound on the MIC:
            MICub = 2^9;
            MMSEM = MCMCupdatePointEstimates(concExcise,zoiExcise,Fmodels{ml},Coefficients,...
                weights(concExcise),'MIC',MIClb,MICub);
        end
        
    end

    legend('location','northwest')
    axis tight
    xlim(xL)
    ylim([0 28])

    xlabel('Levofloxacin concentration (\mug/mL)')
    ylabel('zone of inhibition diameter (mm)')        

    text(0.1,10,['MIC\approx',num2str(MIC,3),'\mug/mL']);

    exportgraphics(gcf,['./figures/ThreeModelExcisedData',num2str(ml),'.pdf'])
    
end


