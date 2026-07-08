clearvars
close all
clc

gr = [1,1,1]/2;
JACK = [0.2 0.5 0.2];

%%

clc
%warning('off')

%dataset 1 from https://europepmc.org/article/PMC/4576963, Table 1

conc = [0.625 , 0.312 , 0.156, 0.08, 0.04];
zoi = [3.44 , 2.4 , 1.65, 1.2, 0.35];
R = 0:0.001:4;

b0bon = [0.093048      -0.1587];
b0expint = [0.037079     -0.11332      0.12774];
b0log = [0   0.599      -4.15];
Iguesses = {b0log,b0expint,b0bon};

% this model is no longer used:
%logF = @(p,r)abs(p(3)) + p(1)*exp((-1 + (1+p(2)*r.^2).^(1/2))/2 + ...
%        log(((-1 + (1+p(2)*r.^2).^(1/2))))/2 + ...
%        (p(2)/2)*r.^2./(-1 + (1+p(2)*r.^2).^(1/2)));

% Here (due to difficulties getting convergence otherwise) ...
% we constrain the frequentist 'fitnlm'approach to positive MIC values
% using a smooth nonlinear function f(p(1))
% MUST remember to take inverse if the parameter value p(1) is needed
% variation/stats are therefore stated in terms of f(p(1)), not p(1) itself

MICfunc = @(p)log10(1+p^2);
MICinversefunc = @(q)sqrt(10.^q-1);

% These impose EUCAST bounds on MICs (but makes no practical difference)
EUCASTMIClb = MICinversefunc(0); %this is zero
EUCASTMICub = MICinversefunc(2^9); %this is massive, overflow? Yes: this is 'inf' in Matlab

% Thus, this MIC nonlinear mapping approach cannot apply EUCAST bounds and
% merely implies the MIC is positive, rather than bounded

logFmod = @(p,r)MICfunc(p(1)) + exp(p(3) + sqrt(1+abs(p(2)).*r.^2 )) .* ( -1 + sqrt(1+abs(p(2))*r.^2) );
expintFmod = @(p,r)(MICfunc(p(1)) + p(3) ./ expint(abs(p(2))*r.^2));
bonevFmod = @(p,r)MICfunc(p(1)).*exp(-r.^2*p(2));
FmodelsMod = {logFmod,expintFmod,bonevFmod};

logF = @(p,r)abs(p(1)) + exp(p(3) + sqrt(1+abs(p(2)).*r.^2 )) .* ( -1 + sqrt(1+abs(p(2))*r.^2) );
expintF = @(p,r)(abs(p(1)) + p(3) ./ expint(abs(p(2))*r.^2));
bonevF = @(p,r)abs(p(1)).*exp(-r.^2*p(2));
Fmodels = {logF,expintF,bonevF};

for pow = [0,1,2]
    
    close all

    weights = @(yhat) 1./(abs(yhat).^pow);
    %weights = @(yhat) ones(size(yhat));
    
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
    
        for D = 1:6
    
            if D == 6
                concExcise = conc;
                zoiExcise = zoi;
                lw = 3;
                lab = labels{ml};
                col = cols{ml};
            else
                concExcise = setdiff(conc,conc(D));
                zoiExcise = setdiff(zoi,zoi(D));
                lw = 2;
                lab = ['Jacknife ',numbers{D}];
                col = JACK;
            end    
            
            converged = 1;
            try
                %as convergence can be delicate, bail out if it fails:
                fitModel = fitnlm(zoiExcise,concExcise,FmodelsMod{ml},Iguesses{ml},'Weights',weights,'Options',opt)
                fitModel.Coefficients.Estimate'

                semilogx(fitModel.feval(R),R,'-','linewidth',2,'color',col,'DisplayName',...
                    [lab,' (adj R^2\approx',num2str(fitModel.Rsquared.Adjusted,2),')'],...
                    'linewidth',lw)                
            catch
                converged = 0;
            end
            
        end
    
        if D == 6
            MIC = fitModel.feval(0);
            text(0.25,0.5,['MIC\approx',num2str(MIC,2),'\mug/mL']);
            Coefficients = fitModel.Coefficients.Estimate;

            % use the MIC-transformed models in the MCMC:
            %MMSEM = MCMCupdatePointEstimates(concExcise,zoiExcise,FmodelsMod{ml},...
            %    Coefficients,weights,'f(MIC)',EUCASTMIClb,EUCASTMICub);

            % do NOT use the MIC-transformed models in the MCMC:
            Coefficients(1) = MICinversefunc(Coefficients(1));
            MMSEM = MCMCupdatePointEstimates(concExcise,zoiExcise,Fmodels{ml},...
                Coefficients,weights,'MIC',0,2^(9));
        end
    
        if ml == 3
            fitLM = fitlm(zoi,log(conc),'Weights',weights(conc));
            semilogx(exp(fitLM.feval(R)),R,'-','linewidth',2,'color',[1,1,1]*0.5,...
                'DisplayName',['Log-linear (adj R^2\approx',num2str(fitLM.Rsquared.Adjusted,2),')'])
            MICllr = exp(fitLM.feval(0));
            text(0.25,0.2,['LLR-MIC\approx',num2str(MICllr,2),'\mug/mL']);
        end

        legend('location','northwest')
        axis tight
        xlim([0.01 1.1])
        ylim([0 4])        
    
        xlabel('erythromycin concentration (\mug/mL)');
        ylabel('zone of inhibition radius (mm)');
    
        exportgraphics(gcf,['./figures/ThreeModelCompareConvex',num2str(ml),num2str(pow),'.pdf'])

    end

end

%warning('on')
