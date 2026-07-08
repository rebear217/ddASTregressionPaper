function MMSEM = MCMCupdatePointEstimates(conc,zoi,func,iGuess,weights,label,MIClb,MICub)

    if nargin < 6
        label = 'MIC';
    end
    
    if nargin < 7
        % EUCAST prior on the MIC bound:
        MICub = 2^9;
        MIClb = 0;
    end
    
    %{
    % MIC prior bounds from data (assume MIC is below all concenrations where zoi data are +ve):
    Z = (zoi > 0);
    cZ = conc(Z);
    [~,ZeroJ] = min(zoi(Z));
    MICub = cZ(ZeroJ);
    % this is an assumption from EUCAST regulations (MICs are always above this value)
    MIClb = 2^(-9);
    %}
    
    data.ydata = conc;
    data.xdata = zoi;
    W = weights(data.ydata);
    
    modelfun = @(x,theta) func(theta,x);
    ssfun = @(theta,data) sum( (data.ydata - modelfun(data.xdata,theta)).^2 .* W );
    
    p = length(iGuess);
    n = length(data.xdata);
    
    model.ssfun  = ssfun;
    options.nsimu = 20000;
    options.updatesigma = 1;
    options.burnintime = floor(options.nsimu/5);
    
    params = {};
    parameterList = {label};
    
    for j = 1:p
        parameterList{j+1} = ['p_',num2str(j)];
    
        if iGuess(j) > 0
            lb = iGuess(j)/10;
            ub = iGuess(j)*10;
        else
            ub = iGuess(j)/10;
            lb = iGuess(j)*10;
        end
    
        if j > 1
            vardetails = {parameterList{j},iGuess(j), lb, ub, iGuess(j)};
        else
            vardetails = {parameterList{j},iGuess(j), MIClb, MICub, iGuess(j)};
        end
        params{j} = vardetails;
    end
    
    [res,chain,s2chain] = mcmcrun(model,data,params,options);
    
    chain = chain(options.burnintime:end,:);
    s2chain = s2chain(options.burnintime:end,:);
    x = (0:0.1:(max(data.xdata)+2))';
    
    %{
    figure(2); clf
    mcmcplot(chain,[],res,'chainpanel');
    
    figure(3); clf
    mcmcplot(chain(options.burnintime:end,:),[],res,'pairs');
    
    figure(4); clf
    %mcmcplot(chain,[],[],'dens')
    %title('Error std posterior')
    mcmcplot(chain,1,res,'hist','kernel');
    xL = xlim;
    xlim([0 min([xL(2),4])]);
    xlabel('estimated MICs (\mu g/mL)')
    title('')
    
    figure(1); clf
    plot(data.xdata,data.ydata,'o');
    hold on
    plot(x,modelfun(x,mean(chain)),'-k')
    legend({'data','model'},'Location','best')
    %}
    
    %semilogx(data.ydata,data.xdata,'.k','markersize',30,'HandleVisibility','off'); 
    %hold on
    out = mcmcpred(res,chain,[],x,modelfun);
    %mcmcpredplot(out);
    
    Olims = out.predlims{1}{1};
    labs = {'99% CIs','95%','90%','50%','MMSEM','50%','90%','95%','99%'};
    
    %for k = [1 9]
    %    semilogx(Olims(k,:),out.data{1},'--','LineWidth',1,'color',[0,0,0],'DisplayName',labs{k})
    %end
    
    k = 5;
    semilogx(Olims(k,:),out.data{1},'-','LineWidth',2,'color',[1,0,1],'DisplayName',labs{k})
    MMSEM.Yresponse = out.data{1};
    MMSEM.Xdose = Olims(k,:);
    
    P = plotshadedYX([Olims(2,:) ; Olims(8,:)],(out.data{1})','k');
    alpha(0.1);
    P.DisplayName = '95% Credible Intervals';
    hold on
    
    axis tight
    
    chainstats(chain,res)

end