function MCMCanalysis(conc,zoi,func,iGuess,weights)

data.ydata = conc;
data.xdata = zoi;
W = sqrt(weights(data.ydata));
%W = ones(size(data.ydata));

modelfun = @(x,theta) func(theta,x);
ssfun = @(theta,data) sum( (data.ydata - modelfun(data.xdata,theta)).^2 .* W );

p = length(iGuess);
n = length(data.xdata);

model.ssfun  = ssfun;
options.nsimu = 200000;
options.updatesigma = 1;
options.burnintime = floor(options.nsimu/5);

params = {};
parameterList = {'MIC'};
for j = 1:p
    parameterList{j+1} = ['p_',num2str(j)];
    if j > 1
        vardetails = {parameterList{j},iGuess(j), -10*abs(iGuess(j)), 10*abs(iGuess(j)), iGuess(j)};
    else
        %vardetails = {parameterList{j},iGuess(j), 1.25, 2.5, iGuess(j)};
        %replace 0 in here with 2^(-9)? Not worth it?
        vardetails = {parameterList{j},iGuess(j), 0, 2^9, iGuess(j)};
    end
    params{j} = vardetails;
end

[res,chain,s2chain] = mcmcrun(model,data,params,options);

chain = chain(options.burnintime:end,:);
s2chain = s2chain(options.burnintime:end,:);

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
x = (0:0.1:11)';
plot(data.xdata,data.ydata,'o');
hold on
plot(x,modelfun(x,mean(chain)),'-k')
legend({'data','model'},'Location','best')

figure(5); clf
semilogy(data.xdata,data.ydata,'o','HandleVisibility','off'); % add data points to the plot
hold on
out = mcmcpred(res,chain,[],x,modelfun);
mcmcpredplot(out);
hold on
xlabel('r (mm)');
ylabel('antibiotic dose (\mug/mL)');
title('Predictive envelopes of the model')
semilogy(data.xdata,data.ydata,'.k','markersize',30); % add data points to the plot
axis tight
xL = xlim;
plot([xL(1),xL(1)],[0.625,1.25],'-k','linewidth',4)

Olims = out.predlims{1}{1};
for k = 1:9
    if ~(k == 5)
        semilogy(out.data{1},Olims(k,:),'-','LineWidth',0.5,'color',[1,1,1]*0.5,'Handlevisibility','off')
    end
end

legend({'99% CIs','95%','90%','50%','MMSEM','data','W interval'},'location','northwest')

chainstats(chain,res)

end