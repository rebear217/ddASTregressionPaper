function MCMCshadedPlot(out,conc,zoi,coreColour)

if nargin < 4
    coreColour = [0,0,0];
end
fg = gcf;

Olims = out.predlims{1}{1};
labs = {'99% CIs','95%','90%','50%','MMSEM','50%','90%','95%','99%'};

for k = [1 2 3 4 6 7 8 9]
    semilogx(Olims(k,:),out.data{1},'-','LineWidth',1,'color',[1,1,1]/3,'DisplayName',labs{k},...
        'HandleVisibility','off')
    hold on
end

for k = 1:4
    col = [1,1,1]*(1-k/5);
    P{k} = plotshadedYX([Olims(k,:) ; Olims(10-k,:)],(out.data{1})','k');
    P{k}.FaceColor = col;
    alpha(0.3);    
    P{k}.DisplayName = labs{k};
end

k = 5;
semilogx(Olims(k,:),out.data{1},'-','LineWidth',4,'color',coreColour,'DisplayName',labs{k})

plot(conc,zoi,'.k','markersize',40,'DisplayName','data')
legend('location','northwest')
ylabel('r (mm)')
xlabel('antibiotic dose (\mug/mL)')
axis tight

end