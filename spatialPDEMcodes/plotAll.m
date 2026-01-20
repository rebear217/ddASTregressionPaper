clear all
close all
clc

load('./mats/allOutputs.mat')

%%

%e.g. output10 = ZOIdoseReponse();

ZOIdoseReponse(output10);
ZOIdoseReponse(output13);
ZOIdoseReponse(output14);
ZOIdoseReponse(output15);
ZOIdoseReponse(output16);
ZOIdoseReponse(output17);
ZOIdoseReponse(output18);
ZOIdoseReponse(output19);
ZOIdoseReponse(output20);

%%

close all
for T = [2,6,8,10,13,14,15,16,20]
    outputT = timesolvePanel(T);
    close all
end

%%

close all
processOutputs({output10,output13,output15,output20});
exportgraphics(gcf,'./figures/transition.PDF')

%%

close all
sl2 = processOutputs({output15part2},0);
sl2.Color = [0        0.447        0.741];
hold on
[sl11,dash1] = processOutputs({output15},2);
[sl12,dash2] = processOutputs({output15},3);
sl12.HandleVisibility = 'off';
sl11.HandleVisibility = 'on';

delete(sl12)
sl11.LineWidth = 4;
%delete(sl1);
%sl12.LineWidth = 4;
sl11.Color = 'k';
dash1.Color = 'b';
dash2.Color = 'r';

%sl22.DisplayName = 'constant density regime';
%sl2.Color = [0 0 0];
%semilogx(1e-3*[1.953 1.953],[0.1 0.85],'--k','linewidth',1,'DisplayName','model fit limit (MFL)');
legend('location','southwest')
axis([1e-5       3.6925      0.11507      0.83722])
exportgraphics(gcf,'./figures/fitModelToModel.PDF')

%%

[mySurface,parameters,types] = runTheModel();
close all
plotWave(mySurface,parameters);

