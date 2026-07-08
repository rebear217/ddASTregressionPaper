% This script plots all figures in the paper which can be done
% by running the file one section at a time: (and do be sure to be in the right Matlab directory)

% This must be installed for MCMC codes to work:
% https://mjlaine.github.io/mcmcstat/
% Do beware that this toolbox overwrites the "boxplot" command.

clearvars
close all
clc

%% Figure 1: no codes are needed for this

%% BOTH Figures 2 and 3

% This may be needed:
cd('regressionMcodes/')

%% Figures 2A&3A
plotBonevModel;

%% Figures 2B&3B
plotExpIntModel;

%% Figures 2C&3C
plotNewLogModel;

%% Figures 2D&3D
% the 3 sections above must be run first before this will work:
plotFit3;

%% Figure 4A
% To test this works, open "MCMCanalysis.m" and change the line
% options.nsimu = 200000; to something with a smaller value, like...
% options.nsimu = 2000;

bonevMCMC;

%% Figure 4B

expintMCMC;

%% Figure 4C

newlogMCMC;

%% Figure 5A-C

plotExcisedDataRegressions;

%% Figure 6A-C

plotNonConvexData;

%% Figure 7A-J

% This may be needed:
cd('..')
cd('spatialPDEMcodes/')

for T = [2,6,8,10,13,14,15,16,20]
    outputT = timesolvePanel(T);
    close all
end

%% Figures 8-9

clearvars
close all
clc
load('./mats/allOutputs.mat')

% This loads in the results of lots of calls like "output20 = ZOIdoseReponse();"
% where 20 here denotes the (default) length of the simulated AST assay which can be
% changed by altering 'T = 20;' in ZOIdoseResponse.m:

% These computations have already been done for various values of T and
% it is much to quicker to load the results of those using the above load command.

%%

% Once loaded, run these:

%% Figure 8A-B

close all
clc
ZOIdoseResponse(output20);

%% Figure 8C

close all
processOutputsCell({output8,output9,output10,output11,output12,output13,output14,output15,output16,output17,output20});
axis tight

exportgraphics(gcf,'./figures/transition.PDF')

%% Figure 8D-E-F

close all
clc
plotON = 1;
plotOFF = 0;

M = 50;
Ns = [90 120 190];
for j = 1:3
    N = Ns(j);
    close all
    clc
    [sl,dash1] = processOutputs({output18},2,plotON,'--',M,N);
    [~,dash2] = processOutputs({output18},3,plotOFF,'-',M,N);
    [~,dash3] = processOutputs({output18},4,plotOFF,'-',M,N);
    
    sl.Color = 'k';
    dash1.Color = 'b';
    dash2.Color = 'r';
    dash3.Color = 'k';
    ylim([0 1.2])
    
    exportgraphics(gcf,['./figures/regressionTransition',num2str(j),'.PDF'])
end

%% Figure 9A

close all
outputCell = {output20};
sampleA0 = [0.0125 0.025 0.05 0.1 0.2 0.4 0.8 1.6 3.2];

figure(1)
thresholdedZOIdoseResponse(outputCell,0.35,sampleA0);
ylim([0 1])
axis tight

exportgraphics(gcf,'./figures/CARSthresholdRegressions.PDF')

%% Figure 9B

close all
figure(1)
outputCell = {output11,output12,output13,output14,output15,output16,output17,output18,output19,output20};
plotOFF = 0;
thresholdedZOIdoseResponse(outputCell,0.35,[],plotOFF);

exportgraphics(gcf,'./figures/CARSthresholdHysteresis.PDF')

