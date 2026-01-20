% This script plots all figures in the paper which can be done
% by un-commenting a line associated with that figure number and
% running that section: (be sure to be in the right directory)

% This must be installed for MCMC codes to work:
% https://mjlaine.github.io/mcmcstat/
% Beware that this toolbox overwrites the "boxplot" command.

clearvars
close all
clc

%% Figure 1: no codes are needed for this

%% BOTH Figures 2 and 3

% This may be needed:
cd('regressionMcodes/')

%%

% Figures 2&3A
plotBonevModel;

%%

% Figures 2&3B
plotExpIntModel;

%%

% Figures 2&3C
plotNewLogModel;

%%

% Figures 2&3D
plotFit3;

%% Figure 4A-D

plotNonConvexData;

%% Figure 5A-F: run all of these separately

close all
bonevMCMC;

%OR:
%close all
%expintMCMC;

%close all
%newlogMCMC;

% To test they work, open "MCMCanalysis.m" and change the line
% options.nsimu = 200000; to something with a smaller value, like...
% options.nsimu = 2000;

%% Figure 6A-J

% This may be needed:
cd('..')
cd('spatialPDEMcodes/')

for T = [2,6,8,10,13,14,15,16,20]
    outputT = timesolvePanel(T);
    close all
end

%% Figure 7

clearvars
close all
clc
load('./mats/allOutputs.mat')

% This loads in the results of lots of calls like "output20 = ZOIdoseReponse();"
% where 20 here denotes the (default) length of the simulated AST assay which can be
% changed by altering 'T = 20;' in ZOIdoseResponse.m:

% These computations have already been done for various values of T and
% it is much to quicker to load the results of those using the above load
% command.

%%

% Once loaded, run these:

%%

% Figure 7A&B
ZOIdoseReponse(output20);

%%

% Figure 7C
close all
processOutputs({output10,output13,output15,output20});

%%

% Figure 7D
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
sl11.Color = 'k';
dash1.Color = 'b';
dash2.Color = 'r';

legend('location','southwest')
axis([1e-5       3.6925      0.11507      0.83722])

