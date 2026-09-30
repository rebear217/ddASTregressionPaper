% make paper figure names:

%%

if ~exist('paperfigures','dir')
    mkdir('paperfigures');
end

%%

oldNamesREG = {'spatialModelDoseResponse20A','spatialModelDoseResponse20B','transition',...
               'regressionTransition1','regressionTransition2','regressionTransition3',...
               'CARSthresholdHysteresis','CARSthresholdRegressions'};

newNamesREG = {'8A','8B','8C','8D','8E','8F','9A','9B'};

N = length(newNamesREG);
for j = 1:N
    srcfile = ['./figures/',oldNamesREG{j},'.pdf'];
    destfile = ['./paperfigures/Figure',newNamesREG{j},'.pdf'];
    try
        movefile(srcfile,destfile);
        disp(['file moved: ',srcfile])
    catch
        disp(['file not found error: ',srcfile])
    end
end

