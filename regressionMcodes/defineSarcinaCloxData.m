function [conc,zoi,R] = defineSarcinaCloxData()

    % data from literature:

    % data taken from
    % https://pmc.ncbi.nlm.nih.gov/articles/PMC546645/pdf/applmicro00228-0032.pdf
    % PAP is unclear (not stated in paper): Bacillus subtilis or Sarcina lutea?
    % (Table 2)

    conc = [3 6 12 24 48];
    zoi = [11.68 16.62 20.96 24.86 28.21];
    R = 0:0.01:32;

end