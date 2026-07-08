function [conc,zoi,R] = defineMicrococcusNisinData()

    % data from literature:

    % Micrococcus luteus (ATCC 10240) and nisin
    % data taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC4576963/pdf/fsn30003-0394.pdf
    % Figure 7 (see legend of A-J)

    conc = [125 75 50 25 18.75 12.5 6.25 2.5 1.25 0.625];
    zoi = [10.7 10.5 10.3 9.5 9 8 6.5 4 0 0];
    R = 0:0.01:12;

end