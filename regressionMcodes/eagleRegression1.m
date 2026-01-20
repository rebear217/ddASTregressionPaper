clc
close all
clearvars

%%  data from literature

%Micrococcus luteus (ATCC 10240) and nisin
%data taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC4576963/pdf/fsn30003-0394.pdf
%Figure 7 (see legend of A-J)

%2 very small non-zero values have been added to 0 in the last 2 entries to
%simplify the code for the Jacknife fits:

LO = 2;

%Micrococcus luteus (ATCC 10240) and nisin
%data taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC4576963/pdf/fsn30003-0394.pdf
%Figure 7 (see legend of A-J)
Ab = [125 75 50 25 18.75 12.5 6.25 2.5 1.25 0.625]';
radius = [10.7 10.5 10.3 9.5 9 8 6.5 4 0.0001 0.000001]';

%data taken from
%https://pmc.ncbi.nlm.nih.gov/articles/PMC546645/pdf/applmicro00228-0032.pdf
%PAP is unclear (not stated in paper): Bacillus subtilis or Sarcina lutea?
%(Table 2)
%Ab = [3 6 12 24 48];
%radius = [11.68 16.62 20.96 24.86 28.21];

Z = zeros(size(radius(1:end-LO)));
ARdata = [Ab(1:end-LO) , radius(1:end-LO)];

R = 0:0.01:12;

bonevArea = @(p,A)abs(p(1))*log(A) + p(2)*exp(-abs(p(3))*A) - abs(p(4));
bonevF = @(p,A)abs(p(1))*log(A) - abs(p(2));

b0 = [1 1 1 1];
weights = @(yhat) 1./(1 + abs(yhat).^2);

fitArea = fitnlm(Ab(1:end-LO),radius(1:end-LO).^2,...
    bonevArea,b0,'Weights',weights)

fitBonev = fitnlm(Ab(1:end-LO),radius(1:end-LO).^2,...
    bonevF,[1 1],'Weights',weights)

RL = exp(-abs(fitArea.ModelCriterion.AIC-fitBonev.ModelCriterion.AIC)/2)

AbFine = 0.5:0.0125:max(Ab);

figure(1)
plot(log2(Ab),radius,'ok','MarkerSize',14);
hold on
plot(log2(Ab(1:end-LO)),sqrt(fitArea.Fitted),'.b','MarkerSize',40);
plot(log2(AbFine),sqrt(fitArea.feval(AbFine)),'-k');

%plot(log2(AbFine),sqrt(fitBonev.feval(AbFine)),'--r');

xlabel('log2 antibiotic dose)')
ylabel('radius (mm)')

exportgraphics(gcf,'./figures/AreaRegression.PDF')

%%

normalForm = @(Ab,a) (log(Ab) + a*exp(-Ab));
Abx = 10.^(-2:0.0005:1.1);

figure(2)
set(2,'pos',[43         255        1238         402])

for a = -4:4
    lw = 1;
    DN = ['a = ',num2str(a)];
    if a == 0
        lw = 5;
        DN = ['a = ',num2str(a),' (Bonev)'];
    end
    if a == 2
        lw = 5;
    end
    nf = normalForm(Abx,a);
    nfp = nf > 0;
    subplot(1,2,1)
    sl = semilogx(Abx(nfp),pi*nf(nfp),'DisplayName',DN,...
        'linewidth',lw);
    hold on
    subplot(1,2,2)
    pl = plot(Abx(nfp),2*sqrt(nf(nfp)),'DisplayName',DN,...
        'linewidth',lw);
    hold on
    if a == 0
        sl.Color = 'k';
        pl.Color = 'k';
    end
    if a == 2
        sl.Color = [1 1 1]*0.6;
        pl.Color = [1 1 1]*0.6;
    end
end

for sp = 1:2
    subplot(1,2,sp)
    title('$\pi(\log(A) + a\exp(-A))$','Interpreter','latex')        
    xlabel('antibiotic dose')
    ylabel('ZoI area')
    axis tight
    legend('location','northwest')
    %yL = ylim;
    %ylim([0 yL(2)])
    H = 1.5*pi;
    if sp == 2
        ylabel('ZoI diameter')    
        H = 2*sqrt(H/pi);
        legend('location','southeast')
        xlim([0 10])
        title('$2(\log(A) + a\exp(-A))^{1/2}$','Interpreter','latex')
    end
    plot(Abx,H*ones(size(Abx)),'--k','DisplayName','one ZoI, 3 dosages');
end

%%

exportgraphics(gcf,'./figures/AreaRegressionNormalForm.PDF')
