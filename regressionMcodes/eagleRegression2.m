%% 

clearvars
close all
clc

X = [85 159 173 198 213 238 255 292 372 451 528 806];
Y = [38 154 203 394 400 415 406 343 259 244 254 259];

Y = 1+(-Y + min(Y))/(max(Y) - min(Y));
Y = 0.8*Y;
Y = 7*Y-3.5;
X = (X - min(X))/(max(X) - min(X));
X = 18*X-7;

figure(1)
plot(X,Y,'.-k','MarkerSize',30)
axis tight
xlabel('penicillin ($\mu g/cm^{3}$)','Interpreter','latex')
ylabel('proportion of surviving organisms')

%%

figure(2)
A = 2.^X;
PLO = 10.^Y;

BLUE = [0.1 0.6 1];

loglog(A,PLO,'.-k','MarkerSize',35,'linewidth',1,...
    'DisplayName','S. aureus data (Eagle-Musselmen 1948, Fig 12)')
axis tight
xlabel('penicillin ($\mu g/cm^{3}$)','Interpreter','latex')
ylabel('proportion of surviving organisms')

b0 = [0.011 39661 196 0.02];
eagle = @(p,A)p(1)*log(abs(A)) + ...
    p(2)*exp(-abs(p(3))*abs(A)) + p(4);

weights = @(yhat) 1./(1 + abs(yhat).^(1/3));
%weights = @(yhat) abs(yhat).^(1);

fitEagle = fitnlm(A,PLO,eagle,b0,'Weights',weights)
pSol = fitEagle.Coefficients.Estimate;

fineA = 2.^(min(X):0.001:max(X));

[fitEagle2,RR,JJ] = nlinfit(A,PLO,eagle,pSol);
[ypred,delta] = nlpredci(eagle,fineA,fitEagle2,RR,'Jacobian',JJ);

hold on
plot(fineA,eagle(pSol,fineA),'-','color',BLUE,...
    'linewidth',4,'DisplayName','Eagle regresion');
legend('location','northeast')
plot(fineA,ypred+delta,'--','color',BLUE,...
    'linewidth',1,'DisplayName','Eagle regresion',...
    'DisplayName','95% CI');
Ylo = ypred-delta;

plot(fineA(Ylo > 0),Ylo(Ylo > 0),'--','color',BLUE,...
    'linewidth',1,'DisplayName','Eagle regresion',...
    'HandleVisibility','off');

ylim([10^(-3.5) 140])
plot(A,PLO,'.-k','MarkerSize',35,'linewidth',1,...
    'HandleVisibility','off')

%%

exportgraphics(gcf,'./figures/EagleMusselRegression.PDF')

%%

figure(2)
eagleLogistic = @(p,A)p(1)*log(abs(A)) + ...
    p(2)./(1+p(5)*exp(-abs(p(3))*abs(A))) + p(4);

b1 = [fitEagle.Coefficients.Estimate ; 1];

weightsLogistic = @(yhat) 1./(1 + abs(yhat).^(1/2));

fitLogisticEagle = fitnlm(A,PLO,eagleLogistic,...
    b1,'Weights',weightsLogistic)
pSolLog = fitLogisticEagle.Coefficients.Estimate;

%plot(fineA,eagleLogistic(pSolLog,fineA),'-','color',BLUE/2,...
%    'linewidth',2,'DisplayName','Logistic-Eagle regresion');

