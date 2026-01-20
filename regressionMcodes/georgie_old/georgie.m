clc
close all
clear vars

%%

model = @(p,x)(p(1) + p(2)*x + p(3)*x.^2);
guess = [1 1 1];

dosages = [16 32 64 128 256];
Dosages = repmat(dosages,6,1);
lg = log10(Dosages(:));
X = 10.^(lg);

%what is this doing? noise?: it's for plotting raw data nicely
R = rand(size(Dosages(:)));
R = 5*(R - 0.5);


for year = 1:2

    if year == 1
        doxy = [4.11900000000000,4.19100000000000,4.30100000000000,5.22800000000000,8.39800000000000;1.74400000000000,2.57300000000000,4.19100000000000,6.15800000000000,7.94200000000000;5.89600000000000,4.94800000000000,5.47400000000000,5.47400000000000,10.0660000000000;3.36500000000000,4.30100000000000,4.30100000000000,6.70000000000000,7.40200000000000;3.07900000000000,3.26900000000000,3.97600000000000,6.42400000000000,7.94200000000000;3.40000000000000,4.48600000000000,5.86900000000000,6.74300000000000,7.64500000000000];
        penicillin = [2.40500000000000,2.21700000000000,4.19100000000000,8.70900000000000,9.07900000000000;0.622000000000000,2.06100000000000,3.97600000000000,6.69700000000000,10.9860000000000;2.89500000000000,3.66400000000000,4.75300000000000,7.64500000000000,12.5660000000000;3.14200000000000,4.75300000000000,7.99200000000000,7.35400000000000,10.5210000000000;1.76700000000000,4.60000000000000,5.55700000000000,7.84300000000000,11.8850000000000;1.96100000000000,3.66400000000000,3.90600000000000,7.45100000000000,9.07900000000000];
        %don't use area, use radius:
        doxy = sqrt(doxy/pi);
        penicillin = sqrt(penicillin/pi);
        %use diameters:
        doxy = doxy*2;
        penicillin = penicillin*2;
    else
        doxy = [2.80000000000000,2.40000000000000,2.40000000000000,2.10000000000000,1.60000000000000,2.35000000000000;3.40000000000000,2.70000000000000,2.70000000000000,2.40000000000000,2.60000000000000,2.55000000000000;4,3.40000000000000,3.40000000000000,3.70000000000000,2.40000000000000,2.90000000000000;4.80000000000000,4.40000000000000,4,4,3.60000000000000,3.05000000000000;5.60000000000000,4.80000000000000,4.80000000000000,4.50000000000000,4.80000000000000,3.83000000000000]';
        penicillin = [2,2.20000000000000,1.50000000000000,1.50000000000000,2,2.58000000000000;2.40000000000000,2.40000000000000,2.10000000000000,1.80000000000000,2.10000000000000,2.90000000000000;2.70000000000000,3.40000000000000,2.30000000000000,2,2.50000000000000,3.75000000000000;3.60000000000000,4.20000000000000,2.90000000000000,3.30000000000000,2.90000000000000,4.10000000000000;4.20000000000000,3,3.50000000000000,3.80000000000000,3.30000000000000,4.80000000000000]';
    end

    %%

    figure(year)
    set(year,'pos',[105         358        1257         448])


    subplot(1,2,1)
    r=semilogx(Dosages(:)+R,penicillin(:),'.k','markersize',20,'color',[1 1 1]*0.7);
    hold on
    p=semilogx(dosages,mean(penicillin),'.k','markersize',40);
    xlabel('multiple of clinical dose')
    ylabel('size of inhibition zone')
    ylim([0 6])
    
    PenModel = fitnlm(lg,penicillin(:),model,guess);
    disp('pen...')
    PenModel.Coefficients.pValue
    
    CIlow = PenModel.Coefficients.Estimate - 1.96*PenModel.Coefficients.SE;
    CIhigh = PenModel.Coefficients.Estimate + 1.96*PenModel.Coefficients.SE;
    
    q=semilogx(X,PenModel.feval(lg),'--k');
    legend([p q r],{'penicillin mean',['datafit (quadratic coefficient 95% CI (',num2str(CIlow(3),2),...
        ',',num2str(CIhigh(3),2),'))'],'image ZoI data'})
    legend('boxoff')
    
    subplot(1,2,2)
    r=semilogx(Dosages(:)+R,doxy(:),'.k','markersize',20,'color',[1 1 1]*0.7);
    hold on
    p=semilogx(dosages,mean(doxy),'.k','markersize',40);
    xlabel('multiple of clinical dose')
    ylabel('size of inhibition zone')
    ylim([0 6])
    
    disp('--------')
    
    DoxModel = fitnlm(lg,doxy(:),model,guess);
    CIlow = DoxModel.Coefficients.Estimate - 1.96*DoxModel.Coefficients.SE;
    CIhigh = DoxModel.Coefficients.Estimate + 1.96*DoxModel.Coefficients.SE;
    q=semilogx(X,DoxModel.feval(lg),'--k');
    legend([p q r],{'doxycycline mean',['datafit (quadratic coefficient 95% CI (',num2str(CIlow(3),2),...
        ',',num2str(CIhigh(3),2),'))'],'image ZoI data'})
    legend('boxoff')

	disp('doxy...')
    DoxModel.Coefficients.pValue
 
   
    export_fig(['modelFit',num2str(year),'.pdf']);
    
end

