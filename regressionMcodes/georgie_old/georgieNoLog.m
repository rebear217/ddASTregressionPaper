clc
close all
clear vars

%%

gray = 0.7*[1 1 1];
red = [1 0 0]*0.75;
blue = [0 0 1]*0.75;
black = 0.1*[1 1 1];

mdlLabel = {'power law','Bonev2008','diffusion kill theory'};
figNames = {'pen1','dox1','pen2','dox2'};
colours = {blue,red,black,gray};
styles = {'-','-','-','-'};
Q = [0 0 0 0];

myModel = @(p,x)(...
p(1)+...
(-1 + (1+p(3)*x.^2).^(1/2))+...
log(abs((-1 + (1+p(3)*x.^2).^(1/2))))+...
p(2)*x.^2./(-1 + (1+p(3)*x.^2).^(1/2)));

models = {@(p,x)(p(2)*x.^p(1)),...
    @(p,x)p(1).*(log(abs(x./p(2)))).^(1/2),...
    myModel};

for mdl = 1:3

    model = models{mdl};

    guess = [0 1];
    if mdl == 3
        guess = [0.5 1 10];
    end

    dosages = [16 32 64 128 256];
    Dosages = repmat(dosages,6,1);
    %lg = log10(Dosages(:));
    %X = 10.^(lg);

    allDosages = (Dosages(:));
    %allDosages = allDosages;

    %what is this for? for plotting? Yep - for moving raw data apart to be visible
    R = rand(size(Dosages(:)));
    R = 5*(R - 0.5);


    for year = 1:2

        if year == 1
            doxy = [4.119,4.191,4.301,5.228,8.398;1.744,2.573,4.191,6.158,7.942;5.896,4.948,5.474,5.47400000000000,10.0660000000000;3.36500000000000,4.30100000000000,4.30100000000000,6.70000000000000,7.40200000000000;3.07900000000000,3.26900000000000,3.97600000000000,6.42400000000000,7.94200000000000;3.40000000000000,4.48600000000000,5.86900000000000,6.74300000000000,7.64500000000000];
            penicillin = [2.405,2.217,4.191,8.709,9.079;0.622,2.061,3.976,6.697,10.986;2.89500000000000,3.66400000000000,4.75300000000000,7.64500000000000,12.5660000000000;3.14200000000000,4.75300000000000,7.99200000000000,7.35400000000000,10.5210000000000;1.76700000000000,4.60000000000000,5.55700000000000,7.84300000000000,11.8850000000000;1.96100000000000,3.66400000000000,3.90600000000000,7.45100000000000,9.07900000000000];
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

        fineallDosages = 1:0.1:1.1*max(allDosages);

        figure(2*(year-1)+1)
        %set(year,'pos',[105         358        1257         448])

        %subplot(1,2,1)
        r=semilogx(Dosages(:)+R,penicillin(:),'.k','markersize',20,'color',[1 1 1]*0.7);
        hold on
        p=semilogx(dosages,mean(penicillin),'.k','markersize',40);
        xlabel('multiple of clinical dose')
        ylabel('size of inhibition zone')
        ylim([0 6])

        PenModel = fitnlm(allDosages,penicillin(:),model,guess)
        disp('pen AIC...')
        PenModel.ModelCriterion.AIC
        AICpen{mdl} = PenModel.ModelCriterion.AIC;
        adjR2{mdl} = PenModel.Rsquared.Adjusted;
        
        %PenModel.Coefficients.pValue

        CIlow = PenModel.Coefficients.Estimate - 1.96*PenModel.Coefficients.SE;
        CIhigh = PenModel.Coefficients.Estimate + 1.96*PenModel.Coefficients.SE;
        
        Y = PenModel.feval(fineallDosages);
        Y = Y(Y > 0);
        q(mdl) = semilogx(fineallDosages(Y > 0),Y,styles{mdl},'color',colours{mdl});

        figure(2*(year-1)+2)
        %subplot(1,2,2)
        rr=semilogx(Dosages(:)+R,doxy(:),'.k','markersize',20,'color',[1 1 1]*0.7);
        hold on
        pp=semilogx(dosages,mean(doxy),'.k','markersize',40);
        xlabel('multiple of clinical dose')
        ylabel('size of inhibition zone')
        ylim([0 6])

        disp('--------')

        DoxModel = fitnlm(allDosages,doxy(:),model,guess)
        Y = DoxModel.feval(fineallDosages);
        Y = Y(Y > 0);
        CIlow = DoxModel.Coefficients.Estimate - 1.96*DoxModel.Coefficients.SE;
        CIhigh = DoxModel.Coefficients.Estimate + 1.96*DoxModel.Coefficients.SE;
        qq(mdl) = plot(fineallDosages(Y > 0),Y,styles{mdl},'color',colours{mdl});
        %legend([p q r],{'doxycycline mean',['datafit (quadratic coefficient 95% CI (',num2str(CIlow(3),2),...
        %    ',',num2str(CIhigh(3),2),'))'],'image ZoI data'})
        %legend('boxoff')

        disp('doxy AIC...')
        DoxModel.ModelCriterion.AIC
        AICdox{mdl} = DoxModel.ModelCriterion.AIC;

        %DoxModel.Coefficients.pValue

        if mdl == 1
            penPower = PenModel.Coefficients.Estimate(1);
            doxPower = DoxModel.Coefficients.Estimate(1);
        end
        
        
        if mdl == 3
            figure(2*(year-1)+1)
            %subplot(1,2,1)
                axis tight

            legend([p q r],{'penicillin mean',[mdlLabel{1},', p\approx',num2str(penPower,2),', (adj R^2 \approx',num2str(adjR2{1},2),')'],...
                [mdlLabel{2},' (adj R^2 \approx',num2str(adjR2{2},2),')'],...
                [mdlLabel{3},' (adj R^2 \approx',num2str(adjR2{3},2),')'],...
                'ZoI data'},'location','southeast')
            %legend('boxoff')

            figure(2*(year-1)+2)
            %subplot(1,2,2)
                axis tight

            legend([pp qq rr],{'doxycycline mean',[mdlLabel{1},', p\approx',num2str(doxPower,2),', (adj R^2 \approx',num2str(adjR2{1},2),')'],...
                [mdlLabel{2},' (adj R^2 \approx',num2str(adjR2{2},2),')'],...
                [mdlLabel{3},' (adj R^2 \approx',num2str(adjR2{3},2),')'],...
                'ZoI data'},'location','southeast')
            %legend('boxoff')
        end

        figure(2*(year-1)+1)
        %subplot(1,2,1)
        xlim([1 270])
        yl = ylim;
        ylim([0 yl(2)]);
        figure(2*(year-1)+2)
        %subplot(1,2,2)
        xlim([1 270])
        yl = ylim;
        ylim([0 yl(2)]);

    end
    
end

for f = 1:4
    figure(f)
	export_fig([figNames{f},'.pdf']);
end

%%


%plot((allDosages),myFit.feval((allDosages)),'ok');
