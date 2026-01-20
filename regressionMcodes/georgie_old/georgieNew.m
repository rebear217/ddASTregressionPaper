clc
close all
clear vars

%%

gray = 0.7*[1 1 1];
red = [1 0 0]*0.75;
blue = [0 0 1]*0.75;
black = 0.1*[1 1 1];

mdlLabel = {'Bonev2008','EI','Exponential rad'};

figNames = {'pen1','dox1','pen2','dox2'};
colours = {blue,red,black,gray};
styles = {'-','-','-','-'};

Q = [0 0 0 0];

logF = @(p,r)p(3) + exp(p(1) + sqrt(1+p(2).*r.^2 )) .* ( sqrt(1+p(2)*r.^2) );
b0F = [0.03 0.4 1.75];

fexpint = @(b,r)(abs(b(3)) + b(1) ./ expint(abs(b(2))*r.^2));
b0ei = [0.5 0.02 1];

bonevF = @(p,r)p(1).*exp(-r.^2*p(2));
b0b = [1.3 -0.03];

models = {bonevF,fexpint,logF};
b0 = {b0b,b0F,b0ei};

for mdl = 1:3

    model = models{mdl};
    guess = b0{mdl};

    dosages = [16 32 64 128 256];
    Dosages = repmat(dosages,6,1);

    allDosages = Dosages(:);
    fineAllDosages = 0.01:0.01:1.1*max(allDosages);

    %what is this for? for plotting? Yep - for moving raw data apart to be visible
    R = rand(size(Dosages(:)));
    R = 5*(R - 0.5);

    for year = 1:2
        if year == 1
            doxyRad = [4.119,4.191,4.301,5.228,8.398;1.744,2.573,4.191,6.158,7.942;5.896,4.948,5.474,5.47400000000000,10.0660000000000;3.36500000000000,4.30100000000000,4.30100000000000,6.70000000000000,7.40200000000000;3.07900000000000,3.26900000000000,3.97600000000000,6.42400000000000,7.94200000000000;3.40000000000000,4.48600000000000,5.86900000000000,6.74300000000000,7.64500000000000];
            penicillinRad = [2.405,2.217,4.191,8.709,9.079;0.622,2.061,3.976,6.697,10.986;2.89500000000000,3.66400000000000,4.75300000000000,7.64500000000000,12.5660000000000;3.14200000000000,4.75300000000000,7.99200000000000,7.35400000000000,10.5210000000000;1.76700000000000,4.60000000000000,5.55700000000000,7.84300000000000,11.8850000000000;1.96100000000000,3.66400000000000,3.90600000000000,7.45100000000000,9.07900000000000];
            %don't use area, use radius:
            doxyRad = sqrt(doxyRad/pi);
            penicillinRad = sqrt(penicillinRad/pi);
                %use diameters:
            %doxyRad = doxyRad*2;
            %penicillinRad = penicillinRad*2;
        else
            doxyRad = [2.80000000000000,2.40000000000000,2.40000000000000,2.10000000000000,1.60000000000000,2.35000000000000;3.40000000000000,2.70000000000000,2.70000000000000,2.40000000000000,2.60000000000000,2.55000000000000;4,3.40000000000000,3.40000000000000,3.70000000000000,2.40000000000000,2.90000000000000;4.80000000000000,4.40000000000000,4,4,3.60000000000000,3.05000000000000;5.60000000000000,4.80000000000000,4.80000000000000,4.50000000000000,4.80000000000000,3.83000000000000]';
            penicillinRad = [2,2.20000000000000,1.50000000000000,1.50000000000000,2,2.58000000000000;2.40000000000000,2.40000000000000,2.10000000000000,1.80000000000000,2.10000000000000,2.90000000000000;2.70000000000000,3.40000000000000,2.30000000000000,2,2.50000000000000,3.75000000000000;3.60000000000000,4.20000000000000,2.90000000000000,3.30000000000000,2.90000000000000,4.10000000000000;4.20000000000000,3,3.50000000000000,3.80000000000000,3.30000000000000,4.80000000000000]';
        end

        fineRadii = 0:0.01:max(penicillinRad(:));        
        figure(2*(year-1)+1)
    
        r = semilogx(Dosages(:),penicillinRad(:),'.k','markersize',20,'color',[1 1 1]*0.7);
        hold on
        p = semilogx(dosages,mean(penicillinRad),'.k','markersize',40);
        xlabel('multiple of clinical dose')
        ylabel('size of inhibition zone')
        ylim([0 6])
    
        PenModel = fitnlm(mean(penicillinRad),dosages,model,guess)
    
        %PenModel = fitnlm(penicillinRad(:),allDosages,model,guess)

        disp('pen AIC...')
        PenModel.ModelCriterion.AIC
        AICpen{mdl} = PenModel.ModelCriterion.AIC;
        adjR2{mdl} = PenModel.Rsquared.Adjusted;      
    
        plot(PenModel.feval(fineRadii),fineRadii,'-k')

    end
end


%%
