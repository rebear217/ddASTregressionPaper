function [mySurface,parameters,types] = runTheModel(option)
    clc

    if nargin < 1
        option = 1;
    end
	Nframes = 2000;

    close all

    options = {'uniform initial S-data',...
        'initial S-data in no-drug region only',...
        'initially S & R everywhere',...
        'no drug at all, carbon in the centre but not cells'};

    parameters = setGlobalParameters();
    [ic,parameters] = setInitialData(parameters,option);

    outState = solveForward(ic,0.001,parameters);

    figure(1)
    set(1,'color','white')
    set(1,'position',[55   470   2*744   359])

    A0 = parameters.A0;
    S0 = parameters.S0;
    C0 = parameters.C0;
    N = parameters.N;
    X = parameters.X;
    yield = parameters.yield;

    initialC = sum(outState(3*N+1:4*N));
    currentC = initialC;

    j = 0;
    mySurface = zeros(Nframes,N);
    
    while j < Nframes && (currentC > initialC * 0.001)
        j = j + 1;
        outState = solveForward(outState,5,parameters);

        S = outState(1:N);
        R = outState(N+1:2*N);
        A = outState(2*N+1:3*N);
        C = outState(3*N+1:4*N);
        currentC = sum(C);

        subplot(1,2,1);
        if j > 1
            delete(p);
            delete(a0);
            delete(a00);
        end
        a0 = area(X,(S+R)/yield);
        set(a0,'facecolor','g','edgecolor','none');
        hold on
        a00 = area(X,R/yield);
        set(a00,'facecolor','r','edgecolor','none');
        p=plot(X,A/A0,'-k');
        if j == 1
            legend('susceptibles','resistants','antibiotic','location','northwest')
            ylabel('density/concentration (normalised)')
            set(gca,'Xtick',[0 0.5 1])
            set(gca,'Xticklabel',{'no drug','spatial location','drug'})
            axis([0 1 0 1]);
        end

        subplot(1,2,2);
        if j > 1
            delete(a1);
            delete(a2);
            delete(p3);
        end
        f = (1e-4 + R)./(1e-4 + S+R);
        a1 = plot(X,f,'-r');
        hold on
        a2 = area(X,R+S);
        set(a2,'facecolor','b');
        p3 = plot(X,C/C0,'-g');
        if j == 1
            axis([0 1 0 1]);
            ylabel('density/concentration (normalised)')
            set(gca,'Xtick',[0 0.5 1])
            set(gca,'Xticklabel',{'no drug','spatial location','drug'})
            legend('freq. resistance','total density','carbon')
        end

        title([num2str(round(1000*currentC/initialC)/10),' % of sugar remaining'])

        drawnow
        mySurface(j,:) = R+S;
    end
    mySurface = mySurface(1:j,:);
    figure(1)
    export_fig(['figures/RFP-GFP-',date,'.pdf'])
    types.R = R;
    types.S = S;   
end