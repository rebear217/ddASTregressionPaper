function output = timesolvePanel(input)

    calc = 1;
    T = 24;

    if nargin > 0
        if isstruct(input)
            A0 = input.A0;
            B = input.B;
            R = input.R;
            S = input.S;
            A = input.A;
            
            popDensity = input.popDensity;
            T = input.T;
            N = input.N;
            ic = input.ic;
            parameters = input.parameters;
            Nframes = input.Nframes;
            option = input.option;
    
            dh = 1/(N-1);
            rspace = 0:dh:1;
            calc = 0;
            output = input;
        else
            T = input;
        end
    end

    close all
    clc
    
    if calc
        option = 1;
        Nframes = 25;
        parameters = setGlobalParameters();
        N = parameters.N;
        dh = 1/(N-1);
        rspace = 0:dh:1;
        [ic,parameters] = setInitialData(parameters,option);
        
        A0 = 0.68;
        
        rspace = rspace(:);

        disp(A0)
        initialCondition = ic;
        initialCondition(2*N+1:3*N) = initialCondition(2*N+1:3*N) * A0;
        
        outState = runTheModel(initialCondition,parameters,Nframes,T);
        S = outState(1:N);
        R = outState(N+1:2*N);
        A = outState(2*N+1:3*N);
        C = outState(3*N+1:4*N);
        B = S+R;
        popDensity = 2*pi*sum(B(:).*rspace)*dh;
    
        output.A0 = A0;
        output.A = A;
        output.B = B;
        output.R = R;
        output.S = S;
        output.C = C;
        output.popDensity = popDensity;
        output.T = T;
        output.N = N;
        output.ic = ic;
        output.parameters = parameters;
        output.Nframes = Nframes;
        output.option = option;
    end

    figure(2)
    set(2,'pos',[686         639        450         393])

    plot(flipud(S),'-','DisplayName',['S @dose ',num2str(A0,2),'\mug/mL, T = ',num2str(T),'h'],'LineWidth',3);
    hold on
    plot(flipud(R),'-','DisplayName','R','LineWidth',3);            
    plot(flipud(A/max(A)),'--k','DisplayName','relative antibiotic concentration','LineWidth',3);            
    plot(flipud(C/max(C)),'-','color',0.75*[0.5 1 0.5],...
        'DisplayName','relative nutrient concentration','LineWidth',3);            
    
    axis tight
    legend('boxoff')
    legend('location','northeast')
    xlabel('space')
    ylim([0 1])    
    set(gca,'Xtick',[1,N])
    set(gca,'Xticklabels',{'drug source','least drug'})
    ylabel('population density (OD)')       
    exportgraphics(gcf,['./figures/spatialModelDoseResponse2-',num2str(T),'.PDF'])
end

function outState = runTheModel(initialCondition,parameters,Nframes,T)
    outState = solveForward(initialCondition,0.1,parameters);    
    for j = 1:Nframes
        outState = solveForward(outState,T,parameters);
    end
end