function output = ZOIdoseReponse(input)

    calc = 1;
    if nargin > 0
        A0s = input.A0s;
        Bs = input.Bs;
        Rs = input.Rs;
        Ss = input.Ss;
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
    end

    close all
    clc
    
    if calc
        option = 1;
        Nframes = 25;
        T = 20;
        parameters = setGlobalParameters();
        N = parameters.N;
        dh = 1/(N-1);
        rspace = 0:dh:1;
        [ic,parameters] = setInitialData(parameters,option);
        A0s = 2.^(-9:0.05:1.8);
        %A0s = 2.^(-20:0.05:-9);
        
        popDensity = zeros(size(A0s));
        Bs = zeros(length(A0s),N);
        Ss = zeros(length(A0s),N);
        Rs = zeros(length(A0s),N);
        rspace = rspace(:);
        parfor j = 1:length(A0s)
            A0 = A0s(j);
            disp(A0)
            initialCondition = ic;
            initialCondition(2*N+1:3*N) = initialCondition(2*N+1:3*N) * A0;
            
            outState = runTheModel(initialCondition,parameters,Nframes,T);
            S = outState(1:N);
            R = outState(N+1:2*N);
            %A = outState(2*N+1:3*N);
            %C = outState(3*N+1:4*N);
            B = S+R;
            popDensity(j) = 2*pi*sum(B(:).*rspace)*dh;
            Bs(j,:) = B;
            Ss(j,:) = S;
            Rs(j,:) = R;        
        end
    
        output.A0s = A0s;
        output.Bs = Bs;
        output.Rs = Rs;
        output.Ss = Ss;
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

    figure(1)
    set(1,'pos',[686         639        910         393])
    
    subplot(1,2,1)
    semilogx(A0s,popDensity,'-k','LineWidth',3,'DisplayName',['dose response @T=',num2str(T),'h'])
    xlabel('antibiotic dose (\mug/mL)')
    ylabel('total bacterial density (OD)')
    axis tight
    hold on
    C = colororder;
    
    for j = 1:3
        J = 34*(j+3);
        col = C(j,:);
    
        subplot(1,2,2)
        plot(fliplr(Bs(J,:)),'color',col,'DisplayName',[num2str(j),': dose ',num2str(A0s(J),2),'\mug/mL'],'LineWidth',2);        
        hold on

        legend('boxoff');
        legend('location','northwest')
        xlabel('space')
        ylim([0 0.8])    
        set(gca,'Xtick',[1,N])
        set(gca,'Xticklabels',{'drug source','least drug'})
        ylabel('population density (OD)')    
    
        subplot(1,2,1)
        plot(A0s(J),popDensity(J),'.','markersize',44,'color',col,'HandleVisibility','off')
        text(A0s(J),popDensity(J)*1.2,num2str(j))    
    end

    subplot(1,2,1)
    yl = ylim;
    ylim([0 yl(2)])
    set(gca,'Xtick',[0.01,0.1,1,10])
    legend('boxoff')

    %{
    subplot(1,3,2)    
    semilogx(A0s,popDensity,'-k','LineWidth',3,'DisplayName','dose response zoom')
    hold on
    xlim([0.35 1])
    xlabel('antibiotic dose (\mug/mL)')
    ylabel('total bacterial density (OD)')
    j = 2;
    J = 34*(j+3);
    col = C(j,:);
    plot(A0s(J),popDensity(J),'.','markersize',64,'color',col,'HandleVisibility','off')
    text(A0s(J),popDensity(J)+0.008,num2str(j))
    legend('boxoff')
    %}

    exportgraphics(gcf,['./figures/spatialModelDoseResponse',num2str(T),'.PDF'])

    figure(2)
    j = 2;
    J = 34*(j+3);
    plot(fliplr(Ss(J,:)),'-','DisplayName',['S @dose ',num2str(A0s(J),2),'\mug/mL, T = ',num2str(T),'h'],'LineWidth',3);
    hold on
    plot(fliplr(Rs(J,:)),'-','DisplayName','R','LineWidth',3);            
    
    legend('boxoff')
    legend('location','northwest')
    xlabel('space')
    ylim([0 0.8])    
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
