function thresholdedZOIdoseResponse(outputTcell,threshold,sampleA0,fitFlag)
    COL1 = [0,1,1]*0.75;
    COL2 = [1,1,0]*0.75;

    L = length(outputTcell);
    sampleZOI = zeros(size(sampleA0));
    weight = @(y)y.^(-2);
    %weight = @(y)ones(size(y));
    pickFitData = L;

    lw = 2;
    if nargin < 4
        fitFlag = 1;
        lw = 6;
    end

    for j = 1:L
        if L > 1
            s = (j-1)/(L-1);
            COL = s*COL1 + (1-s)*COL2;
        else
            COL = (COL1 + COL2) / 2;
        end

        outputT = outputTcell{j};

        if nargin < 2
            threshold = 0.1;
        end
    
        %outputT.B holds S+R (bacterial population density)
        B = outputT.Bs;
    
        [N,M] = size(B);
        A0 = outputT.A0s;
        ZOI = zeros(size(A0));
    
        % rspace variable may not be defined for older data:
        try
            rspace = outputT.rspace;
        catch
            dh = 1/(outputT.N-1);
            rspace = 0:dh:1;
        end
    
        for n = 1:N
            density = B(n,:);
            Fn = find(density > threshold,1,'last');
            if ~isempty(Fn)
                ZOI(n) = 1-rspace(Fn);
            else
                ZOI(n) = NaN;
            end

            %plot(rspace,density)
            %plot(rspace(Fn),0,'ok')
            %hold on
            %pause

        end
    
        if j == 1
            DN = ['CARS ',num2str(outputT.T),'h dose response'];
        else
            DN = [num2str(outputT.T),'h'];
        end
            
        semilogx(A0,ZOI,'-','color',COL,'linewidth',2,...
            'DisplayName',DN,'linewidth',lw)
        hold on

        if j == pickFitData && fitFlag
            for k = 1:length(sampleA0)
                [~,I] = min(abs(A0-sampleA0(k)));
                sampleZOI(k) = ZOI(I);
            end
            semilogx(sampleA0,sampleZOI,'ok','markersize',14,'DisplayName','sampled CARS data')
    
            p0 = [0.01 -15 1];
            colr = {'b','r','k'};
            labels = {'ExpInt','Radical Exp','Bonev'};

            for Zj = 1:3
                if Zj == 3
                    p0 = p0(1:2);
                end
                fitZ = Zfunction(Zj);
                fit = fitnlm(sampleZOI,sampleA0,fitZ,p0,'Weights',weight)
                rPlot = 0:0.01:max(1.1*sampleZOI);
                adjR2 = fit.Rsquared.Adjusted;
                plot(fit.feval(rPlot),rPlot,'-k','DisplayName',...
                    [labels{Zj},' (adj R^2\approx ',num2str(adjR2,3),')'], ...
                    'color',colr{Zj})
            end
            fitLM = fitlm(sampleZOI,log10(sampleA0),'Weights',weight(sampleA0))
            adjR2 = fitLM.Rsquared.Adjusted;
            plot(10.^fitLM.feval(rPlot),rPlot,'-','DisplayName',...
                ['log-linear (adj R^2\approx ',num2str(adjR2,3),')'], ...
                'color',[1 1 1]*0.7)
        end

    end

    xlabel('antibiotic supply dose (\mug/mL)')
    ylabel('ZoI (sub-threshold radius fraction)')
    axis tight
    legend('location','southeast')

end

function f = Zfunction(flag)

    switch flag 
        case 1
            % expint
            f = @(p,r) p(1) + p(3) ./ expint(abs(p(2))*r.^2);
        case 2
            % radical exp
            f = @(p,r) p(1) + p(3)*exp(sqrt(1+abs(p(2)).*r.^2 )) .* ( -1 + sqrt(1+abs(p(2))*r.^2) );
        case 3
            % Bonev
            f = @(p,r) p(1)*exp(-r.^2*p(2));    
        otherwise
            f = @(p,r) p(1)*exp(-r.^2*p(2));    
            
    end

end

