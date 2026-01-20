function plotWave(mySurface,parameters)
    clc
    X = -(parameters.X);
    [N,~] = size(mySurface); 
    times = 1:N;

    figure(1)
    surf(times,X,mySurface','EdgeColor','none','FaceColor','interp');
    set(gcf,'position',[83          63        1455         916])
    view(3)
    box on
    
    set(gca,'Ytick',[-1 -0.5 0])
    set(gca,'Yticklabel',{'drug source','spatial coord.','least drug'},'fontsize',28)
	set(gca,'Xtick',[1 floor(N/2) N])
    set(gca,'Xticklabel',{'early','time','later'})
    zlabel('bacterial density')
    
    axis tight
    camlight left
    %lighting phong
    lighting gouraud
    %axis off
    grid off
    box on
    alpha(1)
    
    hold on
    view(-70,45);
    
    step = 10;
	thickness = 6;

    J = 1;
    for j = 1:step:N
        c = [1 1 1]*0.9 * ((N-j)/N);
        c = [1 1 1]*0.7;
        plot3(times(j)*ones(size(X)),X,mySurface(j,:),'-','linewidth',4,'color',c);
        J = J + 1;
    end
	plot3(times,X(1)*ones(size(times)),mySurface(:,1),'-k');

    %title({'{\bf Spatial stucture and inverted dose-response curves}',...
    %    'Growth of drug-susceptible bacteria creates a wave of resistant phenotypes'})

    figure(2)
    p = {};
    J = 1;
    for j = 1:step:N
        y = mySurface(j,:);
        c = [1 1 1]*0.9 * ((N-j)/N);
        p{J} = plot(X,y,'-','linewidth',(thickness/2)*(j/N),'color',c);
        if j == 1
            hold on
            set(p{J},'color','b')
            set(p{J},'linewidth',2)
            set(gca,'Xtick',[-1 -0.5 0])
            set(gca,'Xticklabel',{'drug source','spatial coord.','least drug'})
        end
        J = J + 1;
    end
    
	set(p{end},'color','k')
    set(p{end},'linewidth',thickness)
	legend([p{1} p{2} p{floor(J/2)} p{end}], ...
        {'initial density','early transients','later transients','final density'},...
        'location','northwest')
    axis tight
    ax = axis;

    axis([min(X) max(X) 0 ax(4)])
    ylabel('bacterial density')
        
    figure(1)
    export_fig(['figures/wave-',date,'.png'])
    figure(2)
    export_fig(['figures/movingWave-',date,'.pdf'])
    
end