function radialPlot(x0,v,colour)
    if nargin < 3
        colour = [1 1 0];
    end
    v = v(:)';
    v = fliplr(v);
    
    m = length(v);

    space = x0*(1:m)/m;
    nTH = 50;
    [R,TH] = meshgrid(space,0:(2*pi/nTH):2*pi);
    X = R.*cos(TH);
    Y = R.*sin(TH);

    Vr = repmat(v,nTH+1,1);
    surf(X,Y,Vr,'linestyle','none','facecolor',colour);
    box on
    grid off
    view([-60 64]);
    %colormap(autumn)
    lighting phong
    axis tight
    xlabel('spatial coordinate')
    ylabel('spatial coordinate')
    zlabel('density')
    
    light('col',colour,'pos',[10 0 10]);
end