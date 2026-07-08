function P = plotshadedYX(x,y,fstr)
    % x: x coordinates
    % y: either just one y vector, or 2xN or 3xN matrix of y-data
    % fstr: format ('r' or 'b--' etc)
    %
    % example
    % x=[-10:.1:10];plotshaded(x,[sin(x.*1.1)+1;sin(x*.9)-1],'r');

    if size(y,1)>size(y,2)
        y=y';
    end
    
    py = [y,fliplr(y)]; % make closed patch
    px = [x(1,:), fliplr(x(2,:))];
    
    P = patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');

    %alpha(.2); % make patch transparent

end