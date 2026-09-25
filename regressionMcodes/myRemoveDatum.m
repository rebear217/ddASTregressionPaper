function removedData = myRemoveDatum(data,removeDatum)

    % this is similar to setdiff(data,removeDatum) where removeDatum is a
    % number and found in data, but we want to ensure at most 1 datum is removed
    % e.g. myRemoveDatum([1,1,2],1) = NaN
    % e.g. myRemoveDatum([1,1,2],2) = [1,1]
    % e.g. myRemoveDatum([1,2,3],2) = [1,3]
    
    % Unlike setdiff, this does not re-order data:

    % Use a Boolean flag to avoid nested if statements:
    resolved = false;

    locate = data - removeDatum;
    L = locate == 0;

    if ~any(L)
        %we cannot remove removeDatum as it is not inside the data vector:
        removedData = data;
        resolved = true;
    end

    if ~resolved && sum(L) > 1
        %there is more than 1 removeDatum inside the data vector
        %this is not desired behaviour, so we throw a problem:
        removedData = NaN;            
        resolved = true;
    end

    if ~resolved
        %there is exactly 1 datum to remove, so let's remove it:
        resolved = true;
        J = find(L);
        N = length(L);
        if N == 1
            removedData = [];
        else
            switch J
                case 1
                    removedData = data(2:end);
                case N
                    removedData = data(1:end-1);
                otherwise
                    removedData = data([1:J-1 J+1:N]);
            end
        end
    end

end