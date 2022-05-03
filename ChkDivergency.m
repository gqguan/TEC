function outTF = ChkDivergency(inVal,n)
    if ~exist('n','var')
        n = 3;
    end
    outTF = false;
    persistent v
    if isempty(v)
        v = nan(1,n);
        v(1) = inVal;
    else
        idx = arrayfun(@(x)isnan(x),v);
        if sum(idx) > 1
            v(1,find(idx,1)) = inVal;
        else
            if sum(idx) == 1
                v(1,find(idx,1)) = inVal;
            else
                v(1) = []; v(n) = inVal;
            end
            [~,sortedI] = sort(v);
            if isequal(sortedI,1:n)
                outTF = true;
            end
        end
    end
