%% 判定n次调用每次输入的数值是否增加
% 功能说明：比较n次调用每次输入的数值，若数值每次比上一次输入数值相等或更大则判定为发散
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
        % 需要调用n次，采集n个数值；若调用次数小于n，将v数组中有元素为nan
        idx = arrayfun(@(x)isnan(x),v);
        if sum(idx) > 1
            v(1,find(idx,1)) = inVal;
        else
            if sum(idx) == 1
                v(1,find(idx,1)) = inVal;
            else % 已采集n个数据，则删除第一个元素后将输入数值添加到数组最后一个元素中
                v(1) = []; v(n) = inVal;
            end
            [~,sortedI] = sort(v);
            if isequal(sortedI,1:n)
                outTF = true;
            end
        end
    end
