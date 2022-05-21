%% 生成case编号便于大规模计算中识别
function caseNo = GenerateCaseNo(name,planMat,varName)
    % 输入参数检查
    argin2OK = ChkArgType({'char','double'},name,planMat);
    if argin2OK
        caseNo = cell(size(planMat,1),1);
        % 检查实验计划考查的因素数目是否与varargin一致
        n1 = size(planMat,2);
        if ~exist('varName','var')
            varName = arrayfun(@(x)char(x),65:65+n1-1,'UniformOutput',false);
        end
        n2 = length(varName);
        if  ~isequal(n1,n2)
            error('输入因素名称变量数目与实验计划不符！')
        end
    else
        error('GenerateCaseNo()输入参数name和planMat有误')
    end
    % 第i个case编号：name+sum_{j=1}^{n1}(varName{j}+planMat(i,j))
    for i = 1:length(caseNo)
        str = '';
        for j = 1:n1
            str = [str,varName{j},'(',num2str(planMat(i,j)),')'];
        end
        caseNo{i} = [name,'(',num2str(i),')','|',str];
    end 

end