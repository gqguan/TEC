%% 从半导体制冷片的性能测试实验结果中拟合TEC参数
%
% by Dr. GUAN Guoqiang @ SCUT on 2020-03-27
%
%  References
%  [1] Xuan X C, et al. Cryogenics, 2002, 42: 273-278.
%  [2] Huang B J, et al. International Journal of Refrigeration 2000, 23(3): 208-218.
%
%% 初始化
clear;
% 数据结构定义
TEC = struct('NumTC', 190, 'NumRatio', 0.9, 'GeomFactor', 0.7e-3, ...
             'HTCoefficient', 270, 'HTArea', 0.0016, ...
             'SeebeckCoefficient', [], 'ElecConductance', [], ...
             'ThermConductance', [], 'Voltage', [], 'Current', [], ...
             'Parameters', []);
% 半导体制冷片的性能测试实验数据
% 从实验数据文件ExpData.txt中导入，实验数据存于工作空间的表变量ExpData中
TE_ImportExpData
%
%% 优化
% 命令行输入需要执行的优化方法
opt = input(' 0 - Optimize (r g) values according to https://doi.org/10.1016/S0011-2275(02)00035-8\n 1 - Optimize (a R K) values according to https://doi.org/10.1016/S0140-7007(99)00046-8\n Input 0 or 1 to select corresponding method to get the TEC parameters: ');
% 设定优化向量的初值
switch opt
    case(0) % 优化r和g值，见参考文献[1]
        opt1a = input('Input type of TEC (1) one-stage or (2) two-stage: ');
        switch opt1a
            case(1)
                TE_log('Given TEC type is one stage');
                TEC.NumRatio = 0;
                x0 = TEC.GeomFactor;
            case(2)
                TE_log('Given TEC type is two stages');
                x0 = [TEC.NumRatio,TEC.GeomFactor];
            otherwise
                TE_log('Given TEC type is unknown', 1);
                return
        end
    case(1) % 优化(a R K)值，见参考文献[2]
        opt2a = input('Input polynomial order to correlate the (a R K) values: ');
        x0 = ones(3, opt2a+1);
        x0(:,2:(opt2a+1)) = 0;
    otherwise
        prompt = sprintf('Unknown running mode of %d in case_GetTECParams.m', opt);
        TE_log(prompt, 1);
        return
end
% 定义优化参数
options = optimset('PlotFcns', @optimplotfval);
% 定义目标函数
fun = @(x)(TE_RMSE(x, TEC, ExpData, opt));
% 获得优化参数
x = fminsearch(fun, x0, options);
%
%% 输出结果
[~,output] = fun(x);
% 输入TEC部件号
output.pid = input('Input TEC part no.: ', 's');
% 参数优化方法
output.opt = opt;
% 构造表
current_tab = struct2table(output, 'AsArray', 1);
% 当前目录存在TEC参数文件时载入表TEC_Params
if exist('TEC_Params.mat', 'file') == 2
    load('TEC_Params.mat')
    TEC_Params = [TEC_Params;current_tab];
else
    TEC_Params = current_tab;
end
% 结果存盘
save('TEC_Params.mat', 'TEC_Params')