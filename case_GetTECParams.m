%% �Ӱ뵼������Ƭ�����ܲ���ʵ���������TEC����
%
% by Dr. GUAN Guoqiang @ SCUT on 2020-03-27
%
%  References
%  [1] Xuan X C, et al. Cryogenics, 2002, 42: 273-278.
%  [2] Huang B J, et al. International Journal of Refrigeration 2000, 23(3): 208-218.
%
%% ��ʼ��
clear global
% 
if exist('TEC','var')
    keepTEC = input('��[1]/��[0]�ù����ռ��TEC���������Ż���');
    switch keepTEC
        case(0)
            clear
            % ���¶������ݽṹ
            TEC = struct('NumTC', 190, 'NumRatio', 7/12, 'GeomFactor', 2.6e-3, ...
                         'HTCoefficient', 270, 'HTArea', 40*40e-6, ...
                         'SeebeckCoefficient', [], 'ElecConductance', [], ...
                         'ThermConductance', [], 'Voltage', [], 'Current', [], ...
                         'Parameters', []);
        case(1)
            clearvars -except TEC
            TE_log(sprintf('ʹ�ù����ռ�ı���TEC���������Ż�'))
        otherwise
            error('�������ָ������')
    end
else
    clear
end

% �뵼������Ƭ�����ܲ���ʵ������
% ��ʵ�������ļ�ExpData.txt�е��룬ʵ�����ݴ��ڹ����ռ�ı����ExpData��
TE_ImportExpData
%
%% �Ż�
% ������������Ҫִ�е��Ż�����
opt = input(' 0 - Optimize (r g K) values according to https://doi.org/10.1016/S0011-2275(02)00035-8\n 1 - Optimize (a R K) values according to https://doi.org/10.1016/S0140-7007(99)00046-8\n Input 0 or 1 to select corresponding method to get the TEC parameters: ');
TE_log(sprintf('Set opt = %d', opt))
% �趨�Ż������ĳ�ֵ
switch opt
    case(0) % �Ż�r��gֵ�����ο�����[1]
        TE_log('Getting TEC type according to TEC.NumRatio');
        if TEC.NumRatio == 0
            TE_log('Given TEC type is one stage')
            TEC.NumRatio = 0;
            x0 = TEC.GeomFactor;
        else
            TE_log('Given TEC type is two stages')
            x0 = [TEC.NumRatio,TEC.GeomFactor,TEC.HTCoefficient];            
        end
    case(1) % �Ż�(a R K)ֵ�����ο�����[2]
        opt2a = input('Input polynomial order to correlate the (a R K) values: ');
        TE_log(sprintf('Using %d-order polynomial function to correlate (a R K)', opt2a), 0, 1)
        x0 = ones(3, opt2a+1);
        x0(:,2:(opt2a+1)) = 0;
    otherwise
        prompt = sprintf('Unknown running mode of %d in case_GetTECParams.m', opt);
        TE_log(prompt, 1);
        return
end
% �����Ż�����
options = optimset('PlotFcns', @optimplotfval, 'MaxFunEvals', 1000*length(x0));
% ����Ŀ�꺯��
fun = @(x)(TE_RMSE(x, TEC, ExpData, opt));
% ����Ż�����
[x,fval,exitflag,optim_output] = fminsearch(fun, x0, options);
TE_log(sprintf('RMSE(exp vs. sim) = %.3f', fval))
switch exitflag
    case(1)
        TE_log(sprintf('Success to regress parameters after %d iterations', optim_output.iterations))
        switch opt
            case(0)
                TE_log('Check TEC properties for (r g k)')
            case(1)
                TE_log('Check TEC.Parameters to correlate (a R K)')
        end
    case(0)
        TE_log('Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.MaxFunEvals', 2)
    case(-1)
        TE_log('The algorithm was terminated by the output function', 2)
end
%
%% ������
[~,output] = fun(x);
% ����TEC������
output.pid = input('Input TEC part no.: ', 's');
% �����Ż�����
output.opt = opt;
% ��־
global TE_LogData
output.log = TE_LogData;
% �����
current_tab = struct2table(output, 'AsArray', 1);
% ��ǰĿ¼����TEC�����ļ�ʱ�����TEC_Params
if exist('TEC_Params.mat', 'file') == 2
    load('TEC_Params.mat')
    TEC_Params = [TEC_Params;current_tab];
else
    TEC_Params = current_tab;
end
% �������
save('TEC_Params.mat', 'TEC_Params')