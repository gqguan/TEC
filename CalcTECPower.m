%% 采用TEC吸热降温，计算给定传热量和冷热侧温度的TEC所需电功
% 当给定传热量大于TEC能力时，自动采用多个TEC并联传热
%
% by Dr. Guan Guoqiang @ SCUT on 2022-05-04
function [E,Q,diffQ,nTEC,newTEC] = CalcTECPower(opStr,givenQ,Th,Tc,TEC,opts)
    % 日志文件
    logFilename = [cd,'\',datestr(datetime,'yyyymmdd'),'.log'];
    logger = log4m.getLogger(logFilename);
    logger.setFilename(logFilename);
    % 初始化
    % 函数调试用
    if nargin == 0
        opStr = 'cooling';
        givenQ = 23.34; % [W]
        Th = 318.15; % [K]
        Tc = 278.15; % [K]
        load('TEC_Params.mat', 'TEC_Params')
        TEC = TEC_Params.TEC(15,:);
        opts = [1,0];
    end
    nTEC = 1;
    switch opStr
        case('cooling')
            x0 = 2; lb = 0.1; ub = 25;
            solOpts = optimoptions(@lsqnonlin, 'Display', 'none');
            fun = @(x)GetTECHeat(x,opStr,Th,Tc,TEC,opts)-givenQ;
            x1 = lsqnonlin(fun,x0,lb,ub,solOpts);
            [Q,E] = GetTECHeat(x1,opStr,Th,Tc,TEC,opts);
            strIU = {'Current','Voltage'};
            strUnit = {'A','V'};
            diffQ = fun(x1);
            if diffQ < -1e-2 % TEC功率不满足需要
                % 用nTEC个TEC并联传热
                nTEC = ceil(givenQ/Q);
                fun1 = @(x)GetTECHeat(x,opStr,Th,Tc,TEC,opts)-givenQ/nTEC;
                x2 = lsqnonlin(fun1,x0,lb,ub,solOpts);
                [Q,E] = GetTECHeat(x2,opStr,Th,Tc,TEC,opts);
                msg = sprintf('给定传热量%.4g[W]大于TEC能力%.4g[W]！采用%d个TEC并联传热，每个TEC输入%s为%.4g[%s]，%s传热量为%.4g[W]，输入电功%.4g[W]\n', ...
                    givenQ,Q,nTEC,strIU{opts(2)+1},x2,strUnit{opts(2)+1},opStr,Q,E);
                logger.info('CalcTECPower',msg)
                Q = Q*nTEC;
                E = E*nTEC;
                diffQ = Q-givenQ;
%                 TE_ShowPerformance(Th,Tc,TEC,opts)
            else
                x2 = x1;
            end
            % 代入求得的TEC设定参数
            newTEC = TEC;
            newTEC.(strIU{opts(2)+1}) = x2;
        case('heating')
        otherwise
            logger.error('CalcTECPower','输入参数opStr内容与预设不符！')
    end

end