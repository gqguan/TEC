%% 采用TEC吸热降温，计算给定传热量和冷热侧温度的TEC所需电功
%
% by Dr. Guan Guoqiang @ SCUT on 2022-05-04
function [E,Q,diffQ] = CalcTECPower(opStr,Q,Th,Tc,TEC,opts)
    % 函数调试用
    if nargin == 0
        opStr = 'cooling';
        Q = 23.34; % [W]
        Th = 318.15; % [K]
        Tc = 278.15; % [K]
        load('TEC_Params.mat', 'TEC_Params')
        TEC = TEC_Params.TEC(14,:);
        opts = [1,0];
    end
    switch opStr
        case('cooling')
            x0 = 2; lb = 0.1; ub = 25;
            solOpts = optimoptions(@lsqnonlin, 'Display', 'none');
            fun = @(x)GetTECHeat(x,opStr,Th,Tc,TEC,opts)-Q;
            x1 = lsqnonlin(fun,x0,lb,ub,solOpts);
            [Q,E] = GetTECHeat(x1,opStr,Th,Tc,TEC,opts);
            diffQ = fun(x1);
            if diffQ < -1e-2
                strIU = {'电流','电压'};
                strUnit = {'A','V'};
                warning('TEC功率不足，当前TEC输入%s为%.4g[%s]，%s传热量为%.4g[W]', ...
                    strIU{opts(2)+1},x1,strUnit{opts(2)+1},opStr,Q);
                xAxis = linspace(lb,x1*1.5);
                yAxis = arrayfun(@(x)GetTECHeat(x,opStr,Th,Tc,TEC,opts),xAxis);
                plot(xAxis,yAxis)
            end
        case('heating')
        otherwise
    end

    function [Qout,E] = GetTECHeat(var,opStr,Th,Tc,TEC,opts)
        switch opts(2)
            case(0)
                TEC.Current = var;
            case(1)
                TEC.Voltage = var;
            otherwise
                error('CalcTECHeat()输入参数opts有误！')
        end
        TECQ = TE_Heat(Th,Tc,TEC,opts(1),opts(2));
        E = TECQ(1)-TECQ(2);
        switch opStr
            case('cooling')
                Qout = TECQ(2);
            case('heating')
                Qout = TECQ(1);
            otherwise
                error('CalcTECHeat()输入参数opStr有误！')
        end
    end

end