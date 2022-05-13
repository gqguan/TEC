%% 寻找DCMD膜组件中集成半导体热泵的设定参数，使TEC传热量满足指定值Q
% 输入参数config为'feedTEHP'时，料液侧集成TEC从渗透液中吸热Q，加热膜组件中的料液
%
% by Dr. Guan Guoqiang @ SCUT on 2022/5/12

function [TECs,profile,flag] = CalcTEHP(config,Q,sIn,TECs,TEXs,membrane,flowPattern,opts)

    % 添加\Common目录以便调用自定义公用函数
    homePath = cd;
    idxPath = strfind(homePath,'\');
    addpath([homePath(1:idxPath(2)),'Common\'])
    % 日志文件
    logFilename = [homePath,'\',datestr(datetime,'yyyymmdd'),'.log'];
    logger = log4m.getLogger(logFilename);
    logger.setFilename(logFilename);
    % 函数测试用数据
    if nargin == 0
        [~,s1,membrane] = InitStruct();
        T1 = 273.15+50; T2 = 273.15+45; % [K]
        W1 = 1.217e-2; W2 = 6.389e-3; % [kg/s]
        config = 'permTEHP'; 
        s1.Temp = T1;
        s1.MassFlow = W1;
        sIn(1) = DCMD_PackStream(s1); % 膜组件料液侧进料
        s2 = s1;
        s2.Temp = T2;
        s2.MassFlow = W2;
        sIn(2) = DCMD_PackStream(s2); % 膜组件渗透侧进料
        % set properties for all TECs
        load('TEC_Params.mat','TEC_Params') % 载入已有的TEC计算参数
        switch config
            case('feedTEHP')
                Q = 3.5; % 料液侧集成TEC的吸热量
                iTEC1 = 9;
                TECs(1) = TEC_Params.TEC(iTEC1);
                % TEC计算参数
                opts = [TEC_Params.opt1(iTEC1),TEC_Params.opt2(iTEC1)];
                iTEC2 = 1; % 近似绝热
                TECs(2) = TEC_Params.TEC(iTEC2); 
                % 边界温度
                TEXs(1) = T2; % 即TECs(1)的冷侧温度初设为渗透液进膜组件温度
                TEXs(2) = 298.15;
            case('permTEHP')
                Q = 13.5; % 渗透侧集成TEC的放热量
                iTEC1 = 1;
                TECs(1) = TEC_Params.TEC(iTEC1);
                iTEC2 = 9; % 近似绝热
                TECs(2) = TEC_Params.TEC(iTEC2); 
                % TEC计算参数
                opts = [TEC_Params.opt1(iTEC2),TEC_Params.opt2(iTEC2)];
                % 边界温度
                TEXs(1) = 298.15;
                TEXs(2) = T1;
        end
        % 流型
        flowPattern = "countercurrent";
    end
    % 参数
    strIU = {'Current','Voltage'};
    x0 = 1.2; lb = 0.5; ub = 10;
    solOpts = optimoptions(@lsqnonlin,'Display','none');

    switch config
        case 'feedTEHP'
            DiffQ2 = @(x)CalcQ2(x,sIn,TECs,TEXs,membrane,flowPattern,opts)-Q;
            % 求TEC设定参数使min(TEC吸热量-指定值Q)
            x1 = lsqnonlin(DiffQ2,x0,lb,ub,solOpts);
            TECs(1).(strIU{opts(2)+1}) = x1;
            [Q2,profile] = CalcQ2(x1,sIn,TECs,TEXs,membrane,flowPattern,opts);
            if Q-Q2 > 1e-4
                flag = 1; % 给定TEC(1)制冷量低于渗透液冷却的冷量要求
            elseif Q-Q2 < 1e-4
                flag = 2; % 给定TEC(1)制冷量超过渗透液冷却的冷量要求
            else
                flag = 0;
            end
        case 'permTEHP'
            DiffQ1 = @(x)CalcQ1(x,sIn,TECs,TEXs,membrane,flowPattern,opts)-Q;
            % 求TEC设定参数使min(TEC吸热量-指定值Q)
            x2 = lsqnonlin(DiffQ1,x0,lb,ub,solOpts);
            TECs(2).(strIU{opts(2)+1}) = x2;
            [~,profile] = CalcQ1(x2,sIn,TECs,TEXs,membrane,flowPattern,opts);
    end

    function [Q2,profile] = CalcQ2(xVal,sIn,TECs,TEXs,membrane,flowPattern,opts)
        % 根据opts(2)设定TECs(1)
        TECs(1).(strIU{opts(2)+1}) = xVal;
        % 计算集成半导体热泵DCMD膜组件中的侧形
        [profile,~] = TEHPiDCMD(sIn,TECs,TEXs,membrane,flowPattern,opts);
        % 料液侧集成半导体热泵的吸热量
        Q2 = sum(cellfun(@(x)x(1,2),profile.QTEC));
    end

    function [Q1,profile] = CalcQ1(xVal,sIn,TECs,TEXs,membrane,flowPattern,opts)
        % 根据opts(2)设定TECs(2)
        TECs(2).(strIU{opts(2)+1}) = xVal;
        % 计算集成半导体热泵DCMD膜组件中的侧形
        [profile,~] = TEHPiDCMD(sIn,TECs,TEXs,membrane,flowPattern,opts);
        % 料液侧集成半导体热泵的吸热量
        Q1 = sum(cellfun(@(x)x(2,1),profile.QTEC));
    end

end