%% ģ��һά�����ȱ�DCMDĤ����еĴ��Ⱥʹ������󣺴�ͳDCMDϵͳ�ĵ�λ�ܺ�
% ��Ĥ������ɰ뵼���ȱ�TEHPʱ����������Ӧ��TEC������������TEC_Params.mat�е�H28��H05��
% ��δ����TEHPʱ��������TEC_Params.mat�е�H00�����ƾ��ȵı߽��������޵繦�����ҵ���ϵ����С��
%
% by Dr. Guan Guoqiang @ SCUT on 2022/05/02
%
% ����̬���������£���Һ������������Q1ΪʹĤ�����Һ������¶����ߵ�Ĥ��������¶ȣ�
% ��͸Һ��ȴ��������Q2ΪĤ�����͸������¶Ƚ��͵�Ĥ��������¶ȡ�
% �ɴˣ���Һ�������ܺİ��繦100%ת��Ϊ�������㣬��ΪE1 = Q1��
% ��͸Һ��ȴ���ܺİ�TEC����ĵ������㣬��ΪE2 = Q2/COP2
%
% ������͵�λ�ܺĵ�ȫ������̬��������Һ������Q1 = QM+WF*CP*(TMH-T0)����ϸ�Ƶ����ʼ�2022/5/3��
% ��͸Һ��ȴ������Q2 = QM+WP*CP*(TMC-TP2)
%
% Ϊȷ��DCMDϵͳ��������Q1��Q2����Ҫ��ÿ�Ĥ����������͸����QM��WP���Ļ����Ͻ�һ��ȷ��Ĥ���¶ȣ�TMH��TMC��
% ����1���ֱ��г�WP��QM�ļ���ʽ��������TMH��TMC��
% ����2����ģ�����е��¶ȷֲ�����T(3)��T(4)����ƽ��ֵ

%% ��ʼ��
clear
% ���ù��ñ������壬���а���DuctGeom���������γߴ磩��Stream�����϶��壩��MembrProps��Ĥ�������ʣ�
CommonDef
% �趨����TEC�༶SFMDϵͳ�ļ���
NumStage = 1;
% �趨Ĥ������ȡ��������¶Ⱥ�����
T1 = 323.15; T2 = 288.15; % [K]
W1 = 1.5e-4; W2 = 1.5e-5; % [kg/s]
% �趨�����¶�
TEXs = [298.15,298.15]; % [K]
% Ĥ����Ȳ���ϳ�ʼ��
s1 = Stream;
s1.Temp = T1;
s1.MassFlow = W1;
% calculate the rest properties of feed-side influent
sIn(1) = DCMD_PackStream(s1);
% Ĥ��������ϳ�ʼ��
s2 = Stream;
s2.Temp = T2;
s2.MassFlow = W2;
% calculate the rest properties of permeate-side influent
sIn(2) = DCMD_PackStream(s2);
% �趨����Ĥ��������
membrane = MembrProps;
% set properties for all TECs
load('TEC_Params.mat') % �������е�TEC�������
% ע�ⰴopt1=0,opt2=1����TEC����������
% opts = [0,1]; TECs(1:(NumStage+1)) = TEC_Params.TEC(18,1);
% opts = [1,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(14,1);
% opts = [0,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(4,1);
opts = [0,0]; TECs(1:(NumStage+1)) = TEC_Params.TEC(1,1); % �൱��δ���ɰ뵼���ȱõ�DCMDĤ���


%% ���㼯���ȱ�DCMDĤ����е��¶ȷֲ�
[profile1,sOut1] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"countercurrent",opts);
% [profile2,sOut2] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"cocurrent",opts);
opStr = 'cooling';

%% DCMDϵͳ��λ�ܺ�
% ������̬����ʱ��Һ������Q(1)����͸Һ������Q(2)
[Q,WP,~,TP1,TP2] = CalcHeat(profile1);
% ����������
E(1) = Q(1);
% ����Ĥ������ð뵼�����书��
% opts = [0,1]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(18,1);
opts = [1,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(14,1);
% opts = [0,0]; exTECs(1:(NumStage+1)) = TEC_Params.TEC(4,1);
[E(2),QTEC,QDiff] = CalcTECPower(opStr,Q(2),TEXs(1),mean([TP1,TP2]),exTECs(1),opts);
% ����ϵͳ���ܺ�
SEC = sum(E)/WP/3600/1000; % [kWh/kg]

%% ���
fprintf('DCMDϵͳ��Ĥ��������¶ȷֱ�Ϊ%.2f[K]��%.2f[K]�����ʷֱ�Ϊ%.4g[kg/s]��%.4g[kg/s]\n',T1,T2,W1,W2)
fprintf('DCMDϵͳ�ȶ�������Ҫ������%.4g[W]��������%.4g[W]\n',Q(1),Q(2))
fprintf('ϵͳ�������TEC�����%s��%s����������Ϊ%.4g[W]����ƫ��Ϊ%.4g[W]\n',TEC_Params.pid{14},opStr,QTEC,QDiff)
fprintf('DCMDϵͳ��ˮ��Ϊ%.4g[kg/h]��������������%.4g[W]����ȴ��������%.4g[W]����λ�ܺ�Ϊ%.4g[kWh/kg]\n',WP*3600,E(1),E(2),SEC)
figObj1 = DispResults(profile1,DuctGeom,membrane);
% figObj2 = DispResults(profile2,DuctGeom,membrane);

function p = DispResults(profile,DuctGeom,membrane)
    persistent iLine
    if isempty(iLine)
        iLine = 0;
    end
    iLine = iLine+1;
    % ��ˮ��
    WP = sum([profile.SM.MassFlow]); % [kg/s]
    JM = WP/membrane.Area*3600; % [kg/m2-h]
    fprintf('��ˮ��Ϊ%.4e[kg/s]����͸ͨ��Ϊ%.3f[kg/m2-h]��\n', WP, JM)
    % TEC����
    E = sum(cell2mat(cellfun(@(x)x(:,1)-x(:,2),profile.QTEC,'UniformOutput',false)),2);
    Q1 = sum(cell2mat(cellfun(@(x)x(:,1),profile.QTEC,'UniformOutput',false)),2);
    Q2 = sum(cell2mat(cellfun(@(x)x(:,2),profile.QTEC,'UniformOutput',false)),2);
    fprintf('DCMDĤ�����Һ�����TEC�繦��Ϊ%.4f[W]��������%.4f[W]���ӻ���������%.4f[W]\n', ...
        E(1), Q1(1), Q2(1))
    fprintf('DCMDĤ�����͸����ȴTEC�繦��Ϊ%.4f[W]��������%.4f[W]���򻷾�������%.4f[W]\n', ...
        E(2), Q2(2), Q1(2))
    % �����¶Ȳ���
    lineColor = {'#7E2F8E';'#A2142F';'#D95319';'#77AC30';'#4DBEEE';'#0072BD'};
    lineStyle = {'-';'--';':';'-.'};
    x = linspace(0,DuctGeom.Length,length(profile.QM));
    p = plot([x;x;x;x;x;x]',profile.T',lineStyle{iLine});
    for iPlot = 1:length(p)
        p(iPlot).Color = lineColor{iPlot};
    end
    if iLine <= 4
        hold on
    else
        hold off
        iLine = 0;
    end
end

function [Q,WP,QM,TP1,TP2] = CalcHeat(profile)
    T0 = 298.15; % �����¶�Ϊ�����¶�[K]
    QM = sum(profile.QM); % ��Ĥ������[W]
    WP = sum(arrayfun(@(x)x.MassFlow,profile.SM)); % ��Ĥ��͸��[kg/s]
    TM = mean(profile.T(3:4,:),2);
    iStart = strfind(profile.Remarks,'��');
    switch profile.Remarks(iStart+1:end)
        case('cocurrent')
            TP1 = profile.S2(1).Temp;
            TP2 = profile.S2(end).Temp;
        case('countercurrent')
            TP1 = profile.S2(end).Temp;
            TP2 = profile.S2(1).Temp;
        otherwise
            error('CalcHeat()�������profile�ֶ�Remarks������Ч��������Ϣ')
    end
    cp1 = mean([profile.S1.SpecHeat]);
    cp2 = mean([profile.S2.SpecHeat]);
    Q(1) = QM+WP*cp1*(TM(1)-T0);
    Q(2) = QM+WP*cp2*(TM(2)-TP2);
end

function [E,Q,fval] = CalcTECPower(opStr,Q,Th,Tc,TEC,opts)
    switch opStr
        case('cooling')
            x0 = 2; lb = 0.1; ub = 12.5;
            solOpts = optimoptions(@lsqnonlin, 'Display', 'none');
            fun = @(x)GetTECHeat(x,opStr,Th,Tc,TEC,opts)-Q;
            x1 = lsqnonlin(fun,x0,lb,ub,solOpts);
            fval = fun(x1);
            [Q,E] = GetTECHeat(x1,opStr,Th,Tc,TEC,opts);
        case('heating')
        otherwise
    end
end

function [Qout,E] = GetTECHeat(var,opStr,Th,Tc,TEC,opts)
    switch opts(2)
        case(0)
            TEC.Current = var;
        case(1)
            TEC.Voltage = var;
        otherwise
            error('CalcTECHeat()�������opts����')
    end
    Q = TE_Heat(Th,Tc,TEC,opts(1),opts(2));
    E = Q(1)-Q(2);
    switch opStr
        case('cooling')
            Qout = Q(2);
        case('heating')
            Qout = Q(1);
        otherwise
            error('CalcTECHeat()�������opStr����')
    end
end

% %% ����SFMD����
% % �����Ż�Ŀ�꺯��Ϊ��λ�ܺ�
% opfun = @(x)interface_SFMD(x, NumStage, SInFeeds, SInPerms, Membranes, TECs);
% % ��ֵ
% x0 = [W1 T1 I1 W2 T2 I2];
% % �����߽�
% lb = [3.5e-6, 40+273.15, 0, 3.5e-6, 5+273.15, 0];
% ub = [2.1e-3, 80+273.15, 8, 2.1e-3, 30+273.15, 8];
% % ����
% opts = optimoptions('fmincon','Display','iter');
% x = fmincon(opfun,x0,[],[],[],[],lb,ub,[],opts);
% 
% %% ������
% [SEC,out1,out2] = opfun(x);
% format short g
% fprintf('Specific energy consumption of %d-stage DCMD is %.4e J/kg (%.4g kWh/kg).\n', NumStage, SEC, SEC*2.778e-7);
% disp(out1)