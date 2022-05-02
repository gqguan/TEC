%% ģ��һά�����ȱ�DCMDĤ����еĴ��Ⱥʹ�������
%
% by Dr. Guan Guoqiang @ SCUT on 2022/05/02

%% ��ʼ��
clear
% ���ù��ñ������壬���а���DuctGeom���������γߴ磩��Stream�����϶��壩��MembrProps��Ĥ�������ʣ�
CommonDef
% �趨����TEC�༶SFMDϵͳ�ļ���
NumStage = 1;
% �趨Ĥ������ȡ��������¶Ⱥ�����
T1 = 323.15; T2 = 288.15; % [K]
W1 = 1.5e-4; W2 = 1.5e-4; % [kg/s]
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
TECs(1:(NumStage+1)) = TEC_Params.TEC(3,1); % ע�ⰴopt=0����TEC����������

%% ���㼯���ȱ�DCMDĤ����е��¶ȷֲ�
[profile1,sOut1] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"countercurrent");
[profile2,sOut2] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"cocurrent");

%% ���
DispResults(profile1,DuctGeom,membrane)
DispResults(profile2,DuctGeom,membrane)

function DispResults(profile,DuctGeom,membrane)
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