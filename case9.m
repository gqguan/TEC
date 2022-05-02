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
T1 = 323.15; T2 = 288.15;
W1 = 1.5e-4; W2 = 1.5e-4;
% �趨�����¶�
TEXs = [298.15,298.15];
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
% [profile2,sOut2] = TEHPiDCMD(sIn,TECs,TEXs,membrane,"cocurrent");

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