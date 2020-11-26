%% �༶�����ȱ�ƽ��Ĥ����Ԫ�Ĳ����Ż���һ��SFMD��͵�λ�ܺ�
%
% by Dr. Guan Guoqiang @ SCUT on 2020/11/25

%% ��ʼ��
clear
% ���ù��ñ������壬���а���DuctGeom���������γߴ磩��Stream�����϶��壩��MembrProps��Ĥ�������ʣ�
CommonDef
% �趨����TEC�༶SFMDϵͳ�ļ���
NumStage = 1;
% �趨Ĥ������ȡ��������¶Ⱥ�����
T1 = 323.15; T2 = 303.15;
W1 = 1.5e-4; W2 = 1.5e-4;
% Ĥ����Ȳ���ϳ�ʼ��
SInFeed = Stream;
SInFeed.Temp = T1;
SInFeed.MassFlow = W1;
% calculate the rest properties of feed-side influent
SInFeed = DCMD_PackStream(SInFeed);
% �趨����Ĥ����Ȳ����
SInFeeds(1:NumStage) = SInFeed;
% Ĥ��������ϳ�ʼ��
SInPerm = Stream;
SInPerm.Temp = T2;
SInPerm.MassFlow = W2;
% calculate the rest properties of permeate-side influent
SInPerm = DCMD_PackStream(SInPerm);
% �趨����Ĥ���������
SInPerms(1:NumStage) = SInPerm;
% �趨����Ĥ��������
Membranes(1:NumStage) = MembrProps;
% set properties for all TECs
load('TEC_Params.mat') % �������е�TEC�������
TECs(1:(NumStage+1)) = TEC_Params.TEC(3,1); % ע�ⰴopt=0����TEC����������
I1 = TECs(1).Current;
I2 = TECs(2).Current;

%% ����SFMD����
% �����Ż�Ŀ�꺯��Ϊ��λ�ܺ�
opfun = @(x)interface_SFMD(x, NumStage, SInFeeds, SInPerms, Membranes, TECs);
% ��ֵ
x0 = [W1 T1 I1 W2 T2 I2];
% �����߽�
lb = [3.5e-6, 40+273.15, 0, 3.5e-6, 5+273.15, 0];
ub = [2.1e-3, 80+273.15, 8, 2.1e-3, 30+273.15, 0];
% ����
x = fmincon(opfun,x0,[],[],[],[],lb,ub);

%% ������
[SEC,out1,out2] = opfun(x);
format short g
fprintf('Specific energy consumption of %d-stage DCMD is %.4e J/kg (%.4g kWh/kg).\n', NumStage, SEC, SEC*2.778e-7);
disp(out1)