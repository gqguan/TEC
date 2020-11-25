function [SEC, Output_Stages, Output_TEHP] = interface_SFMD(vars, NumStage, SInFeeds, SInPerms, Membranes, TECs)
% vars - double array vars(1): W1
%                     vars(2): T1
%                     vars(3): W2
%                     vars(4): T2
if ~all([NumStage,length(SInFeeds),length(SInPerms),length(Membranes),length(TECs)-1])
    fprintf('[ERROR] Ӧ���뵥��SFMD�Ĳ���!\n')
    return
end
% ���빲�ñ���
CommonDef
% �������ʺ��¶�
SInFeeds(1).MassFlow = vars(1);
SInFeeds(1).Temp = vars(2);
SInFeeds = DCMD_PackStream(SInFeeds,DuctGeom);
SInPerms(1).MassFlow = vars(3);
SInPerms(1).Temp = vars(4);
SInPerms = DCMD_PackStream(SInPerms,DuctGeom);
% ���������
[~, Output_Stages, Output_TEHP] = SFMD(NumStage, SInFeeds, SInPerms, Membranes, TECs);
WP = MembrProps(1).Area*Output_Stages.JM;
% ������Һ���ȵ��趨�¶ȵ��ܺ�
Q1S = SInFeeds(1).MassFlow*SInFeeds(1).SpecHeat*(SInFeeds(1).Temp-298.15);
E1 = Q1S*1; % �����100%Ч�ʵ��ܺ�
% ������͸Һ�ﵽ�趨�¶ȵ��ܺ�
Q2S = SInPerms(1).MassFlow*SInPerms(1).SpecHeat*(SInPerms(1).Temp-298.15);
if Q2S >= 0
    E2 = Q2S; % �����100%Ч�ʵ��ܺ�
else
    eta_c = 1-SInPerms(1).Temp/298.15;
    E2 = -Q2S*eta_c/(1-eta_c); % ����ŵѭ��������ܺ�
end
SEC = (E1+E2)/WP;