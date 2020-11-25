function [SEC, Output_Stages, Output_TEHP] = interface_SFMD(vars, NumStage, SInFeeds, SInPerms, Membranes, TECs)
% vars - double array vars(1): W1
%                     vars(2): T1
%                     vars(3): I1
%                     vars(4): W2
%                     vars(5): T2
%                     vars(6): I2
%% ��ʼ��
if ~all([NumStage,length(SInFeeds),length(SInPerms),length(Membranes),length(TECs)-1])
    fprintf('[ERROR] Ӧ���뵥��SFMD�Ĳ���!\n')
    return
end
E = zeros(1,NumStage);
% ���빲�ñ���
CommonDef
% �������ʺ��¶�
SInFeeds(1).MassFlow = vars(1);
SInFeeds(1).Temp = vars(2);
SInFeeds = DCMD_PackStream(SInFeeds,DuctGeom);
TECs(1).Current = vars(3);
SInPerms(1).MassFlow = vars(4);
SInPerms(1).Temp = vars(5);
SInPerms = DCMD_PackStream(SInPerms,DuctGeom);
TECs(2).Current = vars(6);
% �����ˮ�ʡ�TEC�ܺ�
[~, Output_Stages, Output_TEHP] = SFMD(NumStage, SInFeeds, SInPerms, Membranes, TECs);
WP = MembrProps(1).Area*Output_Stages.JM;
E = Output_TEHP.EC;

% ������Һ���ȵ��趨�¶ȵ��ܺ�
QHS = SInFeeds(1).MassFlow*SInFeeds(1).SpecHeat*(SInFeeds(1).Temp-298.15);
EH = QHS*1; % �����100%Ч�ʵ��ܺ�
% ������͸Һ�ﵽ�趨�¶ȵ��ܺ�
QCS = SInPerms(1).MassFlow*SInPerms(1).SpecHeat*(SInPerms(1).Temp-298.15);
if QCS >= 0
    EC = QCS; % �����100%Ч�ʵ��ܺ�
else
    eta_c = 1-SInPerms(1).Temp/298.15;
    EC = -QCS*eta_c/(1-eta_c); % ����ŵѭ��������ܺ�
end
SEC = (EH+EC+sum(abs(E)))/abs(WP);