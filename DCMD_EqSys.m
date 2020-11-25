%% Define the equation system for coupling DCMD with TEHP
%
%  input arguments
%  TIn- (real array(NumStage*6)) ����ΪNumStage*6���¶�����������ÿ��DCMD����
%  6���¶ȣ����ηֱ�Ϊ�Ȳ���桢�Ȳ����塢�Ȳ�Ĥ�桢���Ĥ�桢���������������¶ȣ� 
%                    T(1)=TS1H, T(2)=T1H, T(3)=TM1H, T(4)=TM1C, T(5)=T1C,
%                    T(6)=TS2C, T(7)=TS2H, T(8)=T2H, T(9)=TM2H, 
%                    T(10)=TM2C, T(11)=T2C, T(12) = TS3C
%                    ...
%  TEXs - (real array) hot- and cold-side environmental temperatures [K]
%  TECs - (struct array) properties of TECs TECs(1) = TEC1
%                                           TECs(2) = TEC2
%                                           TECs(3) = TEC3
%                                           ...
%  SInFeeds - (struct array) properties of feed-side influents
%                       SInFeeds(i): feed-side influent of i-th stage
%  SInPerms - (struct array) properties of permeate-side influents
%                       SInFeeds(i): permeate-side influent of i-th stage
%  Membranes - (struct array) properties of membranes
%                             Membranes(i): membrane of i-th stage
%  output arguments
%  F  - (real array(12)) left hand side of eqs.(1)-(12) 
%  Q  - (real matrix(3,2)) absorbed and released heats of TEC1-3
%  QM - (real array(2)) transmembrane heats
%  SM - (struct array(2)) properties of permeation streams
%  SOutFeeds - (struct array(2)) properties of feed-side effluents
%  SOutPerms - (struct array(2)) properties of permeate-side effluents
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-23
%  
function [F, Q, QM, SM, SOutFeeds, SOutPerms] = ...
         DCMD_EqSys(TIn, TEXs, TECs, SInFeeds, SInPerms, Membranes)
%% Check the dimension size of input arguments
% Get the number of DCMD stack
NumStack = length(Membranes);
if (length(SInPerms) == NumStack) && (length(SInFeeds) == NumStack) && ...
   (length(TECs) == NumStack+1) && (length(TIn)/NumStack == 6)
    InputOK = 1;
else
    InputOK = 0;
    fprintf('Error of input arguments in DCMD_EqSys()! \n');
    return;
end
%% Initialize
fvals     = zeros(size(TIn));
Q         = zeros(length(TECs),2);
JHSH      = zeros(size(TECs));
JHSC      = zeros(size(TECs));
JHM       = zeros(size(Membranes));
hF        = zeros(size(Membranes));
hP        = zeros(size(Membranes));
QM        = zeros(size(Membranes));
SOutFeeds = struct(SInFeeds);
SOutPerms = struct(SInPerms);
SM        = struct(SInPerms);
T = reshape(TIn, [6,NumStack]); % ��������¶���������Ϊ�ߴ���[6,NumStack]�Ķ�ά����
%% Define the equations for energy conservation
for i = 1:NumStack
    % Boundary layer adhered the hot side of TECs(i)
    % ����ü�Ĥ������Ȳ�TEC����������
    if i == 1
        Q(i,:) = TE_Heat(T(1,i), TEXs(1), TECs(i));
    else
        Q(i,:) = TE_Heat(T(1,i), T(6,i-1), TECs(i));
    end
    JHSH(i) = Q(i,1)/TECs(i).HTArea; % �ӱ��洫���Ȳ��������ͨ��
    % �����Ȳ������������Ĵ���ϵ����W/m2-K��
    [~, hF(i)] = DCMD_TM(SInFeeds(i), -JHSH(i)); % negative JH1H indicates heat flowing into the feed
    % �����Ȳ���洫��ͨ���Ĳ�ֵ����������Ȳ�����¶�ʹ��������������ͨ������������������
    fvals(1+(i-1)*6) = DCMD_HeatBalance_BL(T(1,i), T(2,i), JHSH(i), hF(i))*Membranes(i).Area; % ע��ú���������Ϊ��ͨ����ֵ
    % Feed-side bulk flow of stage i
    % ��ʼ���ü�DCMDĤ������������Ȳ�����
    SOutFeeds(i) = SInFeeds(i); SOutPerms(i) = SInPerms(i);
    % �趨���Ȳ�����������¶�
    SOutFeeds(i).Temp = T(2,i); SOutPerms(i).Temp = T(5,i);
    % ���¼�����������
    SOutFeeds(i) = DCMD_PackStream(SOutFeeds(i));
    SOutPerms(i) = DCMD_PackStream(SOutPerms(i));
    % �趨Ĥ�����¶�
    Membranes(i).TMH = T(3,i);
    Membranes(i).TMC = T(4,i);
    % �����Ȳ�����Ĥ���������������
    [QM(i), SM(i)] = DCMD_Permeation(-1, Membranes(i), SOutFeeds(i), SOutPerms(i), 1);
    % �����Ȳ�����CV�е�������ֵ����������Ȳ������¶�ʹ�����������ʱ����CV�ľ���������
    fvals(2+(i-1)*6) = DCMD_HeatBalance_BF(T(2,i), Q(i,1), QM(i), SInFeeds(i), SM(i));
    % Boundary layer adhered the hot side of the membrane in the i-th stage
    % �����Ȳ��������Ĥ��Ĵ���ͨ��
    JHM(i) = QM(i)/Membranes(i).Area;
    % �����Ȳ�������Ĥ�洫��ͨ���Ĳ�ֵ����������Ȳ�Ĥ���¶�ʹ������Ĥ��߽����������������Ĥ��߽���ڵĴ�����
    fvals(3+(i-1)*6) = DCMD_HeatBalance_BL(T(2,i), T(3,i), JHM(i), hF(i))*Membranes(i).Area;
    % Boundary layer adhered the cold side of the membrane in the 1st stage
    % ��������Ĥ������������������������
    [QM(i), SM(i)] = DCMD_Permeation(1, Membranes(i), SOutFeeds(i), SOutPerms(i), 1);
    % ��������Ĥ���������Ĵ���ͨ��
    JHM(i) = QM(i)/Membranes(i).Area;
    % �������Ĥ��߽��Ĵ���ϵ��
    [~, hP(i)] = DCMD_TM(SInPerms(i), JHM(i));
    % �������߽���еľ����������������Ĥ���¶�ʹĤ���������Ĵ��������ڿ�Ĥ������
    fvals(4+(i-1)*6) = DCMD_HeatBalance_BL(T(4,i), T(5,i), JHM(i), hP(i))*Membranes(i).Area;
    % Permeate-side bulk flow of stage i+1
    % �����i��DCMD�����TEC������i+1��DCMD�е��Ȳ�TEC���������������
    if i == NumStack
        Q(i+1,:) = TE_Heat(TEXs(2), T(6,i), TECs(i+1));
    else
        Q(i+1,:) = TE_Heat(T(1,i+1), T(6,i), TECs(i+1));
    end
    % �����i��DCMD���������������
    JHSC(i+1) = Q(i+1,2)/TECs(i+1).HTArea;
    % �����������CV�ľ���������������������¶�ʹ�����������ʱ����CV�ľ���������
    fvals(5+(i-1)*6) = DCMD_HeatBalance_BF(T(5,i), Q(i+1,2), QM(i), SInPerms(i), SM(i));
    % Boundary layer adhered the cold side of TEC2
    % ����������߽��ľ��������������������¶�ʹ�߽���ڴ��������ڵ�i+1��DCMD�Ȳ�TEC�����������
    fvals(6+(i-1)*6) = DCMD_HeatBalance_BL(T(5,i), T(6,i), JHSC(i+1), hP(i))*Membranes(i).Area;
end

%% Output results
F = fvals;
end