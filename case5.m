%% Calculate temperatures of stacked DCMD module with TEHPs
%  For the scenario of N-stage DCMD integrated with N+1 TEHPs
%
%  by Dr. GUAN, Guoqiang @ SCUT on 2019-08-30
%
%% Initialize
clear;
NumStage = 2;
W1 = 0.015; W2 = 0.015;
T1 = 323.15; T2 = 303.15;
%  DuctGeom - geometric parameters of flowing channel in MD module
%           .Length (real array(2)) length along the flowing direction in 
%                                   both sides of MD module [m]
%           .Height (real array(2)) height of rectanglarly wetted perimeter
%                                   in both sides of MD module [m]
%           .Width (real array(2))  width of rectanglarly wetted perimeter
%                                   in both sides of MD module [m]
DuctGeom = struct('Length', 0.04,  ...
                 'Height', 0.006, ...
                 'Width',  0.04 );
% initial properties of feed-side influent
SInFeed = struct('Temp', T1,     ...
                 'MassFlow', W1,  ...
                 'MassFraction', 0,  ...
                 'Density', 1e3,     ...
                 'Viscosity', 1e-3,  ...
                 'SpecHeat', 4.18e3, ...
                 'ThermCond', 0.6);
% calculate the rest properties of feed-side influent
SInFeed = DCMD_PackStream(SInFeed);
% get all feed-side influents for each stages
SInFeeds(1:NumStage) = SInFeed;
% initial properties of permeate-side influent
SInPerm = struct('Temp', T2,         ...
                 'MassFlow', W2,     ...
                 'MassFraction', 0,  ...
                 'Density', 1e3,     ...
                 'Viscosity', 1e-3,  ...
                 'SpecHeat', 4.18e3, ...
                 'ThermCond', 0.6);
% calculate the rest properties of permeate-side influent
SInPerm = DCMD_PackStream(SInPerm);
% get all permeate-side influents for each stage
SInPerms(1:NumStage) = SInPerm;
%  MembrProps.TMH: hot-side temperature of membrane [K]
%            .TMC: cold-side temperature of membrane [K]
%            .Area: effective area of membrane [m2]
%            .Thickness: thickness of membrane [m]
%            .MDCoefficient: MD coefficient [kg/s-m2-Pa]
%            .ThermConductivity: thermal conductivity of membrane
MembrProps = struct('TMH', [], 'TMC', [], 'Area', 0.0016, ...
                    'Thickness', 1.5e-4, 'MDCoefficient', 3.2e-7, ...
                    'ThermConductivity', (0.18*0.3+0.025*0.7));
Membranes(1:NumStage) = MembrProps;
%  TEC.NumTC     : Number of thermocouples in TEC
%     .NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                  the 2-stage TEC
%     .GeomFactor: geometry factor of thermcouples in TEC [m]
%     .HTCoefficient     : Overall heat transfer coefficient [W/m2-K]
%     .HTArea            : heat transfer area [m2]
%     .SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%     .ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%     .ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%     .Voltage           : electrical voltage [V]
%     .Current           : electrical current [A]
TEC = struct('NumTC', 190, 'NumRatio', 7/12, 'GeomFactor', 2.6e-3, ...
             'HTCoefficient', 270, 'HTArea', 40*40e-6, ...
             'SeebeckCoefficient', [], 'ElecConductance', [], ...
             'ThermConductance', [], 'Voltage', 12, 'Current', 0.8, ...
             'Parameters', []);
% set properties for all TECs
TECs(1:(NumStage+1)) = TEC;
% set the initial temperatures for all stages
for i=1:NumStage
    T0((1+(i-1)*6):6*i) = [T1+1; T1; T1-1; T2+1; T2; T2-1];
end
% set room temperature as the environmental temperature of both heat source
% and sink
TEXs = [298.15; 298.15];
%% Solve temperatures
opts = optimoptions('fsolve', 'Display', 'Iter', 'MaxFunEvals', 15000, 'MaxIter', 1000);
fun = @(T)DCMD_EqSys(T, TEXs, TECs, SInFeeds, SInPerms, Membranes);
[T, fvals, exitflag] = fsolve(fun, T0, opts);
[~, Q, QM, SM, SOutFeeds, SOutPerms] = fun(T);
%% Calculate the energy efficiency
% Energy consumption of TECs
EC = Q(:,1)-Q(:,2);
% Specific energy consumption
WP_Sum = sum([SM.MassFlow]);
SEC = sum(EC)/WP_Sum;
% Coefficient of performance
COP_H = Q(:,1)./EC; % heating COPs of TEHP
COP_C = Q(:,2)./EC; % refrigrating COPs of TEHP
%% Output results
format short g
fprintf('Specific energy consumption of %d-stage DCMD is %.4e W/kg.\n', NumStage, SEC);
% Performances of membrane separation in each stage
StageNames = cell(NumStage,1);
JM = zeros(NumStage,1);
JH = zeros(NumStage,1);
for i = 1:NumStage
    StageNames{i} = sprintf('Stage-%d', i);
    JM(i) = SM(i).MassFlow/Membranes(i).Area;
    JH(i) = QM(i)/Membranes(i).Area;
end
TOut = reshape(T, [6,length(SM)]); % Temperature profiles of each stage
TSH = TOut(1,:)';
TH  = TOut(2,:)';
TMH = TOut(3,:)';
TMC = TOut(4,:)';
TC  = TOut(5,:)';
TSC = TOut(6,:)';
Output_Stages = table(TSH, TH, TMH, TMC, TC, TSC, JH, JM, ...
                      'RowNames', StageNames);
disp(Output_Stages)
% Performances of TEHPs
TECNames = cell(NumStage+1,1);
for i = 1:(NumStage+1)
    TECNames{i} = sprintf('TEHP-%d', i);
end
TS_Hots(1:NumStage,1)  = TOut(1,:)';
TS_Hots(NumStage+1,1) = TEXs(2);
TS_Colds(1,1) = TEXs(1);
TS_Colds(2:NumStage+1,1) = TOut(6,:)';
Q1 = Q(:,1);
Q2 = Q(:,2);
Output_TEHP = table(COP_H, COP_C, TS_Hots, TS_Colds, Q1, Q2, EC, ...
                    'RowNames', TECNames);
disp(Output_TEHP)