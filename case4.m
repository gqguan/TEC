%% Calculate temperatures of 0-dimension DCMD module with TEHPs
%
% by Dr. GUAN, Guoqiang @ SCUT on 2019-08-17
%
%% Initialize
clear;
global profile;
profile = struct('WP', [], 'Q1H',  [], 'Q1C', [], 'Q2H', [], 'Q2C', [], ...
                 'QM', [], 'TS1H', [], 'TH',  [], 'TMH', [], 'TMC', [], ...
                 'TC', [], 'TS2C', []);
%  SInFeed - properties of feed-side influent
%           .Temp: temperature [K]
%           .MassFlow: mass flowrate [kg/s]
%           .Velocity: velocity [m/s]
%           .MassFraction: mass fraction of NaCl
%           .Density: density [kg/m3]
%           .Viscosity: dynamic viscosity [Pa-s]
%           .SpecHeat: specific heat [J/kg-K]
%           .ThermCond: thermal conductivity [W/m-K]
%           .Enthalpy: enthalpy [W]
SInFeed = struct('Temp', 323.15,     ...
                 'MassFlow', 0.015,  ...
                 'Velocity', 0.015*1e3/(0.04*0.006),  ...
                 'MassFraction', 0,  ...
                 'Density', 1e3,     ...
                 'Viscosity', 1e-3,  ...
                 'SpecHeat', 4.18e3, ...
                 'ThermCond', 0.6,   ...
                 'Enthalpy', 0.015*4180*323.15);
SInPerm = struct('Temp', 303.15,     ...
                 'MassFlow', 0.015,  ...
                 'Velocity', 0.015*1e3/(0.04*0.006),  ...
                 'MassFraction', 0,  ...
                 'Density', 1e3,     ...
                 'Viscosity', 1e-3,  ...
                 'SpecHeat', 4.18e3, ...
                 'ThermCond', 0.6,   ...
                 'Enthalpy', 0.015*4180*303.15);
%  MembrProps.TMH: hot-side temperature of membrane [K]
%            .TMC: cold-side temperature of membrane [K]
%            .Area: effective area of membrane [m2]
%            .Thickness: thickness of membrane [m]
%            .MDCoefficient: MD coefficient [kg/s-m2-Pa]
%            .ThermConductivity: thermal conductivity of membrane
MembrProps = struct('TMH', [], 'TMC', [], 'Area', 0.0016, ...
                    'Thickness', 1.5e-4, 'MDCoefficient', 3.2e-7, ...
                    'ThermConductivity', (0.18*0.3+0.025*0.7));
%  TEC.NumTC     : Number of thermocouples in TEC
%     .NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                  the 2-stage TEC
%     .GeomFactor: geometry factor of thermcouples in TEC [m]
%     .HTCoefficient     : Overall heat transfer coefficient [W/m2-K]
%     .HTArea            : heat transfer area [m2]
%     .SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%     .ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%     .ThermConductance  : thermal conductance of 1 and 2 stage of TEC
TEC = struct('NumTC', 190, 'NumRatio', 0, 'GeomFactor', 3.8e-4, ...
             'HTCoefficient', 270, 'HTArea', 0.0016, ...
             'SeebeckCoefficient', [], 'ElecConductance', [], ...
             'ThermConductance', []);
%% Solve TS1H and TS2C
% Temperatures at the hot surface of TEC1 and the cold surface of TEC2
I1  = 3.0; I2 = 3.0; TEC1 = TEC; TEC2 = TEC;
y   = [293.15 346]; % [TS1C,TS2H]
fun_TS = @(x)DCMD_Diff_TS(x, y, I1, TEC1, I2, TEC2, SInFeed, SInPerm, ...
                       MembrProps);
TS0  = [SInFeed.Temp-1, SInPerm.Temp+1]; % initial values of [TS1H,TS2C]
options = optimset('Display', 'off');
TS = fminsearch(fun_TS, TS0, options); % x=[TS1H,TS2C]
profile.TS1H = TS(1); profile.TS2C = TS(2);
%% Calculate heats of TEC1 and TEC2
QS1 = TE_Heat(TS(1), y(1), I1, TEC1);
profile.Q1H = QS1(1); profile.Q1C = QS1(2);
QS2 = TE_Heat(y(2), TS(2), I2, TEC2);
profile.Q2H = QS2(1); profile.Q2C = QS2(2);
%% Calculate temperatures of feed- and permeate-side bulk flows
fun_T = @(T)DCMD_Diff_THTC(T, QS1(1), QS2(2), SInFeed, SInPerm, MembrProps);
T0 = [SInFeed.Temp SInPerm.Temp];
T = fminsearch(fun_T, T0);
profile.TH = T(1); profile.TC = T(2);
%% Calculate permeation flux and membrane temperature
SFeedSide = SInFeed; % UNDER CONSTRUCTION
SPermSide = SInPerm; % UNDER CONSTRUCTION
SFeedSide.Temp = T(1);
SPermSide.Temp = T(2);
SFeedSide = DCMD_PackStream(SFeedSide);
SPermSide = DCMD_PackStream(SPermSide);
% feed-side heat and mass balance
DirectOpt = -1;
[QTransMembr,STransMembr] = DCMD_SPerm(DirectOpt, MembrProps, ...
                                       SFeedSide, SPermSide);
SOutFeed = DCMD_SOut(DirectOpt, SInFeed, STransMembr, QS1(1), QTransMembr);
profile.WP = STransMembr.MassFlow;
profile.QM = QTransMembr;
profile.TH = SFeedSide.Temp;
profile.TMH = STransMembr.Temp;
% permeate-side heat and mass balance
DirectOpt = 1;
[QTransMembr,STransMembr] = DCMD_SPerm(DirectOpt, MembrProps, ...
                                       SFeedSide, SPermSide);
SOutPerm = DCMD_SOut(DirectOpt, SInPerm, STransMembr, QTransMembr, QS2(2));
profile.TC = SPermSide.Temp;
profile.TMC = STransMembr.Temp;
%% Output results
fprintf('For TS2H = %6.2f K ...\n', y(2));
fprintf('Q1C = %5.2f W, Q1H = %5.2f W, Q2C = %5.2f W, Q2H = %5.2f W\n', ...
        profile.Q1C, profile.Q1H, profile.Q2C, profile.Q2H);
fprintf('TS1H = %6.2f K, TH = %6.2f K, TMH = %6.2f K \n ', ...
        profile.TS1H, profile.TH, profile.TMH);
fprintf('TMC = %6.2f K, TC = %6.2f K, TS2C = %6.2f K\n', ...      
        profile.TMC, profile.TC, profile.TS2C);