%% Calculate temperatures of stacked DCMD module with TEHPs
%  For the scenario of 2-stack DCMD integrated with 3 TEHPs
%
%  by Dr. GUAN, Guoqiang @ SCUT on 2019-08-19
%
%% Initialize
clear;
NumStage = 3;
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
                 'Velocity', 0.015*1e-3/(0.04*0.006),  ...
                 'MassFraction', 0,  ...
                 'Density', 1e3,     ...
                 'Viscosity', 1e-3,  ...
                 'SpecHeat', 4.18e3, ...
                 'ThermCond', 0.6,   ...
                 'Enthalpy', 0.015*4180*323.15);
SInFeeds(1:NumStage) = SInFeed;
SInPerm = struct('Temp', 303.15,     ...
                 'MassFlow', 0.015,  ...
                 'Velocity', 0.015*1e-3/(0.04*0.006),  ...
                 'MassFraction', 0,  ...
                 'Density', 1e3,     ...
                 'Viscosity', 1e-3,  ...
                 'SpecHeat', 4.18e3, ...
                 'ThermCond', 0.6,   ...
                 'Enthalpy', 0.015*4180*303.15);
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
TEC = struct('NumTC', 190, 'NumRatio', 0, 'GeomFactor', 3.8e-4, ...
             'HTCoefficient', 270, 'HTArea', 0.0016, ...
             'SeebeckCoefficient', [], 'ElecConductance', [], ...
             'ThermConductance', [], 'Voltage', 12, 'Current', 0.8);
TECs(1:(NumStage+1)) = TEC;
for i=1:NumStage
    T0((1+(i-1)*6):6*i) = [322.1; 320.0; 319; 302; 303.6; 303.8];
end
TEXs = [298.15; 307.7307];
%% Solve temperatures
opts = optimoptions('fsolve', 'Display', 'Iter');
fun = @(T)DCMD_EqSys(T, TEXs, TECs, SInFeeds, SInPerms, Membranes);
[T, fvals, exitflag] = fsolve(fun, T0, opts);
[~, Q, QM, SM, SOutFeeds, SOutPerms] = fun(T);
%% Calculate the energy efficiency
% Energy consumption of TECs
EC = Q(:,1)-Q(:,2);
% Specific energy consumption
WP = sum([SM.MassFlow]);
SEC = sum(EC)/WP;
%% Output results
TOut = reshape(T, [6,length(SM)]);
TNames = {'TSH';'TH';'TMH';'TMC';'TC';'TSC'};
Output = table(TOut, 'RowNames', TNames);
