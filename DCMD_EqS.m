%% Define the energy balance equation for TEC
%
%  I/O arguments
%  TSH - (i real scalar) hot-side temperature of TEC [K]
%  TSC - (i real scalar) cold-side temperature of TEC [K]
%  TEC - (i struc) properties of TEC
%     .NumTC             : Number of thermocouples in TEC
%     .NumRatio          : ratio of thermocouples in the 1-stage TEC to 
%                          those in the 2-stage TEC
%     .GeomFactor        : geometry factor of thermcouples in TEC [m]
%     .SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%     .ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%     .ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%     .Voltage           : electrical voltage [V]
%     .Current           : electrical current [A]
%  opt - (i integer scalar) opt=0: default output the differences of energy
%                                  balance
%                           opt=1: output the released heat from the hot
%                                  side of TEC
%                           opt=2: output the absorbed heat into the cold
%                                  side of TEC
%  out - (o real scalar) for opt=0 differences of energy balance
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-19
%  
function [out] = DCMD_EqS(TSH, TSC, TEC, opt)
%% Get default value of optional input argument (opt)
% default argument of input opt
if nargin < 4
    opt = 0;
end
%% Initialize
I = TEC.Current;
U = TEC.Voltage;
%% Calculate the heats
Q = TE_Heat(TSH, TSC, TEC);
QH = Q(1); QC = Q(2);
E = U*I;
dE = QH-QC-E;
%% Output according to opt
switch opt
    case 0
        out = dE;
    case 1
        out = QH;
    case 2
        out = QC;
end
%     
end