%% Define the energy balance equation for TEC
%
%  I/O arguments
%  TSH - (i real scalar) hot-side temperature of TEC [K]
%  TSC - (i real scalar) cold-side temperature of TEC [K]
%  I   - (i real scalar) input electrical current of TEC [A]
%  U   - (i real scalar) input electrical voltage of TEC [V]
%  TEC - (i struc) properties of TEC
%       NumTC     : Number of thermocouples in TEC
%       NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                     the 2-stage TEC
%       GeomFactor: geometry factor of thermcouples in TEC [m]
%       SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%       ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%       ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%  dE  - (o real scalar) differences of energy balance
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-19
%  
function [dE, QH, QC] = DCMD_EqS(TSH, TSC, I, U, TEC)
%% Calculate the heats
Q = TE_Heat(TSH, TSC, I, TEC);
QH = Q(1); QC = Q(2);
E = U*I;
dE = QH-QC-E;
end