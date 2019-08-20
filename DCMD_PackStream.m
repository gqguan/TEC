function [SOut, flag] = DCMD_PackStream(SIn)
%% Complete the stream's properties according to its T and MF
%  notes of arguments
%  SIn  - (i struct) stream properties
%           .Temp: temperature [K]
%           .MassFlow: mass flowrate [kg/s]
%           .Velocity: velocity [m/s]
%           .MassFraction: mass fraction of NaCl
%           .Density: density [kg/m3]
%           .Viscosity: dynamic viscosity [Pa-s]
%           .SpecHeat: specific heat [J/kg-K]
%           .ThermCond: thermal conductivity [W/m-K]
%           .Enthalpy: enthalpy [W]
%  SOut - (o struct) stream properties
%           .Temp: temperature [K]
%           .MassFlow: mass flowrate [kg/s]
%           .Velocity: velocity [m/s]
%           .MassFraction: mass fraction of NaCl
%           .Density: density [kg/m3]
%           .Viscosity: dynamic viscosity [Pa-s]
%           .SpecHeat: specific heat [J/kg-K]
%           .ThermCond: thermal conductivity [W/m-K]
%           .Enthalpy: enthalpy [W]
%  flag - (o integer scalar) = 0 normal exit
%                             -1 abnormal
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-16
%
% initialize
flag = 0;
SOut = SIn;
T = SIn.Temp;
W = SIn.MassFlow;
MF = SIn.MassFraction;
% correlate the stream's properties
SOut.Density    = 1e3;
SOut.Viscosity  = 1e-3;
SOut.SpecHeat   = 4.18e3;
SOut.ThermCond  = 0.6;
SOut.Enthalpy   = W*SOut.SpecHeat*T;
%
end