function [SOut] = DCMD_SOut(DirectOpt, SIn, SPerm, QIn, QOut)
%% calculate T and massflow of CV effluent
%  notes of I/O arguments
%  DirectOpt  - (i integer scalar) opt = -1: flow out from CV
%                                         1: flow into CV
%  SIn  - (i struct) variable of CV influent
%     .Temp: temperature [K]
%     .Massflow: mass flowrate [kg/s]
%     .SpecHeat: specific heat [J/kg-K]
%     .Enthalpy: enthalpy [W]
%  SPerm  - (i struct) variable of transmembrane stream
%  QIn  - (i real scalar) input heat of CV [W]
%  QOut - (i real scalar) output heat of CV [W]
%  SOut - (o struct) mass flowrate of CV effluent [kg/s]
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-13
%
%%
% intialize
Win = SIn.MassFlow;
cp  = SIn.SpecHeat;
% transmembrane mass flow, which is negative if CV locates at the feed side
Wp  = SPerm.MassFlow; 
% calculate output of CV according to the material and energy conservation
Wout = Win-Wp; 
SOut.Enthalpy = SIn.Enthalpy+QIn+DirectOpt*SPerm.Enthalpy-QOut;
% output
SOut.Temp = SOut.Enthalpy/Wout/cp;
SOut.Massflow = Wout;
SOut.Specheat = cp;
%
end
