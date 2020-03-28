%% Heat balance in bulk flow
%
% input arguments
% TOut  - (real scalar) effluent temperature {K]
% QIn   - (real scalar) heat flow into the CV [W]
% QOut  - (real scalar) heat flow out from the CV [W]
% Sin   - (struct) influent properties
%         .Temp: temperature [K]
%         .MassFlow: mass flowrate [kg/s]
%         .Velocity: velocity [m/s]
%         .MassFraction: mass fraction of NaCl
%         .Density: density [kg/m3]
%         .Viscosity: dynamic viscosity [Pa-s]
%         .SpecHeat: specific heat [J/kg-K]
%         .ThermCond: thermal conductivity [W/m-K]
%         .Enthalpy: enthalpy [W]
% SubSIn- (struct, optional) substream properties
%         .Temp: temperature [K]
%         .MassFlow: mass flowrate [kg/s] negative value indicates the
%                                         substream flows out from the CV
%         .Velocity: velocity [m/s]
%         .MassFraction: mass fraction of NaCl
%         .Density: density [kg/m3]
%         .Viscosity: dynamic viscosity [Pa-s]
%         .SpecHeat: specific heat [J/kg-K]
%         .ThermCond: thermal conductivity [W/m-K]
%         .Enthalpy: enthalpy [W] negative value indicates the substream
%                                 flows out from the CV
% 
% output arguments
% dQ   - (real scalar) differences of heat balance
%
% by Dr. GUAN Guoqiang @ SCUT on 2019-08-23
%
function [dQ, SOut] = DCMD_HeatBalance_BF(TOut, QIn, QOut, SIn, SubSIn)
% default values of the 5th input argment if it is missed in the list
if nargin < 5
    SubSIn = struct('Temp',         0,  ...
                    'MassFlow',     0,  ...
                    'Velocity',     0,  ...
                    'MassFraction', 0,  ...
                    'Density',      0,  ...
                    'Viscosity',    0,  ...
                    'SpecHeat',     0,  ...
                    'ThermCond',    0,  ...
                    'Enthalpy',     0);
end
% Assume the effluent properties being same as influent
SOut = SIn;
% Set the effluent temperature and mass flow
SOut.Temp = TOut;
SOut.MassFlow = SIn.MassFlow+SubSIn.MassFlow;
SOut = DCMD_PackStream(SOut);
dQ = SIn.Enthalpy+SubSIn.Enthalpy-SOut.Enthalpy+QIn-QOut;
end