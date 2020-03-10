function [SOut, flag] = DCMD_PackStream(SIn, DuctGeom)
%% Complete the stream's velocity, Reynolds number and enthalpy
%  notes of arguments
%  input arguments
%  SIn  - (struct) stream properties
%           .Temp: temperature [K]
%           .MassFlow: mass flowrate [kg/s]
%           .Velocity: velocity [m/s]
%           .MassFraction: mass fraction of NaCl
%           .Density: density [kg/m3]
%           .Viscosity: dynamic viscosity [Pa-s]
%           .SpecHeat: specific heat [J/kg-K]
%           .ThermCond: thermal conductivity [W/m-K]
%           .Enthalpy: enthalpy [W]
%  DuctGeom - geometric parameters of duct in flatsheet MD module
%           .Length (real) length along the flowing direction [m]
%           .Height (real) height of rectanglarly wetted perimeter [m]
%           .Width  (real)  width of rectanglarly wetted perimeter [m]
%
%  output arguments
%  SOut - (struct) stream properties
%           .Temp: temperature [K]
%           .MassFlow: mass flowrate [kg/s]
%           .Velocity: velocity [m/s]
%           .MassFraction: mass fraction of NaCl
%           .Density: density [kg/m3]
%           .Viscosity: dynamic viscosity [Pa-s]
%           .SpecHeat: specific heat [J/kg-K]
%           .ThermCond: thermal conductivity [W/m-K]
%           .Enthalpy: enthalpy [W]
%  flag - (integer scalar) = 0 normal exit
%                             -1 abnormal
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-16
%
% 2019-09-02 (GGQ) add 2nd input argument to calculate avg. velocity
%
% default argument of optional input
if nargin < 2
    DuctGeom = struct('Length', 0.04,  'Height', 0.006, 'Width',  0.04 );
end
% initialize
flag = 0;
SOut = SIn;
T = SIn.Temp;
W = SIn.MassFlow;
xw = SIn.MassFraction;
% correlate the stream's physical properties
rho = SIn.Density  ;
mu  = SIn.Viscosity;
cp  = SIn.SpecHeat ;
k   = SIn.ThermCond;
% get geometry of duct
Length = DuctGeom.Length;
Width  = DuctGeom.Width;
Height = DuctGeom.Height;
% calculate 
v  = W/rho/(Width*Height); % average velocity [m/s]
Dh = (Width*Height)/(Width+Height); % hydraulic diameter [m]
De = sqrt((Width*Height))/pi*4; % equivalent diameter [m]
Re = De*v*rho/mu; % Reynolds number according to equivalent diameter
H  = W*cp*T; % fluid enthalpy [W]
% output
SOut.Temp         = T;
SOut.Velocity     = v;
SOut.MassFraction = xw;
SOut.Density      = rho;
SOut.Viscosity    = mu;
SOut.SpecHeat     = cp;
SOut.ThermCond    = k;
SOut.Enthalpy     = H;
SOut.Re           = Re;
%
end