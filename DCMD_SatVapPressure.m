function [ PSat ] = DCMD_SatVapPressure( T, MassFrac )
%% calculate the saturation vapor pressure
%  notes of I/O arguments
%  T        - (i real scalar) temperature [K]
%  MassFrac - (i real scalar, optional) mass fraction of NaCl
%  PSat     - (o real scalar) saturation vapor pressure [Pa]
%
%  references
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-13
%
%% function body
% default argument of input opt
if nargin < 2
    MassFrac = 0;
end
C = [73.649, -7258.2, -7.3037, 4.1653E-6, 2.0];
PSat = exp(C(1)+C(2)/T+C(3)*log(T)+C(4)*T^C(5));
%
end

