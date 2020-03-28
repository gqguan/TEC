function [ LatentHeat ] = DCMD_LatentHeat( T )
%% calculate the latent heat
%  notes of I/O arguments
%  T          - (i real scalar) temperature [K]
%  LatentHeat - (o real scalar) latent heat [J/kg]
%
%  references
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-13
%
% initialize
Mw = 18;
L = [57.075, -4.3856d-2];
LatentHeat = dot(L, [1,T])/Mw*1e6;
%
end

