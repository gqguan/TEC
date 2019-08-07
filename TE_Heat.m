function [Q, a, R, K] = TE_Heat(Th, Tc, I, N, r)
%% calculate the heats of TEC
%  notes of I/O arguments
%  Th - (i double scalar) hot-side temperature [K]
%  Tc - (i double scalar) cold-side temperature [K]
%  I  - (i double scalar) serial electric current flowing through the two
%                         stages [A]
%  N  - (i double scalar) number of thermocouples in the first stage
%  r  - (i double scalar) number ratio of thermocouples in two stages
%  Q  - (o double array(2)) heats flowing out/in the hot/cold side of TEC
%  a  - (o double array(2)) Seebeck coefficient of 1 and 2 stage of TEC
%  R  - (o double array(2)) electrical resistance of 1 and 2 stage of TEC
%  K  - (o double array(2)) thermal conductance of 1 and 2 stage of TEC
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%% function body
% Get Tm
[Tm, a, R, K] = TE_Tm(Th, Tc, I, r);
% Œ¸°¢∑≈»»¡ø
Q(1) = (I*a(1)*Th+I^2*R(1)/2-K(1)*(Th-Tm))*N*r;
Q(2) = (I*a(2)*Tc-I^2*R(2)/2-K(2)*(Tm-Tc))*N;
%
end
