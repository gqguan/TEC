function [Q, a, R, K] = TE_Heat(Th, Tc, I, N, r, gf)
%% calculate the heats of TEC
%  notes of I/O arguments
%  Th - (i double scalar) hot-side temperature [K]
%  Tc - (i double scalar) cold-side temperature [K]
%  I  - (i double scalar) serial electric current flowing through the two
%                         stages [A]
%  N  - (i double scalar) number of thermocouples in the first stage
%  r  - (i double scalar) number ratio of thermocouples in two stages
%                         = 0 indicates one-stage TEC
%  gf - (i double scalar) geometry factor = A/L [m]
%  Q  - (o double array(2)) heats flowing out/in the hot/cold side of TEC
%  a  - (o double array(2)) Seebeck coefficient of 1 and 2 stage of TEC
%       (o double scalar)                          one-stage TEC
%  R  - (o double array(2)) electrical resistance of 1 and 2 stage of TEC
%       (o double scalar)                          one-stage TEC
%  K  - (o double array(2)) thermal conductance of 1 and 2 stage of TEC
%       (o double scalar)                          one-stage TEC
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%  2019-08-07: add case 0 for calculating the heats in one-stage TEC
%
%% function body
switch r
    case 0
        [a, R, K] = TE_MaterialProp((Th+Tc)/2, gf);
        Q(1) = (I*a*Th+I^2*R/2-K*(Th-Tc))*N;
        Q(2) = (I*a*Tc-I^2*R/2-K*(Th-Tc))*N;
    otherwise
        % Get Tm
        [Tm, a, R, K] = TE_Tm(Th, Tc, I, r, gf);
        % Œ¸°¢∑≈»»¡ø
        Q(1) = (I*a(1)*Th+I^2*R(1)/2-K(1)*(Th-Tm))*N*r;
        Q(2) = (I*a(2)*Tc-I^2*R(2)/2-K(2)*(Tm-Tc))*N;
end
%
end
