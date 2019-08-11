function [Q, TEC] = TE_Heat(Th, Tc, I, TEC)
%% calculate the heats of TEC
%  notes of I/O arguments
%  Th  - (i double scalar) hot-side temperature [K]
%  Tc  - (i double scalar) cold-side temperature [K]
%  I   - (i double scalar) serial electric current flowing through the two
%                         stages [A]
%  TEC - (i struc) struc variable
%        NumTC     : Number of thermocouples in TEC
%        NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                    the 2-stage TEC
%        GeomFactor: geometry factor of thermcouples in TEC [m]
%        SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%        ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%        ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%  Q   - (o double array(2)) heats flowing out/in the hot/cold side of TEC
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%  2019-08-07: add case 0 for calculating the heats in one-stage TEC
%  2019-08-10: update according to the new TE_Tm()
%
%% function body
N0 = TEC.NumTC/(TEC.NumRatio+1);
switch TEC.NumRatio
    case 0
        [a, R, K] = TE_MaterialProp((Th+Tc)/2, TEC.GeomFactor);
        Q(1) = (I*a*Th+I^2*R/2-K*(Th-Tc))*N0;
        Q(2) = (I*a*Tc-I^2*R/2-K*(Th-Tc))*N0;
    otherwise
        % Get Tm
        [Tm, TEC] = TE_Tm(Th, Tc, I, TEC);
        a = TEC.SeebeckCoefficient;
        R = TEC.ElecConductance;
        K = TEC.ThermConductance;
        % Œ¸°¢∑≈»»¡ø
        Q(1) = (I*a(1)*Th+I^2*R(1)/2-K(1)*(Th-Tm))*N0*TEC.NumRatio;
        Q(2) = (I*a(2)*Tc-I^2*R(2)/2-K(2)*(Tm-Tc))*N0;
end
%
end
