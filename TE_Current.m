function [Ival, Tval] = TE_Current(Th, Tc, TEC, opt)
%% Calculate current according hot- and cold-side temperatures
%  notes of I/O arguments
%  Th  - (i double scalar) hot-side temperature [K]
%  Tc  - (i double scalar) cold-side temperature [K]
%  TEC - (i struc) struc variable
%        NumTC     : Number of thermocouples in TEC
%        NumRatio  : ratio of thermocouples in the 1-stage TEC to those in
%                     the 2-stage TEC
%        GeomFactor: geometry factor of thermcouples in TEC [m]
%                 ** the unit of GeomFactor shall be [m-1] due to eq.(6)
%                    and (7) in ref
%        SeebeckCoefficient: Seebeck coefficient of 1 and 2 stage of TEC
%        ElecConductance   : electrical conductance of 1 and 2 stage of TEC
%        ThermConductance  : thermal conductance of 1 and 2 stage of TEC
%  opt - (i integer scalar) optional running mode
%                           = 0 (default) 采用参考文献[1]的方法计算TEC性能参数
%                           = 1 采用参考文献[2]的方法计算TEC性能参数
%  Ival- (o double array(2) for opt = 0) currents in one-stage TEC [A]
%        (o double array for opt = 1) currents in two-stage TEC [A]
%  Tval- (o double scalar for opt = 0) junction temperature [K]
%        (o double scalar for opt = 1) junction temperature [K]
%
%  References
%  [1] Xuan X C, et al. Cryogenics, 2002, 42: 273-278.
%  [2] Huang B J, et al. International Journal of Refrigeration 2000, 23(3): 208-218.
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%  2019-08-10: update according to the new TE_Tm()
%
%% function body
% 初始化
switch nargin
    case(3)
        opt = 0;
    case(4)
        if opt ~= 0 && opt ~= 1
            fprintf('[ERROR] Unknown specified running mode of %d for TE_Current()\n', opt)
            return
        end
end
%  TECStage = 1: currents of one-stage TEC to make Qc = 0
%             2: currents of two-stage TEC to make Qc = 0
if TEC.NumRatio == 0
    TECStage = 1;
else
    TECStage = 2;
end
% use temperature-independant properties at T = (Th+Tc)/2
TEC = TE_MaterialProp((Th+Tc)/2, TEC, opt);
a = TEC.SeebeckCoefficient;
R = TEC.ElecConductance;
K = TEC.ThermConductance;
% calculate the thermocouple number in the first stage of 2-stage TEC
N0 = TEC.NumTC/(TEC.NumRatio+1);
%
switch TECStage
    case 1
        syms I;
        eq12 = 0 == I*a*Tc-I^2*R/2-K*(Th-Tc);
        Ival = eval(solve(eq12, I));
        Tval = 0;
    case 2
        syms I Tm;
        a1 = a; a2 = a; R1 = R; R2 = R; K1 = K; K2 = K;
        r = TEC.NumRatio;
        Qc = 0;
        eq3 = (I*a1*Tm-I^2*R1/2-K1*(Th-Tm))*r == ...
              I*a2*Tm+I^2*R2/2-K2*(Tm-Tc);
        eq4 = Qc == (I*a2*Tc-I^2*R2/2-K2*(Tm-Tc))*N0;
        sol = solve([eq3,eq4], [I,Tm]);
        Ival = TE_Complex2Real(vpa(sol.I), 1e-6);
        Tval = TE_Complex2Real(vpa(sol.Tm), 1e-6);
        Tval = Tval(Tval>Tc & Tval<Th);
    otherwise
        fprintf('[ERROR] Invalid input argument!\n');
end
%
end
