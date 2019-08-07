function [I_calc] = TE_Current(Th, Tc, N, r, opt)
%% Calculate current according hot- and cold-side temperatures
%  notes of I/O arguments
%  Th - (i double scalar) hot-side temperature [K]
%  Tc - (i double scalar) cold-side temperature [K]
%  N  - (i double scalar) number of thermocouples in the first stage
%  r  - (i double scalar) number ratio of thermocouples in two stages
%  opt- (i optional integer scalar) running mode
%       0: (default) currents of one-stage TEC to make Qc = 0
%       1:           currents of two-stage TEC to make Qc = 0
%  I_calc - (o double array(2) for opt = 0) currents in one-stage TEC [A]
%           (o double array(2) for opt = 1) currents in two-stage TEC [A]
%  Q  - (o double array(2)) heats flowing out/in the hot/cold side of TEC
%
%  by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
%% function body
% default argument of input opt
if nargin < 4
    opt = 0;
    r = 0;
end
syms I;
switch opt
    case 0
        [a, R, K ] = TE_MaterialProp((Th+Tc)/2, 0.0015);
        eq12 = 0 == I*a*Tc-I^2*R/2-K*(Th-Tc);
        current = eval(solve(eq12, I));
        I_calc = current;
    case 1
        I_calc = [0 0];
        for j = 1:2
            % intialize
            Tm = (Th+Tc)/2.d0;
            dTm = 1e5;
            while dTm > 1e-5
                % 计算热电偶参数
                % calculate a2 R2 K2
                [a2, R2, K2] = TE_MaterialProp((Tc+Tm)/2, 0.0015);
                % define heat balance equations
                eq4 = 0 == (I*a2*Tc-I^2*R2/2-K2*(Tm-Tc))*N;
                % 计算电流
                current = eval(solve(eq4, I));
                dTm_vec = zeros(size(current));
                % 计算Tm
                for i = 1:length(current)
                    dTm_vec(i) = TE_Tm(Th, Tc, current(i), r)-Tm;       
                end
                dTm = dTm_vec(j);
                Tm = Tm+dTm;
            end
            % 输出电流向量
            I_calc(j) = current(j);
        end
    otherwise
        fprintf('[ERROR] Invalid input argument!\n');
        return
end
%
end
