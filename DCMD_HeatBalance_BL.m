%% Heat balance in boundary layer
%
% input arguments
% T1    - (real scalar) temperature at side 1 {K]
% T2    - (real scalar) temperature at side 2 [K]
% JH_12 - (real scalar) heat flow through 1 to 2 [W/m2]
% h     - (real scalar) overall heat transfer coefficient [W/m2-K]
% 
% output arguments
% dJH   - (real scalar) differences of heat balance
%
% by Dr. GUAN Guoqiang @ SCUT on 2019-08-23
%
function [dJH] = DCMD_HeatBalance_BL(T1, T2, JH_12, h)
    dJH = JH_12-h*(T1-T2);
end