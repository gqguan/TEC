%% investigate the COP vs. Qc/N/(r+1) for two-stage serial TEC
% 
% by Dr. Guan Guoqiang @ SCUT on 2019-08-06
%
clear;
% set the waitbar
ui_waitbar = waitbar(0, 'Setting ... ');
% set parameters
I_Set = 3.25;
N = 190;
Th = 50.6+273.15; 
Tc = 16.5+273.15;
r = 1;
gf = 4.15e-4;
N0 = N/(r+1);
% 计算电流上下边界
IBound = TE_Current(Th, Tc, N0, r, gf, 1);
IMax = max(IBound);
IMin = min(IBound);
dI = (IMax-IMin)/30;
%
if (I_Set>IMax || I_Set<IMin)
    fprintf('Given electrical current is out of range!\n');
    return;
end
% initialize
I = IMin:dI:IMax;
COP = zeros(size(I));
XCoordinate = zeros(size(I));
%
waitbar(0, ui_waitbar, 'Calculating ... ')
% fprintf('Calculating ... [%d]: ', length(I));
for i = 1:length(I)
    str_wb = sprintf('Calculating ... %4.1f %%', (i/length(I)*100));
    waitbar(i/length(I), ui_waitbar, str_wb)    
    [Q, ~, ~, ~] = TE_Heat(Th, Tc, I(i), N0, r, gf);
    COP(i) = Q(2)/(Q(1)-Q(2));
    XCoordinate(i) = Q(2)/N0/(r+1);
end
% 
[Q, ~, ~, ~] = TE_Heat(Th, Tc, I_Set, N0, r, gf);
COP_Set = Q(2)/(Q(1)-Q(2));
XCoordinate_Set = Q(2)/N0/(r+1);
%
delete(ui_waitbar);
fprintf(' Complete calculation! \n');
str_legend = sprintf('Th = %5.1f K, Tc = %5.1f K, r = %4.1f', Th, Tc, r);
plot(XCoordinate, COP); legend(str_legend); hold on;
plot(XCoordinate_Set, COP_Set, 'x');