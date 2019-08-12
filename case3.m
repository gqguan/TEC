%% investigate the COP vs. Qc/N/(r+1) for two-stage serial TEC
% 
% by Dr. Guan Guoqiang @ SCUT on 2019-08-10
%
clear;
% set the waitbar
ui_waitbar = waitbar(0, 'Setting ... ');
% set parameters
<<<<<<< HEAD:case2.m
I_Set = 3.12;
Th = 50.6+273.15; 
Tc = 16.5+273.15;
% initialize TEC parameters
TEC = struct('NumTC', 190, 'NumRatio', 1, 'GeomFactor', 8e-4);
=======
I_Set = 1;
Th = 50.6+273.15; 
Tc = 16.5+273.15;
% initialize TEC parameters
TEC = struct('NumTC', 190, 'NumRatio', 0, 'GeomFactor', 3.8e-4, ...
             'HTCoefficient', 270, 'HTArea', 0.0016);
% calculate the hot and cold junction temperatures
T = TE_JunctionT(Th, Tc, I_Set, TEC);
Th = T(1); Tc = T(2);
>>>>>>> master:case3.m
% 计算电流上下边界
IBound = TE_Current(Th, Tc, TEC, 1);
IMax = max(IBound);
IMin = min(IBound);
dI = (IMax-IMin)/100;
%
if (I_Set>IMax || I_Set<IMin)
    fprintf('Given electrical current %5.3f A is out of range ', I_Set)
    fprintf('(%5.3f %5.3f)!\n', IMin, IMax);
    delete(ui_waitbar);
    return;
end
% initialize
I = IMin:dI:IMax;
COP = zeros(size(I));
Tm_test = zeros(size(I));
XCoordinate = zeros(size(I));
%
waitbar(0, ui_waitbar, 'Calculating ... ')
for i = 1:length(I)
    str_wb = sprintf('Calculating ... %4.1f %%', (i/length(I)*100));
    waitbar(i/length(I), ui_waitbar, str_wb)    
    [Q, TEC] = TE_Heat(Th, Tc, I(i), TEC);
    Tm_test(i) = TE_Tm(Th, Tc, I(i), TEC);
    COP(i) = Q(2)/(Q(1)-Q(2));
    XCoordinate(i) = Q(2)/TEC.NumTC;
end
%
delete(ui_waitbar);
fprintf(' Complete calculation! \n');
str_legend = sprintf('Th = %5.1f K, Tc = %5.1f K, r = %4.1f', ...
                     Th, Tc, TEC.NumRatio);
plot(XCoordinate, COP); legend(str_legend); hold on;
% plot(XCoordinate_Set, COP_Set, 'x');