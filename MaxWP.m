%% 计算DCMD系统给定输入能量的最大产水率
% 求满足计算能耗等于给定值（E）时，产水率最大的膜组件进料温度和流率
% 其中需要调用SimDCMD()计算不同膜组件进料温度和流率的能耗和产水率
%
% by Dr. Guan Guoqiang @ SCUT on 2022/5/17
%
function maxWP = MaxWP(E,config)
% 函数调试参数
if nargin == 0
    E = 30; % [W]
    config = 'classical';
end
sn = sprintf('%s(E=%.4g[W]）Max.WP',config,E);
% 膜组件进料参数操作范围
W1bnd = [1.217e-4*5,1.217e-2];
T1bnd = [273.15+50,273.15+65];
W2bnd = [1.217e-4*5,1.217e-2];
T2bnd = [273.15+30,273.15+45];
x0 = [mean(W1bnd),mean(T1bnd),mean(W2bnd),mean(T2bnd)];
A = []; b = []; Aeq = []; beq = [];
solOpt = optimoptions(@fmincon,'Display','iter');
lb = [W1bnd(1),T1bnd(1),W2bnd(1),T2bnd(1)];
ub = [W1bnd(2),T1bnd(2),W2bnd(2),T2bnd(2)];
f = @(x)-WaterProduction(x,config); % 负号表示求最大值
% 求满足给定能耗条件下的最大产水量
[x,maxWP] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,@EnergyConsumption,solOpt);

    % 给定操作条件下的DCMD系统能耗为给定E值
    function [c,ceq] = EnergyConsumption(xs)
        W1 = xs(1);
        T1 = xs(2);
        W2 = xs(3);
        T2 = xs(4);
        switch config
            case 'classical'
                refluxRatio = inf;
                [outTab,~] = SimDCMD(sn,W1,T1,W2,T2,config,refluxRatio);
        end
        c = -1;
        ceq = outTab.E1+outTab.E2-E;
    end

    % 给定操作条件下的DCMD系统产水量
    function fval = WaterProduction(xs,cfg)
        W1 = xs(1);
        T1 = xs(2);
        W2 = xs(3);
        T2 = xs(4);
        switch cfg
            case 'classical'
                refluxRatio = inf;
                [outTab,~] = SimDCMD(sn,W1,T1,W2,T2,cfg,refluxRatio);
        end
        fval = outTab.WP;
    end

end