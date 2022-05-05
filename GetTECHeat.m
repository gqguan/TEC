function [Qout,E] = GetTECHeat(var,opStr,Th,Tc,TEC,opts)
    switch opts(2)
        case(0)
            TEC.Current = var;
        case(1)
            TEC.Voltage = var;
        otherwise
            error('CalcTECHeat()输入参数opts有误！')
    end
    TECQ = TE_Heat(Th,Tc,TEC,opts(1),opts(2));
    E = TECQ(1)-TECQ(2);
    switch opStr
        case('cooling')
            Qout = TECQ(2);
        case('heating')
            Qout = TECQ(1);
        otherwise
            error('CalcTECHeat()输入参数opStr有误！')
    end
end