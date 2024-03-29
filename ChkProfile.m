%% 检查集成半导体热泵的DCMD膜组件模拟结果
% 总物质平衡
function chkList = ChkProfile(profile,flowPattern)
[DuctGeom,Stream,MembrProps] = InitStruct();
SIN(1) = DCMD_PackStream(profile.S1(1),DuctGeom);
SOUT(1) = DCMD_PackStream(profile.S1(end),DuctGeom);
switch flowPattern
    case 'countercurrent'
        SIN(2) = DCMD_PackStream(profile.S2(end),DuctGeom);
        SOUT(2) = DCMD_PackStream(profile.S2(1),DuctGeom);
    case 'cocurrent'
end
Q1TEC1 = sum(cellfun(@(x)x(1,1),profile.QTEC)); % 膜组件热侧TEC1放热量
Q2TEC1 = sum(cellfun(@(x)x(1,2),profile.QTEC)); % 膜组件热侧TEC1吸热量
Q1TEC2 = sum(cellfun(@(x)x(2,1),profile.QTEC)); % 膜组件冷侧TEC2放热量
Q2TEC2 = sum(cellfun(@(x)x(2,2),profile.QTEC)); % 膜组件冷侧TEC2吸热量
QM = sum(profile.QM); % 跨膜传热量[W]
WP = sum(arrayfun(@(x)x.MassFlow,profile.SM)); % 跨膜渗透量[kg/s]
%
i = 1;
Items{i} = '【1】DCMD膜组件料液侧流率减少量与渗透侧流率增加量偏差（kg/s）';
Results{i} = (SIN(1).MassFlow-SOUT(1).MassFlow)-(SOUT(2).MassFlow-SIN(2).MassFlow);
[Criteria{i},Remarks{i}] = ChkIt(Results{i},1e-8);
i = 2;
Items{i} = 'DCMD膜组件料液焓变（W）';
Results{i} = SIN(1).Enthalpy-SOUT(1).Enthalpy;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 3;
Items{i} = 'DCMD膜组件渗透液焓变（W）';
Results{i} = SIN(2).Enthalpy-SOUT(2).Enthalpy;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 4;
Items{i} = 'DCMD膜组件中热侧TEC1放热量（W）';
Results{i} = Q1TEC1;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 5;
Items{i} = 'DCMD膜组件中冷侧TEC2吸热量（W）';
Results{i} = Q2TEC2;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 6;
Items{i} = '【2】DCMD膜组件输入输出能量偏差（W）';
Results{i} = sum([Results{2:5}]);
[Criteria{i},Remarks{i}] = ChkIt(Results{i},1e-4);
i = 7;
Items{i} = 'DCMD膜组件料液侧物料焓变（W）';
Results{i} = profile.S1(1).Enthalpy-profile.S1(end).Enthalpy;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 8;
Items{i} = 'DCMD膜组件料液侧物料吸热量（Q1TEC1，W）';
Results{i} = Q1TEC1;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 9;
Items{i} = 'DCMD膜组件料液侧物料放热量（QM，W）';
Results{i} = -QM;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 10;
Items{i} = '【3】DCMD膜组件料液侧能量平衡偏差（W）';
Results{i} = sum([Results{7:9}]);
[Criteria{i},Remarks{i}] = ChkIt(Results{i},1e-4);
i = 11;
Items{i} = 'DCMD膜组件渗透侧物料焓变（W）';
Results{i} = profile.S2(end).Enthalpy-profile.S2(1).Enthalpy;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 12;
Items{i} = 'DCMD膜组件渗透侧物料吸热量（QM，W）';
Results{i} = QM;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 13;
Items{i} = 'DCMD膜组件渗透侧物料放热量（Q2TEC2，W）';
Results{i} = -Q2TEC2;
[Criteria{i},Remarks{i}] = ChkIt(Results{i});
i = 14;
Items{i} = '【4】DCMD膜组件渗透侧能量平衡偏差（W）';
Results{i} = sum([Results{11:13}]);
[Criteria{i},Remarks{i}] = ChkIt(Results{i},1e-4);
% 输出
Items = reshape(Items,[],1);
Criteria = reshape(Criteria,[],1);
Results = reshape(Results,[],1);
Remarks = reshape(Remarks,[],1);
chkList = table(Items,Results,Criteria,Remarks);

    function [critValue,remarks] = ChkIt(value,critValue)
        remarks = 'undo';
        if ~exist('critValue','var')
            critValue = missing;
            return
        end
        argType = class(value);
        switch argType
            case 'double'
                if critValue > abs(value)
                    remarks = 'okey';
                else
                    remarks = 'failure';
                end
        end
    end

end