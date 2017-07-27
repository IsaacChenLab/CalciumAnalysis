% load TTs for Vector electrode (has only seven tetrodes, TT1-TT7)

if exist ('TT1_cells.NTT') == 2
    
    [Sc1_TT_Timestamps, Sc1_TT_CellNumbers, Sc1_TT_DataPoints] = importTTs_v3 ('TT1_cells.NTT');
        
else
end

if exist ('TT2_cells.NTT') == 2

    [Sc2_TT_Timestamps, Sc2_TT_CellNumbers, Sc2_TT_DataPoints] = importTTs_v3 ('TT2_cells.NTT');
        
else
end

if exist ('TT3_cells.NTT') == 2

    [Sc3_TT_Timestamps, Sc3_TT_CellNumbers, Sc3_TT_DataPoints] = importTTs_v3 ('TT3_cells.NTT');
        
else
end

if exist ('TT4_cells.NTT') == 2

    [Sc4_TT_Timestamps, Sc4_TT_CellNumbers, Sc4_TT_DataPoints] = importTTs_v3 ('TT4_cells.NTT');
        
else
end

if exist ('TT5_cells.NTT') == 2

    [Sc5_TT_Timestamps, Sc5_TT_CellNumbers, Sc5_TT_DataPoints] = importTTs_v3 ('TT5_cells.NTT');
        
else
end

if exist ('TT6_cells.NTT') == 2

    [Sc6_TT_Timestamps, Sc6_TT_CellNumbers, Sc6_TT_DataPoints] = importTTs_v3 ('TT6_cells.NTT');
       
else
end

if exist ('TT7_cells.NTT') == 2

    [Sc7_TT_Timestamps, Sc7_TT_CellNumbers, Sc7_TT_DataPoints] = importTTs_v3 ('TT7_cells.NTT');
        
else
end

% if exist ('TT8_cells.NTT') == 2
% 
%     [Sc8_TT_Timestamps, Sc8_TT_CellNumbers, Sc8_TT_DataPoints] = importTTs_v3 ('TT8_cells.NTT');
%         
% else
% end
% 
% if exist ('TT9_cells.NTT') == 2
% 
%     [Sc9_TT_Timestamps, Sc9_TT_CellNumbers, Sc9_TT_DataPoints] = importTTs_v3 ('TT9_cells.NTT');
%         
% else
% end
% 
% if exist ('TT10_cells.NTT') == 2
% 
%     [Sc10_TT_Timestamps, Sc10_TT_CellNumbers, Sc10_TT_DataPoints] = importTTs_v3 ('TT10_cells.NTT');
%         
% else
% end

if exist ('TT11_cells.NTT') == 2

    [Sc11_TT_Timestamps, Sc11_TT_CellNumbers, Sc11_TT_DataPoints] = importTTs_v3 ('TT11_cells.NTT');
        
else
end

if exist ('TT12_cells.NTT') == 2

    [Sc12_TT_Timestamps, Sc12_TT_CellNumbers, Sc12_TT_DataPoints] = importTTs_v3 ('TT12_cells.NTT');
        
else
end
