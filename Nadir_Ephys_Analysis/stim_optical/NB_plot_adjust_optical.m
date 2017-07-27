% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 18, 2016

for forever = 1:1000;
    userAdjustY = input('Adjust y-axis offset: 1/2. Save: 0. ');
    
    if userAdjustY == 1
        disp ('Increase offset')
        offsetgraph = offsetgraph + 100
        NB_plot32_regraph
    elseif userAdjustY == 2
        disp ('Decrease offset')
        offsetgraph = offsetgraph - 100
        NB_plot32_regraph
    elseif userAdjustY == 0
        NB_saveFigures
        break
    end
    
    axis tight;
    
end