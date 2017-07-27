% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 18, 2016

disp ('Set axes');
[xTemp, yTemp] = ginput(2)

yMax = yTemp(1)
yMin = yTemp(2)
yAxisChange = 1;

axis([offsetms ms+offsetms yMin*yAxisChange yMax*yAxisChange]);

for forever = 1:1000;
    userAdjustY = input('Adjust offset: 1/2. Adjust y-axis: 3/4. Reset axes: 5. Save: 0. ');
    
    if userAdjustY == 1
        disp ('Increase offset')
        offsetgraph = offsetgraph + 100
        NB_plot32_IO_regraph
        axis tight;
    elseif userAdjustY == 2
        disp ('Decrease offset')
        offsetgraph = offsetgraph - 100
        NB_plot32_IO_regraph
        axis tight;
    elseif userAdjustY == 3
        disp ('Increase axes')
        yAxisChange = yAxisChange*1.1
    elseif userAdjustY == 4
        disp ('Decrease axes')
        yAxisChange = yAxisChange*.9
    elseif userAdjustY == 5
        disp ('Reset axes')
        [xTemp, yTemp] = ginput(2)
        yMax = yTemp(1)
        yMin = yTemp(2)
        yAxisChange = 1;
        axis([offsetms ms+offsetms yMin*yAxisChange yMax*yAxisChange]);    
    elseif userAdjustY == 0
        NB_saveFigures
        break
    end
    
    axis([offsetms ms+offsetms yMin*yAxisChange yMax*yAxisChange]);
    
end