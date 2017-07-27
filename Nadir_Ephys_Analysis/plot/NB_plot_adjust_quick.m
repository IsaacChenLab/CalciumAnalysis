% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 16, 2016


% Quickly adjust y axes
disp ('Set axes');
[xTemp, yTemp] = ginput(2)

yMax = yTemp(1)
yMin = yTemp(2)

axis([offsetms ms+offsetms yMin yMax]);