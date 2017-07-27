% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% July 16, 2016

% This function is used to calculate the ratio of slopes for paired pulse
% stimulation. The function graphs a smoothed version of the first
% derivative

targetChannel = 32

% Plot smoothed slope graph of CSC32
x = ((1:length(CSD_Data(targetChannel,:)))/32)+offsetms;
y = smooth(CSD_Data (targetChannel,:), 5);

figure;
plot(x,y);

disp('Mark minimum and maximum x values on the first stimulation');

% Uses user input to calculate the time window for slopes
[xTemp, yTemp] = ginput(2);
x1Values = [xTemp(1), xTemp(2)];
x1Interp = round( (xTemp-offsetms)*32);
y1Interp = CSD_Data(targetChannel,x1Interp);

% Calculates the time window for the second slope
xTime = x1Interp(2) - x1Interp(1);
x2Interp = x1Interp + 49*32 + xTime;
y2Interp = CSD_Data(32,x2Interp);

% Calculates slopes
slopeOne = (y1Interp(2) - y1Interp(1)) /xTime;
slopeTwo = (y2Interp(2) - y2Interp(1)) /xTime;

% Divides slope2/slope 1
ratioSlopes = slopeTwo/slopeOne

save x1Values.mat x1Values
save ratioSlopes.mat ratioSlopes;