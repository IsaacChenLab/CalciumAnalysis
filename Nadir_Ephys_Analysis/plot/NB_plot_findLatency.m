% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 16, 2016

% This function is used to calculate and write the latency of the
% inflection point after stimulation

[xLatency, yLatency] = ginput(1);

txt = ['\leftarrow ', num2str(round(xLatency,0)), 'ms'];

xLatencyPoint = round( (xLatency-offsetms)*32 );
yLatencyPoint = CSD_Data(32,xLatencyPoint)-32*offsetgraph;

%latencyText = text(xLatency, yLatencyPoint,txt);
latencyText = text(xLatency, yLatency,txt);
%latencyText.BackgroundColor = 'w'
latencyText.FontSize = 16