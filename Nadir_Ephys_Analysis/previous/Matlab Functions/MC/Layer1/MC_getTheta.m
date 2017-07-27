function [data,xtime]=MC_getTheta(file,ch,Time)

if nargin<3
    Time=nan;
end

[d,t,Hz]=MC_openFile(file);

[data]=MC_getElectrodes(file,ch,Time);
NY=Hz/2;

xtime=[1:Time(end)]';

% option 1
L = 1; M = 25;                
N = 24*M;
h = fir1(N,[4 10]/NY,'bandpass');
data=resample(data,L,M,h); % %data=upfirdn(data,h,L,M);

% option 2
% L = 1; M = 1e3;                
% N = 24*M;
% h = fir1(N,[4 10]/NY,'bandpass');
% data=resample(data,L,M,h);
% data=interp1(1:(1e6/Hz):xtime(end),data,xtime,'spline');
% 
% option 3
%[B,A]=butter(2,[4 10]/(Hz/2));
%[B,A]=cheby2(2,30,[4 10]/(Hz/2));
%data=filtfilt(B,A,data);
%data=resample(data,1,25);

return;
