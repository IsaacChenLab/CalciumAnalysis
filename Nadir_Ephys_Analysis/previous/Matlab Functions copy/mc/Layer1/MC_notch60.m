function fltData=MC_notch60(data,Hz)

if nargin<2
    Hz=25e3;
end

[B,A] = butter(2,[55 65]/(Hz/2),'stop');
fltData=filtfilt(B,A,data);

return;