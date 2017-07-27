function MC_ploAllStim(fileName,Time,window)

if nargin<2,
    Time=[0 120];
end
if nargin<3
    window=[-750 750];
end

% rh
f1=figure;
f2=figure;
ea=[13:16 9:12 5:8 1:4];
for i=1:16,
    figure(f1);
    h1=subplot(4,4,ea(i)); hold on;
    figure(f2); 
    h2=subplot(4,4,ea(i)); hold on;
    MC_elecOnAnaTrig(i,1,fileName,window,Time);
    drawnow;
end
cdc;
eval(['print -f1 -depsc2 RHstim1_' fileName]);
eval(['print -f2 -depsc2 RHstim2_' fileName]);
close all;

ea=[5:8 1:4];
f1=figure;f2=figure;
for i=1:8, % amygdala
    figure(f1);
    h1=subplot(4,4,ea(i)); hold on;
    figure(f2); 
    h2=subplot(4,4,ea(i)); hold on;

    h1=subplot(4,4,ea(i)); hold on;
    MC_elecOnAnaTrig(16+i,1,fileName,window,Time);
    drawnow;
end
cdc;
eval(['print -f1 -depsc2 AMYstim1_' fileName]);
eval(['print -f2 -depsc2 AMYstim2_' fileName]);
close all;

% mpfc
f1=figure; f2=figure;
ea=[5 3 1 15 9 13 7 11 25 17 23 21];
for i=1:12,
    figure(f1);
    h1=subplot(5,5,ea(i)); hold on;
    figure(f2); 
    h2=subplot(5,5,ea(i)); hold on;
    MC_elecOnAnaTrig(24+i,1,fileName,window,Time);
    drawnow;
end 
cdc;
eval(['print -f1 -depsc2 MPFCstim1_' fileName]);
eval(['print -f2 -depsc2 MPFCstim2_' fileName]);
close all;


return;

    
    
    