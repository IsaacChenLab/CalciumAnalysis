function [Inxs,f1]=MC_spManualDeletion(waveForms)
    
f1=figure; hold on;

    
Inxs=1:size(waveForms,2);  
wf=waveForms;
% manual deletion of spikes
while(1)
    cla;
    MC_spPlotShape(wf,gca);
    [x,y]=ginput(2);
    if isempty(x)
        break;
    end
    h=line(x,y); set(h,'color','r');
    pause(0.5);
    xxx=MC_spFindCross(repmat([1:size(wf,1)]',1,size(wf,2)),wf,x,y);
    wf(:,xxx)=[];
    Inxs(xxx)=[];
end

%close(f1);

return;



