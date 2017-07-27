function [Inxs]=MC_spManualSeperation(waveForms);
    
    f1=figure;
    % choose according to chi-square distribution of squared distances from average form
    wf=waveForms(:,:);

    hold on;
    MC_spPlotShape(wf,gca);
    
    % manual deletion of spikes
    sepSpikeInx=[];
    while(1)
        [x,y]=ginput(2);
        if isempty(x)
            break;
        end
        h=line(x,y); set(h,'color','g');
        xxx=MC_spFindCross(repmat([1:size(wf,1)]',1,size(wf,2)),wf,x,y);
        sepSpikeInx=[sepSpikeInx; xxx];
        plot(wf(:,xxx),'g');
    end
    Inxs{1}=sepSpikeInx;
    Inxs{2}=setdiff(1:size(wf,2),sepSpikeInx);
    close(f1);

return;



