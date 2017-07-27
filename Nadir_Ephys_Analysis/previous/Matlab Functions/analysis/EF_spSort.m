function [SpikeTimes,SpikeShapes,Grades]=EF_spSort(wav, tms, DIG, inxs)
% function [SpikeTimes,SpikeShapes,Mat,Grades]=MC_spikes(file,ch,Time)
% spike sorting routine, gets a filename, channel number, and possible
% start and end time in seconds (otherwise the whole file is taken).
% returns a cell array of event times in seconds, a cell array of all spike shapes, a sparse matrix of length
% of file in ms times number of spikes found with one where an event
% has occured; and the grades for each spike (1-4, determined by the
% user).
% procedure: 
% 1. calls MC_cutSpikes to threshold spikes and get suggested spike forms
% 2. get clusters of spikes, either by sorting (MC_spikeSort) or amplitude
% selection (MC_spikeAmplitude)
% 3. allows merging and deleting of clusters by looking at the shapes.
% 4. removes spikes from each cluster according to chi-square distribution of squared distances from average
% form
% 5. allows individual deletion of single spikes from clusters.
% 6. request for a grade (1-4, 1 is best) for each spike cluster

if nargin<4
    inxs=[];
end
if nargin<3
    DIG=[];
end

[Mat,Grades]=deal(nan);
ccc={'b','r','g','c','m','k','y','b','r','g'};
scrsz = get(0,'ScreenSize');

try,

% get parameters


Total_ms=tms(end); waveForms=wav; spikeTimes=tms; %spikeSamples=[];


if isempty(spikeTimes)
    error('no spikes');
end

sbps={[2 2],[2 2],[2 3],[3 2],[3 3],[3 3],[4 3],[4 3],[4 3]};

not_done=1; 
if isempty(inxs)
    Inxs{1}=1:length(spikeTimes);
else
    cl=unique (inxs);
    for i=1:length(cl)
        Inxs{i}=find (inxs==cl(i));
    end
end
draw_ind=1;
while not_done,
    if draw_ind,
        ngroups=length(Inxs);
        n=ceil(sqrt(ngroups));
        mnn=min(waveForms(:));
        mxx=max(waveForms(:));
        f1=figure('Position',[1*round(scrsz(3)/10) 3*round(scrsz(4)/10) 6*round(scrsz(3)/10) 5*round(scrsz(4)/10)]); % wave forms & amplitudes
        f2=figure('Position',[1*round(scrsz(3)/10) 3*round(scrsz(4)/10) 6*round(scrsz(3)/10) 5*round(scrsz(4)/10)]); % PCs
        mins=min(waveForms(:,:));
        curr_group=0;
        for i=1:ngroups,
            if length(Inxs{i})>2000
                wv=waveForms(:,Inxs{i}(1:10:end));
            else
                wv=waveForms(:,Inxs{i}(1:end));
            end
            figure(f2); 
            MC_spPlotPCs(wv(1:73,:),ccc{i},f2);
            figure(f1);
            cf=subplot(sbps{ngroups}(1),sbps{ngroups}(2),i);
            MC_spPlotShape(wv,cf,mnn,mxx,ccc{i});
            title(['Group #' num2str(i)]);
            subplot(sbps{ngroups}(1),sbps{ngroups}(2),[((sbps{ngroups}(1)-1)*sbps{ngroups}(2))+1 : sbps{ngroups}(1)*sbps{ngroups}(2)]); hold on;
            plot(spikeTimes(Inxs{i}),mins(Inxs{i}),'.','color',ccc{i},'markersize',4);   
        end
        set(gca,'xlim',[min(spikeTimes) max(spikeTimes)]);
%         set(gca,'xtick',FileTimes); 
%         for i=1:length(FileTimes), ss{i}=sprintf('%3.2f',FileTimes(i)/1e3); end; 
%         set(gca,'xticklabel',ss);
    end
    
    [action,curr_group]=MC_spListActionGroups(ngroups,DIG);
    draw_ind=1;
    switch action,
        case 'Cluster'
            new_inxs=MC_spCluster(waveForms(:,Inxs{curr_group}));
            new_inxs=MC_spInxs(Inxs{curr_group},new_inxs);
            Inxs=MC_spSetGroups('add',Inxs,new_inxs,curr_group);
        case 'Cluster channel'
            p=round(size(waveForms,1)/73);
            if p==1
                sel=1;
            else
                s=[];
                for i=1:p
                    s{i}=num2str(i);
                end
                [sel,ok] = listdlg('PromptString','Channel:','SelectionMode','single','ListString',s);
            end
            wav=waveForms(((sel-1)*73+1):(sel*73),Inxs{curr_group});
            new_inxs=MC_spCluster(wav);
            new_inxs=MC_spInxs(Inxs{curr_group},new_inxs);
            Inxs=MC_spSetGroups('add',Inxs,new_inxs,curr_group);
        case 'Amplitude'
            new_inxs=MC_spAmplitude(waveForms(:,Inxs{curr_group}));
            new_inxs=MC_spInxs(Inxs{curr_group},new_inxs);
            Inxs=MC_spSetGroups('add',Inxs,new_inxs,curr_group);
        case 'Manual'
            new_inxs=MC_spManualSeperation(waveForms(:,Inxs{curr_group}));                
            new_inxs=MC_spInxs(Inxs{curr_group},new_inxs);
            Inxs=MC_spSetGroups('add',Inxs,new_inxs,curr_group);
        case 'ManualCluster'
            new_inxs=MC_spManualCluster(waveForms(:,Inxs{curr_group}));
            new_inxs=MC_spInxs(Inxs{curr_group},new_inxs);
            Inxs=MC_spSetGroups('add',Inxs,new_inxs,curr_group);
        case 'Eliminate Channel'
            p=round(size(waveForms,1)/73);
            if p~=1
                s=[];
                for i=1:p-1
                    s{i}=num2str(i+1);
                end
                if p==2
                    sel=1;
                else
                    [sel,ok] = listdlg('PromptString','Channel:','SelectionMode','single','ListString',s);
                end
                switch sel+1
                    case 2
                        wvs=waveForms(1:73,:);
                        if p>2
                            for i=1:p-2
                                wvs(1+73*i:73*(i+1),:)=waveForms(1+73*(i+1):73*(i+2),:);
                            end
                        end
                        waveForms=wvs;
                        clear wvs;
                    case 3
                        wvs=waveForms(1:73*2,:);
                        if p>3
                            wvs(1+73*2:73*3,:)=waveForms(1+73*3:73*4,:);
                        end
                        waveForms=wvs;
                        clear wvs;
                    case 4
                        wvs=waveForms(1:73*3,:);
                        waveForms=wvs;
                        clear wvs;
                    otherwise

                end
            end
        case 'Invert spikes'
            waveForms(:,Inxs{curr_group})=0-waveForms(:,Inxs{curr_group});
        case 'XCorr'
            if length(curr_group)<=2 && length(curr_group)>=1
                if length(curr_group)==1
                    curr_group(2)=curr_group;
                end
                d1=spikeTimes(Inxs{curr_group(1)});
                if curr_group(2)>ngroups
                    d2=DIG;
                else
                    d2=spikeTimes(Inxs{curr_group(2)});
                end
                f3=figure;
                a=inputdlg('Bin size:', 'Enter binning period (ms)',1,{'1'});
                try 
                    bn=str2num(a{1});
                catch
                    bn=1;
                end
                plot (-.1:bn/1000:(.5-bn/1000),EF_XCorr(d2,d1,[-.1 .5],bn,1000));
                [bx, ~]=ginput(2);
                if ~isempty(bx)
                    ni=ones(length(d1),1);
                    for i=1:length(d2)
                        p=find (d1>d2(i)+bx(1)*1000 & d1<d2(i)+bx(2)*1000);
                        if ~isempty(p)
                            ni(p)=2;
                        end
                    end
                    new_inxs=[];
                    new_inxs{1}=find(ni==1);
                    new_inxs{2}=find(ni==2);
                    new_inxs=MC_spInxs(Inxs{curr_group(1)},new_inxs);
                    Inxs=MC_spSetGroups('add',Inxs,new_inxs,curr_group(1));
                end
                close(f3);
            end
            
        case 'Similarity'
            p=round(size(waveForms,1)/73);
            if p==1
                sel=1;
            else
                s=[];
                for i=1:p
                    s{i}=num2str(i);
                end
                [sel,ok] = listdlg('PromptString','Channel:','SelectionMode','single','ListString',s);
            end
            wav=waveForms(((sel-1)*73+1):(sel*73),Inxs{curr_group});
            rez=EF_spSimilarity(wav);
            if ~strcmp(rez,'OK')
                return;
            end
            global INX;
            new_inxs{1}=find(INX==1);
            if ~isempty(INX==2)
                new_inxs{2}=find(INX==2);
            end
            new_inxs=MC_spInxs(Inxs{curr_group},new_inxs);
            Inxs=MC_spSetGroups('add',Inxs,new_inxs,curr_group);
        case 'Slope'
            new_inxs=EF_spSlope(waveForms(:,Inxs{curr_group}));
            new_inxs=MC_spInxs(Inxs{curr_group},new_inxs);
            Inxs=MC_spSetGroups('add',Inxs,new_inxs,curr_group);
        case 'Merge'
            if length(curr_group)>1
                Inxs=MC_spSetGroups('merge',Inxs,curr_group);
            else
                draw_ind=0;
            end
        case 'Delete'
            Inxs=MC_spSetGroups('delete',Inxs,curr_group);
        case 'Done'
            not_done=0;
        case 'Nothing'
            draw_ind=0;
        case 'Cancel'
            error('stop');
       otherwise
            error('no selection');
    end
    if draw_ind
        close(f1);
        close(f2);
    end

end

ngroups=length(Inxs);
for i=1:ngroups
    cInx=MC_spChiSelection(waveForms(:,Inxs{i}));
    Inxs{i}=MC_spInxs(Inxs{i},cInx);
    [cInx,f1]=MC_spManualDeletion(waveForms(:,Inxs{i}));
    Inxs{i}=MC_spInxs(Inxs{i},cInx);
    %MC_spPlotShape(waveForms(:,Inxs{i}),gca);
    Grades(i)=MC_getGrade;
    close(f1);
    SpikeTimes{i}=spikeTimes(Inxs{i});
    SpikeShapes{i}=waveForms(:,Inxs{i});
end

% final sparse matrix of ones where the spike has occured
% Mat=sparse(Total_ms,ngroups);
% for i=1:ngroups
%     Mat(round(SpikeTimes{i}),i)=1;
% end


catch,
    s=questdlg('now what?','','next channel','quit','next channel');
    if isempty(s) | strcmp(s,'quit')
        close all;
%         error('stop');
        return
    else
        %warning(lasterr);
        SpikeTimes=nan;
        SpikeShapes=nan;
        Mat=nan;
        Grades=nan;
        FileTimes=nan;
        return;
    end
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Grade]=MC_getGrade

gg=questdlg('Cell grade','','1','2','3','2');
if isempty(gg)
    error('no grade');
end
Grade=str2num(gg);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [action,curr_group]=MC_spListActionGroups(ngroups, DIG)

curr_group=1;
actions={'Cluster','Cluster channel','ManualCluster','Delete','Merge','Manual','Amplitude','Slope','Nothing','Eliminate Channel','Invert spikes','Similarity','XCorr','Done'};
[act,ok] = listdlg('PromptString','Action?','SelectionMode','single','ListString',actions);
if ~ok 
    act='Cancel';
    return;
end
action=actions{act};
str=[]; for i=1:ngroups, str{i}=num2str(i); end;

if strcmp(action,'XCorr')
    act='XCorr';
    if ~isempty(DIG)
        str{end+1}='Dig channel';
    end
    [curr_group,ok] = listdlg('PromptString','Group:','SelectionMode','multiple','ListString',str);
    return;
end
if strcmp(action,'Quit')
    act='Cancel';
    return;
end
if strcmp(action,'Done')
    return;
end

if ngroups==9 & ~strcmp(action,'Merge') & ~strcmp(action,'Delete')
    warndlg('Please, no more than 9 groups.');
    action='Nothing';
    return;
end
if ngroups==1
    curr_group=1;
    return;
end

if ngroups==2 & strcmp(action,'Merge')
    curr_group=[1 2];
    return;
end 


if strcmp(action,'Merge'),
    [curr_group,ok] = listdlg('PromptString','Group:','SelectionMode','multiple','ListString',str);
else
    [curr_group,ok] = listdlg('PromptString','Group:','SelectionMode','single','ListString',str);
end
if ~ok
    action='Nothing';
    return;
end

curr_group=unique(curr_group);
    
return;    

