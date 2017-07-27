function [SpikeTimes,SpikeShapes,Mat,Grades,FileTimes]=MC_spSort(file,ch,fileView,Time)
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
    Time=nan;
end
if nargin<3
    fileView=file;
end

[Mat,Grades]=deal(nan);
ccc={'b','r','g','c','m','k','y','b','r','g'};
scrsz = get(0,'ScreenSize');

try,

% get parameters

FROM=10e3; TO=90e3;
rawData=[]; xtime=[]; digData=[];
for j=1:length(fileView),
    [rd,txt,Hz,totalMS]=MC_loadElectrodes(fileView(j),ch,[FROM TO]);
    [dd,tt]=MC_loadDigital(fileView(j),[FROM TO]);
    %fileName=MC_fileName(fileView(j));
    %[dataStream,totalMS,Hz]=MC_openFile(fileName);
    if totalMS>TO,
        %[rd,xt]=MC_getElectrodes(fileName,ch,[FROM TO]);
        rawData=[rawData rd];
        digData=[digData dd];
        xt=txt;
    end
    clear rd dd tt txt;
end   
%[rawData,xt,Hz]=MC_spGetElectrodesDisplay(file,ch);
%rawData=rawData(:,fileView);

if isempty(rawData)
    error('no data');
end
[rms,threshold,outliers,whatevents,evinterv]=MC_spCutSpikesGetParam(rawData,xt,Hz,digData);
clear rawData txt xt;

Total_ms=0; waveForms=[]; spikeTimes=[]; %spikeSamples=[];
% unite spikes over files
for i=1:length(file)
    % get raw data
    [rawData,xtime,Hz,totalMs]=MC_loadElectrodes(file(i),ch,Time);
    [digData]=MC_loadDigital(file(i),Time);
    % cut spikes by threshold
    [wf,st]=MC_spCutSpikes(rawData,xtime,digData,Hz,rms,threshold,outliers,whatevents,evinterv);        
    waveForms=[waveForms wf];
    spikeTimes=[spikeTimes; Total_ms+st];
    clear rawData xtime st wf;
    %spikeSamples=[spikeSamples smp];
    Total_ms=Total_ms+totalMs;
    FileTimes(i)=Total_ms;
end

if isempty(spikeTimes)
    error('no spikes');
end

sbps={[2 2],[2 2],[2 3],[3 2],[3 3],[3 3],[4 3],[4 3],[4 3]};

not_done=1; 
Inxs{1}=1:length(spikeTimes);
draw_ind=1;
while not_done,
    if draw_ind,
        ngroups=length(Inxs);
        n=ceil(sqrt(ngroups));
        mnn=min(waveForms(:));
        mxx=max(waveForms(:));
        f1=figure('Position',[1*round(scrsz(3)/10) 3*round(scrsz(4)/10) 8*round(scrsz(3)/10) 7*round(scrsz(4)/10)]); % wave forms & amplitudes
        f2=figure('Position',[1*round(scrsz(3)/10) 3*round(scrsz(4)/10) 8*round(scrsz(3)/10) 7*round(scrsz(4)/10)]); % PCs
        mins=min(waveForms(:,:));
        curr_group=0;
        for i=1:ngroups,
            wv=waveForms(:,Inxs{i});
            figure(f2); 
            MC_spPlotPCs(wv,ccc{i},f2);
            figure(f1);
            cf=subplot(sbps{ngroups}(1),sbps{ngroups}(2),i);
            MC_spPlotShape(wv,cf,mnn,mxx,ccc{i});
            title(['Group #' num2str(i)]);
            subplot(sbps{ngroups}(1),sbps{ngroups}(2),[((sbps{ngroups}(1)-1)*sbps{ngroups}(2))+1 : sbps{ngroups}(1)*sbps{ngroups}(2)]); hold on;
            plot(spikeTimes(Inxs{i}),mins(Inxs{i}),'.','color',ccc{i},'markersize',4);   
        end
        set(gca,'xlim',[min(spikeTimes) max(spikeTimes)]);
        set(gca,'xtick',FileTimes); 
        for i=1:length(FileTimes), ss{i}=sprintf('%3.2f',FileTimes(i)/1e3); end; 
        set(gca,'xticklabel',ss);
    end
    
    [action,curr_group]=MC_spListActionGroups(ngroups);
    draw_ind=1;
    switch action,
        case 'Cluster'
            new_inxs=MC_spCluster(waveForms(:,Inxs{curr_group}));
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
Mat=sparse(Total_ms,ngroups);
for i=1:ngroups
    Mat(round(SpikeTimes{i}),i)=1;
end


catch,
    s=questdlg('now what?','','next channel','quit','next channel');
    if isempty(s) | strcmp(s,'quit')
        close all;
        error('stop');
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
function [action,curr_group]=MC_spListActionGroups(ngroups);

curr_group=1;
actions={'Cluster','ManualCluster','Delete','Merge','Manual','Amplitude','Nothing','Done'};
[act,ok] = listdlg('PromptString','Action?','SelectionMode','single','ListString',actions);
if ~ok 
    act='Cancel';
    return;
end
action=actions{act};
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

str=[]; for i=1:ngroups, str{i}=num2str(i); end;
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

