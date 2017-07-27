function DA (operation, varargin)

global DAQ LH LL LER TRI TRI_COUNTER;
global SWEEP_TIMES SWEEP_DATA TIME0 EVENTS;

switch lower(operation)
    
    case 'initialize'
        if exist('DAQ', 'var') delete (DAQ); end
        if exist('LH', 'var') delete (LH); end
        if exist('LL', 'var') delete (LL); end
        if exist('LER', 'var') delete (LER); end
        SWEEP_TIMES=[];
        SWEEP_DATA=zeros(10000,1);
        EVENTS=[];
        
        DAQ = daq.createSession('ni');
        DAQ.addAnalogOutputChannel ('Dev1','ao0','Voltage');
        DAQ.addAnalogOutputChannel ('Dev1','ao1','Voltage');
        DAQ.addAnalogInputChannel ('Dev1','ai0','Voltage');
        DAQ.Rate=250000; % input / output sampling rate at 125 KHz
        DAQ.IsContinuous=1;
        DAQ.NotifyWhenDataAvailableExceeds=25000; %save data every .1sec
        DAQ.NotifyWhenScansQueuedBelow=250000;
        LH = DAQ.addlistener ('DataAvailable', @storeData);
        LL = DAQ.addlistener ('DataRequired', @queueData);
        LER = DAQ.addlistener('ErrorOccurred', @dispError);
        TRI = [-0.6:1/625:1.6 1.6:-1/625:-0.6]'; %triangular sweep from -0.6V to 1.6V and back at 400V/s
        TRI(end:25000)=-0.6;
        TRI_COUNTER=0;
        TRI_COUNTER(end:25000)=0;
        TRI=[TRI TRI_COUNTER'];
        TIME0=clock;
        
    case 'destroy'
        if exist('DAQ', 'var') delete (DAQ); end
        if exist('LH', 'var') delete (LH); end
        if exist('LL', 'var') delete (LL); end
        if exist('LER', 'var') delete (LER); end
        if exist('TIMER','var')
            try
                stop(TIMER);
            catch end
            delete (TIMER);
        end
        clear DAQ LH TIMER TRI SWEEP_TIMES SWEEP_DATA EVENTS TIME0;
    
    case 'reset'
        TIME0=clock;
        SWEEP_TIMES=[];
        SWEEP_DATA=zeros(10000,1);
        EVENTS=[];
        
    case 'start'
        DAQ.queueOutputData([TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI]);
        DAQ.startBackground;
        
    case 'stop'
        DAQ.stop;
        
    case 'condition'
%         a=questdlg ('This process takes ~ 45min. Are you sure you want to procede?', 'Question', 'Yes', 'No','No');
%         if a(1)=='N'
%             return
%         end
        a=questdlg('Place the electrode in PBS 7.4 and press OK.','Starting...','OK','Cancel','Cancel');
        if a(1)=='C'
            return;
        end
        TRI = [0:0.0017:3 3:-0.0017:0]'; %triangular sweep from 0V to +3V and back at 70Hz
        TRI_COUNTER=zeros(length(TRI),1);
        TRI=[TRI TRI_COUNTER];
        TRI=[TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI];
        f=figure;
        for i=1:20
            DA('reset');
            DA('start');
            pause(1);
            DA('stop');
            figure(f);
            plot (SWEEP_DATA(:,2:end),'k');
            axis ([0 10000 -10 10]);
            title([num2str(20-i) ' seconds remaining.']);
        end
        TRI = ones (50000,1)*1.5; % hold at +1.5V
        TRI_COUNTER=zeros(length(TRI),1);
        TRI=[TRI TRI_COUNTER];
        TRI=[TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI];
        for i=1:20
            DA('reset');
            DA('start');
            pause(1);
            DA('stop');
            figure(f);
            plot (SWEEP_DATA(:,2:end),'k');
            axis ([0 10000 -10 10]);
            title([num2str(20-i) ' seconds remaining.']);
        end
        close (f)
        disp('done');
        DA('initialize');
        
    case 'background'
        DA('reset');
        TRI = [-0.6:1/625:1.6 1.6:-1/625:-0.6]';
        TRI(end:25000)=-0.6;
        TRI_COUNTER=0;
        TRI_COUNTER(end:25000)=0;
        TRI=[TRI TRI_COUNTER'];
        if isempty(varargin)
            del=2;
        else
            del=varargin{1};
        end
        DA('start');
        pause(del);
        DA('stop');
        TRI = [-0.6:1/625:1.6 1.6:-1/625:-0.6]';
        TRI(end:25000)=-0.6;
        ff=0;
        SW=SWEEP_DATA(:,2:end);
        tmp=mean(SW(700:800,:),1);
        while ff==0
            gd=find (abs(tmp-mean(tmp))<std(tmp)*1.5);
            if length(gd)==length(tmp)
                ff=1;
            end
            SW=SW (:,gd);
            tmp=mean(SW(700:800,:),1);
        end
        TRI_COUNTER = mean(SW (:,gd),2);
        TRI_COUNTER (TRI_COUNTER>10)=10;
        TRI_COUNTER (TRI_COUNTER<-10)=-10;
        TRI_COUNTER(end:25000)=0;
        TRI=[TRI TRI_COUNTER];
    
    case 'activate'
        a=questdlg('Place the electrode in 150mM NaCl pH 9.5 and press OK.','Starting...','OK','Cancel','Cancel');
        if a(1)=='C'
            return;
        end
        TRI = 1.2;
        TRI(end:25000)=1.2;
        TRI_COUNTER=0;
        TRI_COUNTER(end:25000)=0;
        TRI=[TRI' TRI_COUNTER'];
        h=waitbar(0,'Activating electrode ...');
        DA('reset');
        DA('start');
        for i=1:300
            pause (1);
            waitbar(i/300,h);
        end
        close(h);
        DA('stop');
        TRI = 0;
        TRI(end:25000)=0;
        TRI_COUNTER=0;
        TRI_COUNTER(end:25000)=0;
        TRI=[TRI' TRI_COUNTER'];
        DA('reset');
        DA('start');
        pause(0.1);
        DA('stop');
        a=questdlg('Remove the electrode from solution and press OK.','Finished','OK', 'Done', 'Done');
        DA('initialize');
        
end

end


function storeData (~,event)

global SWEEP_DATA SWEEP_TIMES TIME0 DAQ TRI;
    
    SWEEP_TIMES(end+1)=etime(clock,TIME0);
    SWEEP_DATA(:,end+1)=event.Data(1:10000);
end

function queueData (src, ~)

global TRI
    src.queueOutputData([TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI; TRI])
end


function dispError (~, event)

disp(event.Error.getReport())
end