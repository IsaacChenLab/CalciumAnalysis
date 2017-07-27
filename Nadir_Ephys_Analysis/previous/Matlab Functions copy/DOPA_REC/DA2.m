function DA2 (operation, varargin)

global SWEEP_DATA SWEEP_TIMES TIME0 EVENTS AI AO TRI FILE TIMER DIO FIG SCALE;

switch lower(operation)
    
    case 'initialize'
        if isempty (varargin)
            d=date;
            d(d=='-')='_';
            
            FILE=['c:\DATA\TMP_' d '.daq'];
        else
            FILE=[varargin{1} '.daq'];
        end
        if length(varargin)<2
            nsetup=1;
        else
            nsetup=varargin{2};
        end
        
        if exist ('AI', 'var')
            try stop(AI);catch, end;
            delete (AI);
        end
        if exist ('AO', 'var')
            try stop(AO);catch, end;
            delete (AO);
        end
        if exist ('TIMER', 'var')
            try stop(TIMER);catch, end;
            delete (TIMER);
        end
        if exist ('DIO', 'var')
            try stop(DIO);catch, end;
            delete (DIO);
        end
        
        AI=analoginput('nidaq','dev1');
        ch=addchannel(AI,0:nsetup-1);
        set(AI,'SampleRate',200000);
        set(AI,'SamplesPerTrigger',5500);
        set(AI,'TriggerType','Software');
        set(AI,'TriggerChannel',ch(1));
        set(AI,'TriggerCondition','Rising');
        set(AI,'TriggerConditionValue',-0.3);
        set(AI,'TriggerDelayUnits','Samples');
        set(AI,'TriggerDelay',-400);
        set(AI,'TriggerRepeat',inf);
        set(AI,'LoggingMode','Disk&Memory');
        set(AI,'LogFileName',FILE);
        set(AI,'LogToDiskMode','Overwrite');
        set(AI,'TriggerFcn',{@DataStore,clock});
        
        AO=analogoutput ('nidaq','Dev1');
        addchannel(AO,0:nsetup*2-1);
        set(AO,'SampleRate',200000);
        set(AO,'TriggerType','Immediate');
        
        TIMER=timer;
        set(TIMER,'ExecutionMode','fixedRate');
        set(TIMER,'Period',0.1);
        set(TIMER,'TasksToExecute',inf);
        set(TIMER,'TimerFcn',{@DataStart});   %{@DataStart} 'disp (''Hi'')'
        set(TIMER,'BusyMode','drop');
        set(TIMER,'StartDelay',0.1);
        
        DIO=digitalio('nidaq','Dev1');
        names={'licks','reward','csp1','csm1','csp2','csm2'};
        addline(DIO,0:5,'in',names);
        set (DIO,'TimerPeriod',0.01);
        set (DIO,'TimerFcn',{@DigStore,DIO, clock});
        
        EVENTS=sparse(zeros(6,1));
        TRI=[-0.6:1/500:1.4 1.4:-1/500:-0.6];
%         TRI=[-0.6:1/400:(2*500/400-.6) (2*500/400-.6):-1/400:-0.6];
        TRI(end+1:5000)=-.6;
        TRI(4000:4500)=-.4;
        TRI=[TRI' zeros(5000,1)];
        TIME0=clock;
        SWEEP_TIMES=[];
        SWEEP_DATA=zeros(5500,1);
        warndlg ('Board is initialized','Done.');
        
    case 'reset'
        EVENTS=sparse(zeros(6,1));
        TIME0=clock;
        SWEEP_TIMES=[];
        SWEEP_DATA=zeros(5500,1);
        
    case 'start'
        DA2('reset');
        FIG=figure;
        start(AI);
        start(TIMER);
        
    case 'stop'
        try stop (TIMER); catch end
        try stop (AI); catch end
        try close (FIG); catch end
        
    case 'activate'
        DA2('reset');
             f=FIG;
            set(AI,'LoggingMode','Memory');
            set(AI,'TriggerFcn',{@DataStoreMem});
       
        if varargin{1}==1
            TRI=[-0.6:1/500:1.4 1.4:-1/500:-0.6];
    % %         TRI=[-0.6:1/250:3.4 3.4:-1/250:-0.6];
            TRI(1:5000)=1.2;
            TRI=[TRI' zeros(5000,1)];
            a=questdlg ('Place the electrode in the activating solution and press OK','Activating','OK','CANCEL','CANCEL');
            if strcmp(a, 'OK')
                DA2('start');
                pause(300);
                DA2('stop');
                beep;
            end
        elseif varargin{1}==2
            a=questdlg ('Place electrode in PBS 1x and press OK','Activating','OK','CANCEL','CANCEL');
            if strcmp(a, 'OK')
                TRI=[0:3/2500:3 2.9:-3/2500:0];
                TRI(end+1:5000)=0;
                TRI=[TRI' zeros(5000,1)];
                DA2('start');
                pause(20);
                DA2('stop');
                beep;
                TRI=[];
                TRI(1:5000)=1.5;
                TRI=[TRI' zeros(5000,1)];
                DA2('start');
                pause(20);
                DA2('stop');
                beep;
            end
        end
        warndlg ('Process completed.', 'Done!');
        
        FIG=f;
        TRI=[-0.6:1/500:1.4 1.4:-1/500:-0.6];
% %         TRI=[-0.6:1/250:3.4 3.4:-1/250:-0.6];
        TRI(end+1:5000)=-.6;
        TRI(4000:4500)=-.4;
        TRI=[TRI' zeros(5000,1)];
        
    case 'background'
        DA2('reset');
        TRI=[-0.6:1/500:1.4 1.4:-1/500:-0.6];
% %         TRI=[-0.6:1/250:3.4 3.4:-1/250:-0.6];
        TRI(end+1:5000)=-.6;
        TRI(4000:4500)=-.4;
        TRI=[TRI' zeros(5000,1)];
        per=varargin{1};
        set(AI,'LoggingMode','Memory');
        f=FIG;
        set(AI,'TriggerFcn',{@DataStoreMem});
        DA2('start');
        pause(per);
        DA2('stop');
        FIG=f;
        set(AI,'LoggingMode','Disk&Memory');
        set(AI,'TriggerFcn',{@DataStore,clock});
        
        
        ff=0;
        data=SWEEP_DATA(:,2:end);
        
        SWEEP_DATA=zeros(5000,1);
        sh=TRI(3500:5000,1);sh=sh-mean(sh);
        for i=1:size(data,2)
            x=xcorr (sh,data(3500:5000,i)-mean(data(3500:5000,i)),500);
            [~, p]=max(x);
            if p>501
                SWEEP_DATA(1:p-501,end+1)=0;
                SWEEP_DATA(p-500:p-500+4499,end)=data(1:5000,i);
            else
                SWEEP_DATA(1:4500,end+1)=data(502-p:502-p+4499,i);
            end
        end
        SW=SWEEP_DATA(:,2:end);
        [~, p]=max(mean(SW(1:3000,:),2));
        SW=SW(p-1000:p+2500,:);
        tmp=mean(SW(700:800,:),1);
        while ff==0
            gd=find (abs(tmp-mean(tmp))<std(tmp)*1.7);
            if length(gd)==length(tmp)
                ff=1;
            end
            SW=SW (:,gd);
            tmp=mean(SW(700:800,:),1);
        end
        shp = mean(SWEEP_DATA (:,gd+1),2);
        SCALE=shp(4499)-shp(3999);
        shp = mean(SW (:,gd),2);
        shp(end+1:5000)=0;
        TRI=[-0.6:1/500:1.4 1.4:-1/500:-0.6];
%         TRI=[-0.6:1/250:3.4 3.4:-1/250:-0.6];
        TRI(end+1:5000)=-.6;
        TRI(4000:4500)=-.4;
        TRI=[TRI' shp];
end

end

function DataStore (~, ~, clk)

global SWEEP_TIMES TIME0 AI FIG

SWEEP_TIMES(end+1)=etime(clk,TIME0);
try
    dat=getdata(AI);
    figure(FIG);
    plot(dat);
catch
end;
end

function DataStoreMem (~, ~)

global SWEEP_DATA AI FIG

data=getdata(AI);
SWEEP_DATA(:,end+1)=data;
figure(FIG);
plot (data);
end

function DataStart (~, ~)

global TRI AO

stop (AO);
putdata (AO,TRI);
start(AO);
end

function DigStore (~, ~, di, clk)

global TIME0 EVENTS

pos=round(1000*etime(clk,TIME0));
dat=getvalue(di);
EVENTS(dat==1,pos)=1;
end