function [ret] = NI2 (operation,varargin)


ret=nan;
global DIO AI AO AI_chan AO_chan;
global LICKS TIME0 EVENTS TIMER;
global VL RWRD NActive EX;


vol=[0.08 0.08 0.08 0.08];  %0.0441; % *****************************


switch lower(operation)
    
    case 'initialize'
        if isempty(varargin)
            nsetup=1;
        else
            nsetup=varargin{1};
        end
        NActive=nsetup;
        if exist('AI','var') && ~isempty(AI)
            for i=1:length(AI)
                if strcmp(AI{i}.Running,'On'), stop(AI{i}); end
                delete(AI{i});
            end
            AI=[];
        end
        if exist('DIO','var'), delete(DIO); end
        if exist('AO','var'), delete(DIO); end
        if exist('TIMER','var')
            for i=1:length(TIMER)
                delete(TIMER{i});
            end
        end
        AI_chan=[];
        AO=[];
        AO_chan=[];
        
        % set the timer for reward delivery
        for i=1:nsetup
            TIMER{i}=timer;
            set(TIMER{i},'ExecutionMode','singleShot');
            set(TIMER{i},'TasksToExecute',1);
            set(TIMER{i},'TimerFcn',['NI2 (''rewardstop'',' num2str(i) ');']);
            set(TIMER{i},'BusyMode','queue');
        end
        
        % digital - pump control
        DIO=digitalio('nidaq','Dev2'); % *********************
        names={'pump1','pump2','pump3','d','e','f','g','h'};
        addline(DIO,0:7,'out',names);
        putvalue(DIO,zeros(1,8));
        
        % digital - lick monitoring
        addline(DIO,8:11,'in',{'lick1','lick2','lick3','lick4'});
        

        % analog in for lickometer
        for i=1:nsetup
            AI{i}=analoginput('nidaq','Dev2');
            AI_chan{i}=addchannel(AI{i},i-1);
            set(AI{i},'SampleRate',50); % 100 hz
            AI_chan{i}.InputRange = [-10 10];
            set(AI{i},'SamplesPerTrigger', 20);
        end
        
        %analog out for light
        AO=analogoutput('nidaq','Dev2');
        AO_chan=addchannel(AO,0,'light');
        set(AO,'SampleRate',100); % 100 hz
        AO_chan(1).OutputRange = [0 5];
        
        %initialize licking
        VL=zeros(nsetup,1);
        questdlg('Turn on lickometer and press OK.','LICKOMETER','OK','OK');
        VL(:)=1;
        
        TIME0=clock;
        EVENTS=[];
        EVENTS{nsetup,1}=[];
        LICKS=[];
        LICKS{nsetup}=[];
        RWRD=0;
        EX=[];

    case 'destroy'
        try
            putvalue(DIO.Line(1:8),0);
        catch
        end
        try
            for i=1:length(AI)
                stop(AI{i});
                stop(TIMER{i});
                delete (AI{i});
                delete(TIMER{i});
            end
        catch
        end
        AI=[];
        delete(DIO);
        delete (AO);
        AI_chan=[];
        AO=[];
        AO_chan=[];
        VL=[];
        EVENTS=[];
        LICKS=[];
        TIMER=[];
        NActive=0;
        clear DIO AI AO;
        
    case 'testlick'
        if isempty(varargin)
            setup=1;
        else
            setup=varargin{1};
        end
        f=figure;
        NI2 ('startlick',setup,f);
        for i=1:100
            figure(f)
            title ([num2str((100-i)/10) ' seconds remaining...']);
            pause (0.1);
        end
        NI2 ('stoplick',setup);
        bt= questdlg(['Were licks detected on setup ' num2str(setup) '?'],'LICKING','Yes','No','No');
        close (f);
        if strcmp(bt,'Yes')
            ret=1;
        else
            ret=0;
        end
        return;

        
    case 'startlick'
        if isempty(VL) || sum(VL)~=length(VL)
            disp('The lickometer has not been initialized!');
            return;
        end
        if isempty(varargin)
            setup=1;
        else
            setup=varargin{1};
        end
        try
            stop(DIO);
        catch
        end
        
        set (DIO,'TimerPeriod',0.01);
        if length(varargin)<2
            set (DIO,'TimerFcn',{@Lick,DIO,setup});
        else
            set (DIO,'TimerFcn',{@LickFig,DIO,setup,varargin{2}});
        end
        EX=ones(length(setup),1);
        start(DIO);
        
    case 'stoplick'
        stop(DIO);
        
    case 'reward'
        pumps=varargin{1};
        duration=varargin{2}.*vol(pumps);
        for i=1:length(pumps)
            pump=pumps(i);
            f=get(TIMER{pump},'Running');
            if strcmp(f,'on')
                return;
            end
            set(TIMER{pump},'StartDelay',duration(i));
            putvalue(DIO.Line(pump),1);
            es=etime(clock,TIME0);
            start(TIMER{pump});
            EVENTS{pump,1}(end+1)=es;
        end
    case 'rewardstop'
        putvalue(DIO.Line(varargin{1}),0);
        
    case 'light'
        putdata(AO,varargin{1});
        start(AO);
        
    case 'resetall'
        TIME0=clock;
        EVENTS=[];
        s=size(LICKS);
        EVENTS{s(2),1}=[];
        LICKS=[];
        LICKS{s(2)}=[];
        
    case 'resettime'
        TIME0=clock;
        
    case 'event'
        es=etime(clock,TIME0);
        setup=varargin{1};
        if varargin{2}==1
            disp ('Channel restricted for rewards');
            return
        end
        s=size(EVENTS);
        if varargin{2}>s(2)
            EVENTS{setup,varargin{2}}=[];
        end
        EVENTS{setup,varargin{2}}(end+1)=es;
        
    otherwise
        error('unfamiliar operation');
end
    

return


        
function LickFig(~, ~, DIO, setup, fig)
    global TIME0 EX;

    
    dat=getvalue(DIO);
    dat=dat(8+setup);
    if dat || dat==EX
        EX=dat;
        return;
    end
    es=etime(clock,TIME0);
    figure(fig);
    hold on;
    plot (es,0,'kx');
    EX=dat;
return;


function Lick(~, ~, DIO, setup)
    global LICKS TIME0 RWRD EX;

    dat=getvalue(DIO);
    dat=dat(8+setup);
    if sum(dat)==length(dat) || isequal(dat,EX)
        EX=dat;
        return;
    end
    EX=dat;
    es=etime(clock,TIME0);
    s=find (dat==0);
    l=length(s);
    for i=1:l
        LICKS{s(i)}(end+1)=es;
    end
    if RWRD~=0
        NI2('reward',setup,RWRD);
        RWRD=0;
    end
   
return;