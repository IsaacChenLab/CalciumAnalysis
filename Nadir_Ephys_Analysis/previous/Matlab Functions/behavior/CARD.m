function [ret] = CARD(operation,varargin)
% function CARD(operation,varargin)

% panel |    das1200    | lines | def as | purpose
%------------------------------------------
% 51:58 - A0:A7 (port0) -  1:8  - output - board A
% 59:66 - B0:B7 (port1) -  9:16 - output - board B
% 67:70 - C0:C3 (port2) - 17:20 - output - stimulator, pump, to-lever,
% to-acq-bit-1
% 71:74 - C4:C7 (port3) - 21:24 - output - to-acq-bit-2,3,4,5

ret=nan;
global DIO AI ai_chan;
global LICKS;

switch lower(operation)
    
    case 'initialize'
        if exist('AI','var') & isstruct('AI')
            if strcmp(AI.Running,'On'), stop(AI); end
            delete(AI); 
        end
        if exist('DIO','var'), delete(DIO); end
        
        
        % digital
        DIO=digitalio('mcc',2);

        % port 0 dedicated to board A stimulation box, and port 1 to board
        % B, beside the last bit (16)
        names={'a','b','c','d','strobe1','strobe2','polarity','enable'}; % enable is only for port 0
        hA=addline(DIO,0:7,0,'out',names);
        hB=addline(DIO,0:7,1,'out',names);
        putvalue(DIO.Line(1:8),0); % notice that this enables the stimulation board (because 0 goes to 8)
        putvalue(DIO.Line(9:16),0);

        % constant high for the lickometer
        putvalue(DIO.Line(16),1);

        % port 2 is for activating the stimulator, the pump, constant
        % high to the lever, and 1 bit for acquisition
        hStim=addline(DIO,0:3,2,'out');
        putvalue(DIO.Line(17:20),[1 1 0 0]); %pump is activated high, and so is airpuff (stimulator)

        % port 3 is 4 bits to acquisition
        hStim=addline(DIO,0:3,3,'out');
        putvalue(DIO.Line(21:24),[0 0 0 0]);
        
        % analog for lever and lickometer
        AI=analoginput('mcc',2);
        ai_chan=addchannel(AI,0:1,{'lever','lickometer'});
        set(AI,'SampleRate',100); % 100 hz
        ai_chan(1).InputRange = [-10 10]; % lever
        ai_chan(2).InputRange = [-10 10]; % lickometer

        % make sure the reward is off
        putvalue(DIO.Line(18),1);
        % constant high for the lickometer
        putvalue(DIO.Line(16),1);
        
    case 'destroy'
        putvalue(DIO.Line(1:8),0);
        putvalue(DIO.Line(9:16),0);
        putvalue(DIO.Line(17:20),[1 1 0 0]); 
        delete(DIO);
        delete(AI);
        clear DIO AI;
        
    case 'acquire_and_wait'
        % varargin{1} is time in ms
        set(AI,'SamplesPerTrigger', round( varargin{1} / 10));
        set(AI,'TriggerRepeat',1);
        start(AI);
        while strcmp(AI.Running,'On') ; end
        [data_vec tm]=getdata(AI);
        ret=data_vec;
    
    case 'send_dio'
        putvalue(DIO.Line(varargin{1}),1);
        %pause(0.1);
        putvalue(DIO.Line(varargin{1}),0);

    case 'start_acquisition'
        putvalue(DIO.Line(20),1);
        
    case 'end_acquisition'
        putvalue(DIO.Line(20),0);        
        
    case 'lever_up'
        putvalue(DIO.Line(19),1);
        LICKS=0;

    case 'lever_down'
        putvalue(DIO.Line(19),0);
        
    case 'set_stimulus'
        
        cs=varargin{1};
        if any(cs>16) | any(cs<1)
            %error('stimulating electrodes must be between 1 and 16');
        end
        polarity=0;
        
        % board A, first electrode
        vec=dec2binvec(cs(1),4);
        v=[vec 1 0 polarity 0]; % strobe1 is activated, and strobe2 not
        putvalue(DIO.Line(1:8),v);
        % board A, second electrode
        vec=dec2binvec(cs(2),4);
        v=[vec 0 1 polarity 0]; % now strobe2 is on
        putvalue(DIO.Line(1:8),v);
        % board B, first electrode
        vec=dec2binvec(cs(3),4);
        v=[vec 1 0 polarity 0]; % now strobe2 is on
        putvalue(DIO.Line(9:16),v);
        % board B, second electrode
        vec=dec2binvec(cs(4),4);
        v=[vec 0 1 polarity 0]; % now strobe2 is on
        putvalue(DIO.Line(9:16),v);

    case 'send_stimulation'
        putvalue(DIO.Line(17),0);
        %pause(0.05);
        putvalue(DIO.Line(17),1);
    
    case 'reward'
        putvalue(DIO.Line(18),0);
        pause(varargin{1});
        putvalue(DIO.Line(18),1);
        
    case 'trial_start'
        putvalue(DIO.Line(20),1);
        %pause(0.01);
        putvalue(DIO.Line(20),0);
        
    case 'trial_end'
        % signal acquisition and reset values
        putvalue(DIO.Line(1:8),0);
        putvalue(DIO.Line(9:16),0);
        putvalue(DIO.Line(17:20),[0 1 0 0]); 

    case 'lever_sound_reward'
        set(AI,'SamplesPerTrigger', 1);
        set(AI,'TriggerRepeat',0); % only one press per trial
        set(AI,'TriggerChannel',ai_chan);
        set(AI,'TriggerFcn',{@Trigger_Sound_Reward,DIO,AI,varargin{1},varargin{2}});
        set(AI,'TriggerType','Software');
        set(AI,'TriggerCondition','Rising');
        set(AI,'TriggerConditionValue',0.5); % threshold of input current
        start(AI);
        while strcmp(AI.Running,'On'), % infinite wait for a press
            ;
        end

    case 'lever_cs'
        set(AI,'SamplesPerTrigger', 1);
        set(AI,'TriggerRepeat',0); % only one press per trial
        set(AI,'TriggerChannel',ai_chan);
        set(AI,'TriggerFcn',{@Trigger_cs,DIO,AI});
        set(AI,'TriggerType','Software');
        set(AI,'TriggerCondition','Rising');
        set(AI,'TriggerConditionValue',0.5); % threshold of input current
        start(AI);
        while strcmp(AI.Running,'On'), % infinite wait for a press
            ;
        end

    case 'lever_bit'
        set(AI,'SamplesPerTrigger', 1);
        set(AI,'TriggerRepeat',1e9); % only one press per trial
        set(AI,'TriggerChannel',ai_chan(2));
        set(AI,'TriggerFcn',{@Send_lick,DIO,AI,varargin{1}});
        set(AI,'TriggerType','Software');
        set(AI,'TriggerCondition','Rising');
        set(AI,'TriggerConditionValue',0.5); % threshold of input current
        start(AI);

    case 'lever_reward'
        set(AI,'SamplesPerTrigger', 1);
        set(AI,'TriggerRepeat',0); % only one press per trial
        set(AI,'TriggerChannel',ai_chan);
        set(AI,'TriggerFcn',{@Trigger_Reward,DIO,AI,varargin{1}});
        set(AI,'TriggerType','Software');
        set(AI,'TriggerCondition','Rising');
        set(AI,'TriggerConditionValue',0.5); % threshold of input current
        start(AI);
        start_time=clock;
        while strcmp(AI.Running,'On'),
            if etime(clock,start_time) > varargin{2}
                stop(AI);
            end
        end

    case 'lever_reward_2'
        set(AI,'SamplesPerTrigger', 1);
        set(AI,'TriggerRepeat',0); % only one press per trial
        set(AI,'TriggerChannel',ai_chan);
        set(AI,'TriggerFcn',{@Trigger_Reward,DIO,AI,varargin{1}});
        set(AI,'TriggerType','Software');
        set(AI,'TriggerCondition','Rising');
        set(AI,'TriggerConditionValue',0.5); % threshold of input current
        start(AI);

    case 'lever_wait'
        set(AI,'SamplesPerTrigger', 1);
        set(AI,'TriggerRepeat',0); % only one press per trial
        set(AI,'TriggerChannel',ai_chan);
        set(AI,'TriggerFcn',{@Send_bit,DIO,AI});
        set(AI,'TriggerType','Software');
        set(AI,'TriggerCondition','Rising');
        set(AI,'TriggerConditionValue',0.5); % threshold of input current
        start(AI);

    case 'lever_stop'
        if strcmp(AI.Running,'On')
            stop(AI);
        end
    
    case 'check_lever'
        if strcmp(AI.Running,'On')
            ret=1;
        else
            ret=0;
        end
    otherwise
        error('unfamiliar operation');
end
    
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Send_bit(obj,event,DIO,AI,bit)

    CARD('send_dio',bit);
        
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Send_lick(obj,event,DIO,AI,bit)
    global LICKS;
    LICKS=LICKS+1;
    CARD('send_dio',bit);
        
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Trigger_cs(obj,event,DIO,AI)

    stop(AI);
        
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Trigger_Sound_Reward(obj,event,DIO,AI,duration,pitch)

    if nargin<6
        pitch=0;
    end

    N1=0.2;
    Fs=8e3;
    tone1=sin(pitch*[0:(1/Fs):N1]);
    tone1=[tone1' tone1']; % stereo

    CARD('send_dio',21);
    sound(tone1,Fs); % sound for debugging
    pause(1);
    p=rand;
    if p<0.8
        CARD('send_dio',20);
        CARD('reward',duration);
    end
    stop(AI);
        
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Trigger_Reward(obj,event,DIO,AI,duration)
    
    CARD('reward',duration);
    stop(AI);
        
return;




