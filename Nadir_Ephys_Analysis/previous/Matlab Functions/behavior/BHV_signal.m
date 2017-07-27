function EVENTS=BHV_signal(event)

% 
global EVENTS;
global EVENT_INX;
global SAVE_CLOCK_START;

if nargin==0
    return;
end


if isempty(EVENTS),
    EVENT_INX=1;
    SAVE_CLOCK_START=clock;
else
    EVENT_INX=EVENT_INX+1;
end
disp(event);


switch event,
    
    case 'save_events'
        s=date;
        s=sprintf('save %s EVENTS;',s);
        eval(s);
                
    case 'cs_plus',
        EVENTS(EVENT_INX,1)=21;
        CARD('send_dio',21);

    case 'cs_minus',
        EVENTS(EVENT_INX,1)=22;
        CARD('send_dio',22);
      
    case 'success',
        EVENTS(EVENT_INX,1)=20;
        CARD('send_dio',20);

    case 'press',
        EVENTS(EVENT_INX,1)=23;
        CARD('send_dio',23);

    case 'ITI',
        EVENTS(EVENT_INX,1)=24;
        CARD('send_dio',24);

    otherwise,
        error('no such event');
            
end
EVENTS(EVENT_INX,2)=etime(clock,SAVE_CLOCK_START);

if mod(EVENT_INX,50)==0,
    BHV_signal('save_events');
end

return;
            
            
            