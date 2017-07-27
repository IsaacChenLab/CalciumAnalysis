function status=BHV_workMem(RWD,ITI,ISI,PRESI,TRACEI,P,NCS)
% status=BHV_workMem(RWD,ITI,ISI,PRESI,TRACEI,P,NCS)

global EVENTS;
global EVENT_INX;
global SAVE_CLOCK_START;
EVENTS=[];

if nargin<1 | isempty(RWD), % length of reward
    RWD=0.1;
end
if nargin<2 | isempty(ITI), % inter trial interval
    ITI=30;
end
if nargin<3 | isempty(ISI), % interval between the stimuli
    ISI=2;
end
if nargin<4 | isempty(PRESI), % presentaion time of the visual stimulus
    PRESI=0.5;
end
if nargin<5 | isempty(TRACEI), % trace interval to the airpuff
    TRACEI=ISI;
end
if nargin<6 | isempty(P), % proportion of cs+
    P=0.5;
end
if nargin<7 | isempty(NCS), % number of cs-
    NCS=2;
end

close all;

% initiate the das1200 card 
CARD('Initialize');
CARD('lever_up');

% plot the screen
scrsz = get(0,'ScreenSize');
f1=figure('Position',round(scrsz),'color',[0 0 0]);
axes('position',[0 0 1 1]);
set(gca,'xlim',[0 1],'ylim',[0 1]); hold on;
set(gca,'visible','off');
drawnow;
% clear screen
area([0 1],[1 1],0,'facecolor',[0 0 0]); drawnow;

happy='n';
    while ~isempty(happy) & happy~='y'
        stimuliSet=BHV_workMem_genStimuliSet(NCS);
        for i=1:size(stimuliSet,1)
            disp(sprintf('%d %d\n',stimuliSet(i,1),stimuliSet(i,2)));
        end
        happy=input('Happy? ','s');
    end


d=0;
run_ind=1;
while run_ind,

    
    if mod(d,1000)==0
        [stimuli,Minus]=BHV_workMem_genStimuli(stimuliSet,1000,P);
        d=0;
    end
   
    d=d+1;
       
   [x,y,level]=BHV_workMem_quarter(stimuli(d,1));
   CARD('send_dio',20+stimuli(d,1)); % send the quarter thru the matching bit
   h=area(x,y,level,'facecolor',[1 1 1]);
   drawnow; 
   pause(PRESI);
   % clear screen
    area([0 1],[1 1],0,'facecolor',[0 0 0]); drawnow;

   pause(ISI);
   
   [x,y,level]=BHV_workMem_quarter(stimuli(d,2));
   CARD('send_dio',20+stimuli(d,2)); % send the quarter thru the matching bit
   h=area(x,y,level,'facecolor',[1 1 1]);
   drawnow;
   pause(PRESI);
   % clear screen
   area([0 1],[1 1],0,'facecolor',[0 0 0]); drawnow;

   pause(TRACEI);
   
    if ~Minus(d),
        CARD('send_dio',20);
        CARD('send_stimulation');
    end
        
    % inter trial interval
    iti=round(ITI+rand*ITI);
    pause(iti);
    
end

close(f1);
CARD('lever_down');
CARD('destroy');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimuli=BHV_workMem_genStimuliSet(NCS)

% all combinations
q=nchoosek(1:4,2);
q=[1 1; 2 2; 3 3; 4 4; q; [q(:,2) q(:,1)]]; 
% scramble them
x=randperm(size(q,1));
q=q(x,:);

% choose the CS+
t=randsample(size(q,1),1);
stimuli=q(t,:);
q(t,:)=[];

% choose the first CS- so that it has the 1st stimuli similar to the CS+
if NCS>0
    x=find(q(:,1)==stimuli(1,1),1);
    stimuli=[stimuli; q(x,:)];
    q(x,:)=[];
end

% choose the 2nd CS- so that it has the 2nd stimuli similar to the CS+
if NCS>1
    x=find(q(:,2)==stimuli(1,2),1);
    stimuli=[stimuli; q(x,:)];
    q(x,:)=[];
end

% add other CS- if asked
if NCS>2
    t=randsample(size(q,1),NCS-2);
    stimuli=[stimuli; q(t,:)];
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimuli,Minus]=BHV_workMem_genStimuli(stimuliSet,N,P)

ncsplus=round(P*N);
ncsminus=N-ncsplus;

% cs+
csplus=repmat(stimuliSet(1,:),ncsplus,1);

% cs-, all have same prob.
csminusset=stimuliSet(2:end,:);
csminus=repmat(csminusset,ceil(ncsminus/size(csminusset,1)),1);
csminus=csminus(1:ncsminus,:);
% indexes of csminus
minuses=repmat([1:size(csminusset,1)]',ceil(ncsminus/size(csminusset,1)),1);
minuses=minuses(1:ncsminus);

cs=[csplus; csminus];
Minus=[zeros(ncsplus,1); minuses];

% random scramble
x=randperm(N);
stimuli=cs(x,:);
Minus=Minus(x);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,level]=BHV_workMem_quarter(which)

switch which,
    case 1
        x=[0.5 1]; y=[1 1]; level=0.5;
    case 2
        x=[0 0.5]; y=[1 1]; level=0.5;
    case 3
        x=[0 0.5]; y=[0.5 0.5]; level=0;
    case 4
        x=[0.5 1]; y=[0.5 0.5]; level=0;
    otherwise
        error;
end

return;


