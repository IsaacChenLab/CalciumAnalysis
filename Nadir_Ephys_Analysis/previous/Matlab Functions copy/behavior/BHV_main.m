function status=BHV_main(CS,n_trl,cs_freq)
% status=BHV_main(CS,n_trl,cs_freq)
% Main function of the behavioral program
% CS is a vector 8 [csplus_boardA csplus_boardB csminus_boardA
% csminus_boardB] where each number stands for an electrode between 1-16
% n_trl is the number of all trials
% csplus_freq is the frequency of cs+ trials, the rest are cs-

% initiate the das1200 card 
CARD('Initialize');

if nargin<3
    [CS,n_trl,cs_freq,trial_params]=BHV_getParam;
end


% generate random list for trials
[cslist reward]=BHV_GenList(CS,n_trl,cs_freq);

%%%%%% 
% behavioral loop, 
run_ind = 1;
h1=waitbar(0,['Running ' num2str(n_trl) ' trials']);
for cur_trial=1:n_trl,
    BHV_ProcessTrial(cur_trial,cslist(cur_trial,:),reward(cur_trial),trial_params);
    waitbar(cur_trial / n_trl);
    if run_ind==0,
        break;
    end
end

close(h1);
CARD('destroy');

return;
%%% end of BHV_main %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status=BHV_ProcessTrial(tr_inx,cs,reward,params)

% times are in seconds

tl=['Trial #' num2str(tr_inx)];

% signal the stimulation interaface box
CARD('Set_Stimulus',cs);

% send a start trial trigger to the acquisition
CARD('Trial_Start');

% wait a sec for recordings
hmsg=msgbox(sprintf('wait after trial start: %1.1f sec',params(2)),tl,'non-modal');
pause(params(2));

% activate stimulus
sound(sin(cs(1)*1000*[0:0.001:2])); % sound for debugging
pause(params(3));
CARD('Send_Stimulation');

% wait a sec, intended delay
msgbox(sprintf('wait after stimulation: %1.1f sec',params(4)),tl,'replace');
pause(params(4));

% check for lever pressing within X ms

CARD('lever_up');
msgbox(sprintf('wait for lever: %1.1f sec',params(5)),tl,'replace');
data = CARD('acquire_and_wait',params(5)*1000);
CARD('lever_down');
lever_pressed=0;
if any(data > 1), % identify that the lever was pressed
    lever_pressed=1;
end

% if pressed and cs+
if reward==1 & lever_pressed==1,
    msgbox(sprintf('reward time: %1.1f sec',params(1)),tl,'replace');
    CARD('Reward',params(1));
    beep;
    msgbox(sprintf('wait after reward: %1.1f sec',params(6)),tl,'replace');
    pause(params(6));
end

CARD('Trial_End');

% ITI
msgbox(sprintf('Inter trial interval: %1.1f sec',params(7)),tl,'replace');
pause(params(7));

close(hmsg);

% t0=clock;
% sec=etime(clock,t0);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cslist reward]=BHV_GenList(CS,n_trl,cs_freq)

cslist=zeros(n_trl,4);
reward=zeros(n_trl,1);

% number of cs plus and minus trials
csplus_num=round(cs_freq*n_trl);
csminus_num=n_trl-csplus_num;

% geberate random indexes for cs plus and minus
t1=rand(1,n_trl);
[Y,I]=sort(t1);
csplus_inx=I(1:csplus_num);
csminus_inx=I(csplus_num+1:end);

% put the cs plus and minus values in the trial order matrices
cslist(csplus_inx,:) = repmat(CS(1:4),csplus_num,1);
reward(csplus_inx)=1;
cslist(csminus_inx,:)= repmat(CS(5:8),csminus_num,1);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CS,n_trl,cs_freq,trial_params]=BHV_getParam()

while (1),
    CSplus=inputdlg({'Board A, elec 1','Board A, elec 2','Board B, elec 1','Board B, elect 2'},...
        'CS + stimilating elec',[1 1 1 1],{'1','2','1','2'});
    CSminus=inputdlg({'Board A, elec 1','Board A, elec 2','Board B, elec 1','Board B, elect 2'},...
        'CS -',[1 1 1 1],{'3','4','3','4'});
    CS=zeros(1,8);
    for i=1:4,
        CS(i)=str2num(CSplus{i});
        CS(i+4)=str2num(CSminus{i});
    end
    if all(CS)>=1 & all(CS)<=16
        break;
    end
end

while (1),
    s=inputdlg({'number of trials','frequency of CS +'},'Behavior Params',1,{'4','0.5'});
    n_trl=str2num(s{1});
    cs_freq=str2num(s{2});
    if n_trl>=1 & n_trl <= 9999 & cs_freq>=0 & cs_freq<=1
        break;
    end
end

while(1),
    s=inputdlg({'Reward time','Wait after start','Wait between sound and stimulation','Wait after stimulation',...
        'Wait for lever','Wait after reward','Inter Trial Interval'},'Trial parameters',1,...
        {'1','0.4','0.05','1','3','0.5','3'});
    for i=1:7,
        trial_params(i)=str2num(s{i});
    end
    if all(trial_params>=0) & all(trial_params<10)
        break;
    end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




