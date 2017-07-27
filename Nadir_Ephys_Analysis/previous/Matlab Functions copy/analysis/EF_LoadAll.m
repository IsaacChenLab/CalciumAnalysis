function []=EF_LoadAll (ID, type, lst)

if nargin<3
    lst=[];
end

if ID(1)=='M'
    ID=['Z:\EPhys\To Matlab\' ID];
end

cd (ID)
if isempty(lst)
    switch lower(type)
        case 'lfps'
            d=dir ('lfps*');
        case 'spk'
            d=dir('*_spk.mat');
        otherwise
            return;
    end

    for i=1:length(d)
        nm=d(i).name;
        
        switch lower(type)
            case 'lfps'
                nr=str2num(nm(end-5:end-4));
                evalin('base', ['load ' nm]);
                evalin('base', ['lfp' num2str(nr) '=lfp;']);
            case 'spk'
                nr=str2num(nm(end-9:end-8));
                evalin('base', ['load ' nm]);
                evalin('base', ['spk' num2str(nr) '=SpikeTimes;']);
        end
    end
    return;
end
for i=1:length(lst)
    if lst(i)<10
        nm=[type '0' num2str(lst(i)) '.mat'];
    else
        nm=[type num2str(lst(i)) '.mat'];
    end
    evalin('base', ['load ' nm]);
    evalin('base', ['lfp' num2str(lst(i)) '=lfp;']);
end