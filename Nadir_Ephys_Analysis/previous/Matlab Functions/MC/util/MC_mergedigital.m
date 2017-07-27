function MC_mergedigital

global CAT;
global DATA_DIR;
global DATE;

[cats,dates,files]=MC_getDataForSpike;
for i=10:length(cats),
    dcat(cats{i},dates{i});
    
    fileName=MC_fileName(files{i});
    file=MC_fileNumber(fileName);  
    nfl=MC_multiFilesName(file);

    filess={sprintf('Sleep_%s.mat',nfl) , sprintf('Events_%s.mat',nfl) , sprintf('Spikes_%s.mat',nfl) };
    variables={'rem* sws*','Bit*','Sp*'};
    
    cd(DATA_DIR);
    cd('DigData');
    
    clear Sp* rem* sws* Bit* FileTimes;
    for d=1:length(filess),
        
        digFile=MC_digFileName(file);
        if exist(filess{d})
            load(filess{d});
            s=sprintf('save %s %s;',digFile,variables{d});
            eval(s);
        end
    end
        
    q=whos('*Mat*'); dd=[];
    for d=1:length(q),
        dd(d)=size(eval(q(d).name),1);
    end
    if ~all(dd==dd(1))
        error('stop');
    end    
        
end

