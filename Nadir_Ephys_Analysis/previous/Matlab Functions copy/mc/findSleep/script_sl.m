[cats,dates,files]=MC_getDataForSpike;
CH=17;
for i=14:length(cats),    
    dcat(cats{i},dates{i});
    MC_slStart(files{i},CH,25);
end