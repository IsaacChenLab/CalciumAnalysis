[cats,dates,files]=MC_getDataForSpike;
for i=10:26,
    dcat(cats{i},dates{i});
    MC_evStart(files{i});
end


