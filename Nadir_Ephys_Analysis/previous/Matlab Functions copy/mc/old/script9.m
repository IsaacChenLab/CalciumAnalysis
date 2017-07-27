[cats,dates,files]=MC_getDataForSpike;
for i=1:26,
    dcat(cats{i},dates{i});
    MC_swStart(files{i});
end


