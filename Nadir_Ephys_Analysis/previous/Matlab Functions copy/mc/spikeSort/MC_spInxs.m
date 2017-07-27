function new_inxs=MC_spInxs(old_inxs,inxs_groups)

if iscell(inxs_groups)
    for k=1:length(inxs_groups),
        new_inxs{k}=old_inxs(inxs_groups{k});
    end
else
    new_inxs=old_inxs(inxs_groups);
end

    

return;
