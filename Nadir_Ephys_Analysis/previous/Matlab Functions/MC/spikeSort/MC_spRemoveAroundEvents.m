function [cData]=MC_spRemoveAroundEvents(cData,dData,Time,whatevents,evinterv)

if ~isempty(dData) & ~isnan(dData)
    for i=1:length(whatevents)
        x=find(dData==whatevents(i));
        re=MC_runs(x);
        for j=1:size(re,2), 
            rex=[re(1,j)+evinterv(i,1) : re(1,j)+re(2,j)+evinterv(i,2)];
            rex(rex<2)=[]; rex(rex>length(cData)-1)=[]; 
            edges=[rex(1)-1 rex(end)+1];
            cData(rex)=interp1(Time(edges),cData(edges),Time(rex),'linear');
        end
    end 
end

return;

