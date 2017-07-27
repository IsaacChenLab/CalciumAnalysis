function [ao,bo,co]=MC_MakeList(ai,bi,ci)

if nargin~=nargout
    error('mismatch in number of arguments');
end

d=0;
for i=1:length(ai),
    for j=1:length(bi{i})
        for k=1:length(ci{i}{j})
            d=d+1;
            ao(d)=ai(i);
            bo(d)=bi{i}(j);
            co(d)=ci{i}{j}(k);
        end
    end
end

return;

            