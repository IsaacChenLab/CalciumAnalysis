function ind=MC_accept(quest)

if nargin<1
    quest='Accept?';
end

s=questdlg(quest,'','Yes','No','Yes');
if isempty(s)
    error('no selection');
end

if strcmp(s,'Yes')
    ind=1;
else
    ind=0;
end

return;
