function h=MC_sigCC(cc,qq,qq2)

if nargin<2
    qq=10;
end
if nargin<3
    qq2=5;
end

l=size(cc,2);
middle=round(l/2);
for_mean=cc(:,[1:qq end-qq:end]);
f_mean=mean(for_mean,2);
mid_inxs=middle-qq2+1 : middle+qq2;
f_means=repmat(f_mean,1,length(mid_inxs));
p_poiss=1-poisscdf(cc(:,mid_inxs)-1,f_means);
pst=p_poiss<0.05/length(mid_inxs);
h=any(pst,2);

return;
