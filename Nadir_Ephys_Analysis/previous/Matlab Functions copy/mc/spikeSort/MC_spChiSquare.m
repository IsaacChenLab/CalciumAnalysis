function [inx,top]=MC_spChiSquare(cwv,cf,ci)
% function [inx,top]=MC_spChiSquare(cwv,cf)
% choose according to chi-square distribution of squared distances from average form
% cwf is an m*n matrix, where m are samples (rows), and n (columns) are
% observations (e.g. spike shapes: m samples and n spikes).
% plots the distribution if an axes handle is provided in cf.
% returns the list of indexes out of the n observations that are below the
% 0.95 cumulutive density function.

if nargin<3
    ci=0.95;
end
if nargin<2
    cf=nan;
end

mn=mean(cwv');
df=cwv'-repmat(mn,size(cwv,2),1);
nm=sum(abs(df').^2).^(1/2);
[nhist,xhist]=hist(nm,15);
[ab,abconf]=gamfit(nm,ci);
top=gaminv(ci,ab(1),ab(2));
inx=find(nm<top);
if ~isnan(cf)
    axes(cf);
    hold on;
    bar(xhist,nhist,'b');
    ylabel('# of spikes'); xlabel('Sum of squares');
    plot([top top],[0 max(nhist)],'r','LineWidth',2);
end

return;

