function [cc,lags]=MC_xcov(v1,p1,p2,p3)
% function [cc,lags]=MC_xcov(v1,v2,Fs,cc_width)
% cc_width in seconds

Default_Fs=25000; % default MC_rack
Default_cc_width=1; %one second

if nargin==1 | (nargin>=2 & length(p1)==1), %auto
    if nargin==1
        cc_width = round(Default_cc_width*Default_Fs);
        Fs=Default_Fs;
    elseif nargin==3
        cc_width = round(p2*p1);
        Fs=p1;
    else %nargin==2
        cc_width = round(Default_cc_width*p1);
        Fs=p1;
    end
    
    [cc,lags]=xcov(v1,cc_width,'coeff');
    
else % cross
    if nargin==2
        cc_width = round(Default_cc_width*Default_Fs);
        Fs=Default_Fs;
    elseif nargin==3
        cc_width = round(p2*Default_cc_width);
        Fs=p2;
    else % nargin==4
        cc_width = round(p3*p2);
        Fs=p2;
    end
    
    [cc,lags]=xcov(v1,p1,cc_width,'coeff');
    
end

lags=lags/Fs;
plot(lags,cc);
set(gca,'XLim',[lags(1) lags(end)]);
set(gca,'XTick',[lags(1) 0 lags(end)]);
set(gca,'YLim',[min(cc)-range(cc)/5 max(cc)+range(cc)/5]);

return;




    
   


