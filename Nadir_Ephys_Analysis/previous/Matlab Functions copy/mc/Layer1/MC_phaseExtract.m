function [phase]=MC_phaseExtract(y,method)

if nargin<2,
    method='extrema';
end
dy=diff(y)'; % signal derivative 
phase=nan*zeros(size(y,1),1);
switch (method),
    case 'minima',
        x=MC_runs(dy>0);
        inx=x(1,:)';
        phase(inx,1)=-pi;
        phase(inx(2:2:end),1)=pi;
    case 'maxima',
        x=MC_runs(dy<0);
        inx=x(1,:)';
        phase(inx,1)=0;
        phase(inx(2:2:end),1)=2*pi;
    case 'extrema',
        x=MC_runs(dy>0);
        inx=x(1,:)';
        phase(inx,1)=-pi;
        phase(inx(2:2:end),1)=pi;
        x=MC_runs(dy<0);
        inx=x(1,:)';
        phase(inx,1)=0;
    case 'up',
        x=MC_runs(y>0);
        inx=x(1,:)';
        phase(inx,1)=-pi/2;
        phase(inx(2:2:end),1)=nan;
    case 'down',
        x=MC_runs(y<0);
        inx=x(1,:)';
        phase(inx,1)=pi/2;
    case 'zero',
        x=MC_runs(y>0);
        inx=x(1,:);
        phase(inx,1)=-pi/2;
        x=MC_runs(y<0);
        inx=x(1,:)';
        phase(inx,1)=pi/2;        
    case 'hilbert'
        error('nothing yet');
        
    otherwise,
        error('method unrecognized');
        
end
phase(1)=0; phase(end)=0;
phase=interp1(find(~isnan(phase)),phase(~isnan(phase)),1:length(phase),'linear')';
switch method,
    case 'maxima'
        phase(find(phase>pi))=phase(find(phase>pi))-2*pi;
end

return;
