function [Inxs]=MC_spChiSelection(waveForms)

f1=figure;
% choose according to chi-square distribution of squared distances from average form

ci=0.99;

str_cf='0.95'; ci=0.99;
while 1,
    
    cf=subplot(3,3,9);
    cla;
    [Inxs,top]=MC_spChiSquare(waveForms,cf,ci);
    outinx=setdiff(1:size(waveForms,2),Inxs);
    wfout=waveForms(:,outinx);

    cf=subplot(3,3,[1 2 4 5]); cla; drawnow; hold on;
    MC_spPlotShape(waveForms,cf);
    plot(wfout,'r'); mxx=max(waveForms(:)); mnn=min(waveForms(:));
    text(3,mxx-round((mxx-mnn)/20),[sprintf('%2.1f',100*size(wfout,2)/size(waveForms,2)) '%'],'color','r');

    accept=questdlg('Accept?','','Yes',str_cf,'No','Yes');
    if isempty(accept)
        error('quit');
    elseif strcmp(accept,'No')
        Inxs=1:size(waveForms,2);
        break;
    elseif strcmp(accept,'Yes')
        break;
    elseif strcmp(accept,'0.95')
        str_cf='0.99'; ci=0.95;
    elseif strcmp(accept,'0.99')
        str_cf='0.95'; ci=0.99;
    end
end

close(f1);

return;

