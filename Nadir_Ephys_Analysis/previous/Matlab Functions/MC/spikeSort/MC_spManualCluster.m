function [Inxs]=MC_spManualCluster(waveForms)


f1=figure;
% use inpolygon
MC_spPlotPCs(waveForms,'k',f1);
s=questdlg('PCs?','','1-2','1-3','2-3','1-2');
if isempty(s)
    error('no selection');
end
close(f1);

switch s
    case '1-2'
        inx1=1; inx2=2;
    case '1-3'
        inx1=1; inx2=3;
    case '2-3'
        inx1=2; inx2=3;
    otherwise
        error('unknown');
end

ind2=1;

f1=figure; hold on;

while ind2,
    
    mat=waveForms';
    [mat_c,pc,exp,lat,coef]=MC_spPCA(mat);
    plot(mat_c(:,inx1),mat_c(:,inx2),'.k','markersize',4);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    axis square;

    ind1=1;
    d=1;
    points=[];
    while ind1,
        [x,y]=ginput(1);
        points(d,:)=[x y];
        if d>1
            h=line([points(d-1,1) points(d,1)],[points(d-1,2) points(d,2)]);
            set(h,'color','b');
        end
        if d>3
            xxx=MC_spFindCross(points(1:d-1,1),points(1:d-1,2),points(d-1:d,1),points(d-1:d,2));
            if ~isempty(xxx)
                ind1=0;
                break;
            end
        end
        d=d+1;
    end

    IN=inpolygon(mat_c(:,inx1),mat_c(:,inx2),points(:,1),points(:,2));
    plot(mat_c(find(IN),inx1),mat_c(find(IN),inx2),'.b','markersize',4);
    s=questdlg('Accept?','','Yes','No','Yes');
    if strcmp(s,'Yes'),
        break;
        ind2=0;
    elseif strcmp(s,'No')
        cla;
    elseif isempty(s)
        error('');
    end
    
end
close(f1);

Inxs{1}=find(IN);
Inxs{2}=find(~IN);

return;
