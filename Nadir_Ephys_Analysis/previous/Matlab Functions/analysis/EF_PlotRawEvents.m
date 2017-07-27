function EF_PlotRawEvents (ev, data, fq, wnd, flt)


if flt
    data=MC_spHighPass(double(data),fq);
end
% data=double(data);
xs=0:1/fq:length(data)/fq;
xs=xs(1:end-1);
ln=round((wnd(2)-wnd(1))*fq);
xx=(wnd(1):1/fq:wnd(2)*2).*1000;
xx=xx(1:ln+1);
md=zeros(ln+1,1);
figure; hold on
for i=1:length(ev)
    f=round(fq*(ev(i)+wnd(1)));
    try
    d=data(f:f+ln);
    d=d-mean(d(0-round(wnd(1)*fq):0-round(wnd(1)*fq)+round(fq*0.01)));
    md=md+d./length(ev);
    plot (xx,d,'k')
    catch
    end
end
plot (xx,md,'r','LineWidth',3)