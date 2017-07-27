function EF_RasterPlot (varargin)

%ref, dat, wnd, fs, fig, char,bn
ref=varargin{1};
dat=varargin{2};
wnd=varargin{3};
fs=varargin{4};
bn=1;
char='k.';
fig=[];
if length(varargin)>4
    for i=5:2:length(varargin)-1
        op=varargin{i};
        vl=varargin{i+1};
        switch lower(op)
            case 'figure'
                fig=vl;
            case 'char'
                char=vl;
            case 'bin'
                bn=vl;
        end
    end
end
wnd=wnd*fs;

if isempty(fig)
    figure;    
else
    figure (fig);
 
end
subplot(4,4,1:12); hold on
for i=1:length(ref)
    d=dat(dat>ref(i)+wnd(1))-ref(i)-wnd(1)+1;
    if ~isempty(d(d<wnd(2)-wnd(1)))
        plot (d(d<wnd(2)-wnd(1))+wnd(1),i,char);
    end
end
subplot(4,4,13:16); hold on
plot (EF_XCorr(ref,dat,wnd./fs,bn,fs),char(1));
