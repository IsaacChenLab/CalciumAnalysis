function MC_plotXcovPairs(matrix,Fs,cc_width,channels)
% function MC_plotXcovPairs(matrix)

[m,nchannels]=size(matrix);
nsubplot=nchannels;

if nargin<4
    channels=1:nchannels;
end

figure;
q=linspace(8,4,12);
fnt=round(q(nchannels));

for i=1:nchannels,
    for j=i:nchannels,
        subplot(nsubplot,nsubplot,(i-1)*nsubplot+j);
        hold on; set(gca,'FontSize',fnt);
        if i==j,
            if nargin==1,
                MC_xcov(matrix(:,i));
            elseif nargin==2
                MC_xcov(matrix(:,i),Fs);
            else
                MC_xcov(matrix(:,i),Fs,cc_width);
            end 
        else
            if nargin==1,
                MC_xcov(matrix(:,i),matrix(:,j));
            elseif nargin==2
                MC_xcov(matrix(:,i),matrix(:,j),Fs);
            else
                MC_xcov(matrix(:,i),matrix(:,j),Fs,cc_width);
            end 
        end
        %axis square;
    end
end

ah= axes('position',[0,0,1,1]);   
set(ah,'visible','off');
x=linspace(0,1,nchannels+3);
y=linspace(0,1,nchannels+2);
for i=1:nchannels,
    text(x(i+2)-0.03,0.96,num2str(channels(i)),'HorizontalAlignment','center');
    text(0.08,1-y(i+1),num2str(channels(i)),'HorizontalAlignment','center');
end


return;
