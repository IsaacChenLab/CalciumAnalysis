function MC_spPlotPCs(wv,c,f)

if nargin<3
    f=figure;
end

figure(f);
mat=wv';
[mat_c,pc,exp,lat,coef]=MC_spPCA(mat);

num_coef=size(mat_c,2);

%nn=nchoosek(1:3,2);
nn=[1 2; 1 3; 2 3];
n=ceil(sqrt(size(nn,1)+1));

% variance explained
subplot(n,n,1); hold on;
cms=cumsum(exp); 
plot(cms(1:25),'*-'); 
xlabel('Coefficients'); ylabel('% of variance'); 
plot([num_coef+0.5 num_coef+0.5],[cms(1) cms(27)],'r');
set(gca,'xlim',[0.5 25.5]);
% set(gca,'xtick',[1 num_coef 20]);
set(gca,'YLim',[10*floor(cms(1)/10) 100]);


% plot first coefficients
for i=1:size(nn,1),
    subplot(n,n,i+1); hold on;
    plot(mat_c(:,nn(i,1)),mat_c(:,nn(i,2)),['.' c],'markersize',3);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    xlabel(num2str(nn(i,1)));
    ylabel(num2str(nn(i,2)));
    axis square;
end

return;

