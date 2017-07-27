function plotMeanWave(waves,channel,convert)

    meanwave = squeeze(mean(waves(channel,:,:),3));
    stdwave = squeeze(std(double(waves(channel,:,:)),0,3));
    
    hold on;
    plot(meanwave*convert,'k','LineWidth',2);
    plot((meanwave+stdwave)*convert,'k:','LineWidth',0.5);
    plot((meanwave-stdwave)*convert,'k:','LineWidth',0.5);
end
