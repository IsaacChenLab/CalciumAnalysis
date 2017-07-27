function [maxamp ch] = maxChannel(waves)

% returns maxamp - amplitude of the peak point across all channels
% ch - channel of that point (first channel = 1)
meanwave = squeeze(mean(waves,3));
[maxamp ch] = max(max(meanwave'));

end
