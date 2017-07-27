function [cc,cc_std]=MC_spikeCrossCorr(sp1,sp2,winLen)

if nargin<3
    winLen=500;
end

