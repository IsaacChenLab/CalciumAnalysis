function time=MC_sToT(sample,Hz)

if nargin<2
    Hz=1e3;
end

time=round(sample/Hz);

return;
