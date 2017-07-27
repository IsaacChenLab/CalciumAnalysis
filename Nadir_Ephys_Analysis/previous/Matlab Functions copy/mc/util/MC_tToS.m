function sample=MC_timeToSamples(time,Hz)

if nargin<2
    Hz=1e3;
end

sample=round(time*Hz);

return;
