function y = importCSC_JAW (samples)

[m,n] = size(samples);

y = reshape (samples, 1, m*n);

end

