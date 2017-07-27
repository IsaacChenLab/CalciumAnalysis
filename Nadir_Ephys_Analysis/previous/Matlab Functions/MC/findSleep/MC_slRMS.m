function y=MC_slRMS(x,bins)

if nargin<2
    bins=50;
end

B=ones(1,bins);
B=B/sum(B);

for i=1:size(x,2)
    y(:,i)=sqrt(filtfilt(B,1,x(:,i).^2));
end 

return;
