function filterdata = CSCBandPass(data,frange,samplerate,order);

temp = zeros(size(data,1),size(data,2));
for i=1:size(data,1)
    temp(i,:) = buttfilt(data(i,:),frange,samplerate,'bandpass',order);
end

filterdata=temp;

end
