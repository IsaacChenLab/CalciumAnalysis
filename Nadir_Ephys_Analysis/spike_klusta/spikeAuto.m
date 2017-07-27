function autoc = spikeAuto(tstamps,zerolag)

    k=1;
    a=zeros(length(tstamps)*(length(tstamps)-1),1);
    
    for i=1:(length(tstamps)-1)
        for j=1:length(tstamps)
            a(k)=tstamps(j)-tstamps(i);
            k=k+1;
        end
    end
    
    if zerolag==0
        a = a(a~=0);
    end
    
    autoc=a;
    
end