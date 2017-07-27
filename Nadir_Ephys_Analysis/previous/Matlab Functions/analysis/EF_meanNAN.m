function rez=EF_meanNAN(m,dim)

s=size(m);
if length(s)==1
    rez=[];
else
    str='zeros(';
    for i=1:length(s)
        if dim~=i
            str=[str num2str(s(i)) ','];
        end
    end
    str=[str(1:end-1) ')'];
    rez=eval(str);
end
