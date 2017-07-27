function mat=EF_Smooth2 (mat, wnd, rep)

s=size(mat);
for z=1:rep
    for i=1:s(1)
        mat(i,:)=MC_smooth(mat(i,:)',wnd(1));
    end

    for i=1:s(2)
        mat(:,i)=MC_smooth(mat(:,i),wnd(2));
    end

end