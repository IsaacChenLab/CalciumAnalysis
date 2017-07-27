function varargout=MC_readArgs(inputs,defaults)

L=length(defaults)/2;

if nargout ~= L
    error('mismatch in number of output arguments and defaults');
end

for i=1:nargout,
    code=defaults{(i-1)*2+1};
    found=0;
    for j=1:length(inputs)
        if strncmp(code,inputs{j},length(code))
            varargout(i)=inputs(j+1);
            found=1;
            break;
        end
    end
    if ~found
        varargout(i)=defaults((i-1)*2+2);
    end
end

return;


