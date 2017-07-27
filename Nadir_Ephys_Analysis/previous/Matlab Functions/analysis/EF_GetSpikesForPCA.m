function EF_GetSpikesForPCA (grades, nmpool, mxpool, FL)

% collects spike shapes from all currently sorted cells of 'grades' [1 2 3]
% It also averages the shapes of nmpool consecutive spikes of the same cell
% and adds up to mxpool final shapes to the result. The final output is
% saved in the FL file, in the form 73 x nnn matrix named SPKS
%
% Default working deroctory is Z:\EPhys\To Matlab
dfdir='Z:\EPhys\To Matlab';
r=dir(dfdir);
SPKS=zeros(73,1);
h=waitbar(0);
for i=3:length(r)
    if r(i).isdir
        cd ([dfdir '\' r(i).name]);
        waitbar(i/length(r));
        rr=dir('*_spk.mat');
        if ~isempty(rr)
            for j=1:length(rr)
                try
                    load (rr(j).name,'SpikeShapes', 'Grades');
                    for z=1:length(Grades)
                        if ~isempty(find(Grades(z)==grades))
                            sp=SpikeShapes{z};
                            sz=size(sp,2);
                            n=0;
                            while n<=mxpool && (n+1)*nmpool<sz
                                SPKS(:,end+1)=mean(sp(1:73,n*nmpool+1:(n+1)*nmpool),2);
                                n=n+1;
                            end
                        end
                    end
                catch
                end
            end
        end
    end
end
SPKS=SPKS(:,2:end);
save ([dfdir '\' FL],'SPKS');
close (h)
