function [A,C] = FC_crosscorr2(s)
% Functional connectivity using cross-correlation

% Pairwise correlations between neurons calculated using the maximum normalized cross-correlation.
% Onset of calcium events are used to generate a fluorescence traces, Fi and Fj for neurons i and j.
% The maximum normalized cross-correlation between i and j at lags 0 to
% 500ms is recorded.
% For statistical testing of functional connection between i and j,
% surrogate onset times are generated for j by random permutation.
% Surrogate fluorescence trace for j is computed. Normalized
% cross-correlation between Fi and Fsur_j is computed. This is done 1000
% times. If original max cross-correlation is >99 percentile of surrogate,
% a functional connection between i and j exists
%
%%
%fprintf('\tUsing original method for FC but with native threshold on \n');

load('params.mat');
thresh = params.FC_thresh;
%fprintf("xcorr threshold = %.2f\n", thresh);

%fprintf("applying threshold\n");
%fprintf("active cells only\n");

q = 0;

[N,T] = size(s.dF_cell);
A = zeros(N);
C = zeros(N);
Nsur = params.FC.CC.Nsur;
maxlag = params.FC.CC.maxlag; % 500 ms
upd = textprogressbar(N);
for i=1:N
    multiWaitbar('Functional connectivity: cross-correlation',i/N);
    upd(i);
    %if(params.parallel)
    if(false)
        parfor j=1:N
            if(i~=j && ~isempty(s.Spikes_cell{i}) && ~isempty(s.Spikes_cell{j}))
                
                Csur = zeros(Nsur,1);
                Fi = Surrogate_Fluorescence(s.Spikes_cell{i},T,s.fps);
                Fj = Surrogate_Fluorescence(s.Spikes_cell{j},T,s.fps);
                r = xcorr(Fi,Fj,ceil(maxlag*s.fps),'coeff');
                C(i,j) = max(r);
                for k=1:Nsur
                    Nspks = length(s.Spikes_cell{j});
                    sur_spks = sort(randsample(T,Nspks));
                    F_sur_j = Surrogate_Fluorescence(sur_spks,T,s.fps);
                    r = xcorr(Fi,F_sur_j,ceil(maxlag*s.fps),'coeff');
                    Csur(k) = max(r);
                    
                end
                if(C(i,j)>prctile(Csur,99))
                    A(i,j) = 1;
                end
            end
        end
    else
        for j=1:N
             if(i~=j && ~isempty(s.Spikes_cell{i}) && ~isempty(s.Spikes_cell{j}))
                
                dF_i = s.dF_cell(i,:);
                dF_j = s.dF_cell(j,:);
                
                r = max(xcorr(dF_i,dF_j,ceil(maxlag*s.fps),'coeff'));
                
                %q = q+1;
                
                 if(r > thresh)
                     A(i,j) = 1;
                     C(i,j) = r;
                     q = q+1;
                     %fprintf("found connection %d\n", q);                     
                 end
                
%                 Csur = zeros(Nsur,1);
%                 Fi = Surrogate_Fluorescence(s.Spikes_cell{i},T,s.fps);
%                 Fj = Surrogate_Fluorescence(s.Spikes_cell{j},T,s.fps);
%                 r = xcorr(Fi,Fj,ceil(maxlag*s.fps),'coeff');
%                 C(i,j) = max(r);
%                 
%                 
%                 for k=1:Nsur
%                     Nspks = length(s.Spikes_cell{j});
%                     sur_spks = sort(randsample(T,Nspks));
%                     F_sur_j = Surrogate_Fluorescence(sur_spks,T,s.fps);
%                     r = xcorr(Fi,F_sur_j,ceil(maxlag*s.fps),'coeff');
%                     Csur(k) = max(r);
%                     
%                 end
%                 
%                 if(C(i,j)>prctile(Csur,99))
%                     A(i,j) = 1;
%                     q = q+1;
%                 else
%                     C(i,j) = 0;
%                 end

                
            end
        end
    end
end

fprintf("number of connections %d\n", q);