% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% June 22, 2016

% Get all CSC files in dir
CSCfiles = dir('*.ncs');

if length(CSCfiles) > 0
    % Loop through each
    for id = 1:length(CSCfiles)
      % Get the file name (minus the extension)
        [~, f] = fileparts(CSCfiles(id).name);
        % Find where the _ suffix is and remove it from the temp name
        suffix = find(ismember(f,'_'),1,'last');
        if f(suffix) == '_'; f(suffix:end) = [];    
            % Rename file and add in the extension
            movefile(CSCfiles(id).name, strcat(f,'.ncs'));
        end
    end
end

% Get all TT files in dir
TTfiles = dir('*.ntt');

if length(TTfiles) > 0
    % Loop through each
    for id = 1:length(TTfiles)
        [~, f] = fileparts(TTfiles(id).name);
        suffix = find(ismember(f,'_'),1,'last');
        if f(suffix) == '_'; f(suffix:end) = [];          
            movefile(TTfiles(id).name, strcat(f,'.ntt'));
        end
    end
end


% Get all EV files in dir
EVfiles = dir('*.nev');

if length(EVfiles) > 0
    for id = 1:length(EVfiles)
        [~, f] = fileparts(EVfiles(id).name);
        suffix = find(ismember(f,'_'),1,'last');
        if f(suffix) == '_'; f(suffix:end) = [];
            movefile(EVfiles(id).name, strcat(f,'.nev'));
        end
    end
end


% Get all RAW files in dir
RAWfiles = dir('*.nrd');

if length(RAWfiles) > 0
    for id = 1:length(RAWfiles)
        [~, f] = fileparts(RAWfiles(id).name);
        suffix = find(ismember(f,'_'),1,'last');
        if f(suffix) == '_'; f(suffix:end) = []; 
            movefile(RAWfiles(id).name, strcat(f,'.nrd'));
        end
    end
end