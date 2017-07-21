% finds the location of every channel

% OUTPUT
%   Column 1 = channel number
%   Columns 2,3 = x and y coordinates of channel location


fprintf("\tchoose .kwik file...");
[data_file, data_path] = uigetfile('*.kwik','Choose .kwik file...');
fprintf("Selected!\n");
filename = strcat(data_path, data_file);

fprintf("\tSelect folder to save figures to...");
save_path = uigetdir('','Select folder to save figures to');
fprintf("Selected!\n");

channels = 0:31;
numChan = length(channels);
locPath = '/channel_groups/0/channels/';

x = zeros(numChan,1);
y = zeros(numChan,1);
 
for c = 1:numChan
    f = strcat(locPath, num2str(channels(c)));
    z = h5readatt(filename, f, 'position');
    x(c) = z(1);
    y(c) = z(2);
end
     
output = [channels' x y];
save(strcat(save_path,'/', 'channelLocations.mat'),'output');