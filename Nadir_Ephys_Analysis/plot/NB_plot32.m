% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% June 22, 2016

% Edited from plot32.m

% Plots the averaged responses from all 32 channels in msec time scale.
% Modify offsetgraph to get a better separation of the responses.

offsetgraph = 500;

figure('units','normalized','outerposition',[0 0 1 1]);

hold on

for i = 1:32;  
   
    % Calculates x coordinate in miliseconds
    x = ((1:length(CSD_Data(i,:)))/32)+offsetms;
    
    % Calculates y coordinate based on CSD_Data (which is based on averaged 
    % CSC_wave data
    y = CSD_Data (i,:) - i*offsetgraph;
    
    h = plot(x,y);
end

% Holds axes 
axis tight;



% Set title based on folder
%currentDirectory = pwd;
%[upperPath, deepestFolder, ~] = fileparts(currentDirectory)
%figureTitle = strcat('Average CSC Responses for:', 32, deepestFolder);
%title(figureTitle,'interpreter','none');
%title('Average CSC Responses for: SPVAT03, TS 10477-10592, 470nm LED, 20ms, 0.5Hz, 50 stim');

%figure('KeyPressFcn',@KeyPress);
%function KeyPress(ObjH, EventData);
%Key = get(ObjH,'CurrentCharacter');
%switch Key
%    case char


%set(handles.mainWindow.'keypressfcn',{@currentWindowKeypress,handles});
%fucntion currentWindowKeyPress(hObject, eventdata, handles);

%if strcmp(currentKey, '