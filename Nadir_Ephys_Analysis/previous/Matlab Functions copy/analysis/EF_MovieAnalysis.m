function varargout = EF_MovieAnalysis(varargin)
% EF_MOVIEANALYSIS MATLAB code for EF_MovieAnalysis.fig
%      EF_MOVIEANALYSIS, by itself, creates a new EF_MOVIEANALYSIS or raises the existing
%      singleton*.
%
%      H = EF_MOVIEANALYSIS returns the handle to a new EF_MOVIEANALYSIS or the handle to
%      the existing singleton*.
%
%      EF_MOVIEANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EF_MOVIEANALYSIS.M with the given input arguments.
%
%      EF_MOVIEANALYSIS('Property','Value',...) creates a new EF_MOVIEANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EF_MovieAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EF_MovieAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EF_MovieAnalysis

% Last Modified by GUIDE v2.5 23-Oct-2013 14:25:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EF_MovieAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @EF_MovieAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EF_MovieAnalysis is made visible.
function EF_MovieAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EF_MovieAnalysis (see VARARGIN)

% Choose default command line output for EF_MovieAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EF_MovieAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EF_MovieAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global MOV POS

a=get (hObject,'Value');
POS=round(MOV.NumberOfFrames*a);
POS=min(max(POS,1),MOV.NumberOfFrames);
UpdateDisplay;


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sl_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to sl_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
UpdateDisplay;

% --- Executes during object creation, after setting all properties.
function sl_zoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in but_loadmovie.
function but_loadmovie_Callback(hObject, eventdata, handles)
% hObject    handle to but_loadmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MOVIE_FL MOV POS HAND

[f,p,~] = uigetfile('*.avi');
MOVIE_FL=[p f];
set (handles.movie_name,'String',f);
MOV = VideoReader(MOVIE_FL);
f=read(MOV,1);
axes(handles.axes_movie);
image (f);
set(gca,'xtick',[]); set(gca,'ytick',[]);
POS=1;
set (handles.TIMEP,'String','0.001');
set (handles.slider1,'SliderStep',[2/MOV.NumberOfFrames 0.01]);
HAND=handles;



% --- Executes on button press in but_setfolderep.
function but_setfolderep_Callback(hObject, eventdata, handles)
% hObject    handle to but_setfolderep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CELLS FLD EVS COM COR
[~,FLD,~] = uigetfile('*_spk.mat');
f=find(FLD=='\');
set (handles.ep_dir,'String',FLD(f(end-1)+1:end-1));
if exist ([FLD '\events.mat'],'file')
    load ([FLD '\events.mat']);
    EVS=EVENTS;
    COM=COMMENTS;
    try
        COM{end+1}=' + (new)';
    catch
        c{1}=COM;
        COM=c;
        COM{end+1}=' + (new)';
    end
    COR=CORRECTION;
else
    q=questdlg ('Do you want to synchronize movie and recording?','Recording is NOT synchronized!','Yes','No','Yes');
    if strcmp(q,'No')
        figure1_CloseRequestFcn(handles.figure1, eventdata, handles);
    end
    q=questdlg ('Choose method:','Method','Function','Points','Function');
    if strcmp(q,'Function')
        i=inputdlg ('Enter function : ', 'Provide conversion function',1,{'DIG=DIG-DIG(.)+.'});
        load ([FLD 'digi01.mat'],'DIG');
        DIG=find(DIG);
        d=DIG;
        eval ([i{1} ';']);
        EVS{1}=d;
        EVS{1}(2,:)=DIG;
        COR=1;
    else
        load ([FLD 'digi01.mat'],'DIG');
        DIG=find(DIG)./1000;
        rep=1;
        while rep
            ans=inputdlg ('Light ON at (sec.):','Enter known timepoints',1,{'0 11.3'});
            try
                pnts=str2num(ans{1});rep=0;
            catch
                warndlg ('Please enter numbers separated by space.','Wrong format');
            end
        end
        ddig=diff(DIG);
        for i=2:length(pnts)
            if pnts(i)-pnts(i-1)>1.5*min(ddig)
                xx=floor((pnts(i)-pnts(i-1))/min(ddig));
                intt=(pnts(i)-pnts(i-1))/xx;
                pnts(end+1:end+xx-1)=pnts(i-1)+intt:intt:pnts(i)-1;
            end
        end
        pnts=unique(pnts);
        md=mean(ddig(1:end-4));
        mp=mean(diff(pnts));
        COR=mp/md;
        d=DIG;
        DIG=DIG-DIG(1)+pnts(1)/COR;
        DIG=1000*DIG;
        EVS{1}=d*1000;
        EVS{1}(2,:)=DIG;
    end
    COM{1}='Sync.';
    COM{end+1}=' + (new)';
end
d=dir([FLD '*_spk.mat']);
CELLS=[];

lst=[];
n=0;
for i=1:length(d)
    load ([FLD d(i).name],'SpikeTimes','Grades');
    for j=1:length(Grades)
        n=n+1;
        CELLS{n}=SpikeTimes{j};
        f=find(d(i).name=='_');
        lst{n}=['    ' d(i).name(f-2:f-1) '            ' num2str(j) '            ' num2str(Grades(j))];
    end
end

UpdateFiringRates(str2num(get(handles.edit1,'String')));
set(handles.CELLS,'String',lst);
set(handles.EVENTS,'String',COM);

% --- Executes on selection change in CELLS.
function CELLS_Callback(hObject, eventdata, handles)
% hObject    handle to CELLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CELLS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CELLS


% --- Executes during object creation, after setting all properties.
function CELLS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CELLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in but_selectcells.
function but_selectcells_Callback(hObject, eventdata, handles)
% hObject    handle to but_selectcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SEL

g=get(handles.CELLS,'Value');
if length(g)<1 || length(g)>5
    warndlg ('Please select a number of cells between 1 and 5.', 'Selection out of range','modal');
    return;
end
SEL=g;
UpdateDisplay;

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MOV POS  FRS CELLS FLD EVS COM SEL

q=questdlg ('Are you sure?', 'Exit program','Yes','No','No');
if strcmp(q,'Yes')
    MOV=[];
    POS=[];
    FRS=[];
    CELLS=[];
    FLD=[];
    EVS=[];
    COM=[];
    SEL=[];
    figure1_CloseRequestFcn(handles.figure1, eventdata, handles);
end

% --- Executes on button press in but_stepb.
function but_stepb_Callback(hObject, eventdata, handles)
% hObject    handle to but_stepb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MOV POS

if POS-1<1
    return;
end
POS=POS-1;
UpdateDisplay;

% --- Executes on button press in but_stepf.
function but_stepf_Callback(hObject, eventdata, handles)
% hObject    handle to but_stepf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MOV POS

if POS+1>MOV.NumberOfFrames
    return;
end
POS=POS+1;
UpdateDisplay;

% --- Executes on button press in but_play.
function but_play_Callback(hObject, eventdata, handles)
% hObject    handle to but_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TIMER MOV HAND

a=get(handles.but_play,'String');
if strcmp(a,'PLAY')
    TIMER=timer;
    set(TIMER,'ExecutionMode','fixedRate');
    set(TIMER,'TasksToExecute',1000000);
    set(TIMER,'TimerFcn',@PlayTM);
    set(TIMER,'BusyMode','drop');
    set(TIMER,'Period',1);
    set (handles.but_play,'String','STOP');
    start(TIMER);
else
    set (handles.but_play,'String','PLAY');
    stop(TIMER);
end

% --- Executes on selection change in EVENTS.
function EVENTS_Callback(hObject, eventdata, handles)
% hObject    handle to EVENTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EVENTS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EVENTS


% --- Executes during object creation, after setting all properties.
function EVENTS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EVENTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in but_selectevents.
function but_selectevents_Callback(hObject, eventdata, handles)
% hObject    handle to but_selectevents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SELD COM EVS

g=get(handles.EVENTS,'Value');
if length(g)<1 || length(g)>5
    warndlg ('Please select a number of events between 1 and 5.', 'Selection out of range','modal');
    return;
end
if g(end)==length(COM)
    an=inputdlg ('Name: ','Create a new event?',1,{''});
    if isempty(an)
        return;
    end
    COM{end}=an{1};
    COM{end+1}=' + (new)';
    EVS{end+1}=[];
    set (handles.EVENTS,'String',COM);
    set (handles.EVENTS,'Value',length(COM));
end
SELD=g;
UpdateDisplay;

% --- Executes on button press in but_addevent.
function but_addevent_Callback(hObject, eventdata, handles)
% hObject    handle to but_addevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global POS EVS COM MOV COR

s=get(handles.EVENTS,'Value');
if length(s)>1 || s==length(COM)
    warndlg ('Please select only one valid event to mark!','Incorrect selection','modal');
    return;
end
p=POS/COR;
tm=p*1000/MOV.FrameRate-EVS{1}(2,1)+EVS{1}(1,1);
if isempty(find(EVS{s}==tm))
    EVS{s}(end+1)=tm;
    EVS{s}=sort(EVS{s});
end
UpdateDisplay;

% --- Executes on button press in but_deleteevents.
function but_deleteevents_Callback(hObject, eventdata, handles)
% hObject    handle to but_deleteevents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global POS EVS COM MOV COR

s=get(handles.EVENTS,'Value');
if length(s)>1 || s==length(COM) || s(1)==1
    warndlg ('Please select only one valid event to delete!','Incorrect selection','modal');
    return;
end
p=POS/COR;
tm=p*1000/MOV.FrameRate-EVS{1}(2,1)+EVS{1}(1,1);
if ~isempty(find(EVS{s}==tm))
    EVS{s}=EVS{2}(EVS{2}~=tm);
    EVS{s}=sort(EVS{s});
else
    warndlg ('No event at this timepoint.', 'Action not possible.','modal');
end
UpdateDisplay;



% --- Executes on button press in but_saveevents.
function but_saveevents_Callback(hObject, eventdata, handles)
% hObject    handle to but_saveevents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLD EVS COM COR


EVENTS=EVS;
for i=1:length(COM)-1
    COMMENTS{i}=COM{i};
end
CORRECTION=COR;
save ([FLD '\events.mat'],'EVENTS','COMMENTS','CORRECTION');


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in but_setresolution.
function but_setresolution_Callback(hObject, eventdata, handles)
% hObject    handle to but_setresolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateFiringRates(str2num(get(handles.edit1,'String')));
UpdateDisplay;


function UpdateDisplay
global MOV POS  FRS CELLS FLD EVS COM SEL SELD HAND COR

handles=HAND;
cols_cells={'k' 'b' 'r' 'g' 'm'};
cols_dig = {'r' 'c' 'g' 'y' 'm'};
if ~isempty(MOV)
    f=read(MOV,POS);
    axes(handles.axes_movie);
    image (f);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    set(handles.slider1,'Value',POS/MOV.NumberOfFrames);
    set (handles.TIMEP,'String',num2str(POS/MOV.FrameRate,'%f'));
end
p=POS/COR;
if ~isempty(SEL)
    axes(handles.axes_spikes);
    zm=get(handles.sl_zoom,'Value');
    ofs=(EVS{1}(2,1)-EVS{1}(1,1))/MOV.FrameRate;
    per=MOV.NumberOfFrames-zm*(MOV.NumberOfFrames-MOV.FrameRate);
    x1=round(max(1,p-(per)/2));
    x2=min(x1+round(per),size(FRS,2));
    if x2-x1<per
        x1=x2-round(per);
        x1=max(1,x1);
    end
    hold off;
    plot (x1,0,'w.');hold on;
    mx=max(max(FRS(SEL,:)));
    for i=1:length(SEL)
        plot (FRS(SEL(i),x1:x2),cols_cells{i});
    end
    if ~isempty(SELD) && SELD(end)~=length(COM)
        for i=1:length(SELD)
            if ~isempty(EVS{SELD(i)})
                ev=EVS{SELD(i)}(1,:);
                ev=ev-EVS{1}(1,1)+EVS{1}(2,1);
                ev=ev(1,find(ev(1,:)>=x1*1000/MOV.FrameRate & ev(1,:)<=x2*1000/MOV.FrameRate));
                if ~isempty(ev)
                    ev=ev*MOV.FrameRate/1000-x1+1; %%%
                    plot (ev,mx/2,['*' cols_dig{i}]);
                end
            end
        end
    end
    line ([p-x1+1 p-x1+1],[0 mx],'Color','r'); %%%
    axis ([0 round(per) 0 mx]);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
end


function UpdateFiringRates (per)
global FRS MOV CELLS

FRS=zeros(1,MOV.NumberOfFrames);
for n=1:length(CELLS)
    fr=EF_FiringRate (CELLS{n}, per,1/MOV.FrameRate,1000);
    FRS(n,1:length(fr))=fr(:,1);
end

function PlayTM (~,~,~)
global TIMER POS MOV

if POS==MOV.NumberOfFrames
    stop(TIMER)
    set (handles.but_play,'String','PLAY');
    return;
end
POS=POS+2;
UpdateDisplay;

% --- Executes during object creation, after setting all properties.
function but_loadmovie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to but_loadmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global MOV POS  FRS CELLS FLD EVS COM SEL HAND
MOV=[];
POS=[];
FRS=[];
CELLS=[];
FLD=[];
EVS=[];
COM=[];
SEL=[];
HAND=handles;


% --- Executes on mouse press over axes background.
function axes_spikes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global MOV EVS POS
% 
% C = get (handles.axes_spikes, 'CurrentPoint');
% zm=get(handles.sl_zoom,'Value');
% ofs=(EVS{1}(2,1)-EVS{1}(1,1))/MOV.FrameRate;
% per=MOV.NumberOfFrames-zm*(MOV.NumberOfFrames-MOV.FrameRate);
% x1=round(max(1,POS-(per)/2));
% x2=min(x1+round(per),size(FRS,2));
% if x2-x1<per
%     x1=x2-round(per);
% end
