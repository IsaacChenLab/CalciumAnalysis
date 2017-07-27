function varargout = EF_EventSelection(varargin)
% EF_EVENTSELECTION MATLAB code for EF_EventSelection.fig
%      Function takes as parameter a cell string of events (timepoints),
%      where each cell contains the timestamps of one event. Timescale can
%      be any format, but consitent across events. The rezult is a string of timepoints
%      that correspond to the selected main event, filtered by the selected
%      rules.
%      Example: rez=EF_EventSelection({[1 3 7...]; [0.2 0.8 1.4...]}); you
%      can select rules to get type 1 events that occur more than 4 sec
%      appart and are not followed by type 2 events for 0.2 sec.
%
%
%      EF_EVENTSELECTION, by itself, creates a new EF_EVENTSELECTION or raises the existing
%      singleton*.
%
%      H = EF_EVENTSELECTION returns the handle to a new EF_EVENTSELECTION or the handle to
%      the existing singleton*.
%
%      EF_EVENTSELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EF_EVENTSELECTION.M with the given input arguments.
%
%      EF_EVENTSELECTION('Property','Value',...) creates a new EF_EVENTSELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EF_EventSelection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EF_EventSelection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EF_EventSelection

% Last Modified by GUIDE v2.5 02-Apr-2013 13:03:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EF_EventSelection_OpeningFcn, ...
                   'gui_OutputFcn',  @EF_EventSelection_OutputFcn, ...
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


% --- Executes just before EF_EventSelection is made visible.
function EF_EventSelection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EF_EventSelection (see VARARGIN)

% Choose default command line output for EF_EventSelection
handles.output = hObject;
global EVS REF REZ

EVS=varargin{1};
REF=2;
REZ=EVS{REF};
set (handles.popupmenu7,'Value',2);
set (handles.text9,'String',['# : ' num2str(length(REZ))]);
PlotEvs (handles,REF);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EF_EventSelection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EF_EventSelection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global REZ
varargout{1} = REZ;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
if get(hObject,'Value')
    set (handles.popupmenu2,'Enable','on');
    set (handles.present2,'Enable','on');
    set (handles.absent2,'Enable','on');
    set (handles.edit2,'Enable','on');
else
    set (handles.popupmenu2,'Enable','off');
    set (handles.present2,'Enable','off');
    set (handles.absent2,'Enable','off');
    set (handles.edit2,'Enable','off');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
if get(hObject,'Value')
    set (handles.popupmenu4,'Enable','on');
    set (handles.present4,'Enable','on');
    set (handles.absent4,'Enable','on');
    set (handles.edit4,'Enable','on');
else
    set (handles.popupmenu4,'Enable','off');
    set (handles.present4,'Enable','off');
    set (handles.absent4,'Enable','off');
    set (handles.edit4,'Enable','off');
end


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
if get(hObject,'Value')
    set (handles.popupmenu3,'Enable','on');
    set (handles.present3,'Enable','on');
    set (handles.absent3,'Enable','on');
    set (handles.edit3,'Enable','on');
else
    set (handles.popupmenu3,'Enable','off');
    set (handles.present3,'Enable','off');
    set (handles.absent3,'Enable','off');
    set (handles.edit3,'Enable','off');
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global REF REZ EVS

s=get (handles.pushbutton1,'String');
if s(1)=='S'
    set (handles.pushbutton1,'String','RE-SET');
    set(handles.checkbox1,'Enable','on');
    set(handles.checkbox2,'Enable','on');
    set(handles.checkbox3,'Enable','on');
    set(handles.checkbox4,'Enable','on');
    set(handles.popupmenu7,'Enable','off');
else
    set (handles.pushbutton1,'String','SET');
    set(handles.checkbox1,'Enable','off');
    set(handles.checkbox2,'Enable','off');
    set(handles.checkbox3,'Enable','off');
    set(handles.checkbox4,'Enable','off');
    set(handles.popupmenu7,'Enable','on');
    
end
REF=get(handles.popupmenu7,'Value');
REZ=EVS{REF};
set (handles.text9,'String',['# : ' num2str(length(REZ))]);
PlotEvs (handles,REF)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin ('base','global REZ;');
delete (handles.figure1);



% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PlotEvs (handles, rf)

global EVS REZ REF
cold={'k.'; 'b.'; 'g.'; 'r.'; 'm.'; 'c.'};
col={'k'; 'b'; 'g'; 'r'; 'm'; 'c'};
if nargin<2 || isempty(rf)
    RZ=REZ;
    rf=get(handles.popupmenu7,'Value');
else
    RZ=EVS{rf};
end
s=get(handles.edit7,'String');
f=find(s==':');
try
    p1=str2num(s(1:f-1));
    p2=str2num(s(f+1:end));
catch
    warndlg ('Incorrect values for Peri-stim. time.');
    return;
end
try
    hold (handles.axes1, 'off');
    hold (handles.axes2, 'off');
    for i=1:length(EVS)
        plot (handles.axes1,EVS{i},i,cold{i});
        hold (handles.axes1, 'on');
        t=EF_XCorr(RZ*1000,EVS{i}*1000,[p1 p2],100,1000);
        plot (handles.axes2,p1:.1:(p2-.1),t/max(t),col{i});
        hold (handles.axes2, 'on');
    end
    plot (handles.axes1,RZ,rf,'ro');
    
catch
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateRez (handles);
s=get(handles.popupmenu7,'Value');
ss=get(handles.pushbutton1,'String');
if ss(1)=='S'
    PlotEvs (handles,s)
else
    PlotEvs(handles)
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

if get(hObject,'Value')
    set (handles.popupmenu1,'Enable','on');
    set (handles.present1,'Enable','on');
    set (handles.absent1,'Enable','on');
    set (handles.edit1,'Enable','on');
else
    set (handles.popupmenu1,'Enable','off');
    set (handles.present1,'Enable','off');
    set (handles.absent1,'Enable','off');
    set (handles.edit1,'Enable','off');
end


function UpdateRez (handles)

global EVS REZ REF

oREZ=REZ;
REZ=EVS{REF};
if get (handles.checkbox1,'Value')
    rf=get (handles.popupmenu1,'Value');
    pr=get (handles.present1,'Value');
    s=get(handles.edit1,'String');
    f=find(s==':');
    try
        p1=str2num(s(1:f-1));
        p2=str2num(s(f+1:end));
    catch
        warndlg ('Incorrect values for Rule 1, time interval');
        REZ=oREZ;
        return;
    end
    nREZ=[];
    for i=1:length(REZ)
        df=EVS{rf}-REZ(i);
        if xor (isempty (find (df>= p1 & df<=p2)), pr)
            nREZ(end+1)=REZ(i);
        end
    end
    REZ=nREZ;
end

if get (handles.checkbox2,'Value')
    rf=get (handles.popupmenu2,'Value');
    pr=get (handles.present2,'Value');
    s=get(handles.edit2,'String');
    f=find(s==':');
    try
        p1=str2num(s(1:f-1));
        p2=str2num(s(f+1:end));
    catch
        warndlg ('Incorrect values for Rule 2, time interval');
        REZ=oREZ;
        return;
    end
    nREZ=[];
    for i=1:length(REZ)
        df=EVS{rf}-REZ(i);
        if xor (isempty (find (df>= p1 & df<=p2)), pr)
            nREZ(end+1)=REZ(i);
        end
    end
    REZ=nREZ;
end

if get (handles.checkbox3,'Value')
    rf=get (handles.popupmenu3,'Value');
    pr=get (handles.present3,'Value');
    s=get(handles.edit3,'String');
    f=find(s==':');
    try
        p1=str2num(s(1:f-1));
        p2=str2num(s(f+1:end));
    catch
        warndlg ('Incorrect values for Rule 3, time interval');
        REZ=oREZ;
        return;
    end
    nREZ=[];
    for i=1:length(REZ)
        df=EVS{rf}-REZ(i);
        if xor (isempty (find (df>= p1 & df<=p2)), pr)
            nREZ(end+1)=REZ(i);
        end
    end
    REZ=nREZ;
end

if get (handles.checkbox4,'Value')
    rf=get (handles.popupmenu4,'Value');
    pr=get (handles.present4,'Value');
    s=get(handles.edit4,'String');
    f=find(s==':');
    try
        p1=str2num(s(1:f-1));
        p2=str2num(s(f+1:end));
    catch
        warndlg ('Incorrect values for Rule 4, time interval');
        REZ=oREZ;
        return;
    end
    nREZ=[];
    for i=1:length(REZ)
        df=EVS{rf}-REZ(i);
        if xor (isempty (find (df>= p1 & df<=p2)), pr)
            nREZ(end+1)=REZ(i);
        end
    end
    REZ=nREZ;
end
set (handles.text9,'String',['# : ' num2str(length(REZ))]);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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
