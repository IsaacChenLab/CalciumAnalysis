function varargout = CPP(varargin)
% CPP M-file for CPP.fig
%      CPP, by itself, creates a new CPP or raises the existing
%      singleton*.
%
%      H = CPP returns the handle to a new CPP or the handle to
%      the existing singleton*.
%
%      CPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CPP.M with the given input arguments.
%
%      CPP('Property','Value',...) creates a new CPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CPP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CPP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CPP

% Last Modified by GUIDE v2.5 01-Apr-2013 11:09:55

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_EVENTS: 1 - go in BOX1; 2 - go out BOX1; 3 - go in BOX2;
%           4 - go out BOX2; -1 - end of session; 0 - start session


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CPP_OpeningFcn, ...
                   'gui_OutputFcn',  @CPP_OutputFcn, ...
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


% --- Executes just before CPP is made visible.
function CPP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CPP (see VARARGIN)

% Choose default command line output for CPP
handles.output = hObject;

% Update handles structure
global T_TOTAL T_BOX1 T_BOX2 T_BEGIN S_BOX1 S_BOX2 B1_BEGIN B2_BEGIN T_EVENTS;

T_TOTAL=0;
T_BOX1=0;
T_BOX2=0;
T_BEGIN=[];
S_BOX1=0;
S_BOX2=0;
B1_BEGIN=[];
B2_BEGIN=[];
T_EVENTS=[];

guidata(hObject, handles);

% UIWAIT makes CPP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CPP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global T_TOTAL T_BOX1 T_BOX2 T_BEGIN S_BOX1 S_BOX2 B1_BEGIN B2_BEGIN T_EVENTS;

if ~isempty(T_BEGIN)
    %save data
    
    %fix last box in
    es=clock;
    stat=get (handles.pushbutton1,'String');
    if stat(2)~='a'
        T_TOTAL=etime(es, T_BEGIN);
        T_EVENTS(end+1,2)=-1;
        T_EVENTS(end,1)=T_TOTAL;
        if S_BOX1
            T_BOX1=T_BOX1+etime (es,B1_BEGIN);
        elseif S_BOX2
            T_BOX2=T_BOX2+etime (es,B2_BEGIN);
        end
    end
    set (handles.text5, 'String', [num2str(T_BOX1*100/T_TOTAL) ' %']);
    set (handles.text6, 'String', [num2str(T_BOX2*100/T_TOTAL) ' %']);
    set (handles.text4, 'String', 'Idle');
    set (handles.text3, 'String', 'Time (min): 0.00');
    set (handles.pushbutton1, 'String', 'Save');
    
    id=get(handles.edit3,'String');
    ss=get(handles.edit4,'String');
    s=get(handles.edit5,'String');
    fld=get(handles.edit6,'String');
    fname=[fld id '_' ss '.mat'];
    
    if exist(fname,'file')
        ans=questdlg ('File exists. Overwrite?', 'Conflict detected', 'Yes', 'No', 'No');
        if ans(1)=='N';
            return;
        end
    end
    
    eval (['save ' fname ' T_TOTAL T_BOX1 T_BOX2 T_EVENTS s']);
    set (handles.pushbutton1, 'String', 'Start');
    T_TOTAL=0;
    T_BOX1=0;
    T_BOX2=0;
    T_BEGIN=[];
    set (handles.edit1, 'String', 'OUT');
    set (handles.edit2, 'String', 'OUT');
    set (handles.edit5, 'String', '');
    set (handles.text5, 'String', '  ');
    set (handles.text6, 'String', '  ');
    S_BOX1=0;
    S_BOX2=0;
    B1_BEGIN=[];
    B2_BEGIN=[];
    
    T_EVENTS=[];
else
    T_BEGIN=clock;
    T_TOTAL=0;
    set (handles.edit1, 'String', 'OUT');
    set (handles.edit2, 'String', 'OUT');
    S_BOX1=0;
    S_BOX2=0;
    B1_BEGIN=[];
    B2_BEGIN=[];
    T_EVENTS(1,1:2)=0;
    set (handles.text4, 'String', 'Running');
    set (handles.text5, 'String', '  ');
    set (handles.text6, 'String', '  ');
    set (handles.pushbutton1, 'String', 'Stop&Save');
end;
    



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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global T_TOTAL T_BOX1 T_BOX2 T_BEGIN S_BOX1 S_BOX2 B1_BEGIN B2_BEGIN T_EVENTS;

es=clock;
if isempty(T_BEGIN)
    return
end
T_TOTAL=etime (es,T_BEGIN);
set (handles.text3,'String', ['Time (min): ' num2str(T_TOTAL/60)]);
if S_BOX2
    S_BOX2=0;
    T_BOX2=T_BOX2+etime (es,B2_BEGIN);
    B2_BEGIN=[];
    T_EVENTS(end+1,1)=etime (es,T_BEGIN);
    T_EVENTS(end,2)=4;
    set(handles.edit2,'String','OUT');
else
    S_BOX2=1;
    B2_BEGIN=clock;
    set(handles.edit2,'String','IN');
    if S_BOX1
        S_BOX1=0;
        T_BOX1=T_BOX1+etime (es,B1_BEGIN);
        B1_BEGIN=[];
        set(handles.edit1,'String','OUT');
        T_EVENTS(end+1,1)=etime (es,T_BEGIN);
        T_EVENTS(end,2)=2;
    end
    T_EVENTS(end+1,1)=etime (es,T_BEGIN);
    T_EVENTS(end,2)=3;
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global T_TOTAL T_BOX1 T_BOX2 T_BEGIN S_BOX1 S_BOX2 B1_BEGIN B2_BEGIN T_EVENTS;

es=clock;
if isempty(T_BEGIN)
    return
end
T_TOTAL=etime (es,T_BEGIN);
set (handles.text3,'String', ['Time (min): ' num2str(T_TOTAL/60)]);
if S_BOX1
    S_BOX1=0;
    T_BOX1=T_BOX1+etime (es,B1_BEGIN);
    B1_BEGIN=[];
    set(handles.edit1,'String','OUT');
    T_EVENTS(end+1,1)=etime (es,T_BEGIN);
    T_EVENTS(end,2)=2;
else
    S_BOX1=1;
    B1_BEGIN=clock;
    set(handles.edit1,'String','IN');
    if S_BOX2
        S_BOX2=0;
        T_BOX2=T_BOX2+etime (es,B2_BEGIN);
        B2_BEGIN=[];
        set(handles.edit2,'String','OUT');
        T_EVENTS(end+1,1)=etime (es,T_BEGIN);
        T_EVENTS(end,2)=4;
    end
    T_EVENTS(end+1,1)=etime (es,T_BEGIN);
    T_EVENTS(end,2)=1;
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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
