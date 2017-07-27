function varargout = CPP2(varargin)
% CPP2 MATLAB code for CPP2.fig
%      CPP2, by itself, creates a new CPP2 or raises the existing
%      singleton*.
%
%      H = CPP2 returns the handle to a new CPP2 or the handle to
%      the existing singleton*.
%
%      CPP2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CPP2.M with the given input arguments.
%
%      CPP2('Property','Value',...) creates a new CPP2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CPP2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CPP2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CPP2

% Last Modified by GUIDE v2.5 01-Apr-2013 13:28:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CPP2_OpeningFcn, ...
                   'gui_OutputFcn',  @CPP2_OutputFcn, ...
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


% --- Executes just before CPP2 is made visible.
function CPP2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CPP2 (see VARARGIN)

% Choose default command line output for CPP2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global TB TTB CUR TIMER TS BIN Running STPW

TB=[]; TB{3}(1,1)=0;
TB{2}(1,1)=0;
TB{1}(1,1)=0;
TTB=0;TTB(3)=0;
TS=[];
BIN=20;
TIMER=timer ('ExecutionMode','fixedRate','Period',1,'TimerFcn',{@UpdateTime, handles});
Running=0;
STPW=[];
CUR=2;
set(handles.BOX1,'Enable','off');
set(handles.BOX2,'Enable','off');
set(handles.BOX3,'Enable','off');

% UIWAIT makes CPP2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CPP2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BOX1.
function BOX1_Callback(hObject, eventdata, handles)
% hObject    handle to BOX1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TB TTB CUR  TS  STPW
tm=toc(STPW);
te=toc(TS);
STPW=tic;
TTB(CUR)=TTB(CUR)+tm;
TB{CUR}(end+1,1)=0;
TB{CUR}(end,2)=te;
CUR=1;
TB{1}(end+1,1)=1;
TB{1}(end,2)=te;

% --- Executes on button press in BOX2.
function BOX2_Callback(hObject, eventdata, handles)
% hObject    handle to BOX2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TB TTB CUR  TS  STPW
tm=toc(STPW);
te=toc(TS);
STPW=tic;
TTB(CUR)=TTB(CUR)+tm;
TB{CUR}(end+1,1)=0;
TB{CUR}(end,2)=te;
CUR=2;
TB{2}(end+1,1)=1;
TB{2}(end,2)=te;

% --- Executes on button press in BOX3.
function BOX3_Callback(hObject, eventdata, handles)
% hObject    handle to BOX3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TB TTB CUR  TS  STPW
tm=toc(STPW);
te=toc(TS);
STPW=tic;
TTB(CUR)=TTB(CUR)+tm;
TB{CUR}(end+1,1)=0;
TB{CUR}(end,2)=te;
CUR=3;
TB{3}(end+1,1)=1;
TB{3}(end,2)=te;

% --- Executes on button press in START.
function START_Callback(hObject, eventdata, handles)
% hObject    handle to START (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TB TTB CUR TIMER TS BIN Running STPW

if Running==0
    try
        BIN=str2num(get(handles.edit2,'String'));
    catch
        warndlg('Please enter a valid bin size in seconds.');
        return
    end
    set (handles.START,'String','STOP');
    Running=1;
    set(handles.BOX1,'Enable','on');
    set(handles.BOX2,'Enable','on');
    set(handles.BOX3,'Enable','on');
    TB=[]; TB{3}(1,1)=0;
    TB{2}(1,1)=0;
    TB{1}(1,1)=0;
    TTB=0;TTB(3)=0;
    CUR=2;
    TS=tic;
    STPW=TS;
    TB{2}(end+1,1)=1;
    TB{2}(end,2)=toc(TS);
    start(TIMER)
elseif Running==1
    stop(TIMER);
    tm=toc(STPW);
    te=toc(TS);
    TTB(CUR)=TTB(CUR)+tm;
    TB{CUR}(end+1,1)=0;
    TB{CUR}(end,2)=te;
    set (handles.START,'String','DONE');
    Running=3;
    set(handles.BOX1,'Enable','off');
    set(handles.BOX2,'Enable','off');
    set(handles.BOX3,'Enable','off');

    bns=floor(te/BIN)+2;
    pos=0:BIN:BIN*(bns-1);
    BOX=zeros(bns,3);
    BOX(1,1)=TTB(1);
    BOX(1,2)=TTB(2);
    BOX(1,3)=TTB(3);
    for j=1:3
        ins = TB{j}(TB{j}(:,1)==1,2);
        outs= TB{j}(TB{j}(:,1)==0,2);
        outs=outs(2:end);
        for i=1:length(ins)
            bn1=floor(ins(i)/BIN);
            bn2=floor(outs(i)/BIN);
            if bn1==bn2
                BOX(bn1+2,j)=BOX(bn1+2,j)+outs(i)-ins(i);
            else
                inl=ins(i);
                for z=bn1:bn2
                    tm=min((z+1)*BIN-inl, outs(i)-z*BIN);
                    BOX(z+2,j)=BOX(z+2,j)+tm;
                    inl=(z+1)*BIN;
                end
            end
        end
    end
    fl=get(handles.edit1,'String');
    save(fl,'BOX');
else
    delete(handles.figure1);
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


function UpdateTime (~, ~, hd)
global  TTB CUR STPW

tmp=toc(STPW);
t=TTB;
t(CUR)=t(CUR)+tmp;
set (hd.tb1,'String',num2str(t(1)));
set (hd.tb2,'String',num2str(t(2)));
set (hd.tb3,'String',num2str(t(3)));


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    t=1;
catch
end