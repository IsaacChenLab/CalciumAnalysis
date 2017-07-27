function varargout = EF_spSimilarity(varargin)
% EF_SPSIMILARITY MATLAB code for EF_spSimilarity.fig
%      EF_SPSIMILARITY by itself, creates a new EF_SPSIMILARITY or raises the
%      existing singleton*.
%
%      H = EF_SPSIMILARITY returns the handle to a new EF_SPSIMILARITY or the handle to
%      the existing singleton*.
%
%      EF_SPSIMILARITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EF_SPSIMILARITY.M with the given input arguments.
%
%      EF_SPSIMILARITY('Property','Value',...) creates a new EF_SPSIMILARITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EF_spSimilarity_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EF_spSimilarity_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EF_spSimilarity

% Last Modified by GUIDE v2.5 23-Sep-2013 11:53:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EF_spSimilarity_OpeningFcn, ...
                   'gui_OutputFcn',  @EF_spSimilarity_OutputFcn, ...
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

% --- Executes just before EF_spSimilarity is made visible.
function EF_spSimilarity_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EF_spSimilarity (see VARARGIN)

% Choose default command line output for EF_spSimilarity
handles.output = 'Yes';

% Update handles structure
guidata(hObject, handles);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.
if isempty(varargin)
    return
end
global WAV SHPS PCAs SH1 SH2 RDBT INX

% avg extr1 extr2 flat noise pca pre1 pre2 pre3 pre4 pre5 pre6
WAV=varargin{1};
WAV=WAV-repmat(mean(WAV,1),size(WAV,1),1);
SHPS=mean(WAV,2);
v=SimComp (SHPS(:,1), WAV);
[~, p]=max(v);
SHPS(:,2)=WAV(:,p);
[~, p]=min(v);
SHPS(:,3)=WAV(:,p);
SHPS(:,4)=rand(size(SHPS,1),1);
SHPS(:,4)=SHPS(:,4)-mean(SHPS(:,4));
SHPS(:,5)=SHPS(:,4);
[~,p]=min(SHPS(:,1));
SHPS(p-1:p+1,5)=-4;
if exist('PCA_spike_shapes.mat','file')
    load('PCA_spike_shapes.mat');
else
    [pc, exp]=EF_ExtractPCA (WAV');
end
PCAs=pc;
SHPS(:,6)=PCAs(:,1);
if exist ('Preset_spike_shapes.mat','file')
    load ('Preset_spike_shapes.mat');
    SHPS(:,7:12)=SPKS;
else
    s=[71 72 81 82 91 92 101 102 111 112 121 122];
    for i=1:6
        st=eval(['handles.radiobutton' num2str(s(i))]);
        set(st,'Enable','Off');
    end
    SHPS(:,7:9)=0;
end
vm=max(SHPS,[],1)-min(SHPS,[],1);
SHPS=SHPS./repmat(vm,size(SHPS,1),1);
SHPS(:,4)=SHPS(:,4)./4;
for i=3:14
    st=eval(['handles.axes' num2str(i)]);
    plot (st,SHPS(:,i-2),'k');
    axis (st,[0 74 -1 1]);
end
RDBT=[11 21 31 41 51 61 71 81 91; 12 22 32 42 52 62 72 82 92];

for i=1:18
    st=eval(['handles.radiobutton' num2str(RDBT(i))]);
    set (st,'Value',0);
end

SH1=1;
set (handles.radiobutton11,'Value',1);
SH2=4;
set (handles.radiobutton42,'Value',1);
INX=ones (size(WAV,2),1);

SimPlot (handles);

% axis (handles.axes2,[0 1 0 1]);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes EF_spSimilarity wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = EF_spSimilarity_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.figure1);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global WAV SHPS SH1 SH2 INX

val1= SimComp (SHPS(:,SH1), WAV);
val2= SimComp (SHPS(:,SH2), WAV);
ind2=1;

f1=figure; hold on;

while ind2,
    plot(val1,val2,'.k','markersize',4);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    axis square;

    ind1=1;
    d=1;
    points=[];
    while ind1,
        [x,y]=ginput(1);
        points(d,:)=[x y];
        if d>1
            h=line([points(d-1,1) points(d,1)],[points(d-1,2) points(d,2)]);
            set(h,'color','b');
        end
        if d>3
            xxx=MC_spFindCross(points(1:d-1,1),points(1:d-1,2),points(d-1:d,1),points(d-1:d,2));
            if ~isempty(xxx)
                ind1=0;
                break;
            end
        end
        d=d+1;
    end

    IN=inpolygon(val1,val2,points(:,1),points(:,2));
    plot(val1(IN),val2(IN),'.b','markersize',4);
    s=questdlg('Accept?','','Yes','No','Yes');
    if strcmp(s,'Yes'),
        break;
        ind2=0;
    elseif strcmp(s,'No')
        cla;
    elseif isempty(s)
        error('');
    end
    
end
close(f1);
INX(IN)=1;
INX(~IN)=2;
SimPlot (handles);



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global WAV SHPS SH1 SH2 INX

val1= SimComp (SHPS(:,SH1), WAV);
val2= SimComp (SHPS(:,SH2), WAV);
INX(:)=1;
INX(val2>val1)=2;
SimPlot (handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global WAV SHPS SH1 SH2 INX

val1= SimComp (SHPS(:,SH1), WAV);
val2= SimComp (SHPS(:,SH2), WAV);

[INX,~,~] = kmeans([val1; val2]',2,'display','notify','replicates',20,'distance','sqEuclidean');
SimPlot (handles);



% --- Executes on button press in radiobutton91.
function radiobutton91_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton91 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton91
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=9;
SimPlot (handles);

% --- Executes on button press in radiobutton92.
function radiobutton92_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton92 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton92
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=9;
SimPlot (handles);


% --- Executes on button press in radiobutton81.
function radiobutton81_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton81
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=8;
SimPlot (handles);



% --- Executes on button press in radiobutton82.
function radiobutton82_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton82
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=8;
SimPlot (handles);


% --- Executes on button press in radiobutton72.
function radiobutton72_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton72
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=7;
SimPlot (handles);

% --- Executes on button press in radiobutton71.
function radiobutton71_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton71
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=7;
SimPlot (handles);

% --- Executes on button press in radiobutton61.
function radiobutton61_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton61
global  SHPS  SH1 SH2 PCAs

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
n=inputdlg('PCA: ','Choose the component number',1,{'1'});
try
    n=round(str2num(n{1}));
    if n<1 || n>50
        return
    end
catch
    return;
end
SHPS(:,6)=PCAs(:,n);
SH1=6;

st=handles.axes8;
plot (st,SHPS(:,6),'k');
axis (st,[0 74 -1 1]);
SimPlot (handles);

% --- Executes on button press in radiobutton62.
function radiobutton62_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton62
global  SHPS  SH1 SH2 PCAs

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
n=inputdlg('PCA: ','Choose the component number',1,{'1'});
try
    n=round(str2num(n{1}));
    if n<1 || n>50
        return
    end
catch
    return;
end
SHPS(:,6)=PCAs(:,n);
SH2=6;

st=handles.axes8;
plot (st,SHPS(:,6),'k');
axis (st,[0 74 -1 1]);

SimPlot (handles);


% --- Executes on button press in radiobutton51.
function radiobutton51_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton51
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=5;
SimPlot (handles);



% --- Executes on button press in radiobutton52.
function radiobutton52_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton52
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=5;
SimPlot (handles);



% --- Executes on button press in radiobutton42.
function radiobutton42_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton42
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=4;
SimPlot (handles);

% --- Executes on button press in radiobutton41.
function radiobutton41_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton41
global  SHPS  SH1 SH2 

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=4;
SimPlot (handles);

% --- Executes on button press in radiobutton111.
function radiobutton111_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton111 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton111
global  SHPS  SH1 SH2 

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=11;
SimPlot (handles);

% --- Executes on button press in radiobutton121.
function radiobutton121_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton121 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton121
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=12;
SimPlot (handles);

% --- Executes on button press in radiobutton112.
function radiobutton112_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton112 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton112
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=11;
SimPlot (handles);

% --- Executes on button press in radiobutton22.
function radiobutton22_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton22
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=2;
SimPlot (handles);


% --- Executes on button press in radiobutton21.
function radiobutton21_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton21
global  SHPS  SH1 SH2 

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=2;
SimPlot (handles);

% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton11
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=1;
SimPlot (handles);

% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton12
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=1;
SimPlot (handles);


function SimPlot (handles)
global WAV SHPS SH1 SH2 INX

hold (handles.axes2,'off');
plot (handles.axes2,SimComp (SHPS(:,SH1),WAV(:,INX==1)),SimComp (SHPS(:,SH2),WAV(:,INX==1)),'k.');
hold (handles.axes2,'on');
hold (handles.axes_sh,'off');
plot (handles.axes_sh,mean(WAV(:,INX==1),2),'k');
hold (handles.axes_sh,'on');
if ~isempty (INX==2)
    plot (handles.axes2,SimComp (SHPS(:,SH1),WAV(:,INX==2)),SimComp (SHPS(:,SH2),WAV(:,INX==2)),'r.');
    plot (handles.axes_sh,mean(WAV(:,INX==2),2),'r');
end
xs=get (handles.axes2,'XLim');
ys=get (handles.axes2,'YLim');
line ([-1 1],[-1 1],'Color','r', 'Parent', handles.axes2);
axis(handles.axes2,[xs ys]);

function rez = SimComp (sh, wav)

mat=repmat(sh,1,size(wav,2));
rez=dot (wav,mat,1)./sqrt(dot (wav,wav,1).*dot(mat,mat,1));


% --- Executes on button press in radiobutton101.
function radiobutton101_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton101 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton101
global  SHPS  SH1 SH2 

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=10;
SimPlot (handles);

% --- Executes on button press in radiobutton102.
function radiobutton102_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton102 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton102
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=10;
SimPlot (handles);

% --- Executes on button press in radiobutton122.
function radiobutton122_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton122
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=12;
SimPlot (handles);

% --- Executes on button press in radiobutton32.
function radiobutton32_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton32
global  SHPS  SH1 SH2

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH2) '2']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH2=3;
SimPlot (handles);

% --- Executes on button press in radiobutton31.
function radiobutton31_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton31
global  SHPS  SH1 SH2 

if ~get(hObject,'Value')
    return
end
oldb=eval(['handles.radiobutton' num2str(SH1) '1']);
set (oldb,'Value',0);
set (hObject,'Value',1);
SH1=3;
SimPlot (handles);
