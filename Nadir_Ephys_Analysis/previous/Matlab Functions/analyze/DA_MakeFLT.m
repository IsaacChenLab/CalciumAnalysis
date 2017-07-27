function varargout = DA_MakeFLT(varargin)
% DA_MAKEFLT MATLAB code for DA_MakeFLT.fig
%      DA_MAKEFLT, by itself, creates a new DA_MAKEFLT or raises the existing
%      singleton*.
%
%      H = DA_MAKEFLT returns the handle to a new DA_MAKEFLT or the handle to
%      the existing singleton*.
%
%      DA_MAKEFLT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DA_MAKEFLT.M with the given input arguments.
%
%      DA_MAKEFLT('Property','Value',...) creates a new DA_MAKEFLT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DA_MakeFLT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DA_MakeFLT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DA_MakeFLT

% Last Modified by GUIDE v2.5 10-Jul-2012 12:18:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DA_MakeFLT_OpeningFcn, ...
                   'gui_OutputFcn',  @DA_MakeFLT_OutputFcn, ...
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


% --- Executes just before DA_MakeFLT is made visible.
function DA_MakeFLT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DA_MakeFLT (see VARARGIN)

global VISIBLE AX_pos AX_len FLT
% Choose default command line output for DA_MakeFLT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DA_MakeFLT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DA_MakeFLT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BT_Done.
function BT_Done_Callback(hObject, eventdata, handles)
% hObject    handle to BT_Done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
