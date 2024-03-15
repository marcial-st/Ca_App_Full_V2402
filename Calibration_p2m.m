function varargout = Calibration_p2m(varargin)
% CALIBRATION_P2M MATLAB code for Calibration_p2m.fig
%      CALIBRATION_P2M, by itself, creates a new CALIBRATION_P2M or raises the existing
%      singleton*.
%
%      H = CALIBRATION_P2M returns the handle to a new CALIBRATION_P2M or the handle to
%      the existing singleton*.
%
%      CALIBRATION_P2M('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATION_P2M.M with the given input arguments.
%
%      CALIBRATION_P2M('Property','Value',...) creates a new CALIBRATION_P2M or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Calibration_p2m_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Calibration_p2m_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Calibration_p2m

% Last Modified by GUIDE v2.5 17-Jun-2020 09:19:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Calibration_p2m_OpeningFcn, ...
                   'gui_OutputFcn',  @Calibration_p2m_OutputFcn, ...
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


% --- Executes just before Calibration_p2m is made visible.
function Calibration_p2m_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Calibration_p2m (see VARARGIN)

% Choose default command line output for Calibration_p2m
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Calibration_p2m wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Calibration_p2m_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
global m_ref
global n_ref
global m_ser
global n_ser

set(handles.axes1,'Visible','Off')
set(handles.edit1,'Visible','Off')
set(handles.text2,'Visible','Off')
set(handles.text3,'Visible','Off')
set(handles.text4,'Visible','Off')
set(handles.uitoggletool1,'Enable','Off')
set(handles.uitoggletool2,'Enable','Off')
set(handles.uitoggletool3,'Enable','Off')
set(handles.pb_undo,'Enable','Off')
set(handles.pb_refpoints,'Enable','Off')
set(handles.pb_done,'Enable','Off')
[filename, pathname] = uigetfile({'*.tif;*.jpg;*.png','Image Files (*.tif,*.jpg,*.png)';...                                  
                                  '*.*',  'All Files (*.*)'},...
                                  'Select reference image for calibration');
if (filename ~= 0)
    img_ref = rgb2gray(imread([pathname,filename]));
    [m_ref,n_ref] = size(img_ref);    
    img_ser = getappdata(0,'im340_1');
    [m_ser,n_ser] = size(img_ser);
    set(handles.axes1,'Visible','On')
    set(handles.uitoggletool1,'Enable','On')
    set(handles.uitoggletool2,'Enable','On')
    set(handles.uitoggletool3,'Enable','On')
    set(handles.pb_undo,'Enable','Off')
    set(handles.pb_refpoints,'Enable','On')
    imshow(img_ref,'Parent',handles.axes1)    
end


% --------------------------------------------------------------------
function pb_refpoints_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_refpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global m_ref
global n_ref
global m_ser
global n_ser
global scale_ref
global dist

axes(handles.axes1);
points = ginput(2);
% scale_ref = inputScaleReference;
line(points(:,1),points(:,2),'Color','magenta','LineWidth',2.5)
set(handles.pb_done,'Enable','On')

if (m_ref == m_ser)&&(n_ref == n_ser)
    delta_x = 1;
    delta_y = 1;
else
    rescale = wantToRescale;
    if rescale
        delta_x = m_ser/m_ref;
        delta_y = n_ser/n_ref;
    else
        delta_x = 1;
        delta_y = 1;
    end    
end
dist = sqrt(delta_x*(points(1,1)-points(2,1))^2+delta_y*(points(1,2)-points(2,2))^2);
set(handles.edit1,'Visible','On')
set(handles.text2,'Visible','On')
set(handles.pb_undo,'Enable','On')
set(handles.pb_done,'Enable','On')

function rescale = wantToRescale(~)
% 'Reference and series images have different sizes. Would you like to rescale reference to series?'
answer = questdlg('Reference image size is different to series image size. Would you like to rescale?', ...
	'Size missmatch','Yes','No','Yes');
% Handle response
switch answer
    case 'Yes'
        rescale = 1;        
    case 'No'
        rescale = 0;
    otherwise
        rescale = 0;
end

function scale_ref = inputScaleReference(~)
% 'Reference and series images have different sizes. Would you like to rescale reference to series?'
prompt = {'Scale_ref = '};
dlgtitle = 'Please input scale reference in meters';
dims = [1 35];
definput = {'0'};
scale_ref = str2double(inputdlg(prompt,dlgtitle,dims,definput))';


% --------------------------------------------------------------------
function pb_undo_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout = Calibration_p2m_OutputFcn(hObject, eventdata, handles)
% Calibration_p2m_OutputFcn



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global scale_ref
global dist
global p2m

scale_ref = get(handles.edit1,'String');
[scale_ref,~] = eng2num(scale_ref)
if (dist~=0)
    dist
    p2m = scale_ref/dist;
    p2m_str = num2eng(p2m);
    p2m_str = [p2m_str,'m/pixel'];
    set(handles.text4,'String',p2m_str)
    set(handles.text3,'Visible','On')
    set(handles.text4,'Visible','On')
end
set(handles.pb_done,'Enable','On')

 

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


% --------------------------------------------------------------------
function pb_done_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global p2m
p2m
if (p2m == 0)
    warndlg("Please provide a scale reference")
elseif isempty(p2m)
    warndlg("Please provide a scale reference")
else    
    setappdata(0,'p2m',p2m);
    close();
end
