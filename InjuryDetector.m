function varargout = InjuryDetector(varargin)
% INJURYDETECTOR MATLAB code for InjuryDetector.fig
%      INJURYDETECTOR, by itself, creates a new INJURYDETECTOR or raises the existing
%      singleton*.
%
%      H = INJURYDETECTOR returns the handle to a new INJURYDETECTOR or the handle to
%      the existing singleton*.
%
%      INJURYDETECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INJURYDETECTOR.M with the given input arguments.
%
%      INJURYDETECTOR('Property','Value',...) creates a new INJURYDETECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before InjuryDetector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to InjuryDetector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InjuryDetector

% Last Modified by GUIDE v2.5 24-Jun-2023 09:15:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InjuryDetector_OpeningFcn, ...
                   'gui_OutputFcn',  @InjuryDetector_OutputFcn, ...
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


% --- Executes just before InjuryDetector is made visible.
function InjuryDetector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InjuryDetector (see VARARGIN)

% Choose default command line output for InjuryDetector
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes InjuryDetector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = InjuryDetector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
global inj_skel
global inj_bin
global inj_img
global process_list

process_list = getappdata(0,'process_list');
img = getappdata(0,'img_inj');

process_list = processListOnlyEnhancement(process_list);
img = processListApply(process_list,img);
% img = imread('img0-218.tif');

set(handles.uipanel1,'Visible','off')
set(handles.text2,'Visible','off')
set(handles.edit1,'Visible','off')
set(handles.pb_go,'Visible','off')
set(handles.slider1,'Visible','off')
axes(handles.axes1)
imshow(img)
[inj_bin,inj_img] = injurySegmentation(img);
imshow(inj_img)

set(handles.pb_go,'Visible','off')
set(handles.slider1,'Visible','off')

inj_skel = bwmorph(inj_bin,'skel',Inf);
inj_skel_props = regionprops(inj_skel,'Area','Perimeter');

set(handles.slider1, 'Min', min([inj_skel_props.Area]));
set(handles.slider1, 'Max', max([inj_skel_props.Area]));
set(handles.slider1, 'Value', mean([inj_skel_props.Area]));
set(handles.uipanel1,'Visible','on')
set(handles.edit1,'Visible','on')
set(handles.slider1,'Visible','on')


% --- Executes on button press in pb_go.
function pb_go_Callback(hObject, eventdata, handles)
% hObject    handle to pb_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global inj_band
global inj_bin

if isempty(inj_band)
    warndlg("Please define the injury band.")
elseif isempty(inj_bin) 
    warndlg("Please define the injury band.")
else
    setappdata(0,'inj_band',inj_band)
    setappdata(0,'inj_bin',inj_bin)
    close(gcf)
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global area_thr
global inj_override

area_thr = get(handles.slider1,'Value')
set(handles.edit1,'String',num2str(area_thr))
inj_override=0;
plot_injury(hObject, eventdata, handles)
set(handles.pb_go,'Visible','on')

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined nhhin a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global area_thr
global inj_override
area_thr = str2double(get(handles.edit1,'String'))
set(handles.slider1,'Value',area_thr)
inj_override=0;
plot_injury(hObject, eventdata, handles)
set(handles.pb_go,'Visible','on')

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

function plot_injury(hObject, eventdata, handles)
global inj_skel
global inj_bin
global inj_img
global inj_band
global area_thr
global inj_filtered
global inj_override

if ~inj_override
    inj_band = find(sum(inj_filtered) >= 1);
end
inj_filtered = bwareafilt(inj_skel,[area_thr 1000]);
[h_i,~] = size(inj_filtered);
inj_skel = bwmorph(inj_bin,'skel',Inf);
inj_skel_props = regionprops(inj_skel,'Area','Perimeter');
img_class_inj_filt = inj_img;
img_class_inj_filt(:,:,2) = img_class_inj_filt(:,:,2)+uint8(inj_filtered)*255;
axes(handles.axes1)
imshow(img_class_inj_filt)
rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)


% --- Executes on button press in pb_band_start.
function pb_band_start_Callback(hObject, eventdata, handles)
% hObject    handle to pb_band_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global inj_band

[x_1,~]=ginput(1);
inj_band(1)=floor(x_1);


% --- Executes on button press in pb_band_end.
function pb_band_end_Callback(hObject, eventdata, handles)
% hObject    handle to pb_band_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global inj_band
global inj_override

[x_end,~]=ginput(1);
if length(inj_band) == 1
    inj_band(2)=floor(x_end);
else
    inj_band(length(inj_band))=floor(x_end);
end
inj_override = 1;
plot_injury(hObject, eventdata, handles)
set(handles.pb_go,'Visible','on')