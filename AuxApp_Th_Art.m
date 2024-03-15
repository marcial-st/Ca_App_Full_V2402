function varargout = AuxApp_Th_Art(varargin)
% AUXAPP_TH_ART MATLAB code for AuxApp_Th_Art.fig
%      AUXAPP_TH_ART, by itself, creates a new AUXAPP_TH_ART or raises the existing
%      singleton*.
%
%      H = AUXAPP_TH_ART returns the handle to a new AUXAPP_TH_ART or the handle to
%      the existing singleton*.
%
%      AUXAPP_TH_ART('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUXAPP_TH_ART.M with the given input arguments.
%
%      AUXAPP_TH_ART('Property','Value',...) creates a new AUXAPP_TH_ART or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AuxApp_Th_Art_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AuxApp_Th_Art_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AuxApp_Th_Art

% Last Modified by GUIDE v2.5 10-Jan-2023 09:41:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AuxApp_Th_Art_OpeningFcn, ...
                   'gui_OutputFcn',  @AuxApp_Th_Art_OutputFcn, ...
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


% --- Executes just before AuxApp_Th_Art is made visible.
function AuxApp_Th_Art_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AuxApp_Th_Art (see VARARGIN)
global ca_prof_smooth
global t
global th_artifact
global prof_sel
global n_prof
global n_frames
global th_prom
global th_diff

prof_sel = 1;
[n_prof,n_frames] = size(ca_prof_smooth);
ca_prof_smooth=getappdata(0,'ca_prof_smooth');
t=getappdata(0,'t');
th_prom=getappdata(0,'th_prom');
th_artifact = getappdata(0,'th_artifact');
th_diff = getappdata(0,'th_diff');
if ~isempty(th_prom)
    set(handles.edit1,'String',num2str(th_artifact))
end
% Choose default command line output for AuxApp_Th_Art
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
update_plot(hObject, handles)
% UIWAIT makes AuxApp_Th_Art wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AuxApp_Th_Art_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global th_artifact
th_artifact = str2double(get(handles.edit1,'String'));
update_plot(hObject, handles)


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
function tb_th_inc_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_th_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global th_artifact
th_artifact = th_artifact+0.5;
set(handles.edit1,'String',num2str(th_artifact))
update_plot(hObject, handles)

% --------------------------------------------------------------------
function tb_th_dec_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_th_dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global th_artifact
th_artifact = th_artifact-0.5;
set(handles.edit1,'String',num2str(th_artifact))
update_plot(hObject, handles)


% --------------------------------------------------------------------
function tb_go_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global th_artifact
setappdata(0,'th_artifact',th_artifact)
close();

% --------------------------------------------------------------------
function update_plot(hObject, handles)
global ca_prof_smooth
global prof_sel
global t
global n_prof
global n_frames
global th_artifact
global th_diff
global th_prom

disp("Funciona")
global_events = zeros(n_prof,n_frames);
global_events_w = zeros(n_prof,n_frames);

[std_prof, std_idx]=sort(std(ca_prof_smooth,0,2),'ascend');
n_comp = 20;
l_fit = zeros(n_comp,2);
for i_comp = 1:n_comp            
    l_fit(i_comp,:) = polyfit(t,ca_prof_smooth(std_idx(i_comp),:),1);
    yfit = l_fit(i_comp,1)*t+l_fit(i_comp,2);
end
l_fit = mean(l_fit);
y_comp = (l_fit(1)*t);
profiles_smooth_comp = ca_prof_smooth-(ones(n_prof,1)*y_comp);
profiles_smooth_bu = ca_prof_smooth;
ca_prof_smooth = profiles_smooth_comp;

profiles_diff = diff(ca_prof_smooth,1,2)>th_diff;
[~,inj_frame] = max(sum(profiles_diff,1));

for i_prof=1:n_prof
    [pks_all,pk_max_loc_f,pks_w_all,pks_prom_all,tri_vector] = findpeaks_wrapper(ca_prof_smooth(i_prof,:),0,i_prof);        
    global_events(i_prof,pk_max_loc_f) = pks_prom_all;
    global_events_w(i_prof,pk_max_loc_f) = pks_w_all;
end 
global_events_bin = global_events>th_prom;
global_events_vector = sum(global_events_bin,1);
%TODO: zero the region global_events_vector(inj_frame:inj_frame+peak_width_f) = 0;
global_events_vector(inj_frame:inj_frame+30) = 0;
plot(handles.axes1,t,global_events_vector)
th_artifact_str = strcat("Th_art=",num2str(th_artifact));
yline(handles.axes1,th_artifact,'r','Label',th_artifact_str,'Interpreter','none','LabelHorizontalAlignment','left')
xlim([t(1) t(end)])
if th_artifact>max(global_events_vector)
    ylim([0 th_artifact+2])
else
    ylim([0 max(global_events_vector)+2])
end
xlabel(handles.axes1,"t[s]")
ylabel(handles.axes1,"F2F Activity Freq")
