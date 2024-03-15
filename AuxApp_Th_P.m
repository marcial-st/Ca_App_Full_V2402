function varargout = AuxApp_Th_P(varargin)
% AUXAPP_TH_P MATLAB code for AuxApp_Th_P.fig
%      AUXAPP_TH_P, by itself, creates a new AUXAPP_TH_P or raises the existing
%      singleton*.
%
%      H = AUXAPP_TH_P returns the handle to a new AUXAPP_TH_P or the handle to
%      the existing singleton*.
%
%      AUXAPP_TH_P('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUXAPP_TH_P.M with the given input arguments.
%
%      AUXAPP_TH_P('Property','Value',...) creates a new AUXAPP_TH_P or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AuxApp_Th_P_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AuxApp_Th_P_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AuxApp_Th_P

% Last Modified by GUIDE v2.5 09-Jan-2023 10:42:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AuxApp_Th_P_OpeningFcn, ...
                   'gui_OutputFcn',  @AuxApp_Th_P_OutputFcn, ...
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


% --- Executes just before AuxApp_Th_P is made visible.
function AuxApp_Th_P_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AuxApp_Th_P (see VARARGIN)
global ca_prof_smooth
global t
global th_prom
global prof_sel
global n_prof
global n_frames

prof_sel = 1;
[n_prof,n_frames] = size(ca_prof_smooth);
ca_prof_smooth=getappdata(0,'ca_prof_smooth');
t=getappdata(0,'t');
th_prom=getappdata(0,'th_prom');
if ~isempty(th_prom)
    set(handles.edit1,'String',num2str(th_prom))
end
% plot(handles.axes1,t,ca_prof_smooth,'Color',[0.85 0.85 0.85])
% hold(handles.axes1,"on")
% plot(handles.axes1,t,ca_prof_smooth(prof_sel,:),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
% xlabel("t[s]")
% ylabel("Ratio")
% ylim([min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
% xlim([t(1) t(end)])
% Choose default command line output for AuxApp_Th_P
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
update_plot(hObject, handles)

% UIWAIT makes AuxApp_Th_P wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AuxApp_Th_P_OutputFcn(hObject, eventdata, handles) 
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
global th_prom
th_prom = str2double(get(handles.edit1,'String'));
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
function tb_thp_inc_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_thp_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global th_prom
th_prom = th_prom+0.001;
set(handles.edit1,'String',num2str(th_prom))
update_plot(hObject, handles)

% --------------------------------------------------------------------
function tb_thp_dec_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_thp_dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global th_prom
th_prom = th_prom-0.001;
set(handles.edit1,'String',num2str(th_prom))
update_plot(hObject, handles)

% --------------------------------------------------------------------
function tb_go_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global th_prom
setappdata(0,'th_prom',th_prom)
close();

% --------------------------------------------------------------------
function tb_prf_next_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_prf_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global prof_sel
global n_prof

if prof_sel== n_prof
    prof_sel = 1;
else
    prof_sel = prof_sel+1;
end
update_plot(hObject, handles)

% --------------------------------------------------------------------
function tb_prof_prev_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_prof_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global prof_sel
global n_prof

if prof_sel==1
    prof_sel = n_prof;
else
    prof_sel = prof_sel-1;
end
update_plot(hObject, handles)

function update_plot(hObject, handles)
global ca_prof_smooth
global prof_sel
global t
global n_prof
global n_frames
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



for i_prof=1:n_prof
    [pks_all,pk_max_loc_f,pks_w_all,pks_prom_all,tri_vector] = findpeaks_wrapper(ca_prof_smooth(i_prof,:),0,i_prof);        
    global_events(i_prof,pk_max_loc_f) = pks_prom_all;
    global_events_w(i_prof,pk_max_loc_f) = pks_w_all;
end 
global_events_bin = global_events>th_prom;
all_prominences = nonzeros(global_events);

cla(handles.axes1)
plot(handles.axes1,t,ca_prof_smooth,'Color',[0.9 0.9 0.9])
hold(handles.axes1,"on")
plot(handles.axes1,t,ca_prof_smooth(prof_sel,:),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
plot(handles.axes1,t(global_events_bin(prof_sel,:)),0.02+ca_prof_smooth(prof_sel,global_events_bin(prof_sel,:)),'rv','MarkerSize',9)
xlabel(handles.axes1,"t[s]")
ylabel(handles.axes1,"Ratio")
ylim(handles.axes1,[min(min(ca_prof_smooth)) max(max(ca_prof_smooth))])
xlim(handles.axes1,[t(1) t(end)])
text(handles.axes1,0.72*t(end),0.97*max(max(ca_prof_smooth)),strcat("Profile #",num2str(prof_sel)))


nbins=20;
hbins = linspace(min(all_prominences),max(all_prominences),nbins);
hlim = floor(nbins/4);
cla(handles.axes2)
hist(handles.axes2,all_prominences,hbins)
xlabel(handles.axes2,"Left prominence")
ylabel(handles.axes2,"Freq")
% xlim(handles.axes2,[hbins(hlim) hbins(3*hlim)])
ylim(handles.axes2,[0 150])
th_prom_str = strcat("Th_prom=",num2str(th_prom));
xline(handles.axes2,th_prom,'r','Label',th_prom_str,'Interpreter','none','LabelVerticalAlignment','top')
