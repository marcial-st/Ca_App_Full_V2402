function varargout = PeakGT(varargin)
% PEAKGT MATLAB code for PeakGT.fig
%      PEAKGT, by itself, creates a new PEAKGT or raises the existing
%      singleton*.
%
%      H = PEAKGT returns the handle to a new PEAKGT or the handle to
%      the existing singleton*.
%
%      PEAKGT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEAKGT.M with the given input arguments.
%
%      PEAKGT('Property','Value',...) creates a new PEAKGT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PeakGT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PeakGT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PeakGT

% Last Modified by GUIDE v2.5 08-Apr-2021 00:23:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PeakGT_OpeningFcn, ...
                   'gui_OutputFcn',  @PeakGT_OutputFcn, ...
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


% --- Executes just before PeakGT is made visible.
function PeakGT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PeakGT (see VARARGIN)

% Choose default command line output for PeakGT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PeakGT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PeakGT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
global profiles_matrix
global profiles_smooth
global t
global nframes
global pidx
global spks

pidx = get(handles.listbox1,'Value');
cla(handles.axes1)
[pks,pks_loc,pks_w,pks_prom] = findpeaks(profiles_smooth(pidx,:),'WidthReference','halfprom');
plot(t,profiles_matrix(pidx,:),'LineWidth',2,'Color',[0.8 0.8 0.8],'Parent',handles.axes1)
hold on
plot(t,profiles_smooth(pidx,:),'LineWidth',1.5,'Color',[0 0 0.8],'Parent',handles.axes1)

if not(isempty(spks(pidx).pks))
    pks = spks(pidx).pks;
    plot(pks(:,1),pks(:,2),'r*','MarkerSize',6,'Parent',handles.axes1)
    set(handles.pushbutton2,'Enable','on')
else
    set(handles.pushbutton2,'Enable','off')
end

xlim([1 t(end)])
ylim([min(min(profiles_matrix)) max(max(profiles_matrix))])

set(handles.p_peak,'Enable','On')




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


% --------------------------------------------------------------------
function m_file_Callback(hObject, eventdata, handles)
% hObject    handle to m_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to m_loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global profiles_matrix;
global profiles_smooth;
global t;
global prof_labels;
global nframes;
global nprof;
global data;
global spks
global spksflag;

spks = struct();
Ts = 3;
wfilter = 11;
[filename, pathname] = uigetfile('*.mat','Select the 340nm first image');
if (filename ~= 0)
    data = load(fullfile(pathname,filename));
    profiles_matrix = data.profiles_matrix;
    spksflag = isfield(data,'spks');
 
    [nprof,nframes]=size(profiles_matrix);
    t = (0:1:nframes-1).*Ts;
    for i_smooth = 1:1:nprof    
        s_prof = sgolayfilt(profiles_matrix(i_smooth,:),3,wfilter);
        profiles_smooth(i_smooth,:) = s_prof;        
    end
 
    for ind=1:nprof
         prof_labels{ind} = strcat("Prof ",num2str(ind,'%03d'));
         if spksflag == 0
             spks(ind).pks = [];
         end
    end
    set(handles.listbox1,'string',prof_labels);
    if spksflag == 0
       msgbox('<-- Please select a profile');
    else
       spks = data.spks;
       msgbox('<--Please select a profile, last session detected and loaded :)');
    end
%     plot(t,profiles_matrix(1,:),'Parent',handles.axes1)
%     xlim([1 nframes])
%     ylim([min(min(profiles_matrix)) max(max(profiles_matrix))])
end
disp('Hola')


% --- Executes on button press in p_peak.
function p_peak_Callback(hObject, eventdata, handles)
% hObject    handle to p_peak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pidx
global spks

pks = spks(pidx).pks;
[x,y] = ginput(1);
x = floor(x);
pks = [pks;[x,y]];
plot(pks(:,1),pks(:,2),'r*','MarkerSize',6,'Parent',handles.axes1)
spks(pidx).pks = pks;
set(handles.pushbutton2,'Enable','on')
disp('Hola')


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pidx
global spks
global nframes
global nprof
global profiles_matrix
global profiles_smooth
global t


pks = spks(pidx).pks;

if not(isempty(pks))
    [x,y] = ginput(1);
    x = floor(x);
    del_idx = find(vecnorm(pks-[x,y],2,2)==min(vecnorm(pks-[x,y],2,2)));
    pks(del_idx,:)=[];
    spks(pidx).pks = pks;
    if sum(size(pks))~= 0
        cla(handles.axes1)
        plot(t,profiles_matrix(pidx,:),'LineWidth',2,'Color',[0.8 0.8 0.8],'Parent',handles.axes1)
        hold on
        plot(t,profiles_smooth(pidx,:),'LineWidth',1.5,'Color',[0 0 0.8],'Parent',handles.axes1)
        plot(pks(:,1),pks(:,2),'r*','MarkerSize',6,'Parent',handles.axes1)
        xlim([1 t(end)])
        ylim([min(min(profiles_matrix)) max(max(profiles_matrix))])
    end
else
    set(handles.pushbutton2,'Enable','off')
end

% --------------------------------------------------------------------
function m_save_Callback(hObject, eventdata, handles)
% hObject    handle to m_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global spks
global profiles_matrix
global profiles_smooth

[filename,path,uaccept] = uiputfile('*.mat','Save File','test.mat');
if uaccept ==1    
    save(fullfile(path,filename),'profiles_matrix','profiles_smooth','spks')
    msgbox('Session saved, thanks a lot!');
end