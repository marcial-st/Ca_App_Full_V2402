function varargout = PeakGT(varargin)
% PEAKGT MATLAB code for PeakGT.fig
%      PEAKGT, by itself, creates a new PEAKGT or raises the existinggrid o
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

% Last Modified by GUIDE v2.5 18-Aug-2021 16:36:20

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
global data
global pidx
global spks
global pk_gt_en
global pk_all_en
global Ts
global th_p
global th_w
global t_gray

pidx = get(handles.listbox1,'Value');
cla(handles.axes1)
[pks_all,pks_loc_all,pks_w_all,pks_prom_all] = findpeaks(profiles_smooth(pidx,:),'WidthReference','halfprom');
plot(t,profiles_matrix(pidx,:),'LineWidth',2,'Color',[0.8 0.8 0.8],'Parent',handles.axes1)
hold on
plot(t,profiles_smooth(pidx,:),'LineWidth',1.5,'Color',[0 0 0.8],'Parent',handles.axes1)

if pk_all_en
    if not(isempty(pks_all))
        
        fwp = (pks_w_all*Ts)>=th_w;
        fvp = pks_prom_all>=th_p;       
        
        fp = fwp.*fvp;
        fn = not(fp);
        
        plot(nonzeros(pks_loc_all*Ts.*fp),nonzeros(pks_all.*fp),'kv','MarkerSize',8,'MarkerFaceColor','green','Parent',handles.axes1)
        if t_gray
            plot(nonzeros(pks_loc_all*Ts.*fn),nonzeros(pks_all.*fn),'v','MarkerSize',7,'MarkerFaceColor',[0.85 0.8 0.8],'Parent',handles.axes1)
        end
        set(handles.pushbutton2,'Enable','on')
    else
        set(handles.pushbutton2,'Enable','off')
        disp("pks_all empty")
    end
end

if pk_gt_en     
    if not(isempty(spks(pidx).pks))
        pks_gt = spks(pidx).pks;
        [prom_gt,~,~,wact_gt,~] = getFeatures(pks_gt,profiles_smooth(pidx,:));
        
        fwp = (wact_gt*Ts)>=th_w;
        fvp = prom_gt>=th_p;       
        
        fp = fwp.*fvp;
        fn = not(fp);
        
        plot(nonzeros(pks_gt(:,1).*fp),nonzeros(pks_gt(:,2).*fp),'rp','MarkerSize',10,'MarkerFaceColor','magenta','Parent',handles.axes1)
        if t_gray
            plot(nonzeros(pks_gt(:,1).*fn),nonzeros(pks_gt(:,2).*fn),'p','MarkerSize',9,'MarkerFaceColor',[0.85 0.8 0.8],'Parent',handles.axes1)
        end
        set(handles.pushbutton2,'Enable','on')
    else
        set(handles.pushbutton2,'Enable','off')
    end
end

xlim([1 t(end)])
% ylim([min(min(profiles_matrix)) max(max(profiles_matrix))])
ylim([0.75 1.55])
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
global Ts;
global th_p
global th_w

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
    
     if isfield(data,'params')
        if isfield(data.params,'th_p')
            th_p = data.params.th_p;            
        else
            th_p = 0;
        end
        set(handles.edit_p,'String',num2str(data.params.th_p))
            
        if isfield(data.params,'th_w')
            th_w = data.params.th_w;
        else
            th_w = 0;
        end
        set(handles.edit_w,'String',num2str(data.params.th_w))
    end
    
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
listbox1_Callback(hObject, eventdata, handles)
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
%         ylim([min(min(profiles_matrix)) max(max(profiles_matrix))])
        ylim([0.75 1.55])
    end
else
    set(handles.pushbutton2,'Enable','off')
end
listbox1_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function m_save_Callback(hObject, eventdata, handles)
% hObject    handle to m_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global spks
global profiles_matrix
global profiles_smooth
global data
global th_p
global th_w

data.params.th_p = th_p;
data.params.th_w = th_w;
data.profiles_matrix = profiles_matrix;
data.profiles_smooth = profiles_smooth;
data.spks = spks;

[filename,path,uaccept] = uiputfile('*.mat','Save File','test.mat');
if uaccept ==1    
%     save(fullfile(path,filename),'profiles_matrix','profiles_smooth','spks')
    save(fullfile(path,filename),'-struct','data')
%     save(fullfile(path,filename),'data')
    msgbox('Session saved, thanks a lot!');
end


% --------------------------------------------------------------------
function m_tools_Callback(hObject, eventdata, handles)
% hObject    handle to m_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text4,'Enable','on');
set(handles.text5,'Enable','on');
set(handles.edit_p,'Enable','on');
set(handles.edit_w,'Enable','on');
set(handles.t_gtpeaks,'Enable','on');
set(handles.t_allpeaks,'Enable','on');

% --------------------------------------------------------------------
function m_peakstats_Callback(hObject, eventdata, handles)
% hObject    handle to m_peakstats (see GCBO)
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
global Ts

[n_cells n_frames] = size(profiles_smooth);

prom_gt_global = [];
wact_gt_global = [];

prom_all_global = [];
wact_all_global = [];

for i=1:1:n_cells
    
    %% Ground Truth
    [prom_gt,~,~,wact_gt,~] = getFeatures(spks(i).pks,profiles_smooth(i,:));
    prom_gt_global = [prom_gt_global;prom_gt];
    wact_gt_global = [wact_gt_global;wact_gt];
    %% Plain Data
    [pks_all,pks_loc_all,pks_w_all,pks_prom_all] = findpeaks(profiles_smooth(i,:),'WidthReference','halfprom');
    [prom_all,~,~,wact_all,~] = getFeatures([pks_loc_all*Ts;pks_all]',profiles_smooth(i,:));
    prom_all_global = [prom_all_global;prom_all];
    wact_all_global = [wact_all_global;wact_all];
%     prom_all_global = [prom_all_global;pks_prom_all];
%     wact_all_global = [wact_all_global;pks_w_all];
end

hedges = 0:0.02:0.4;
ghist = figure;
hpl = histogram(prom_gt_global,'BinEdges',hedges,'FaceColor','r')
title('Ground Truth Global Prominence')
xlabel('Prom Value')
ylabel('Freq')
x = hpl.BinEdges ;
Prom_L = hpl.Values ;
text(x(1:end-1),Prom_L,num2str(Prom_L'),'vert','bottom','horiz','center');

ahist = figure;
hpl = histogram(prom_all_global,'BinEdges',hedges,'FaceColor',[0.75 0.75 0.75])
title('Automatic Global Prominence')
xlabel('Prom Value')
ylabel('Freq')
x = hpl.BinEdges ;
Prom_L = hpl.Values ;
text(x(1:end-1),Prom_L,num2str(Prom_L'),'vert','bottom','horiz','center');

fwidth=figure;
wedges = [0:5:40].*Ts;
hwact = histogram(wact_gt_global.*Ts,'BinEdges',wedges,'FaceColor','m')
title('Ground Truth Global Activity Width')
xlabel('Time')
ylabel('Freq')
x = hwact.BinEdges ;
y = hwact.Values ;
text(x(1:end-1),y,num2str(y'),'vert','bottom','horiz','center');

fwidtha=figure;
hwact = histogram(wact_all_global.*Ts,'BinEdges',wedges,'FaceColor','g')
title('Automatic Global Activity Width')
xlabel('Time')
ylabel('Freq')
x = hwact.BinEdges ;
y = hwact.Values ;
text(x(1:end-1),y,num2str(y'),'vert','bottom','horiz','center');

function edit_w_Callback(hObject, eventdata, handles)
% hObject    handle to edit_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_w as text
%        str2double(get(hObject,'String')) returns contents of edit_w as a double
global th_w
th_w = str2num(get(handles.edit_w,'String'))
listbox1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_w_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p as text
%        str2double(get(hObject,'String')) returns contents of edit_p as a double
global th_p
th_p = str2num(get(handles.edit_p,'String'))
listbox1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function t_gtpeaks_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to t_gtpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pk_gt_en
pk_gt_en = get(handles.t_gtpeaks,'State');
listbox1_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function t_allpeaks_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to t_allpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pk_all_en
pk_all_en = get(handles.t_allpeaks,'State');
listbox1_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function t_gray_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to t_gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t_gray
t_gray = get(handles.t_gray,'State');
listbox1_Callback(hObject, eventdata, handles)
