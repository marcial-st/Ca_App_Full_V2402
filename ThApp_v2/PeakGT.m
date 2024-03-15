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

% Last Modified by GUIDE v2.5 20-Sep-2021 08:32:03

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
global load_control
set(handles.popupmenu1,'Enable','Off')
load_control = 0;


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
global th_k
global t_gray
global TP
global TN
global FP
global FN
global total
global params
global fullrates
global epsilon
global tstart
global tstim
global tstop
global tstop_f
global approach_sel

global sensitivity
global specificity
global accuracy 

pidx = get(handles.listbox1,'Value');
tstop_f = floor(tstop/Ts);
cla(handles.axes1)

[TP,TN,FP,FN,total,fullrates] = peakassessment(spks,profiles_smooth,profiles_matrix,th_p,th_w,th_k,epsilon,params,approach_sel);
% disp("TN = ")
% TN
% disp("FP = ")
% FP
[pks_all,pks_loc_all,pks_w_all,pks_prom_all] = findpeaks_wrapper(profiles_smooth(pidx,tstart:tstop_f));

plot(handles.axes1,t(tstart:tstop_f),profiles_matrix(pidx,tstart:tstop_f),'LineWidth',2,'Color',[0.8 0.8 0.8],'Parent',handles.axes1)
hold on
plot(handles.axes1,t(tstart:tstop_f),profiles_smooth(pidx,tstart:tstop_f),'LineWidth',1.5,'Color',[0 0 0.8],'Parent',handles.axes1)
tstim_f = floor(tstim/Ts);
xline(handles.axes1,t(tstim_f),'--',{strcat('t_{stim}=',num2str(tstim))},'LabelVerticalAlignment','bottom')

if pk_all_en
    if not(isempty(pks_all))
        
%         fwp = (pks_w_all*Ts)>=th_w;
%         fvp = pks_prom_all>=th_p;       
%         
%         fp_a = and(fwp,fvp)';
%         fn = not(fp_a);        
%         plot(pks_loc_all.*fp_a*Ts,pks_all.*fp_a,'kv','MarkerSize',8,'MarkerFaceColor','green','Parent',handles.axes1)
        filt_a_p = [fullrates.filt_a_p{pidx}];
        filt_a_n = [fullrates.filt_a_n{pidx}];
        pks_all = [fullrates.pks_all{pidx}];
        pks_loc_all = [fullrates.pks_loc_all{pidx}];
        plot(handles.axes1,(pks_loc_all(filt_a_p)*Ts)-Ts,pks_all(filt_a_p),'kv','MarkerSize',8,'MarkerFaceColor','green','Parent',handles.axes1)
        if t_gray
            plot(handles.axes1,(pks_loc_all(filt_a_n)*Ts)-Ts,pks_all(filt_a_n),'v','MarkerSize',7,'MarkerFaceColor',[0.85 0.8 0.8],'Parent',handles.axes1)
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
        [txt_pk,~]= size(pks_gt);
        [prom_gt,~,~,wact_gt,~] = getFeatures(pks_gt,profiles_smooth(pidx,:));
        
        fwp = (wact_gt*Ts)>=th_w;
        fvp = prom_gt>=th_p;       
        
        fp = fwp.*fvp;
        fn = not(fp);
        plot(handles.axes1,pks_gt(:,1),pks_gt(:,2),'rp','MarkerSize',10,'MarkerFaceColor','magenta','Parent',handles.axes1)        
%         plot(nonzeros(pks_gt(:,1).*fp),nonzeros(pks_gt(:,2).*fp),'rp','MarkerSize',10,'MarkerFaceColor','magenta','Parent',handles.axes1)
%         if t_gray
%             plot(nonzeros(pks_gt(:,1).*fn),nonzeros(pks_gt(:,2).*fn),'p','MarkerSize',9,'MarkerFaceColor',[0.85 0.8 0.8],'Parent',handles.axes1)
%         end
        set(handles.pushbutton2,'Enable','on')
    else
        set(handles.pushbutton2,'Enable','off')
    end
end

    yrange = max(max(profiles_smooth))-min(min(profiles_smooth));
 
%            strcat("  TP = ",num2str(TP/sum([fullrates.count_pos])*100),"%"),...
%            strcat("  FP = ",num2str(FP/sum([fullrates.count_pos])*100),"%"),...
%            strcat("  TN = ",num2str(TN/sum([fullrates.count_neg])*100),"%"),...
%            strcat("  FN = ",num2str(FN/sum([fullrates.count_neg])*100),"%"),...    
   
    if sum([fullrates.count_pos])~=0
        tpr = TP/sum([fullrates.count_pos])*100;
        fpr = FP/sum([fullrates.count_pos])*100;
    else
        tpr = 0;
        fpr = 0;        
    end
           
   if sum([fullrates.count_neg])~=0
        tnr = TN/sum([fullrates.count_neg])*100;
        fnr = FN/sum([fullrates.count_neg])*100;
   else
       tnr = 0;
       fnr = 0;
   end
   
   acc = (TP+TN)/(sum([fullrates.count_pos])+sum([fullrates.count_neg])); 
%             " ",...
%            strcat("  TP = ",num2str(tpr),"%"),...
%            strcat("  FP = ",num2str(fpr),"%"),...
%            strcat("  TN = ",num2str(tnr),"%"),...
%            strcat("  FN = ",num2str(fnr),"%"),...    
    txt = {"Global rates:",...
           strcat("  Sensitivity = ",num2str(tpr/(tpr+fnr))),...
           strcat("  Specificity = ",num2str(tnr/(tnr+fpr))),...           
           strcat("  Accuracy = ",num2str(acc)),...           
           " ",...
           "Local rates:",...fullrates
           strcat("  TP = ",num2str(fullrates.TP(pidx)),"     ","FP = ",num2str(fullrates.FP(pidx))),...           
           strcat("  TN = ",num2str(fullrates.TN(pidx)),"     ","FN = ",num2str(fullrates.FN(pidx))),...          
           strcat("  total = ",num2str(fullrates.total(pidx))) 
        };
   
    sensitivity = tpr/(tpr+fnr);
    specificity = tnr/(tnr+fpr);
    accuracy = acc;
  
    text(handles.axes1,t(tstop_f)-0.25*t(tstop_f),min(min(profiles_smooth))+0.8*yrange,txt)

xlim(handles.axes1,[1 t(tstop_f)])
ylim(handles.axes1,[min(min(profiles_smooth))-0.01 max(max(profiles_smooth))+0.01])
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
global th_k 
global TP
global TN
global FP
global FN
global total
global params
global fullrates
global epsilon
global tstart
global tstim
global tstop
global approach_sel
global filename
global pathname

global load_control

if load_control 
    profiles_matrix = []; profiles_smooth = []; t = []; prof_labels = []; nframes = []; 
    nprof = []; data = []; spks = []; spksflag = []; th_p = []; th_w = []; th_k  = [];
    TP = []; TN = []; FP = []; FN = []; total = []; params = []; fullrates = [];
    tstart = []; tstim = []; tstop = []; approach_sel = [];
end

if ~load_control 
    [filename, pathname] = uigetfile('*.mat','Select the 340nm first image');    
    epsilon = 3;
end

Ts = 3;
spks = struct();
wfilter = 11;

% if (filename ~= 0)
if not(isempty(filename))
    
    data = load(fullfile(pathname,filename));
    profiles_matrix = data.profiles_matrix;
    spksflag = isfield(data,'spks');
    
    params = data.params;
    
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
        
%          params = data.params;
        
        if isfield(data.params,'th_p')
            th_p = data.params.th_p;            
        else
            th_p = 0;
        end
        set(handles.edit_p,'String',num2str(th_p))
            
        if isfield(data.params,'th_w')
            th_w = data.params.th_w;
        else
            th_w = 0;
        end
        set(handles.edit_w,'String',num2str(th_w))
        
        if isfield(data.params,'th_k')
            th_k = data.params.th_k;            
        else
            th_k = 0;
        end
        set(handles.e_k,'String',num2str(th_k))
        
        if isfield(data.params,'tstart')
            if data.params.tstart == 0
                tstart = 1;
            else
                tstart = data.params.tstart;            
            end
        else
            tstart = 1;
        end                        
        
        if isfield(data.params,'tstim')
            tstim = data.params.tstim;            
        else
            tstim = 20;
        end                        
      
        if isfield(data.params,'tstop')
            tstop = data.params.tstop;            
        else
            tstop = nframes*Ts;
        end                        
        
     else
         th_p = 0;
         th_w = 0;
         th_k = 0;
         epsilon = 5;
         params = 0;
         tstart = 1;
         tstim = 20;
         tstop = nframes*Ts;
     end
     
     if isempty(approach_sel)
         approach_sel = 1;
     end

    set(handles.edit_epsilon,'String',num2str(epsilon))
     
    if spksflag == 0
       msgbox('<-- Please select a profile');
    else
       spks = data.spks;       
       if load_control
           disp('Last session detected and loaded :)');
       else
           [TP,TN,FP,FN,total,fullrates] = peakassessment(spks,profiles_smooth,profiles_matrix,th_p,th_w,th_k,epsilon,params,approach_sel);
           msgbox('<--Please select a profile, last session detected and loaded :)');
       end       
    end
    set(handles.t_path,'String',fullfile(pathname,filename))
    set(handles.popupmenu1,'Enable','On')
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
    [pks_all,pks_loc_all,pks_w_all,pks_prom_all] = findpeaks_wrapper(profiles_smooth(i,:));
%     [pks_all,pks_loc_all,pks_w_all,pks_prom_all] = findpeaks(profiles_smooth(i,:),'WidthReference','halfprom');
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



function e_k_Callback(hObject, eventdata, handles)
% hObject    handle to e_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_k as text
%        str2double(get(hObject,'String')) returns contents of e_k as a double
global th_k
th_k = str2num(get(handles.e_k,'String'))
listbox1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function e_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_epsilon_Callback(hObject, eventdata, handles)
% hObject    handle to edit_epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_epsilon as text
%        str2double(get(hObject,'String')) returns contents of edit_epsilon as a double
global epsilon
epsilon = str2num(get(handles.edit_epsilon,'String'))
listbox1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_epsilon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function m_settime_Callback(hObject, eventdata, handles)
% hObject    handle to m_settime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global params
global tstart
global tstim
global tstop
global nframes
global Ts

prompt = {'t_start:','t_stim:','t_stop'};
dlgtitle = 'Please input times';
dims = [1 35];
definput = {'','',''};
tans = inputdlg(prompt,dlgtitle,dims,definput)
if str2double(tans{1})==0
    tstart = 1;
else
    tstart = str2double(tans{1});
end
tstim = str2double(tans{2});

if str2double(tans{3})>(nframes*Ts)
    tstop = nframes*Ts;
    disp(strcat('tstop exceeded limits, max allowed is: ',tstop))
else
    tstop = str2double(tans{3});
end

params.tstart = tstart;
params.tstim = tstim;
params.tstop = tstop;
listbox1_Callback(hObject, eventdata, handles)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global approach_sel

approach_sel = get(handles.popupmenu1,'Value');
listbox1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function tb_mytask_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_mytask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
global th_k
global t_gray
global TP
global TN
global FP
global FN
global total
global params
global fullrates
global epsilon
global tstart
global tstim
global filename
global approach_sel
global sensitivity
global specificity
global accuracy 
global load_control
global pathname
global filename
global epsilonv

spath = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\ReportGen";
pathname = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Report_reduced_classification\Unified\220502"; %090721
% pathname = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Report_reduced_classification\Unified\111721";
% pathname = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\Report_reduced_classification\Unified\111721_full_basal";
filenames = ["100528_1.mat","100528_4.mat","100528_5.mat","100528_6.mat","100528_8.mat","100601_1.mat","100604_9.mat"];
exp_aux = 1:1:length(filenames);

a_id = [1 2];
a_name = {'Default','Fine_v1'};
a_dict = containers.Map(a_id,a_name);

disp("Executing myTask section ...")

th_p_vector = 0.020:0.001:0.026;
% th_p_vector = 0.024;
np_th = length(th_p_vector);

th_w_vector = 0;
nw_th = length(th_w_vector);

% th_k_vector = [4 6 8 10 12 14];
% th_k_vector = [1 3 4 7 9 10 13];
th_k_vector = 0;
nk_th = length(th_k_vector);

epsilonv = [10];

mat_sen_def = zeros(np_th,length(epsilonv),length(filenames));
mat_spe_def = zeros(np_th,length(epsilonv),length(filenames));
mat_acc_def = zeros(np_th,length(epsilonv),length(filenames));
mat_sen_fine = zeros(np_th,length(epsilonv),length(filenames));
mat_spe_fine = zeros(np_th,length(epsilonv),length(filenames));
mat_acc_fine = zeros(np_th,length(epsilonv),length(filenames));


load_control = 1;

% for kf = 1:1:length(filenames)
%     %% Parameter scanning
%     filename = filenames(kf);
%     m_loadfile_Callback(hObject, eventdata, handles)
%     for ee = 1:1:length(epsilonv)        
%         epsilon = epsilonv(ee);
%         set(handles.edit_epsilon,'String',num2str(epsilon))
%         for ii = 1:1:2
%             approach_sel = ii;    
%             set(handles.popupmenu1,'Value',approach_sel); 
%             for j=1:np_th    
%                 th_p = th_p_vector(j);
%                 disp(strcat("   Approach= ",a_dict(approach_sel)))
%                 disp(strcat("   test th= ",num2str(th_p)))
%                 set(handles.edit_p,'String',num2str(th_p))
%                 pause(0.1)
%                 listbox1_Callback(hObject, eventdata, handles)
%                 pause(0.1)
%                 if approach_sel == 1
%                     mat_sen_def(j,ee,kf) = sensitivity;
%                     mat_spe_def(j,ee,kf) = specificity;
%                     mat_acc_def(j,ee,kf) = accuracy;
%                 else
%                     mat_sen_fine(j,ee,kf) = sensitivity;
%                     mat_spe_fine(j,ee,kf) = specificity;
%                     mat_acc_fine(j,ee,kf) = accuracy;
%                 end               
%             end
%         end
%         
%     end
% par_var_approach = [1,2]; %% 1 for th_p,th_w,epsilon      2 for th_k,th_w,epsilon
par_var_approach = 1;
for kf = 1:1:length(filenames)
% for kf = 6
    filename = filenames(kf)
    m_loadfile_Callback(hObject, eventdata, handles)
	for pva=1:1:length(par_var_approach)
		for ee = 1:1:length(epsilonv)        
			epsilon = epsilonv(ee);
			set(handles.edit_epsilon,'String',num2str(epsilon))
			for ii = 1:1:2
				approach_sel = ii;    
				set(handles.popupmenu1,'Value',approach_sel);
				if pva==1
					th_k = 0;
					set(handles.edit_p,'String',num2str(th_k))
					for j=1:np_th    
						th_p = th_p_vector(j);
						set(handles.edit_p,'String',num2str(th_p))
						for jw = 1:1:nw_th
							th_w = th_w_vector(jw);
							set(handles.edit_w,'String',num2str(th_w))
							pause(0.1)
							disp(strcat("   Approach= ",a_dict(approach_sel)))
							disp(strcat("   test th= ",num2str(th_p)))
							listbox1_Callback(hObject, eventdata, handles) %% Asssessment here (:
							pause(0.1)
							if approach_sel == 1
								mat_pw_sen_def(j,jw,ee,kf) = sensitivity;
								mat_pw_spe_def(j,jw,ee,kf) = specificity;
								mat_pw_acc_def(j,jw,ee,kf) = accuracy;                                
                                if ((j==1)&&(jw==1)&&(ee==1)) 
                                    exp_data(kf).gt_match_global = [fullrates.gt_match_rate];
                                end
							else
								mat_pw_sen_fine(j,jw,ee,kf) = sensitivity;
								mat_pw_spe_fine(j,jw,ee,kf) = specificity;
								mat_pw_acc_fine(j,jw,ee,kf) = accuracy;
							end               
						end
					end						
				elseif pva == 2
					th_p = 0;
					set(handles.edit_p,'String',num2str(th_p))
					for j=1:nk_th    
						th_k = th_k_vector(j);
						set(handles.edit_p,'String',num2str(th_k))
						for jw = 1:1:nw_th
							th_w = th_w_vector(jw);
							set(handles.edit_w,'String',num2str(th_w))
							pause(0.1)
							disp(strcat("   Approach= ",a_dict(approach_sel)))
							disp(strcat("   test th= ",num2str(th_p)))
							listbox1_Callback(hObject, eventdata, handles) %% Asssessment here (:
							pause(0.1)
							if approach_sel == 1
								mat_kw_sen_def(j,jw,ee,kf) = sensitivity;
								mat_kw_spe_def(j,jw,ee,kf) = specificity;
								mat_kw_acc_def(j,jw,ee,kf) = accuracy;
							else
								mat_kw_sen_fine(j,jw,ee,kf) = sensitivity;
								mat_kw_spe_fine(j,jw,ee,kf) = specificity;
								mat_kw_acc_fine(j,jw,ee,kf) = accuracy;
							end  
						end						            
					end
				end
			end
			
		end
	end    
% %     %% Data Visualization
% %     fig1 = figure;    
% %     cvect=['r','g','b'];
% %     for fidx = 1:1:length(epsilonv)
% %         subplot(3,1,1)
% %         hold on;
% %         plot(th_p_vector,mat_sen_def(:,fidx,kf),'-','Color',cvect(fidx));
% %         plot(th_p_vector,mat_sen_fine(:,fidx,kf),'--','Color',cvect(fidx));
% %         xlabel('th_p','Interpreter','none')
% %         ylabel('Sensitivity','Interpreter','none')
% %         subplot(3,1,2)
% %         hold on;
% %         plot(th_p_vector,mat_spe_def(:,fidx,kf),'-','Color',cvect(fidx));
% %         plot(th_p_vector,mat_spe_fine(:,fidx,kf),'--','Color',cvect(fidx));
% %         xlabel('th_p','Interpreter','none')
% %         ylabel('Specificity','Interpreter','none')
% %         subplot(3,1,3)
% %         hold on;
% %         plot(th_p_vector,mat_acc_def(:,fidx,kf),'-','Color',cvect(fidx));
% %         plot(th_p_vector,mat_acc_fine(:,fidx,kf),'--','Color',cvect(fidx));
% %         xlabel('th_p','Interpreter','none')
% %         ylabel('Accuracy','Interpreter','none')
% %         leg_array{fidx*2-1}=strcat("Default e=",num2str(epsilonv(fidx)));
% %         leg_array{fidx*2}=strcat("Fine e=",num2str(epsilonv(fidx)));
% %     end  
% %     legend(leg_array,'Location','best');
% %     sgtitle(strcat("Experiment ",filenames(kf)),'Interpreter','none')    
% %     saveas(gcf,fullfile(spath,strcat("im_b",num2str(kf),".png")))
end

 fig2 = figure;
[X,Y] = meshgrid(exp_aux,th_p_vector);
hold on
for fidx = 1:1:length(epsilonv)
    surf(X,Y,reshape(mat_sen_def(:,fidx,:),[51,7]),'FaceAlpha',0.45,'FaceLighting','gouraud','FaceColor',[1,0,0],'LineStyle',':')        
    surf(X,Y,reshape(mat_sen_fine(:,fidx,:),[51,7]),'FaceAlpha',0.45,'FaceLighting','gouraud','FaceColor',[0,0,1],'LineStyle',':')
    title("Sensitivity")
    xlabel("Exp")
    ylabel("th_p",'Interpreter','none')
    zlabel("Sensitivity")
end
hold off
view(-145,8)
grid on
legend('Default','Fine','Location','best')
saveas(gcf,fullfile(spath,strcat("im_a",num2str(01),".png")))

fig3 = figure;
[X,Y] = meshgrid(exp_aux,th_p_vector);
hold on
for fidx = 1:1:length(epsilonv)
    surf(X,Y,reshape(mat_spe_def(:,fidx,:),[51,7]),'FaceAlpha',0.45,'FaceLighting','gouraud','FaceColor',[1,0,0],'LineStyle',':')        
    surf(X,Y,reshape(mat_spe_fine(:,fidx,:),[51,7]),'FaceAlpha',0.45,'FaceLighting','gouraud','FaceColor',[0,0,1],'LineStyle',':')
    title("Specificity")
    xlabel("Exp")
    ylabel("th_p",'Interpreter','none')
    zlabel("Specificity")
    table1 = table(th_p_vector',mat_spe_def(:,fidx,1),mat_spe_def(:,fidx,2),mat_spe_def(:,fidx,3),mat_spe_def(:,fidx,4),mat_spe_def(:,fidx,5),mat_spe_def(:,fidx,6),mat_spe_def(:,fidx,7),mat_spe_fine(:,fidx,1),mat_spe_fine(:,fidx,2),mat_spe_fine(:,fidx,3),mat_spe_fine(:,fidx,4),mat_spe_fine(:,fidx,5),mat_spe_fine(:,fidx,6),mat_spe_fine(:,fidx,7),mat_sen_def(:,fidx,1),mat_sen_def(:,fidx,2),mat_sen_def(:,fidx,3),mat_sen_def(:,fidx,4),mat_sen_def(:,fidx,5),mat_sen_def(:,fidx,6),mat_sen_def(:,fidx,7),mat_sen_fine(:,fidx,1),mat_sen_fine(:,fidx,2),mat_sen_fine(:,fidx,3),mat_sen_fine(:,fidx,4),mat_sen_fine(:,fidx,5),mat_sen_fine(:,fidx,6),mat_sen_fine(:,fidx,7),mat_acc_def(:,fidx,1),mat_acc_def(:,fidx,2),mat_acc_def(:,fidx,3),mat_acc_def(:,fidx,4),mat_acc_def(:,fidx,5),mat_acc_def(:,fidx,6),mat_acc_def(:,fidx,7),mat_acc_fine(:,fidx,1),mat_acc_fine(:,fidx,2),mat_acc_fine(:,fidx,3),mat_acc_fine(:,fidx,4),mat_acc_fine(:,fidx,5),mat_acc_fine(:,fidx,6),mat_acc_fine(:,fidx,7),'VariableNames',{'th_p','def_spec1','def_spec2','def_spec3','def_spec4','def_spec5','def_spec6','def_spec7','fine_spec1','fine_spec2','fine_spec3','fine_spec4','fine_spec5','fine_spec6','fine_spec7','def_sen1','def_sen2','def_sen3','def_sen4','def_sen5','def_sen6','def_sen7','fine_sen1','fine_sen2','fine_sen3','fine_sen4','fine_sen5','fine_sen6','fine_sen7','def_acc1','def_acc2','def_acc3','def_acc4','def_acc5','def_acc6','def_acc7','fine_acc1','fine_acc2','fine_acc3','fine_acc4','fine_acc5','fine_acc6','fine_acc7'});
    writetable(table1,'Threshold_assessment.xls','Sheet',strcat('epsilon',num2str(epsilonv(fidx))),'WriteVariableNames',false);
end
hold off
view(-145,8)
grid on
legend('Default','Fine','Location','best')
saveas(gcf,fullfile(spath,strcat("im_a",num2str(02),".png")))

fig4 = figure;
[X,Y] = meshgrid(exp_aux,th_p_vector);
hold on
for fidx = 1:1:length(epsilonv)
    surf(X,Y,reshape(mat_acc_def(:,fidx,:),[51,7]),'FaceAlpha',0.45,'FaceLighting','flat','FaceColor',[1,0,0],'LineStyle',':')        
    surf(X,Y,reshape(mat_acc_fine(:,fidx,:),[51,7]),'FaceAlpha',0.45,'FaceLighting','flat','FaceColor',[0,0,1],'LineStyle',':')
    title("Accuracy")
    xlabel("Exp")
    ylabel("th_p",'Interpreter','none')
    zlabel("Accuracy")
end
hold off
view(-145,8)
grid on
legend('Default','Fine','Location','best')    
saveas(gcf,fullfile(spath,strcat("im_a",num2str(03),".png")))



load_control = 0;

% fig1 = figure;
% hold on;
% plot
% plot(th_p_vector,mat_sen_def,'r-');
% plot(th_p_vector,mat_sen_fine,'r--');
% 
% fig1=figure;
% hold on;
% plot(th_p_vector,mat_sen_def,'r+-');
% plot(th_p_vector,mat_spe_def,'bo-');
% plot(th_p_vector,mat_acc_def,'gv-');
% xlabel('th_p','Interpreter','none')
% title(strcat("Experiment=",filename,",  Approach= ",a_dict(approach_sel)),'Interpreter','none')
% legend('Sensitivity','Specificity','Accuracy')
% hold off;
% 
% fig2=figure;
% plot3(mat_sen_def,mat_spe_def,th_p_vector)
% xlabel('Sensitivity')
% ylabel('Specificity')
% zlabel('th_p','Interpreter','none')
% title(strcat("Experiment=",filename,",  Approach= ",a_dict(approach_sel)),'Interpreter','none')
% grid on;

% mat_sen_def
% mat_spe_def
% mat_acc_def
genReport
disp("End myTask section ...")


function genReport

global epsilonv

spath = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\ReportGen";

import mlreportgen.report.*
import mlreportgen.dom.*;

filePattern = fullfile(spath, '*.png'); % Change to whatever pattern you need.
theFiles = struct2table(dir(filePattern));
theFiles = sortrows(theFiles,"name","ascend");

d = Document(fullfile(spath,"myImageReport"),"pdf");
% open(d);

layout = d.CurrentPageLayout;
layout.PageSize.Height = '8.5in';
layout.PageSize.Width = '11in';
layout.PageSize.Orientation  ='landscape';
layout.PageMargins.Left  = '1cm';
layout.PageMargins.Right = '1cm';
layout.PageMargins.Top = '1cm';
layout.PageMargins.Bottom = '1cm';

pagenumber = PageNumber();
layout.Style = [{pagenumber}];

% Create the footer and add a page number to it
myfooter = PDFPageFooter();
para = Paragraph();
para.HAlign = 'center';
append(para,Page());

% Add the page number to the footer
append(myfooter,para);
layout.PageFooters = myfooter;

% f = waitbar(0,"Generating pdf report ...");
% ch = Chapter(strcat("Chapter Report of Exp: ",exp_name));
% append(d,ch)

p = Paragraph(strcat("Peak Detection, Scan th_p = [0,0.05], step=0.001, epsilon = ",num2str(epsilonv)));
p.Bold = true;
append(d,p);
disp(strcat("Generating report, please wait ..."));

for k=1:2
   plot1 = Image(fullfile(theFiles.folder{k},theFiles.name{k}));
   plot1.Style = { HAlign('center') };
   plot1.Height = "3.8in";      
   append(d,plot1); 
end

% if ~isempty(mtable)
    % table1 = MATLABTable(mtable);
    % append(d,table1);
% end

% if ~isempty(data)
    % par=Paragraph(' ');
    % par.Style = {OuterMargin('0.5in','0in','0in','12pt')};
    % append(d,par);
    % fnames = fieldnames(data);
    % for kp=1:1:length(fieldnames(data))
        % par = Paragraph(strcat(fnames{kp}," = ",num2str(data.(fnames{kp}),'%f04')));
        % append(d,par);
    % end
% end

disp(strcat("     ... 0 %"));
for k=3:length(theFiles.name)
   plot1 = Image(fullfile(theFiles.folder{k},theFiles.name{k}));
   plot1.Style = { HAlign('center') };
   plot1.Height = "3.8in";      
%    plot1.Style = {ScaleToFit};
%    plot1.Width = "4in";
   append(d,plot1); 
%    waitbar(k/length(theFiles.name),f)
   disp(strcat("     ... ",num2str(round(k*100/length(theFiles.name)))," %"));
end
disp(strcat("Opening pdf ..."));
rptview(d);
disp(strcat("Report Finished (:"));