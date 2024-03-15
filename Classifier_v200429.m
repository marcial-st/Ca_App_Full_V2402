  function varargout = Classifier_v200429(varargin)
% CLASSIFIER_V200429 MATLAB code for Classifier_v200429.fig
%      CLASSIFIER_V200429, by itself, creates a new CLASSIFIER_V200429 or raises the existing
%      singleton*.
%
%      H = CLASSIFIER_V200429 returns the handle to a new CLASSIFIER_V200429 or the handle to
%      the existing singleton*.
%
%      CLASSIFIER_V200429('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLASSIFIER_V200429.M with the given input arguments.
%
%      CLASSIFIER_V200429('Property','Value',...) creates a new CLASSIFIER_V200429 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Classifier_v200429_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Classifier_v200429_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools m  enu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Classifier_v200429

% Last Modified by GUIDE v2.5 12-Oct-2023 11:48:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Classifier_v200429_OpeningFcn, ...
                   'gui_OutputFcn',  @Classifier_v200429_OutputFcn, ...
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


% --- Executes just before Classifier_v200429 is made visible.
function Classifier_v200429_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Classifier_v200429 (see VARARGIN)

% Choose default command line output for Classifier_v200429
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global CELLS
global im340_1
global im340_n
global process_list
global spath
global exp_name
global params
global net_profiles
global en_ignore_low_resp
global en_show_clusters
global en_classlevel2
global en_classlevel3

% spath = "C:\Users\MST\Google Drive\INAOE Admisión - Marcial\Research\ReportGen";
% 
% theFiles = dir(spath);
% for k = 1 : length(theFiles)
%   baseFileName = theFiles(k).name;
%   fullFileName = fullfile(spath, baseFileName);
%   fprintf(1, 'Now deleting %s\n', fullFileName);
%   delete(fullFileName);
% end

% disp(handles)
% set(handles.lb_classes,'Visible','off')
% load  prof_v190921.mat;
% im340_1=imread('img0-003.tif');
% im340_n=imread('img0-218.tif');
en_ignore_low_resp = 0;
en_show_clusters = 0;
en_classlevel2 = 0;
en_classlevel3 = 0;

CELLS=getappdata(0,'CELLS');
im340_1=getappdata(0,'im340_1');
im340_n=getappdata(0,'im340_n');
process_list = getappdata(0,'process_list');
exp_name = getappdata(0,'exp_name');
params = getappdata(0,'params');
net_profiles = load('net_profiles.mat');
net_profiles = net_profiles.trainedNetwork_1;
set(handles.lb_classes,'Enable','off')
set(handles.pb_injury,'Enable','off')
set(handles.edit_th_p,'Enable','off')
set(handles.pb_aux_th_lp,'Enable','off')
set(handles.edit5,'Enable','off')
set(handles.pb_aux_art,'Enable','off')
set(handles.pb_refresh,'Enable','off')
set(handles.text8,'Enable','off')
set(handles.t_artifact,'Enable','off')
set(handles.chb_plot,'Enable','off')
set(handles.ax_classim,'Visible','off')
set(handles.pb_injury,'Enable','off')
set(handles.pb_emulate,'Enable','off')
set(handles.pb_classify,'Enable','off')
set(handles.chb_plot,'Enable','off')
set(handles.cb_artifacts2,'Enable','off')
set(handles.cb_ignore_low_resp,'Enable','off')
set(handles.cb_point_out_cell,'Enable','off')
set(handles.cb_show_original,'Enable','off')
set(handles.cb_show_clusters,'Enable','off')


% UIWAIT makes Classifier_v200429 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Classifier_v200429_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
getOnlineCells(hObject, handles)

function getOnlineCells(hObject, handles)
    global CELLS
    global cells_online
    global online_features
    global clustTreeEuc
    global online_profiles
    global t
    global spath
    global pclass
    global params
    global ca_prof_smooth
    global t_sample
    global n_frames
    global exp_data
    global im340_1
    global th_prom
    global th_diff
    global th_artifact
    global global_events_bin
    global global_events
    global en_artifact_detector_2
    global Ts
    global en_ignore_low_resp
    global en_show_clusters
    

    if (isempty(CELLS) && isempty(exp_data))
        disp("Empty input data")
    else
        if ~isempty(exp_data)
            CELLS = exp_data.CELLS;
            cells_online = exp_data.cells_online;
            pclass = exp_data.class_vector;
            params = exp_data.params;  
            im340_1 = exp_data.current_img;
        else
        %     t_start = params(1);   % start of experiment
        %     w_mainpk = params(2);
        %     th_plateu = params(3);
        %     e_plateu = params(4);
        %     shft_l = params(5);
        %     shft_r = params(6);
        %     wfilter = params(7);
        %     kmain = params(8);
        %     ksec = params(9);
        %     t_basal_0 = params(10);
        %     t_basal_n = params(11);        
            params =[18,40,0.01,0.25,16,15,wfilter,5,2.75,1,12; %1
                 30,40,0.01,0.25,15,15,wfilter,7,2.25,15,25; %2 1.4 std ratio threshold
                 50,40,0.01,0.25,10,15,wfilter,7,2.75,15,40; %3 
                 40,40,0.01,0.25,10,15,wfilter,7,2.75,15,35; %4 
                 65,40,0.01,0.25,25,15,wfilter,7,2.75,20,40; %5 bueno
                 3*65,40,0.01,0.25,25,15,wfilter,7,2.75,10,55; %6 
                 60,40,0.01,0.25,15,15,wfilter,7,2.75,5,55; %7 bueno
                 45,40,0.01,0.25,10,15,wfilter,7,2.75,10,40; %8 muy pocas muestras, análisis no viable
                 60,40,0.01,0.25,10,15,wfilter,7,2.75,20,50]; %9 
       end
        inj_frame = [33,67,76,69,67,65]; % exacto
        wfilter = 11;

        
        t_sample =3;
        %     guidata(hObject, handles);
        cells_online = CELLS(contains({CELLS.status},'ONLINE'));            % Get only ONLINE cells
        cells_offline = CELLS(contains({CELLS.status},'OFFLINE'));

        Ts = 3;
        en_add_offline_cells = 0;
        delta_frames = 0;
        
        if en_add_offline_cells
            i_2map = 1;
            f_min = floor(delta_frames+exp_data.params.tstim/Ts)
            n_frames = length(cells_online(1).ca_profile);
            disp(strcat("looking for offline profiles beyond ",num2str(f_min)," frames"))
            for i = 1:1:length(cells_offline)
                if length(cells_offline(i).ca_profile)>f_min
                    padval = cells_offline(i).ca_profile(end);
                    i_frame = length(cells_offline(i).ca_profile);
                    cells_offline(i).ca_profile(i_frame+1:n_frames) = padval;
                    cells_offline2map(i_2map) = cells_offline(i); 
                    i_2map = i_2map+1;
                end
            end
            global_cells = [cells_online,cells_offline2map];
            cells_online = global_cells;
        end
        th_artifact =  str2double(get(handles.edit5,'String'));
        if isempty(th_artifact)
            th_artifact = 7;
            disp("Warning. using th_artifact=7 default value")
        end
        th_diff =  str2double(get(handles.edit_th_diff,'String'));
        if isempty(th_diff)
            th_diff = 0.04;
            disp("Warning. using th_diff=0.04 default value")
        end
        th_prom =  str2double(get(handles.edit_th_p,'String'));
        if isempty(th_prom)
            th_prom = 0.024;
            disp("Warning. using th_prom=0.024 default value")
        end
        
        Ts = 3;
        
        
        peak_width_f = 30;

        n_online = max(size(cells_online))        
        n_frames = max(size(cells_online(1).ca_profile));            % only valid since Class2 are online cells    
        t = 0:1:n_frames-1;
        t = t.*t_sample;
        if isempty(exp_data)
            cells_online(end)=[];
        end
        online_profiles = reshape([cells_online(:).ca_profile],n_frames,[]); % profiles array Profile-Column
        [ca_prof_smooth,~] = smoothProfiles(online_profiles',11);
        %     pclass = class_expert_v201007(online_profiles',params(7,:),0);
        %     pclass = class_expert_reduced_v2105(online_profiles',params(6,:),0);
        %     pclass = class_net(cells_online,net_profiles);
        if ~isempty(exp_data)
            pclass = [];
            global_events_bin = [];
            global_events= [];
            [pclass,global_events_bin,global_events] = classifier_leftprom_basic4(ca_prof_smooth,th_artifact,n_frames,Ts,th_diff,params,th_prom,peak_width_f,en_artifact_detector_2);
            
        end
    end
   
        
function edit_nclass_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nclass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nclass as text
%        str2double(get(hObject,'String')) returns contents of edit_nclass as a double


% --- Executes during object creation, after setting all properties.
function edit_nclass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nclass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_classify.
function pb_classify_Callback(hObject, eventdata, handles)
% hObject    handle to pb_classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NCLS
global CELLS
global clustTreeEuc
global im340_1
global cells_online
global spath
global pclass
global clabels
global imgClass
global t_sample
global exp_data
global filename
global inj_band
global inj_bin
global pix2um
global exp_name
global t
global online_profiles
global ca_prof_smooth
global global_events
global global_events_bin
global en_artifact_detector_2
global th_prom
global th_diff
global th_artifact
global n_frames
global Ts
global params
global nNCLS
global imgLabel
global en_ignore_low_resp
global en_show_clusters
global en_classlevel2
global en_classlevel3
global hist_stack_normalized
global bbox_wh_mean
% Best colors for future
% colors = [0.2588,0.7098,0.8588; 
%     0.5765,0.8588,0.2588;
%     0.8588,0.5647,0.2588;
%     0.7765,0.2667,0.8588;
%     0.5,0.5,0;
%     0,1,0;
%     0,0.5,0;
%     0,1,1;
%     0,0.5,0.5;
%     0,0,1;
%     0,0,0.5;
%     1,0,1;
%     0.5,0,0.5;
%     0.7569,0.4706,0.0471;
%     0.9451,0.7686,0.0588;
%     0.8235,0.7059,0.8706;
%     ];


colors = [0.2431,0.4784,0.9216;
    1,0,0;
    0.5,0,0;
    1,1,0;
    0.5,0.5,0;
    0,1,0;
    0,0.5,0;
    0,1,1;
    0,0.5,0.5;
    0,0,1;
    0,0,0.5;
    1,0,1;
    0.5,0,0.5;
    0.7569,0.4706,0.0471;
    0.9451,0.7686,0.0588;
    0.8235,0.7059,0.8706;
    ];

% pix2um = 0.4277;
set(handles.cb_point_out_cell,'Visible','on')
set(handles.cb_point_out_cell,'Enable','on')

% cindex = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
% clabels = ["Unknown";
%     "Non response";
%     "Peak";
%     "Peak+plateau";
%     "Peak+oscillations";
%     "Peak+plateu+oscillations";
%     "Oscillations";
%     "Oscillations+plateau";
%     "Basal activity+Peak";
%     "Basal activity+Peak+plateau";
%     "Basal activity+Peak+oscillations";
%     "Basal activity+ Peak+plateu+oscillations";
%     "Basal activity+Oscillations";
%     "Basal activity+Oscillations+plateau";
%     "Basal activity";
%     "Basal activity +plateu";
%     ]

% cindex = [0,1,2,3,4];
% clabels = ["Unknown"; %0
%     "Non response"; %1
%     "Basal"; %2
%     "Response"; %3
%     "Basal+Response";%4
%     ];
cindex = [0,1,2,3];
clabels = ["BR=0,RR=0"; %0
    "BR=1,RR=0"; %1
    "BR=0,RR=1"; %2
    "BR=1,RR=1"; %3
    ];
cdict = containers.Map(cindex,clabels);

if ~isempty(exp_data)
    bbox_wh_mean = mean(reshape([CELLS.bbox_wh],[2 length([CELLS.bbox_wh])/2])',1);
end

inj_band = getappdata(0,'inj_band');
inj_bin = getappdata(0,'inj_bin');
if not(isempty(inj_band))
    inj_band_center = inj_band(1)+((inj_band(end)-inj_band(1))/2);
    inj_band_bu = inj_band;
    inj_band(1) = inj_band_center-round(bbox_wh_mean(1)/2); % Update inj_band to shrink the injury region
    inj_band(end) = inj_band_center+round(bbox_wh_mean(1)/2);
    % inj_bin(:,1:inj_band(1)) = 0;
    % inj_bin(:,inj_band(end):end) = 0;
    inj_bin(:,1:inj_band_bu(1)) = 0;
    inj_bin(:,inj_band_bu(end):end) = 0;
end

if en_artifact_detector_2
   pclass = [];
   global_events_bin = [];
   global_events= [];
   [pclass,global_events_bin,global_events] = classifier_leftprom_basic4(ca_prof_smooth,th_artifact,n_frames,Ts,th_diff,params,th_prom,20,en_artifact_detector_2);
end

NCLS = [];
if ~isempty(exp_data)
    current_img_p = processListExecute(exp_data.current_img,exp_data.process_list);
    LO=bwconncomp(current_img_p);
    imgLabel=labelmatrix(LO);
    disp('imgLabel computed')
else
    imgLabel=getappdata(0,'LBase');
    disp('imgLabel got from exp_data file')
end

pclass_u = unique(pclass);

n_classes = length(pclass_u);

set(handles.lb_classes,'Visible','on')
set(handles.pb_injury,'Visible','on')
string_class = cell(1,n_classes);

for i=1:1:length(pclass_u)
    NCLS(i).cVector = find(pclass==pclass_u(i))';
    disp(strcat("size(Class ",num2str(pclass_u(i))," ) = ",num2str(size(NCLS(i).cVector))))
    string_class(i) = {strcat("Class ",num2str(pclass_u(i)),": ",values(cdict,{pclass_u(i)}))};
%     NCLS(i).label = strcat("Class",num2str(pclass_u(i),'%02u'));
    NCLS(i).label = string_class(i);
    NCLS(i).color = 255.*colors(pclass_u(i)+1,1:3);
    NCLS(i).title = strcat("Class ",num2str(pclass_u(i)),": ",values(cdict,{pclass_u(i)}));
end 

if ~isempty(exp_data)
    NCLS = classAnalizer(NCLS,cells_online,0);
else
    NCLS = classAnalizer(NCLS,CELLS,0);
end

% cla(handles.ax_classim) 
[h_i,~] = size(im340_1);
if isempty(inj_band)
    en_ind_class = 0;
else
    en_ind_class = 1;
end
en_gen_class = 0;
en_global_class = 1;
en_plot_clusters_in_global = en_show_clusters; 
en_visual_per_cluster = en_show_clusters;
en_y_histograms = 0;

pix2um = get_px2num(im340_1,40,0);

if en_global_class
    fig_all_classes = figure
    fig_all_classes = gca
    % axes('Position',[0.3 0.3 0.7 0.7])
    [imgClassj,~] = imshowCellClasses(NCLS,im340_1,cells_online,imgLabel);
    
    newimg = imgClassj;
      % imshow(newimg)
    if not(isempty(inj_bin))
        newimg(:,:,1) = newimg(:,:,1)+uint8(inj_bin(:,:,1)).*255;
        imshow(newimg)
        rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
    end
    colors = reshape([NCLS.color],3,[])';
    imlegend(colors/255,[NCLS.label])
    newimg_global = newimg;

    
        hold on
        for j=1:1:length(NCLS)
            cells_ci = cells_online(NCLS(j).cVector);
            cells_ci_loc= zeros(length(cells_ci),2);
            for i_cell = 1:1:length(cells_ci);
                cells_ci_loc(i_cell,:) = cells_ci(i_cell).xy_hist(1,:);
            end
            [IDX,C,SUMD,K] = best_kmeans(cells_ci_loc);            
            if en_plot_clusters_in_global
                newimg = insertMarker(newimg,floor(C),'+','color',NCLS(j).color,'size',12);
                for jC=1:1:K
                    cent_v = (ones(length(cells_ci),2).*floor(C(jC,:)));
                    dist_v = sqrt(diag((cent_v-cells_ci_loc)*(cent_v-cells_ci_loc)'));
                    cent_rad_max = floor(max(dist_v(IDX==jC)));
                    cent_rad_mean = floor(mean(dist_v(IDX==jC)));
                    drawellipse('Center',floor(C(jC,:)),'SemiAxes',[cent_rad_max,cent_rad_max],'StripeColor',(NCLS(j).color/255));
                    drawellipse('Center',floor(C(jC,:)),'SemiAxes',[cent_rad_mean,cent_rad_mean],'StripeColor',(NCLS(j).color/255));
                end            
            end
        end        
        % title("All clusters")
        title(filename,'Interpreter','none')
        hold off    
end


if en_ind_class
    % cells_ci_loc_y = []; cells_ci_loc_y_r = []; cells_ci_loc_y_l = []; % Y
    for j=1:1:length(NCLS)
        cells_ci = cells_online(NCLS(j).cVector);
        % if not(isempty(inj_band))
            % cells_ci_loc = zeros(length(cells_ci),3);
        % else
            cells_ci_loc= zeros(length(cells_ci),2);
        % end

        i_cell_r = 0; i_cell_l = 0;
        cells_ci_loc_x_l = [];cells_ci_loc_x_r = []; 
        cells_ci_d2inj_l = []; cells_ci_d2inj_r = []; cells_ci_d2inj = [];  % X
        cells_ci_loc_y_l = []; cells_ci_loc_y_r = [];
        
        for i_cell = 1:1:length(cells_ci)
            cells_ci_loc(i_cell,1:2) = cells_ci(i_cell).xy_hist(1,:);
            if not(isempty(inj_band))
                % cells_ci_loc_y(j,i_cell) = cells_ci_loc(i_cell,2); 
                if inj_band(1)>=cells_ci_loc(i_cell,1)
                    i_cell_l = i_cell_l+1;
                    cells_ci_loc_x_l(i_cell_l,:) = cells_ci_loc(i_cell,:); 
                    cells_ci_loc_y_l(i_cell_l,:) = cells_ci_loc(i_cell,2); 
                    cells_ci_d2inj(i_cell) = abs(inj_band(1)-cells_ci_loc(i_cell,1))*pix2um;
                elseif inj_band(end)<=cells_ci_loc(i_cell,1)
                    i_cell_r = i_cell_r+1;
                    cells_ci_loc_x_r(i_cell_r,:) = cells_ci_loc(i_cell,:); 
                    cells_ci_loc_y_r(i_cell_r,:) = cells_ci_loc(i_cell,2); 
                    cells_ci_d2inj(i_cell) = abs(inj_band(end)-cells_ci_loc(i_cell,1))*pix2um;
                else
                    cells_ci_loc(i_cell,:) = [];
                    cells_ci_d2inj(i_cell) = 9999;
                end
            end
        end
        band_offset =  floor((inj_band(end)-inj_band(1))/2);
        % if ~isempty(cells_ci_loc_x_l) cells_ci_d2inj_l = abs(inj_band(1)-cells_ci_loc_x_l(:,1))*pix2um; end   % Left side distance
        % if ~isempty(cells_ci_loc_x_r) cells_ci_d2inj_r = abs(inj_band(end)-cells_ci_loc_x_r(:,1))*pix2um; end % Rigth side distance
        if ~isempty(cells_ci_loc_x_l) cells_ci_d2inj_l = abs((inj_band(1)+band_offset)-cells_ci_loc_x_l(:,1))*pix2um; end   % Left side distance
        if ~isempty(cells_ci_loc_x_r) cells_ci_d2inj_r = abs((inj_band(end)-band_offset)-cells_ci_loc_x_r(:,1))*pix2um; end % Rigth side distance

        figure;
        nbins = 20;
        subplot(2,2,1:2);histogram(cells_ci_d2inj,2*nbins,"FaceColor",NCLS(j).color/255); ylabel("Frequency"); xlabel("Distance [um]"); text(300,2.5,strcat("n=",num2str(length(cells_ci_d2inj))));title(strcat(filename,": All cells"));
        if ~isempty(cells_ci_loc_x_l)
            subplot(2,2,3);  histogram(cells_ci_d2inj_l,nbins,"FaceColor",NCLS(j).color/255); ylabel("Frequency"); xlabel("Distance [um]"); text(300,2.5,strcat("n=",num2str(length(cells_ci_d2inj_l))));title(strcat(filename,": Left side cells"));
        end
        if ~isempty(cells_ci_loc_x_r)
            subplot(2,2,4);  histogram(cells_ci_d2inj_r,nbins,"FaceColor",NCLS(j).color/255); ylabel("Frequency"); xlabel("Distance [um]"); text(300,2.5,strcat("n=",num2str(length(cells_ci_d2inj_r))));title(strcat(filename,": Rigth side cells"));
        end
        sgtitle(strcat(filename,": ",NCLS(j).title,", Distance to Injury, X"))

        cells_ci_loc_symetric = [-cells_ci_d2inj_l;cells_ci_d2inj_r];
        figure;
        if (~isempty(cells_ci_loc_symetric))
            [y_max,x_max]=size(im340_1);
            % h_limit = (max([inj_band(1),x_max-inj_band(end)])+10)*pix2um;
            % h_limit = max(abs(cells_ci_loc_symetric))+10;
            h_lim_x_l = -(inj_band(1)+band_offset)*pix2um;  % this is more accurate
            h_lim_x_r= (x_max-(inj_band(end)-band_offset))*pix2um;
            v_edges_x = linspace(h_lim_x_l,h_lim_x_r,2*nbins);
            h_lim_y = y_max*pix2um;
            v_edges_y = linspace(0,h_lim_y,2*nbins);
            % h_limit = floor(x_max*pix2um/2);
            % v_edges = linspace(-h_limit,h_limit,2*nbins);           
            histogram(cells_ci_loc_symetric,'BinEdges',v_edges_x,"FaceColor",NCLS(j).color/255); ylabel("Frequency"); xlabel("Distance [um]"); text(300,2.5,strcat("n=",num2str(length(cells_ci_d2inj_r))));title("Reference Center");
            xlim([h_lim_x_l,h_lim_x_r]);

            hist_stack(:,j)=histcounts(cells_ci_loc_symetric,'BinEdges',v_edges_x);
            hist_stack_normalized = hist_stack./sum(hist_stack,2)*100;
            hist_stack_normalized(isnan(hist_stack_normalized))=0;

            hist_stack_y_l(:,j)=histcounts(cells_ci_loc_y_l*pix2um,'BinEdges',v_edges_y);
            hist_stack_y_l_normalized = hist_stack_y_l./sum(hist_stack_y_l,2)*100;
            hist_stack_y_l_normalized(isnan(hist_stack_y_l_normalized))=0;

            hist_stack_y_r(:,j)=histcounts(cells_ci_loc_y_r*pix2um,'BinEdges',v_edges_y);
            hist_stack_y_r_normalized = hist_stack_y_r./sum(hist_stack_y_r,2)*100;
            hist_stack_y_r_normalized(isnan(hist_stack_y_r_normalized))=0;

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sub Class Analysis START
                        
            switch j
                case 1
                    j
                    disp("Sub class analysis: 0 begin")                
                    disp("Sub class analysis: 0 end")               
                case 2
                    j
                    disp("Sub class analysis: 1 begin")
                    disp("Sub class analysis: 1 end")
                case 3
                    j
                    disp("Sub class analysis: 2 begin")                    

                    [NCLS16,NCLS4,th_vals]=subclass2_classifier(exp_data.spks,[NCLS(j).cVector],3,exp_data.params.tstim,ca_prof_smooth,global_events,global_events_bin);
                    if en_classlevel2
                        fig_sc2_prof_all = figure;
                        plot(t,ca_prof_smooth(NCLS(j).cVector,:));title("All Profiles Class 2");xlim([t(1) t(end)]);xlabel("Time [s]");ylabel("Ratio");

                        fig_sc2_sclass4 = figure;                    
                        [imgSClass2,~] = imshowCellClasses(NCLS4,im340_1,cells_online,imgLabel);
                        if not(isempty(inj_bin))
                            imgSClass2(:,:,1) = imgSClass2(:,:,1)+uint8(inj_bin(:,:,1)).*255;
                            imshow(imgSClass2)
                            rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
                        end
                        title("Class2: Sub class analysis - 4 Class" )
                        colors = reshape([NCLS4.color],3,[])';
                        imlegend(colors/255,[NCLS4.label])
                    end

                    % fig_sc2_sclass16 = figure;                    
                    % [imgSClass3,~] = imshowCellClasses(NCLS16,im340_1,cells_online,imgLabel);
                    % if not(isempty(inj_bin))
                    %     imgSClass3(:,:,1) = imgSClass3(:,:,1)+uint8(inj_bin(:,:,1)).*255;
                    %     imshow(imgSClass3)
                    %     rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
                    % end
                    % title("Class2: Sub class analysis - 16 Class" )
                    % colors = reshape([NCLS16.color],3,[])';
                    % imlegend(colors/255,[NCLS16.label])
                    
                    % figure;imshow(imgSClass2)
                    [NCSL2_0] = subclass2_s0_classifier(exp_data.spks,[NCLS4(1).cVector],t_sample,exp_data.params.tstim,th_vals,ca_prof_smooth,filename,global_events,global_events_bin,en_ignore_low_resp);
                    if en_classlevel3
                        fig_sc2_0_x = figure;                    
                        [imgSClass2,~] = imshowCellClasses(NCSL2_0,im340_1,cells_online,imgLabel);
                        if not(isempty(inj_bin))
                            imgSClass2(:,:,1) = imgSClass2(:,:,1)+uint8(inj_bin(:,:,1)).*255;
                            imshow(imgSClass2)
                            rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
                        end
                        title("Class2.0.X" )
                        colors = reshape([NCSL2_0.color],3,[])';
                        imlegend(colors/255,[NCSL2_0.label])
                    end

 
                    if en_classlevel2
                        [nNCLS,string_class] = meltNCLS(NCLS,NCLS4,string_class);
                    end
                    if en_classlevel3
                        [nNCLS,string_class] = meltNCLS(nNCLS,NCSL2_0,string_class);
                    end
                    if (~en_classlevel2&&~en_classlevel3)
                        nNCLS = NCLS;
                    end


                    if (en_classlevel2||en_classlevel3)
                        fig_cdiagram = figure;
                        imshow(imread('class_diagram.png'))

                        fig_melt = figure;                    
                        [imgMelt,~] = imshowCellClasses(nNCLS,im340_1,cells_online,imgLabel);
                        if not(isempty(inj_bin))
                            imgMelt(:,:,1) = imgMelt(:,:,1)+uint8(inj_bin(:,:,1)).*255;
                            imshow(imgMelt)
                            rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
                        end
                        title("Melt Classes" )
                        colors = reshape([nNCLS.color],3,[])';
                        imlegend(colors/255,[nNCLS.label])
                    end

                    disp("Sub class analysis: 2 end")
                case 4
                    j
                    disp("Sub class analysis: 3 begin")
                    disp("Sub class analysis: 3 end")
                otherwise
                    disp("sub class: custom")
            end            
            
            % Sub Class Analysis END
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            newimg_global = insertMarker(newimg_global,cells_ci_loc(:,1:2),'*','color','green','size',8);
            if j==length(NCLS)
                figure;
                subplot(2,1,1)
                barabs = bar(v_edges_x(2:end),hist_stack,'stacked');title("Absolute");xlabel("Distance [um]");ylabel("Freq")
                legend('Class_0','Class_1','Class_2','Class_3','Interpreter','none')
                for j_bar =1:1:length(NCLS)
                    barabs(j_bar).FaceColor = NCLS(j_bar).color./255;
                end

                subplot(2,1,2)
                barnorm = bar(v_edges_x(2:end),hist_stack_normalized,'stacked');title("Normalized %");xlabel("Distance [um]");ylabel("Freq")
                legend('Class_0','Class_1','Class_2','Class_3','Interpreter','none')
                for j_bar =1:1:length(NCLS)
                    barnorm(j_bar).FaceColor = NCLS(j_bar).color./255;
                end

                sgtitle(strcat(filename,": Class-distance Stats"))
                %================================
                %---- Plot figure + histogram----
                %================================
                fig_all_classes_clone = figure
                % fig_all_classes_clone.Units='normalized';
                % fig_all_classes_clone.Position=[0 0 .6 .8];                
                mesh_x = floor((v_edges_x+abs(min(v_edges_x)))/pix2um)+1;
                mesh_y = floor((v_edges_y)/pix2um);mesh_y(1)=1;
                newimg_mesh = uint8(zeros(size(newimg_global)));
                % newimg_mesh(:,mesh_x,:)=255;
                % newimg_mesh(mesh_y,:,:)=255;
                newimg_global(:,mesh_x,:)=0;
                newimg_global(mesh_y,:,:)=0;
                % newimg_global = newimg_global+newimg_mesh;



                ax_fluo = axes('Position',[0.2 0.25 0.6 0.6])                
                imshow(newimg_global,'Parent',ax_fluo)
                if not(isempty(inj_bin))
                    rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
                end
                colors = reshape([NCLS.color],3,[])';
                imlegend(colors/255,[NCLS.label])
                % imlegend([],[])
                %================================             
                axes('Position',[0.19 0.1 0.62 0.15])
                barnorm = bar(v_edges_x(2:end),hist_stack_normalized,'stacked');xlabel("Distance [um]");ylabel("Frequency");%title("Normalized %")                
                for j_bar =1:1:length(NCLS)
                    barnorm(j_bar).FaceColor = NCLS(j_bar).color./255;
                end
                xbarCnt = vertcat(barnorm.XEndPoints); 
                ybarTop = vertcat(barnorm.YEndPoints);
                ybarCnt = ybarTop - hist_stack_normalized'/2; 
                txt = compose('%d',hist_stack');
                hist_stack_bin = hist_stack'>0;
                th = text(xbarCnt(hist_stack_bin), ybarCnt(hist_stack_bin), txt(hist_stack_bin),'HorizontalAlignment', 'center','VerticalAlignment', 'middle','Color', 'k','FontSize', 7);
                % hist_stack_normalized = hist_stack./sum(hist_stack,2);
                % hist_stack_normalized(isnan(hist_stack_normalized))=0;
                %================================
if en_y_histograms                
                axes('Position',[0.1 0.2375 0.09 0.625])
                barnorm = barh(v_edges_y(2:end),hist_stack_y_l_normalized,'stacked');ylabel("Distance [um]");xlabel("Freq");%title("Normalized %")
                set(gca, 'YDir','reverse')                
                for j_bar =1:1:length(NCLS)
                    barnorm(j_bar).FaceColor = NCLS(j_bar).color./255;
                end
                xbarCnt = vertcat(barnorm.XEndPoints); 
                ybarTop = vertcat(barnorm.YEndPoints);
                ybarCnt = ybarTop - hist_stack_y_l_normalized'/2; 
                txt = compose('%d',hist_stack_y_l');
                hist_stack_y_l_bin = hist_stack_y_l'>0;
                th = text(ybarCnt(hist_stack_y_l_bin), xbarCnt(hist_stack_y_l_bin), txt(hist_stack_y_l_bin),'HorizontalAlignment', 'center','VerticalAlignment', 'middle','Color', 'k','FontSize', 7);            
                %================================
                axes('Position',[0.81 0.2375 0.09 0.625])
                axy2 = gca;
                barnorm = barh(v_edges_y(2:end),hist_stack_y_r_normalized,'stacked');
                set(gca, 'YDir','reverse')
                xlabel("Freq");
                yyaxis right                
                set(gca, 'YDir','reverse')
                ylabel("Distance [um]");%title("Normalized %")               
                for j_bar =1:1:length(NCLS)
                    barnorm(j_bar).FaceColor = NCLS(j_bar).color./255;
                end
                xbarCnt = vertcat(barnorm.XEndPoints); 
                ybarTop = vertcat(barnorm.YEndPoints);
                ybarCnt = ybarTop - hist_stack_y_r_normalized'/2; 
                txt = compose('%d',hist_stack_y_r');
                hist_stack_y_r_bin = hist_stack_y_r'>0;
                th = text(ybarCnt(hist_stack_y_r_bin), xbarCnt(hist_stack_y_r_bin), txt(hist_stack_y_r_bin),'HorizontalAlignment', 'center','VerticalAlignment', 'middle','Color', 'k','FontSize', 7);
end

                sgtitle(strcat(filename))
                %-------------------------------
            end
        end

            if not(isempty(inj_band))
            aux_cells_ciloc=cells_ci_loc==[0,0];
            aux_idx_cells_ciloc=find(aux_cells_ciloc(:,1));
            if not(isempty(aux_idx_cells_ciloc))
                cells_ci_loc(aux_idx_cells_ciloc,:)=[];
            end            
        end

        if not(isempty(inj_band))

            IDX = [];  C=[];  K=0;
            IDXc = []; Cc=[]; Kc=0; cells_ci_loc_c=[];

            n_kmeans = 2;

            for i_kmeans=1:1:n_kmeans                
                if not(isempty(cells_ci_loc_x_l))
                    [IDXl,Cl,~,Kl] = best_kmeans(cells_ci_loc_x_l);   %IDX: cluster index, C: cluster centroid, K: number of clusters
                else
                    IDXl=[];Cl=[];Kl=0;
                end
                
                if not(isempty(cells_ci_loc_x_r))
                    [IDXr,Cr,~,Kr] = best_kmeans(cells_ci_loc_x_r);   %IDX: cluster index, C: cluster centroid, K: number of clusters
                else
                    IDXr=[];Cr=[];Kr=0;
                end
                
                if ~isempty(IDXl)
                    IDX = [IDXl;IDXr+max(IDXl)];
                else
                    IDX = [IDXr];
                end                
                C   = [Cl;Cr];
                K   = Kl+Kr;

                disp(strcat("   k iteration ",num2str(i_kmeans),"/",num2str(n_kmeans)))

                if i_kmeans==1
                    IDXc = [IDX];
                else
                    IDXc = [IDXc;IDX+max(IDX)];
                end                
                Cc = [Cc;C];
                Kc = Kc+K;
                cells_ci_loc_c=[cells_ci_loc_c;cells_ci_loc];

            end
        else
            [IDX,C,~,K] = best_kmeans(cells_ci_loc);   %IDX: cluster index, C: cluster centroid, K: number of clusters
        end
        if en_visual_per_cluster
            figure
            [imgClassj,~] = imshowCellClasses(NCLS(j),im340_1,cells_online,imgLabel);
            newimg = insertMarker(imgClassj,floor(C(:,1:2)),'+','color',NCLS(j).color,'size',12);
            newimg = insertMarker(newimg,cells_ci_loc(:,1:2),'*','color','black','size',10);
            if not(isempty(inj_bin))
                newimg(:,:,1) = newimg(:,:,1)+uint8(inj_bin(:,:,1)).*255;
            end
            imshow(newimg)
            if not(isempty(inj_bin))
                rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
            end
            title(NCLS(j).title)
            C=floor(C);
            for jC=1:1:K
                cent_v = ones(length(cells_ci_loc),2).*C(jC,1:2);
                dist_v = sqrt(diag((cent_v-cells_ci_loc(:,1:2))*(cent_v-cells_ci_loc(:,1:2))'));
                cent_rad_max = floor(max(dist_v(IDX==jC)));
                cent_rad_mean = floor(mean(dist_v(IDX==jC)));
                drawellipse('Center',floor(C(jC,1:2)),'SemiAxes',[cent_rad_max,cent_rad_max],'StripeColor',(NCLS(j).color/255));
                drawellipse('Center',floor(C(jC,1:2)),'SemiAxes',[cent_rad_mean,cent_rad_mean],'StripeColor',(NCLS(j).color/255));
            end
    
            if not(isempty(inj_band))
                % class_data(j).cluster_features = cluster_feature_extraction(IDX,C(:,1:2),K,cells_ci_loc(:,1:2),inj_band);
                class_data(j).cluster_features = cluster_feature_extraction(IDXc,Cc(:,1:2),Kc,cells_ci_loc_c(:,1:2),inj_band,pix2um);
                plot_features(class_data(j),NCLS(j).color/255,NCLS(j).title)
                sgtitle(strcat(NCLS(j).title,", Features"))
            end
        end
    end
end
if en_gen_class
figure;
[imgClass,NCLS] = imshowCellClasses(NCLS,im340_1,cells_online,imgLabel);
h_gca = gca;
if ~isempty(filename)
    title(strcat(filename,": Spatial classification"),'Interpreter','none')
else
    title('Spatial classification')
end
if not(isempty(inj_bin))
    imgClass(:,:,1) = imgClass(:,:,1)+uint8(inj_bin(:,:,1)).*255;  
    imshow(imgClass)
    rectangle('Position',[inj_band(1) 1 inj_band(end)-inj_band(1) h_i-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
end
end

% saveas(h_gca,fullfile(spath,"Class.png"));    
set(handles.lb_classes,'String',string_class)
set(handles.pb_report,'Enable','on')
set(handles.cb_point_out_cell,'Enable','on')
set(handles.cb_show_original,'Enable','on')
set(handles.menu_tsample,'Enable','on')
set(handles.menu_save,'Enable','On')
 if not(isempty(inj_bin))
    set(handles.pb_emulate,'Enable','on')
 end



% --- Executes on selection change in lb_classes.
function lb_classes_Callback(hObject, eventdata, handles)
% hObject    handle to lb_classes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_classes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_classes
global NCLS
global online_profiles
global ca_prof_smooth
global t
global spath
global selclass
global params
global t_sample
global en_plot
global nNCLS

if not(isempty(nNCLS))
    NCLS_backup = NCLS;
    NCLS = [];
    NCLS = nNCLS;
end

selclass = get(handles.lb_classes,'Value')
% figure;
if en_plot
    % plotClasses(online_profiles,NCLS(selclass),selclass,t,params,t_sample)
    plotClasses(ca_prof_smooth',NCLS(selclass),selclass,t,params,t_sample)
    h_gca = gca;
    stitle = strcat("class",num2str(selclass,'%02d'));
end
% title(stitle)
% saveas(h_gca,fullfile(spath,strcat(NCLS(selclass).label,".png")));

set(handles.text_info,'String',strcat("T=",num2str(t_sample)," s"))
set(handles.text3,'Visible','on')
set(handles.slider1,'Visible','on')
set(handles.slider1, 'Max',length(NCLS(selclass).cVector));
set(handles.slider1, 'Min',1);
set(handles.slider1, 'Value',1);
set(handles.slider1, 'SliderStep', [1/length(NCLS(selclass).cVector) 0.1]);
slider1_Callback(hObject, eventdata, handles)
set(handles.menu_save,'Enable','On')





% --- Executes during object creation, after setting all properties.
function lb_classes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_classes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_emulate.
function pb_emulate_Callback(hObject, eventdata, handles)
% hObject    handle to pb_emulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NCLS
global inj_band
global CELLS
global en_classlevel3
global en_classlevel2
global nNCLS
global cells_online
global params
global ca_prof_smooth
global pix2um
global hist_stack_normalized
global fitmodel
global t
global global_events_bin
global global_events
global en_ignore_low_resp
global th_vals

% inj_band = getappdata(0,'inj_band');

video_en=1;
current_path = pwd;
video_pathv = fullfile(pwd,'videogen');
% video_pathv="C:\Users\marci\OneDrive\Documentos\MATLAB\synthimage";


if not(isempty(inj_band))
    [nNCLS] = classGetDistance(nNCLS,CELLS,inj_band);
    % [fis]=simFisGen(NCLS);
end

probabilities_x = fliplr(hist_stack_normalized(1:floor(end/2),:)');
probabilities_x(:,not(sum(probabilities_x) > 0))=[]; % clean up columns zeros



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% este es el chido
%%%      Profiles Class characterization
%% ----------------------------------------------
if not(isempty(cells_online))
    % nNCLS = classAnalizer(nNCLS,cells_online,1);
    %%%      Class 0 - BR=0, RR=0
    fitmodel_c0 = 'poly3';
    [nNCLS] = classAnalizer_c0(nNCLS,cells_online,ca_prof_smooth,params,inj_band,pix2um,fitmodel_c0);
    %%%      Class 2 - BR=0, RR=1
    fitmodel_c2_basal = 'poly3';
    fitmodel_c2_response = 'sin8';
    [nNCLS] = classAnalizer_c2(nNCLS,cells_online,ca_prof_smooth,params,inj_band,pix2um,global_events,global_events_bin,fitmodel_c2_basal,fitmodel_c2_response,th_vals,en_ignore_low_resp);

    % fit_line = coeffs(1)*cells_ci_d2inj+coeffs(2);
%% ----------------------------------------------
%%%      Class 1 - BR=0, RR=1
end
prompt = {'Rows of cells:','Columns of cells:'};
dlgtitle = 'Matrix size';
dims = [1 35];
definput = {'15',num2str(length(probabilities_x))};

msize = inputdlg(prompt,dlgtitle,dims,definput);
% msize = 1;
grid_en = 0;
if not(isempty(msize))
    sim_row = floor(str2double(msize(1)));
    sim_col = floor(str2double(msize(2)));
    if sim_col>length(probabilities_x)
        warning("User selected more columns than available form experimental data, cells out of bounds will repeat the boundary probabilities")
    end
    %%%% For Thesis figures begin
    % sim_row = 12;   
    % sim_col = 24;
    %%%% For Thesis figures end
    [sim_cell,sim_canv,sim_cw,sim_ch] = sytheticFrameGen(sim_row,sim_col,grid_en);     % Generate simulation cells & canvas
    synt_fig = figure;
    figure(synt_fig)
    sytheticFramePlot(sim_cell,sim_canv)
    % [inj_x,~] = ginput(1);    % Get injury and plot
    % inj_x = round(inj_x);
    %%%% For Thesis figures begin
    [~,inj_x] = size(sim_canv);
    inj_x = round(inj_x/2);
    %%%% For Thesis figures end
    [m_canv,n_canv]=size(sim_canv); 
    color_matrix = reshape([NCLS(:).color],3,[])'/255;
    [sim_cell] = syntheticClassAssign(sim_cell,inj_x,color_matrix,probabilities_x,sim_col,sim_row,m_canv,n_canv,sim_cw,sim_ch,pix2um);
    % [sim_cell] = simClassify(sim_cell,inj_x,fis,color_matrix); % Fuzzy
    % approach
    sytheticFramePlot(sim_cell,sim_canv)
    rectangle('Position',[inj_x 1 1 m_canv-1],'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3)
    [sim_cell] = applyCellModelpClass(sim_cell,nNCLS,t);
    sim_cell = scaleRandomClasses(sim_cell);
    sytheticFramePlotSeq(sim_cell,sim_canv,video_en,video_pathv,inj_x)
end



% --- Executes on button press in pb_injury.
function pb_injury_Callback(hObject, eventdata, handles)
% hObject    handle to pb_injury (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im340_n
global process_list
global exp_data

if ~isempty(exp_data)
    img_inj = exp_data.inj_img;
    setappdata(0,'img_inj',img_inj) 
else
    setappdata(0,'img_inj',im340_n)
end
setappdata(0,'process_list',process_list);

InjuryDetector
set(handles.pb_classify,'Enable','on')
set(handles.pb_emulate,'Visible','on')


% --- Executes on button press in pb_report.
function pb_report_Callback(hObject, eventdata, handles)
% hObject    handle to pb_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global exp_name
global spath
import mlreportgen.dom.*

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


p = Paragraph(strcat("Report of Exp: ",exp_name));
p.Bold = true;
append(d,p);
for k=1:length(theFiles.name)
   plot1 = Image(fullfile(theFiles.folder{k},theFiles.name{k}));
   plot1.Style = { HAlign('center') };
   plot1.Height = "3.8in";      
%    plot1.Style = {ScaleToFit};
%    plot1.Width = "4in";
   append(d,plot1); 
end

close(d);
rptview(d);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global CELLS
global cells_online
global NCLS
global online_profiles
global t
global selclass
global params
global point_out
global imgClass
global ca_prof_smooth
global show_original
global t_sample
global nNCLS
global im340_1
global imgLabel
global en_ignore_low_resp
global en_show_clusters

selp = get(handles.slider1,'Value');
if length(NCLS(selclass).cVector)<2
    set(handles.slider1,'Enable','Off')
else
    set(handles.slider1,'Enable','On')
end
set(handles.text3,'String',num2str(floor(selp)));
set(handles.text4,'String',strcat("/ ",num2str(length(NCLS(selclass).cVector))));
set(handles.text5,'String',strcat("P:",num2str(NCLS(selclass).cVector(floor(selp)))));
axes(handles.ax_classim)
cla(handles.ax_classim)

if show_original
    plot(t,online_profiles(:,NCLS(selclass).cVector),'Color',[0.8,0.8,0.85]); hold on
    plot(t,online_profiles(:,NCLS(selclass).cVector(floor(selp))),'LineWidth',1.5,'Color',NCLS(selclass).color./255)
else
    plot(t,ca_prof_smooth(NCLS(selclass).cVector,:),'Color',[0.8,0.8,0.85]); hold on
    plot(t,ca_prof_smooth(NCLS(selclass).cVector(floor(selp)),:),'LineWidth',1.5,'Color',NCLS(selclass).color./255)    
end

if isstruct(params)
%     xline(params.tstart,'--',{'t_{start}'},'LabelVerticalAlignment','bottom','Color',[0.45,0.45,0.45])
    xline(params.tstim,'--',{strcat("t_{inj}=",num2str(params.tstim),"s")},'LabelVerticalAlignment','top','Color',[0.45,0.45,0.45])
    xlim([t(1) params.tend]);
else
    xline(params(1)*t_sample,'--',{'t_{start}'},'LabelVerticalAlignment','bottom','Color',[0.45,0.45,0.45])
    xlim([t(1) t(end)]);
end
% set(gca,'Color',[0.1392,0.1471,0.2020])
ylim([min(min(online_profiles)) max(max(online_profiles))])
xlabel('t [s]');ylabel('Ratio')
title(NCLS(selclass).title);
if point_out
    imgc_local = figure;
    figure(imgc_local)
    [imgClassInd,~] = imshowCellClasses(nNCLS(selclass),im340_1,cells_online,imgLabel);    
    idx = NCLS(selclass).cVector(floor(selp));
    rectangle('Position',[cells_online(idx).xy_hist(1,1)-10 cells_online(idx).xy_hist(1,2)-10 cells_online(idx).bbox_wh(1)+10 cells_online(idx).bbox_wh(2)+10],'LineWidth',2,'EdgeColor','r','LineStyle',':')
    title(strcat("Profile ",num2str(idx)))
end





% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in cb_point_out_cell.
function cb_point_out_cell_Callback(hObject, eventdata, handles)
% hObject    handle to cb_point_out_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_point_out_cell
global point_out;

point_out = get(handles.cb_point_out_cell,'Value');
slider1_Callback(hObject, eventdata, handles)

% --- Executes on button press in cb_show_original.
function cb_show_original_Callback(hObject, eventdata, handles)
% hObject    handle to cb_show_original (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_show_original
global show_original;

show_original = get(handles.cb_show_original,'Value');
slider1_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_settings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tsample_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
global t_sample
global n_frames

prompt = {'Enter T sample please'};
dlgtitle = 'Input T sample';
dims = [1 32];
definput = {'3'};
t_sample = inputdlg(prompt,dlgtitle,dims,definput);
if not(isempty(t_sample))
    t_sample = str2double(t_sample);
end
if (isempty(t_sample)||(t_sample == 0))    
    t_sample = 3;
    fwr = warndlg('T sample empty or 0, default value 3s set instead');
end
set(handles.text_info,'String',strcat("T=",num2str(t_sample)," s"))
t = 0:1:n_frames-1;
t = t.*t_sample;
slider1_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NCLS
global cells_online
global pclass
global params
global global_events_bin
global global_events
global th_prom
global th_artifact

calcium_true_events = global_events.*double(global_events_bin);
calcium_all_candidates = global_events;
parameters = params;
parameters.tinj = params.tstim;
parameters.th_prom = th_prom;
parameters.th_art = th_artifact;
parameters = rmfield(parameters,{'tstim','inj_w','inj_x'});
parameters = struct2table(parameters);

defname = strcat("analysis_",char(datetime('now','Format','yMMd')),".xls");
[filename,path,uaccept] = uiputfile('*.xls','Save File',defname);
if uaccept ==1    
%     save(fullfile(path,filename),'NCLS','cells_online','pclass','calcium_true_events','calcium_all_candidates')
%     filename_xls = strrep(filename,'.mat','.xls');
    writetable(parameters,fullfile(path,filename),'Sheet','parameters')
    writematrix(pclass,fullfile(path,filename),'Sheet','class_vectors')
    writematrix(calcium_true_events,fullfile(path,filename),'Sheet','calcium_true_events')
    writematrix(calcium_all_candidates,fullfile(path,filename),'Sheet','calcium_all_candidates')
    msgbox('Data structures saved, thanks a lot!');
end


% --------------------------------------------------------------------
function menu_load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global exp_data
global filename
global selclass
[filename, pathname] = uigetfile('*.mat','Select the data structure');
if filename~=0
    load(fullfile(pathname,filename));
    getOnlineCells(hObject, handles)
    set(handles.t_file,'String',fullfile(pathname,filename))
    set(handles.lb_classes,'Enable','on')
    set(handles.pb_injury,'Enable','on')
    set(handles.edit_th_p,'Enable','on')
    set(handles.pb_aux_th_lp,'Enable','on')
    set(handles.edit5,'Enable','on')
    set(handles.pb_aux_art,'Enable','on')
    set(handles.pb_refresh,'Enable','on')
    set(handles.text8,'Enable','on')
    set(handles.t_artifact,'Enable','on')
    set(handles.chb_plot,'Enable','on')
    set(handles.pb_classify,'Enable','on')
    set(handles.ax_classim,'Visible','on')    
     set(handles.pb_emulate,'Enable','off')
    set(handles.pb_classify,'Enable','off')
    set(handles.chb_plot,'Enable','on')
    set(handles.cb_artifacts2,'Enable','on')
    set(handles.cb_ignore_low_resp,'Enable','on')
    set(handles.cb_point_out_cell,'Enable','on')
    set(handles.cb_show_original,'Enable','on')
    set(handles.cb_show_clusters,'Enable','on')
    pb_classify_Callback(hObject, eventdata, handles)
    selclass = 0;
    lb_classes_Callback(hObject, eventdata, handles)    
end



function edit_th_p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_th_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_th_p as text
%        str2double(get(hObject,'String')) returns contents of edit_th_p as a double
global th_prom
th_prom = str2double(get(handles.edit_th_p,'String'));
pb_classify_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_th_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_th_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_th_diff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_th_diff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_th_diff as text
%        str2double(get(hObject,'String')) returns contents of edit_th_diff as a double
global th_diff
th_diff = str2double(get(handles.edit_th_diff,'String'));
pb_classify_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_th_diff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_th_diff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
global th_artifact
th_artifact =  str2double(get(handles.edit5,'String'));
pb_classify_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_aux_th_lp.
function pb_aux_th_lp_Callback(hObject, eventdata, handles)
% hObject    handle to pb_aux_th_lp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global CELLS
    global cells_online
    global online_features
    global clustTreeEuc
    global online_profiles
    global t
    global spath
    global pclass
    global params
    global ca_prof_smooth
    global t_sample
    global n_frames
    global exp_data
    global th_prom
setappdata(0,'ca_prof_smooth',ca_prof_smooth) 
setappdata(0,'t',t) 
setappdata(0,'th_prom',th_prom)
AuxApp_Th_P


% --- Executes on button press in pb_refresh.
function pb_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to pb_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global th_prom
global th_artifact

th_prom=getappdata(0,'th_prom');
if ~isempty(th_prom)
    set(handles.edit_th_p,'String',num2str(th_prom))
end
th_artifact=getappdata(0,'th_artifact');
if ~isempty(th_artifact)
    set(handles.edit5,'String',num2str(th_artifact))
end
getOnlineCells(hObject, handles)
pb_classify_Callback(hObject, eventdata, handles)


% --- Executes on button press in pb_aux_art.
function pb_aux_art_Callback(hObject, eventdata, handles)
% hObject    handle to pb_aux_art (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global CELLS
    global cells_online
    global online_features
    global clustTreeEuc
    global online_profiles
    global t
    global spath
    global pclass
    global params
    global ca_prof_smooth
    global t_sample
    global n_frames
    global exp_data
    global th_prom
    global th_artifact
    global th_diff
th_artifact = str2double(get(handles.edit5,'String'));
setappdata(0,'ca_prof_smooth',ca_prof_smooth) 
setappdata(0,'t',t) 
setappdata(0,'th_prom',th_prom) 
setappdata(0,'th_diff',th_diff) 
setappdata(0,'th_artifact',th_artifact) 
AuxApp_Th_Art


% --- Executes on button press in chb_plot.
function chb_plot_Callback(hObject, eventdata, handles)
% hObject    handle to chb_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global en_plot
en_plot = get(handles.chb_plot,'Value')
% Hint: get(hObject,'Value') returns toggle state of chb_plot


% --- Executes on button press in cb_artifacts2.
function cb_artifacts2_Callback(hObject, eventdata, handles)
% hObject    handle to cb_artifacts2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global en_artifact_detector_2
en_artifact_detector_2 = get(handles.cb_artifacts2,'Value')
getOnlineCells(hObject, handles)
pb_classify_Callback(hObject, eventdata, handles)
slider1_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of cb_artifacts2


% --- Executes on button press in cb_ignore_low_resp.
function cb_ignore_low_resp_Callback(hObject, eventdata, handles)
% hObject    handle to cb_ignore_low_resp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global en_ignore_low_resp
en_ignore_low_resp = get(handles.cb_ignore_low_resp,'Value')
getOnlineCells(hObject, handles)
pb_classify_Callback(hObject, eventdata, handles)
slider1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of cb_ignore_low_resp


% --- Executes on button press in cb_show_clusters.
function cb_show_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to cb_show_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global en_show_clusters
en_show_clusters = get(handles.cb_show_clusters,'Value')
getOnlineCells(hObject, handles)
pb_classify_Callback(hObject, eventdata, handles)
slider1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of cb_show_clusters


% --- Executes on button press in cb_classlevel2.
function cb_classlevel2_Callback(hObject, eventdata, handles)
% hObject    handle to cb_classlevel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global en_classlevel2
global en_classlevel3
en_classlevel2 = get(handles.cb_classlevel2,'Value')
if ~en_classlevel2
    en_classlevel3 = 0;
    set(handles.cb_classlevel3,'Value',0)    
end
getOnlineCells(hObject, handles)
pb_classify_Callback(hObject, eventdata, handles)
slider1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of cb_classlevel2


% --- Executes on button press in cb_classlevel3.
function cb_classlevel3_Callback(hObject, eventdata, handles)
% hObject    handle to cb_classlevel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global en_classlevel3
global en_classlevel2
en_classlevel3 = get(handles.cb_classlevel3,'Value')
if en_classlevel3
    en_classlevel2 = 1;
    set(handles.cb_classlevel2,'Value',1)    
end
getOnlineCells(hObject, handles)
pb_classify_Callback(hObject, eventdata, handles)
slider1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of cb_classlevel3
