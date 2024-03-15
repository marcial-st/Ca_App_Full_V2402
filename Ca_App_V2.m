function varargout = Ca_App_V2(varargin)
% CA_APP_V2 MATLAB code for Ca_App_V2.fig
%      CA_APP_V2, by itself, creates a new CA_APP_V2 or raises the existing
%      singleton*.
%
%      H = CA_APP_V2 returns the handle to a new CA_APP_V2 or the handle to
%      the existing singleton*.
%
%      CA_APP_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CA_APP_V2.M with the given input arguments.
%
%      CA_APP_V2('Property','Value',...) creates a new CA_APP_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ca_App_V2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ca_App_V2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ca_App_V2

% Last Modified by GUIDE v2.5 22-Jul-2020 10:55:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ca_App_V2_OpeningFcn, ...
                   'gui_OutputFcn',  @Ca_App_V2_OutputFcn, ...
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


% --- Executes just before Ca_App_V2 is made visible.
function Ca_App_V2_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ca_App_V2 (see VARARGIN)

% Choose default command line output for Ca_App_V2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global tb_auxEdgesCtrl;
global tb_manualSelectionCtrl;
global A;       % System matrix, Kalman estimator
global H;       % Observation matrix, Kalman estimator
global Q;       % System randomness, Kalman estimator
global R;       % Observation error covariance matrix, Kalman estimator
global nroi;

nroi=1;
tb_auxEdgesCtrl=0;
tb_manualSelectionCtrl=0;
set(handles.pb_newExp,'Enable','off')
set(handles.pb_calibrate,'Enable','off')
set(handles.tb_manualSelection,'Enable','off')
set(handles.tb_auxEdges,'Enable','off')
set(handles.pb_zoomOut,'Enable','off')
set(handles.pb_zoomIn,'Enable','off')
set(handles.pb_roi,'Enable','off')
set(handles.pb_allRois,'Enable','off')
set(handles.pb_start,'Enable','off')
set(handles.pb_save,'Enable','off')
logo=imread('logo_new.jpg');
imshow(logo,'Parent', handles.ax_im)
set(handles.ax_prof,'Visible','off')
sx=0.001;sy=0.001;ra=0.0269;rb=0.0993;T=1;
[A,H,Q,R]=k_eval(T,sx,sy,ra,rb);
% UIWAIT makes Ca_App_V2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Ca_App_V2_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in lb_status.
function lb_status_Callback(hObject, eventdata, handles)
% hObject    handle to lb_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_status contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_status


% --- Executes during object creation, after setting all properties.
function lb_status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function pb_loadSeries_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_loadSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im340;
global im380;
global D340;
global D380;
global numCycles;
global frameCounter;
global fc;
global fpc;

fc=1;
[D340,D380,numCycles,frameCounter,fpc]=get_images;
idcs   = strfind(D340(1).folder,filesep);
set(handles.st_path,'String',D340(1).folder(1:idcs(end)))
im340=imread(strcat(D340(1).folder,'\',D340(1).name));
im380=imread(strcat(D380(1).folder,'\',D380(1).name));
imshow(im340(:,:,1),'Parent', handles.ax_im)
if im340~=0
    set(handles.ax_prof,'Visible','on')
    hold(handles.ax_prof)
    grid(handles.ax_prof)
    xlim(handles.ax_prof,[1 frameCounter])
    xlabel(handles.ax_prof,'Frame','FontSize',9)
    ylabel(handles.ax_prof,'Ratio','FontSize',9)
    set(handles.pb_calibrate,'Enable','on')
    set(handles.pb_newExp,'Enable','on')
    set(handles.menu_tools,'Enable','on')
else
    logo=imread('logo_new.jpg');
    imshow(logo,'Parent', handles.ax_im)
    set(handles.ax_prof,'Visible','off')
end


% --------------------------------------------------------------------
function pb_newExp_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_newExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global
global A;       % System matrix, Kalman estimator
global H;       % Observation matrix, Kalman estimator
global Q;       % System randomness, Kalman estimator
global R;       % Observation error covariance matrix, Kalman estimator
global nroi;
global tb_auxEdgesCtrl;
global tb_manualSelectionCtrl;

set(handles.pb_newExp,'Enable','on')
set(handles.pb_loadSeries,'Enable','on')
set(handles.pb_calibrate,'Enable','off')
set(handles.tb_manualSelection,'Enable','off')
set(handles.tb_auxEdges,'Enable','off')
set(handles.pb_zoomOut,'Enable','off')
set(handles.pb_zoomIn,'Enable','off')
set(handles.pb_roi,'Enable','off')
set(handles.pb_allRois,'Enable','off')
set(handles.pb_start,'Enable','off')
set(handles.pb_save,'Enable','off')
set(handles.lb_status,'string','','Fontsize',11)
set(handles.edit1,'string','','Fontsize',11)
set(handles.textw,'string','','Fontsize',11)
cla(handles.ax_prof,'reset')
set(handles.ax_prof,'Visible','off')
logo=imread('logo_new.jpg');
imshow(logo,'Parent', handles.ax_im)
set(handles.st_path,'String','')
nroi=1;
tb_auxEdgesCtrl=0;
tb_manualSelectionCtrl=0;
sx=0.001;sy=0.001;ra=0.0269;rb=0.0993;T=1;
[A,H,Q,R]=k_eval(T,sx,sy,ra,rb);

% --------------------------------------------------------------------
function pb_calibrate_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global D340;
global D380;
global im340;
global im380;
global frameCounter;
global w;
setappdata(0,'w',100);
setappdata(0,'im340',im340(:,:,1));
setappdata(0,'im380',im380(:,:,1));
setappdata(0,'D340',D340);
setappdata(0,'D380',D380);
setappdata(0,'frame_count',frameCounter);
SegmentationTool
% setappdata(0,'fc',frameCounter);
% BinSettings
set(handles.pb_newExp,'Enable','on')
set(handles.pb_calibrate,'Enable','on')
set(handles.tb_manualSelection,'Enable','on')
set(handles.tb_auxEdges,'Enable','on')
set(handles.pb_zoomOut,'Enable','on')
set(handles.pb_zoomIn,'Enable','on')
set(handles.pb_roi,'Enable','on')
set(handles.pb_allRois,'Enable','on')
set(handles.pb_start,'Enable','off')
set(handles.pb_save,'Enable','off')


% --------------------------------------------------------------------
function tb_auxEdges_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_auxEdges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bin_340;
global n_ratio; 
global frameCounter; 
global nroi; 
global CELLS; 
global eb;
global fc;
global tb_auxEdgesCtrl;
global w;
global im340;
global process_list
% global nroi;

% [bin_340,~,~,~,~]=bin_sequence(im340,w);
bin_340=getappdata(0,'bin_340');
process_list = getappdata(0,'process_list');
process_list_enh = processListOnlyEnhancement(process_list);
im340_enh = processListApply(process_list_enh,im340(:,:,fc));
tb_auxEdgesCtrl=tb_auxEdgesCtrl+1;

if(mod(tb_auxEdgesCtrl,2)==0)
%     imshow(im340(:,:,fc), 'Parent', handles.ax_im)
    imshow(im340_enh, 'Parent', handles.ax_im)
else
    % visualization image in axes_cells
    eb=edge(bin_340(:,:,fc),'canny'); % edges from first frame b340
    eb=imdilate(eb,strel('disk', 1, 4));
    mask1=im340_enh.*uint8(not(eb));
%     mask1=im340(:,:,fc).*uint8(not(eb));
    im0(:,:,1)=mask1+uint8(eb)*0; % mask R
    im0(:,:,2)=mask1+uint8(eb)*0; % mask G
    im0(:,:,3)=mask1+uint8(eb)*0; % mask B
    imshow(im0, 'Parent', handles.ax_im)
end
	if (isempty(CELLS)==0)
        for ip=1:1:nroi-1
            rectangle('Position',CELLS(ip).bbox,'EdgeColor',CELLS(ip).color,'LineWidth',1) % draw ROI, bounding box
            text(CELLS(ip).bbox(1),CELLS(ip).bbox(2)-6,CELLS(ip).user_label,'Color',CELLS(ip).color) % draw user defined labe
        end
	end
% if nroi>1
%   for i=1:1:nroi-1
%      rectangle('Position',CELLS(i).bbox,'EdgeColor',CELLS(i).color,'LineWidth',1.5) % draw ROI, bounding box
%      text(CELLS(i).bbox(1),CELLS(i).bbox(2)-6,CELLS(i).user_label,'Color',CELLS(i).color) % draw user defined labels
%   end
% end

function tb_auxEdges_OnCallback(hObject, eventdata, handles)
% ---------------------------------------------------------
function tb_auxEdges_OffCallback(hObject, eventdata, handles)
% ---------------------------------------------------------


% --------------------------------------------------------------------
function tb_manualSelection_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tb_manualSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tb_manualSelectionCtrl;
tb_manualSelectionCtrl=tb_manualSelectionCtrl+1;
if(mod(tb_manualSelectionCtrl,2)==0)
    set(handles.pb_allRois,'Enable','on')
else
    set(handles.pb_allRois,'Enable','off')
end


% --------------------------------------------------------------------
function pb_roi_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global LBase;
global Base;
global nroi;
global CELLS;
global fc;
global A;
global H;
global Q;
global R;
global tb_manualSelectionCtrl;
global bin_340;
global ratio;
global process_list;

if tb_manualSelectionCtrl==0
    process_list = getappdata(0,'process_list');
%     w=getappdata(0,'w');
%     k=getappdata(0,'k');
    bin_340=getappdata(0,'bin_340');
    LBase=getappdata(0,'LBase');
	ratio=getappdata(0,'ratio');
%     set(handles.textw,'string',num2str(w),'Fontsize',11) 
end
set(handles.pb_start,'Enable','on') 
set(handles.pb_allRois,'Enable','off') 
%% Cell selection mode functions
if(mod(tb_manualSelectionCtrl,2)==0)
    [cell,nroi]=sel_auto(LBase,nroi,fc,A,H,Q,R,ratio(:,:,fc));
else
   [cell,nroi]=sel_manual(Base,LBase,nroi,fc,A,H,Q,R,ratio(:,:,fc));
end
%% Plot ROIs and Labels
if isstruct(cell)
   rectangle('Position',cell.bbox,'EdgeColor',cell.color,'LineWidth',1.5) % draw ROI, bounding box
   text(cell.bbox(1),cell.bbox(2)-6,cell.user_label,'Color',cell.color) % draw user defined labels
   %% Build CELLS structure
%    if nroi<=2
   if nroi<2       
       CELLS=cell;
   else
       CELLS=[CELLS cell];
   end
end


% --------------------------------------------------------------------
function pb_allRois_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_allRois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LBase;
global Base;
global nroi;
global CELLS;
global fc;
global A;
global H;
global Q;
global R;
global ratio;

% Cell selection mode functions
LBase=getappdata(0,'LBase');
ratio=getappdata(0,'ratio');
n_cells=nonzeros(unique(LBase));

for nroi=1:1:length(n_cells)
   [cell]=sel_auto_all(LBase,nroi,fc,A,H,Q,R,ratio(:,:,fc));
    %% Plot ROIs and Labels
    if isstruct(cell)
        rectangle('Position',cell.bbox,'EdgeColor',cell.color,'LineWidth',1) % draw ROI, bounding box
        text(cell.bbox(1),cell.bbox(2)-6,cell.user_label,'Color',cell.color) % draw user defined labels
%         set(handles.lb_status,'string',c_t,'Fontsize',12)
        %% Build CELLS structure
%         if nroi<=2
        if nroi<2            
            CELLS=cell;
        else
            CELLS=[CELLS cell];
        end
    end
end    
set(handles.pb_start,'Enable','on')


% --------------------------------------------------------------------
function pb_start_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fc;
global frameCounter;
global D340;
global D380;
global fpc;
global numCycles;
global CELLS;
global bin_340;
global im340;
global LBase;
% global im380;
global ratio;
global w;
global nroi;
global A;
global H;
global Q;
global R;
global masterCount;
global segnet_medium
global segnet_rough
global process_list
global reliability_limit
global entropy_frame
global process_time

segnet_medium = load('SegNet_medium_v200714.mat');
segnet_rough = load('SegNet_rough_v200527.mat');
process_list_enh = processListOnlyEnhancement(process_list);

   def_limit ={'3'};      %   default label
   dlg_title='RLimit';
   prompt='Please input the reliability limit';
   num_lines=1;
   RL= inputdlg(prompt,dlg_title,num_lines,def_limit);
   if isempty(RL)                        %   if user doesn't give one
        reliability_limit = str2num(cell2mat(def_limit));   %   assign default limit
   else
        reliability_limit = str2num(cell2mat(RL));           %  user limit
   end
  set(handles.edit1,'string',RL,'Fontsize',11)
  set(handles.pb_loadSeries,'Enable','off')
  set(handles.pb_calibrate,'Enable','off')
  set(handles.tb_manualSelection,'Enable','off')
  set(handles.tb_auxEdges,'Enable','off')
  set(handles.pb_zoomOut,'Enable','off')
  set(handles.pb_zoomIn,'Enable','off')
  set(handles.pb_roi,'Enable','off')
  set(handles.pb_allRois,'Enable','off')
  imDim=size(bin_340);
  packIndex=zeros(numCycles+1,1);
  packIndex(2:numCycles+1)=fpc;
  if(mod(frameCounter,fpc)~=0)
      packIndex(numCycles+1)=mod(frameCounter,fpc);
  end
  lbase_props = regionprops(LBase,'MajoraxisLength');
  vecinity = floor(mean([lbase_props.MajorAxisLength]));
  im2s = zeros(imDim(1),imDim(2),3);
  start_time = tic; 
  masterCount=2;
  for i1=1:1:numCycles
% for i1=1:1:1
      limInf=sum(packIndex(1:i1))+1;
      limSup=sum(packIndex(1:i1+1));
%       [im340,bin_340,ratio]=img_subset(D340,D380,limInf,limSup,imDim,w,segnet_medium,segnet_rough);
      [im340,bin_340,ratio] = img_subset_plist(D340,D380,limInf,limSup,imDim,process_list,segnet_medium,segnet_rough);
      uint8(im340);
%***************************************       
%     Insert Here the Analysis Kernel 
%   Begin Kernel
      if i1==1 % Initial conditions
          i2_start=2;
          i2_end=packIndex(i1+1);
          cc1=bwconncomp(bin_340(:,:,2));
          LM_fc1=labelmatrix(cc1);          
          LM_fc2=LBase;
      else    % Next blocks
          i2_start=1;
          i2_end=packIndex(i1+1);
          cc1=bwconncomp(bin_340(:,:,1));
          LM_fc1=labelmatrix(cc1);
      end
      
      for i2=i2_start:1:i2_end % Cycle per Img Subset
          entropy_frame((i1-1)*100+i2) = entropy(im340(:,:,i2));
          im340_i_enh = processListApply(process_list_enh,im340(:,:,i2));
          im2s(:,:,1)=im340_i_enh;
          im2s(:,:,2)=im340_i_enh;
          im2s(:,:,3)=im340_i_enh+128*uint8(bin_340(:,:,i2));
          imshow(uint8(im2s), 'Parent', handles.ax_im) % Update image to show
          
          for i3=1:1:nroi-1 % Cycle per ROI
              if reliability_limit == CELLS(i3).reliability_counter
                  CELLS(i3).status='OFFLINE';
              else
                  kernel_ps=and(LM_fc2==CELLS(i3).track_label,LM_fc1);  %pseudo kernel
                  cck=bwconncomp(kernel_ps);
                  n_cck=cck.NumObjects;
                  if n_cck == 1
                      [new_t_l]=int_tl(cck,LM_fc1);
                      stats_new = regionprops('table',LM_fc1==new_t_l,'Centroid');
                      stats_old = regionprops('table',LM_fc2==CELLS(i3).track_label,'Centroid');
                      dist = pdist([stats_new.Centroid;stats_old.Centroid]);
                      if dist > 1*vecinity
                          [new_t_l]=9999;                          
                      end                      
%                   elseif n_cck > 1
%                       [kernel_sub] = kernel_vecinity(imDim,vecinity,CELLS(i3).bbox,LM_fc1);
%                       stats_old = struct2array(regionprops(LM_fc2==CELLS(i3).track_label,'Centroid','Area','Perimeter','BoundingBox','MajoraxisLength','MinoraxisLength','Orientation'));
%                       stats_new = reshape(struct2array(regionprops(kernel_sub,'Centroid','Area','Perimeter','BoundingBox','MajoraxisLength','MinoraxisLength','Orientation')),11,[])';
%                       [new_t_l] = knnsearch(stats_new,stats_old);
                  else
                      [new_t_l]=9999;                          
                  end
%                       [x_e_aux,~]=kalman_v1(A,H,Q,R,CELLS(i3).x_meas(1:2)',CELLS(i3).x_est,CELLS(i3).ecov_matrix,1);
%                       [new_t_l]=LM_fc1(floor(x_e_aux(4)),floor(x_e_aux(1)));
                  [cell_i,bbox_i]=store_data(CELLS(i3),new_t_l,masterCount,LM_fc1,A,H,Q,R,ratio(:,:,i2));
                  rectangle('Position',bbox_i,'EdgeColor',cell_i.color,'LineWidth',1.5) % draw ROI, bounding box
                  text(cell_i.bbox(1),cell_i.bbox(2)-6,cell_i.user_label,'Color',cell_i.color) % draw user defined labels
                  plot(handles.ax_prof,cell_i.ca_profile,'Color',cell_i.color,'LineWidth',2)
                  CELLS(i3)=cell_i;
%                    pause
              end
          end
          pause(0.01)
          if i2 < i2_end
              cc1=bwconncomp(bin_340(:,:,i2+1));     % Labeled Frames Update 
              LM_fc1=labelmatrix(cc1);
          end
          cc2=bwconncomp(bin_340(:,:,i2));
          LM_fc2=labelmatrix(cc2);
          
          c_t=cell_text(CELLS);    % Update status list
          set(handles.lb_status,'string',c_t,'Fontsize',11)
          masterCount=masterCount+1;
      end
%       c_t=cell_text(CELLS);
%       set(handles.list_cell_status,'string',c_t,'Fontsize',12)
%       count=0;
%   End Kernel      
%***************************************
  end
 [CELLS]=cellsCleanOverlap(CELLS);
 process_time = toc(start_time); 
 msgbox('Analysis completed!');
 set(handles.pb_save,'Enable','on')
 set(handles.pb_start,'Enable','off')
 set(handles.pb_class,'Enable','on')
 
% for fc=2:1:frameCounter;
%   c_t=cell_text(CELLS);
%   set(handles.lb_status,'string',c_t,'Fontsize',11)
% end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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
function pb_save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CELLS
global DATA_table
global masterCount
global D340
global reliability_limit
global process_list
global entropy_frame
global process_time

javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');

% strcat(D340(1).folder,'\',D340(1).name)
% [DATA_table]=save_data(CELLS,masterCount);
idcs   = strfind(D340(1).folder,filesep);
folder = D340(1).folder(1:idcs(end))
file = ['CELLS_RL',num2str(reliability_limit),'_',datestr(now,'yymmdd_hhMM'),'.mat'];
[FileName,PathName] = uiputfile([folder,file],'Save Data')
save([PathName,FileName],'CELLS','process_list','reliability_limit','entropy_frame','process_time')

% xlwrite(fileName, xlsData, sheetName, startRange);
% xlswrite(strcat(PathName,FileName),DATA_table); 
% xlwrite(strcat(PathName,FileName),DATA_table); 
% writecell(DATA_table,strcat(PathName,FileName));

h = msgbox('File succesfully saved !')


% --------------------------------------------------------------------
function pb_class_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CELLS
global im340
global process_list

setappdata(0,'CELLS',CELLS);
setappdata(0,'CELLS',CELLS);
setappdata(0,'im340_1',im340(:,:,1));
setappdata(0,'im340_n',im340(:,:,end));
setappdata(0,'process_list',process_list);

Classifier_v200429


% --------------------------------------------------------------------
function menu_tools_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_p2m_Callback(hObject, eventdata, handles)
% hObject    handle to menu_p2m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im340

setappdata(0,'im340_1',im340(:,:,1));
Calibration_p2m


% --------------------------------------------------------------------
function importRois_Callback(hObject, eventdata, handles)
% hObject    handle to importRois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LBase;
global Base;
global nroi;
global CELLS;
global fc;
global A;
global H;
global Q;
global R;
global ratio;

% Cell selection mode functions
LBase=getappdata(0,'LBase');
ratio=getappdata(0,'ratio');
n_cells=nonzeros(unique(LBase));

[roifilename, roipathname] = uigetfile({'*.xls;*.xlsx','Spreadsheets files (*.xls,*.xlsx)';...                                  
                                  '*.*',  'All Files (*.*)'},...
                                  'Select ROIs spreadsheet');

if roifilename ~= 0
    
    T = readtable([roipathname,roifilename]);
    pos = table2array(T);
    [n_pos,~] = size(pos);       
    pos_updated = updateROIS(pos,LBase);
    
    axes(handles.ax_im);
    for nroi=1:1:n_pos
       [cell]=sel_auto_import(LBase,nroi,fc,A,H,Q,R,ratio(:,:,fc),pos_updated);
        %% Plot ROIs and Labels
        if isstruct(cell)            
            rectangle('Position',cell.bbox,'EdgeColor',cell.color,'LineWidth',1) % draw ROI, bounding box
            text(cell.bbox(1),cell.bbox(2)-6,cell.user_label,'Color',cell.color) % draw user defined labels
    %         set(handles.lb_status,'string',c_t,'Fontsize',12)
            %% Build CELLS structure
%             if nroi<=2            
            if nroi<2
                CELLS=cell;
            else
                CELLS=[CELLS cell];
            end
        end
    end
    nroi = max(size(CELLS));
    set(handles.pb_start,'Enable','on')    
end
