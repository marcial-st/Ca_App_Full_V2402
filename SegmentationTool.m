function varargout = SegmentationTool(varargin)
% SEGMENTATIONTOOL MATLAB code for SegmentationTool.fig
%      SEGMENTATIONTOOL, by itself, creates a new SEGMENTATIONTOOL or raises the existing
%      singleton*.
%
%      H = SEGMENTATIONTOOL returns the handle to a new SEGMENTATIONTOOL or the handle to
%      the existing singleton*.
%
%      SEGMENTATIONTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENTATIONTOOL.M with the given input arguments.
%
%      SEGMENTATIONTOOL('Property','Value',...) creates a new SEGMENTATIONTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SegmentationTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SegmentationTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SegmentationTool

% Last Modified by GUIDE v2.5 24-Jul-2020 00:21:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SegmentationTool_OpeningFcn, ...
                   'gui_OutputFcn',  @SegmentationTool_OutputFcn, ...
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


% --- Executes just before SegmentationTool is made visible.
function SegmentationTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SegmentationTool (see VARARGIN)

% Choose default command line output for SegmentationTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SegmentationTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SegmentationTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
global process_dir
global process_list
global segnet_medium
global segnet_rough
global current_img

segnet_medium = load('SegNet_medium_v200714.mat');
segnet_rough = load('SegNet_rough_v200527.mat');
process_dir = {'--- Enhancement ---',...
               'MedianFilter','Contrast Eq','Adaptive Contrast Eq',...
               '--- Binarization ---',...
               'Adaptative Bin','SegNet Medium','SegNet Rough','Kmeans',...
               '--- Post Bin ---',...
               'Erode','Dilate','OpenFilt','CloseFilt','AreaThreshold','FillHoles'};
set(handles.menu_process,'String',process_dir)
process_list = struct('name',{},'params',{});

D340=getappdata(0,'D340');
D380=getappdata(0,'D380');
frame_count = getappdata(0,'frame_count');
current_frame = 1;

set(handles.axes1, 'Visible','On');
current_img = imread(strcat(D340(current_frame).folder,'\',D340(current_frame).name));
imshow(current_img,'Parent', handles.axes1)

set(handles.listbox1, 'Visible','On');
set(handles.pb_apply, 'Visible','On');
set(handles.pb_add_process, 'Visible','On');
set(handles.menu_process, 'Visible','On');
set(handles.sliderImg, 'Visible','On');
set(handles.text_current_frame, 'Visible','On');
set(handles.uipanel1, 'Visible','On');
set(handles.sliderImg, 'Max',frame_count);
set(handles.sliderImg, 'Min',1);
set(handles.sliderImg, 'Value',1);
set(handles.sliderImg, 'SliderStep', [1/frame_count 0.1]);
set(handles.pb_savelist,'Enable','Off')

% set(handles.listbox1, 'Visible','Off');
% set(handles.pb_apply, 'Visible','Off');
% set(handles.pb_add_process, 'Visible','Off');
% set(handles.menu_process, 'Visible','Off');
% set(handles.sliderImg, 'Visible','Off');
% set(handles.text_current_frame, 'Visible','Off');




% --------------------------------------------------------------------
function pb_load_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global D340
global D380
global current_frame
global current_img
global frame_count

current_frame = 1;

[D340,D380,~,frame_count,~]=get_images;

set(handles.axes1, 'Visible','On');
current_img = imread(strcat(D340(current_frame).folder,'\',D340(current_frame).name));
imshow(current_img,'Parent', handles.axes1)

set(handles.listbox1, 'Visible','On');
set(handles.pb_apply, 'Visible','On');
set(handles.pb_add_process, 'Visible','On');
set(handles.menu_process, 'Visible','On');
set(handles.sliderImg, 'Visible','On');
set(handles.text_current_frame, 'Visible','On');
set(handles.uipanel1, 'Visible','On');
set(handles.sliderImg, 'Max',frame_count);
set(handles.sliderImg, 'Min',1);
set(handles.sliderImg, 'Value',1);
set(handles.sliderImg, 'SliderStep', [1/frame_count 0.1]);


% --- Executes on slider movement.
function sliderImg_Callback(hObject, eventdata, handles)
% hObject    handle to sliderImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global current_frame
global current_img
global frame_count
global D340

current_frame = round(get(handles.sliderImg,'Value'));
set(handles.text_current_frame,'String',strcat(num2str(current_frame),"/",num2str(frame_count)));
strcat(D340(current_frame).folder,'\',D340(current_frame).name)
current_img = imread(strcat(D340(current_frame).folder,'\',D340(current_frame).name));
% imshow(current_img,'Parent', handles.axes1)
processListExecute(handles)


% --- Executes during object creation, after setting all properties.
function sliderImg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in menu_process.
function menu_process_Callback(hObject, eventdata, handles)
% hObject    handle to menu_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_process contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_process
global process_new

process_new = get(handles.menu_process,'Value');


% --- Executes during object creation, after setting all properties.
function menu_process_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_add_process.
function pb_add_process_Callback(hObject, eventdata, handles)
% hObject    handle to pb_add_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global process_dir
global process_list
global process_new
global idx

process_new = get(handles.menu_process,'Value');
% process_dir = {'--- Enhancement ---',...
%                'MedianFilter','Contrast Eq','Adaptive Contrast Eq',...
%                '--- Binarization ---',...
%                'Adaptative Bin','SegNet Medium','SegNet Rough',...
%                '--- Post Bin ---',...
%                'Erode','Dilate','OpenFilt','CloseFilt','FillHoles','AreaThreshold'};
selection = string(process_dir(process_new));
switch selection
    case '--- Enhancement ---'
        update_en = 0;
    case '--- Binarization ---'
        update_en = 0;
    case '--- Post Bin ---'
        update_en = 0;
    case 'MedianFilter'
        params = processDefMedianFilter(handles);
        update_en = 1;
    case 'Adaptative Bin'
        params = processDefAdaptFilter(handles);
        update_en = 1;
    case 'Kmeans'
        params = processDefKmeans(handles);
        update_en = 1;
    case 'Erode'
        params = processDefStrel(handles);
        update_en = 1;
    case 'Dilate'
        params = processDefStrel(handles);
        update_en = 1;
    case 'OpenFilt'
        params = processDefStrel(handles);
        update_en = 1;
    case 'CloseFilt'
        params = processDefStrel(handles);
        update_en = 1;
    otherwise
        params = [];
        update_en = 1;
end

if update_en == 1
    idx = max(size(process_list));
    process_list(idx+1).name = process_dir(process_new);
    process_list(idx+1).params = params;
    processListUpdate(handles)
    processListExecute(handles)
end
set(handles.pb_savelist,'Enable','On')



% text_m = cell(idx+1,1);
% for i=1:1:idx+1
%     text_m(i) = strcat(process_list(i).name,{' '},num2str(process_list(i).params));
% end
% set(handles.listbox1,'string',text_m,'Fontsize',8)


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
% global idx
% idx = get(handles.listbox1,'Value');

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


% --- Executes on button press in pb_apply.
function pb_apply_Callback(hObject, eventdata, handles)
% hObject    handle to pb_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global process_list
global current_img
global current_img_p
global D340
global D380

current_img = imread(strcat(D340(1).folder,'\',D340(1).name));
processListExecute(handles)
LO=bwconncomp(current_img_p);
LBase=labelmatrix(LO);

im340=imread(strcat(D340(1).folder,'\',D340(1).name));
im380=imread(strcat(D380(1).folder,'\',D380(1).name));
ratio=(double(im340)./double(im380));

setappdata(0,'process_list',process_list);   
setappdata(0,'bin_340',current_img_p);
setappdata(0,'LBase',LBase);
setappdata(0,'ratio',ratio);

close();

function params = processDefMedianFilter(~)
prompt = {'Rows:','Cols:'};
dlgtitle = 'Filter Size (pixels)';
dims = [1 35];
definput = {'7','7'};
params = str2double(inputdlg(prompt,dlgtitle,dims,definput))';

function params = processDefAdaptFilter(~)
prompt = {'W:'};
dlgtitle = 'Window Size (pixels)';
dims = [1 35];
definput = {'40'};
params = str2double(inputdlg(prompt,dlgtitle,dims,definput))';

function params = processDefStrel(~)
prompt = {'S:'};
dlgtitle = 'Strel-Disk Size (pixels)';
dims = [1 35];
definput = {'2'};
params = str2double(inputdlg(prompt,dlgtitle,dims,definput))';

function params = processDefKmeans(~)
prompt = {'K:'};
dlgtitle = 'K classes';
dims = [1 35];
definput = {'7'};
params = str2double(inputdlg(prompt,dlgtitle,dims,definput))';

function params = processDefAreaFilter(~)
prompt = {'Amin:','Amax:'};
dlgtitle = 'Area Thresholds (pixels)';
dims = [1 35];
definput = {'100','400'};
params = str2double(inputdlg(prompt,dlgtitle,dims,definput))';

function processListUpdate(handles)
global process_list
idx = max(size(process_list));
text_m = cell(idx,1);
if idx ~= 0
    for i=1:1:idx
        text_m(i) = strcat(process_list(i).name,{' '},num2str(process_list(i).params));
    end
end
set(handles.listbox1,'string',text_m,'Fontsize',8,'Value',idx)


function processListExecute(handles)
global process_list
global current_img
global current_img_p
global segnet_medium
global segnet_rough
% process_dir = {'--- Enhancement ---',...
%                'MedianFilter','Contrast Eq','Adaptive Contrast Eq',...
%                '--- Binarization ---',...
%                'Adaptative Bin','SegNet Medium','SegNet Rough',...
%                '--- Post Bin ---',...
%                'Erode','Dilate','OpenFilt','CloseFilt','AreaThreshold'};
current_img_p = current_img;
idx = max(size(process_list));
if idx ~= 0    
    for i =1:1:idx
        selection = string(process_list(i).name);
        switch selection
            case 'MedianFilter'
                current_img_p = medfilt2(current_img_p,process_list(i).params);
                disp(selection)
            case 'Contrast Eq'
                current_img_p = histeq(current_img_p);
                disp(selection)
            case 'Adaptive Contrast Eq'
                current_img_p = adapthisteq(current_img_p);
                disp(selection)
            case 'Adaptative Bin'
                w = process_list(i).params;
                padding='replicate';
                img=double(current_img_p);
                imgs=img.*img;
                meani = averagefilter(img, [w w], padding);
                meani2 = averagefilter(imgs, [w w], padding);
                va2=meani2-meani.*meani;
                stdi2=sqrt(va2);
                current_img_p= true(size(img));
                current_img_p(img < (meani+stdi2)) = 0;
                logical(current_img_p);
                disp(selection)
            case 'SegNet Medium'
                C = semanticseg(current_img_p,segnet_medium.net);
                current_img_p = C=='cell';
                current_img_p(1:10,:) = 0;
                current_img_p(:,1:10) = 0;
                current_img_p(end-10:end,:) = 0;
                current_img_p(:,end-10:end) = 0;
                logical(current_img_p);                             
                disp(selection)
            case 'SegNet Rough'
                C = semanticseg(current_img_p,segnet_rough.net);
                current_img_p = C=='cell';
                current_img_p(1:10,:) = 0;
                current_img_p(:,1:10) = 0;
                current_img_p(end-10:end,:) = 0;
                current_img_p(:,end-10:end) = 0;
                logical(current_img_p);
                disp(selection)
            case 'Kmeans'
                
                chLevel=0:255;
%                 intervals=7.0;
                intervals = process_list(i).params;
                winSz=9;
                nlevels=3;
    
                nrows=size(current_img_p,1);
                ncols=size(current_img_p,2);
%                 bands=[0,floor(ncols*0.10), floor(ncols*0.825), ncols];
%                 b0=bands(2)+1;
%                 b1=bands(3);
%                 imRes=uint8(zeros(nrows,b1-b0+1,3));
                imRes=uint8(zeros(nrows,ncols,3));
    
                his=imhist(current_img_p(:,:));
%                 his=imhist(current_img_p(:,b0:b1));
                [vMx,idmx]=max(his);
                hlast=find(his(idmx:end) ==0,1,'first')+idmx-1;
                delta= floor(double(hlast-idmx)/intervals);
                if delta <1
                    delta = 1;
                end
                dIntensity=floor(250.0/intervals);
                centroids(1,1)=10;
                for ii=2:intervals+4
                    centroids(ii,1)=idmx + (ii-6)*delta;
                end
                ab = reshape(double(current_img_p),nrows*ncols,1);
%                 ab = reshape(double(current_img_p(:,b0:b1)),nrows*(b1-b0+1),1);

                nClus=length(centroids);
                [cluster_idx, C]=kmeans(ab,nClus, 'Distance','cityblock', ...
                    'Start',centroids);
                [Cs, idC]=sort(C);
%                 pixel_labels = reshape(cluster_idx,nrows,(b1-b0+1));
                pixel_labels = reshape(cluster_idx,nrows,ncols);
                rFac=floor(255/nClus);
                for k=nClus-(nlevels-1):nClus
                    imRes(:,:,1)=imRes(:,:,1) + uint8((pixel_labels==idC(k)))*k*rFac;
                end
                gFac=floor(255/(nClus-nlevels));
                for k=5:nClus-(nlevels)
                    imRes(:,:,2)=imRes(:,:,2) + uint8((pixel_labels==idC(k)))*k*gFac;
                end
                bFac=floor(255/4);
                for k=2:4
                    imRes(:,:,3)=imRes(:,:,3) + uint8((pixel_labels==idC(k)))*k*bFac;
                end
                current_img_p = imRes(:,:,1);
                disp(selection)
            case 'Erode'
                se = strel('disk',process_list(i).params);
                current_img_p = imerode(current_img_p,se);
                disp(selection)
            case 'Dilate'
                se = strel('disk',process_list(i).params);
                current_img_p = imdilate(current_img_p,se);
                disp(selection)
            case 'OpenFilt'
                se = strel('disk',process_list(i).params);
                current_img_p = imopen(current_img_p,se);
                disp(selection)
            case 'CloseFilt'
                se = strel('disk',process_list(i).params);
                current_img_p = imclose(current_img_p,se);
                disp(selection)
            case 'FillHoles'
                current_img_p = imfill(current_img_p,'holes');
                disp(selection)
            case 'AreaThreshold'
                if islogical(current_img_p)
                    if isempty(process_list(i).params)
                        stats = regionprops(current_img_p,'Area');
                        figure;histogram([stats.Area]);title('Area histogram');xlabel('Area');ylabel('Frequency');
                        range = processDefAreaFilter;
                        process_list(i).params = range;                    
                    end                    
                    current_img_p = bwareafilt(current_img_p,process_list(i).params);
                    processListUpdate(handles)
                else
                    warndlg('Area filter only for binary images','Warning');
                end
                disp(selection)
            otherwise
                disp('holis')
        end
    end
    if islogical(current_img_p)
        current_img_p(1:10,:) = 0;
        current_img_p(:,1:10) = 0;
        current_img_p(end-10:end,:) = 0;
        current_img_p(:,end-10:end) = 0;
        [m,n]=size(current_img_p);
        imgrgb=uint8(zeros(m,n,3));
        imgrgb(:,:,1)=uint8(current_img); imgrgb(:,:,2)=imgrgb(:,:,1); imgrgb(:,:,3)=imgrgb(:,:,1)+uint8(current_img_p).*255;
        axes(handles.axes1)%pongo la imagen en la ventana uno
        imshow(imgrgb) 
        imshow(imgrgb,'Parent', handles.axes1)
    else
        imshow(current_img_p,'Parent', handles.axes1)
    end    
else
    imshow(current_img,'Parent', handles.axes1)
end


% --- Executes on key press with focus on listbox1 and none of its controls.
function listbox1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
global process_list

process_bye = get(handles.listbox1,'Value');
switch eventdata.Key
    case 'delete'
        process_list(process_bye) = [];
    case 'backspace'
        process_list(process_bye) = [];
end
processListUpdate(handles)
processListExecute(handles)


% --------------------------------------------------------------------
function pb_loadlist_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_loadlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global process_list

[filename, pathname] = uigetfile('*.mat','Select process list file');
if filename~=0
    pl = load([pathname,filename]);
    process_list = pl.process_list;
    processListUpdate(handles)
    processListExecute(handles)
end


% --------------------------------------------------------------------
function pb_savelist_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pb_savelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global process_list

[FileName,PathName] = uiputfile(strcat('processList_v',datestr(now,'yymmdd'),'.mat'),'Save process list');
if FileName~=0
    process_list
    save([PathName,FileName],'process_list')
    h = msgbox('File succesfully saved !')
end