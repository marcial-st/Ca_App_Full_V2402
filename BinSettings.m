    function varargout = BinSettings(varargin)
% BINSETTINGS MATLAB code for BinSettings.fig
%      BINSETTINGS, by itself, creates a new BINSETTINGS or raises the existing
%      singleton*.
%
%      H = BINSETTINGS returns the handle to a new BINSETTINGS or the handle to
%      the existing singleton*.
%
%      BINSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BINSETTINGS.M with the given input arguments.
%
%      BINSETTINGS('Property','Value',...) creates a new BINSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BinSettings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BinSettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BinSettings

% Last Modified by GUIDE v2.5 28-May-2020 18:30:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BinSettings_OpeningFcn, ...
                   'gui_OutputFcn',  @BinSettings_OutputFcn, ...
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


% --- Executes just before BinSettings is made visible.
function BinSettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BinSettings (see VARARGIN)

% Choose default command line output for BinSettings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global img
global img380
global t
global w
global roiMask
global roiMask_i
global Irgb
global applyB_count
global segnet_medium
global segnet_rough

segnet_medium = load('SegNet_medium_v200527.mat');
segnet_rough = load('SegNet_rough_v200527.mat');

applyB_count=0;
t=5;
w=getappdata(0,'w');
img=getappdata(0,'im340'); %would be im340
img380=getappdata(0,'im380');
set(handles.tText,'string',num2str(t))
set(handles.wText,'string',num2str(w))
set(handles.autoradio,'visible','off');
set(handles.axes1,'Visible','on');
axes(handles.axes1)%pongo la imagen en la ventana uno
imshow(img)
axis off;
% Code
    [m,n]=size(img);
    roiMask=zeros(m,n);
    roiMask_i=roiMask;
    Irgb=uint8(zeros(m,n,3));
    Irgb(:,:,1)=img;
    Irgb(:,:,2)=img;
    Irgb(:,:,3)=img; %+(uint8(noiseMask)*255);
    axes(handles.axes1)%pongo la imagen en la ventana uno
    imshow(Irgb)
    set(handles.pb_open,'visible','off');
    set(handles.undoB,'visible','on');
    set(handles.undoB,'enable','off');
    set(handles.newroiB,'visible','on');
    set(handles.applyB,'visible','on');
    undoEnable = 0;
    msg = msgbox('Thanks! Now, please use the top toolbar & select some cells, then "Preview" and "Let´s Go" once ready.');
    uiwait(msg)
% Start Visible/Enable On/Off
% set(handles.pb_open,'visible','off');
% set(handles.applyB,'enable','off');
% set(handles.undoB,'enable','off');
% set(handles.newroiB,'enable','off');
% set(handles.upT,'enable','off');
% set(handles.tText,'enable','off');
% set(handles.newroiB,'enable','off');
% set(handles.letsgoB,'enable','off');
% UIWAIT makes BinSettings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BinSettings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function tText_Callback(hObject, eventdata, handles)
% hObject    handle to tText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of tText as text
%        str2double(get(hObject,'String')) returns contents of tText as a double
global img
global w
global t
t = str2double(get(handles.tText,'string'));
binImg =bin_and_clean(img,w);
axes(handles.axes1)%pongo la imagen en la ventana uno
imshow(binImg)

% --- Executes during object creation, after setting all properties.
function tText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usua<lly have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wText_Callback(hObject, eventdata, handles)
% hObject    handle to wText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of wText as text
%        str2double(get(hObject,'String')) returns contents of wText as a double
global img
global w
global t
wMax=min(size(img))/2;
if wMax ~= 0
    w = str2double(get(handles.wText,'string'));
    binImg =bin_and_clean(adapthisteq(img),w);
    [m,n]=size(img);
    imgrgb=uint8(zeros(m,n,3));
    imgrgb(:,:,1)=uint8(img); imgrgb(:,:,2)=imgrgb(:,:,1); imgrgb(:,:,3)=imgrgb(:,:,1)+uint8(binImg).*255;
    axes(handles.axes1)%pongo la imagen en la ventana uno
    imshow(imgrgb)
end
% 
% binImg =bin_and_clean(img,w,t);
% axes(handles.axes1)%pongo la imagen en la ventana uno
% imshow(binImg)


% --- Executes during object creation, after setting all properties.
function wText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in upT.
function upT_Callback(hObject, eventdata, handles)
% hObject    handle to upT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img
global w
global t
t = str2double(get(handles.tText,'string'));
if t>=100 
    t=100;
else
    t = t+1;
end
set(handles.tText,'string',num2str(t))
binImg =bin_and_clean(img,w);
axes(handles.axes1)%pongo la imagen en la ventana uno
imshow(binImg)

% --- Executes on button press in downT.
function downT_Callback(hObject, eventdata, handles)
% hObject    handle to downT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img
global w
global t
t = str2double(get(handles.tText,'string'));
if t==0
    t=0;
else
    t = t-1;
end
set(handles.tText,'string',num2str(t))
binImg =bin_and_clean(img,w);
axes(handles.axes1)%pongo la imagen en la ventana uno
imshow(binImg)

% --- Executes on button press in upW.
function upW_Callback(hObject, eventdata, handles)
% hObject    handle to upW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img
global w
global t
global wMax
wMax=min(size(img))/2;
if wMax ~= 0
    w = str2double(get(handles.wText,'string'));
    if w>=wMax 
        w=wMax;
    else
        w = w+1;
    end
    set(handles.wText,'string',num2str(w))
%     binImg =bin_and_clean(img,w,t);
%     axes(handles.axes1)%pongo la imagen en la ventana uno
%     imshow(binImg)
    binImg =bin_and_clean(adapthisteq(img),w);
    [m,n]=size(img);
    imgrgb=uint8(zeros(m,n,3));
    imgrgb(:,:,1)=uint8(img); imgrgb(:,:,2)=imgrgb(:,:,1); imgrgb(:,:,3)=imgrgb(:,:,1)+uint8(binImg).*255;
    axes(handles.axes1)%pongo la imagen en la ventana uno
    imshow(imgrgb)
else
    set(handles.wText,'string','NoImg')
end

% --- Executes on button press in downW.
function downW_Callback(hObject, eventdata, handles)
% hObject    handle to downW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img
global w
global t
wMin=1;
wMax=min(size(img))/2;
w = str2double(get(handles.wText,'string'));
if wMax ~= 0
    if w==wMin
        w=wMin;
    else
        w = w-1;
    end
    set(handles.wText,'string',num2str(w))
    binImg =bin_and_clean(adapthisteq(img),w);
    [m,n]=size(img);
    imgrgb=uint8(zeros(m,n,3));
    imgrgb(:,:,1)=uint8(img); imgrgb(:,:,2)=imgrgb(:,:,1); imgrgb(:,:,3)=imgrgb(:,:,1)+uint8(binImg).*255;
    axes(handles.axes1)%pongo la imagen en la ventana uno
    imshow(imgrgb)
else
    set(handles.wText,'string','NoImg')
end


% --- Executes on button press in autoradio.
function autoradio_Callback(hObject, eventdata, handles)
% hObject    handle to autoradio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoradio
global img
global noiseMask
global roiMask
global roiMask_i
global Irgb
global undoEnable
aEnable=get(handles.autoradio,'value');
if aEnable == 1
%     msg = msgbox('Need a noisy sample, please select a big one!');
%     uiwait(msg)
%     hnoise=imfreehand(handles.axes1);
%     noiseMask=hnoise.createMask();
    [m,n]=size(img);
    roiMask=zeros(m,n);
    roiMask_i=roiMask;
    Irgb=uint8(zeros(m,n,3));
    Irgb(:,:,1)=img;
    Irgb(:,:,2)=img;
    Irgb(:,:,3)=img; %+(uint8(noiseMask)*255);
    axes(handles.axes1)%pongo la imagen en la ventana uno
    imshow(Irgb)
    msg = msgbox('Thanks! Now, please use the top toolbar & select some cells, then push APPLY.');
    uiwait(msg)
    set(handles.undoB,'visible','on');
    set(handles.undoB,'enable','off');
    set(handles.newroiB,'visible','on');
    set(handles.applyB,'visible','on');
    undoEnable = 0;
else
    set(handles.undoB,'visible','off');
    set(handles.newroiB,'visible','off');
    set(handles.applyB,'visible','off');
    axes(handles.axes1)%pongo la imagen en la ventana uno
    imshow(img)
end

%     Ib=hnoise.createMask();


% --- Executes on button press in applyB.
function applyB_Callback(hObject, eventdata, handles)
% hObject    handle to applyB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img
global roiMask
global noiseMask
global undoEnable
global roiStats
global w
global t
global applyB_count


autoEnable = max(max(roiMask));
if autoEnable == 1
    applyB_count=applyB_count+1;
    if applyB_count==1
%     msg = msgbox('Place here the automatic tune super chido code XD');
%     uiwait(msg)
        roiStats  = regionprops('table',logical(roiMask),'Area');
        w = round(2*sqrt(mean(roiStats.Area)));
%     zNoise=sum(noiseMask(:) == 0);
%     noiseSample=uint8(noiseMask).*img;
%     [co,bl]=imhist(noiseSample);
%     co(1)=co(1)-zNoise;
%     meanNoise = co'*bl/sum(co);
%     sigmaNoise = sqrt(sum(co.*(bl-meanNoise).^2)/sum(co));
%     t=round(sigmaNoise);
        set(handles.wText,'string',num2str(w))
    else
        w=str2num(get(handles.wText,'string'));
    end
%     set(handles.tText,'string',num2str(t))
    binImg =bin_and_clean(adapthisteq(img),w);
    [m,n]=size(img);
    imgrgb=uint8(zeros(m,n,3));
    imgrgb(:,:,1)=uint8(img); imgrgb(:,:,2)=imgrgb(:,:,1); imgrgb(:,:,3)=imgrgb(:,:,1)+uint8(binImg).*255;
    axes(handles.axes1)%pongo la imagen en la ventana uno
    imshow(imgrgb)
    set(handles.letsgoB,'enable','on');
else
    msg = msgbox('Please select some cells, more cells, better tune');
    uiwait(msg)
end

% --------------------------------------------------------------------
function newroiB_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to newroiB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img
global roiMask
global roiMask_i
global Irgb
global undoEnable
undoEnable = 1;
set(handles.undoB,'enable','on');
hroi=imfreehand(handles.axes1);
roiMask_i = hroi.createMask();
roiMask = roiMask+roiMask_i;
Irgb(:,:,1)=img+(uint8(roiMask)*255);
axes(handles.axes1)%pongo la imagen en la ventana uno
imshow(Irgb)


% --------------------------------------------------------------------
function undoB_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to undoB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img
global roiMask
global roiMask_i
global Irgb
global undoEnable
if undoEnable==1
    roiMask = roiMask-roiMask_i;
    Irgb(:,:,1)=img+(uint8(roiMask)*255);
    axes(handles.axes1)%pongo la imagen en la ventana uno
end
imshow(Irgb)    


% --------------------------------------------------------------------
function letsgoB_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to letsgoB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global w
global im340 
global im380
global t
global w
global Base
global LBase
global bin_340
global fc
global k 
global ratio

if (w ~= 999)&&(w ~= 998)
    [bin_340,m,n,k]=bin_sequence(im340,w);
end

Base=bin_340(:,:,fc);
LO=bwconncomp(Base);
LBase=labelmatrix(LO);
ratio=(double(im340(:,:,1))./double(im380(:,:,1)));

setappdata(0,'w',w);
setappdata(0,'k',k);
setappdata(0,'bin_340',bin_340);
setappdata(0,'LBase',LBase);
setappdata(0,'ratio',ratio);


 close();


% --- Executes on button press in pb_segnet1.
function pb_segnet1_Callback(hObject, eventdata, handles)
% hObject    handle to pb_segnet1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global segnet_medium
global im340 
global im380
global w
global LBase
global ratio
global img
global bin_340


w = 999;

% img=imread('img0-003.tif');

C = semanticseg(img,segnet_medium.net);
bin_340 = C=='cell';
bin_340(1:10,:) = 0;
bin_340(:,1:10) = 0;
bin_340(end-10:end,:) = 0;
bin_340(:,end-10:end) = 0;
se = strel('disk',4);
bin_340 = imopen(bin_340,se);
bin_340 = imfill(bin_340,'holes');
% se = strel('disk',1);
% bin_340 = imdilate(bin_340,se);

[m,n]=size(img);
imgrgb=uint8(zeros(m,n,3));
imgrgb(:,:,1)=uint8(img); imgrgb(:,:,2)=imgrgb(:,:,1); imgrgb(:,:,3)=imgrgb(:,:,1)+uint8(bin_340).*255;
axes(handles.axes1)%pongo la imagen en la ventana uno
imshow(imgrgb)

LO=bwconncomp(bin_340);
LBase=labelmatrix(LO);
ratio=(double(im340(:,:,1))./double(im380(:,:,1)));

set(handles.letsgoB,'enable','on');

% --- Executes on button press in pb_segnet2.
function pb_segnet2_Callback(hObject, eventdata, handles)
% hObject    handle to pb_segnet2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global segnet_rough
global img 
global w
global LBase
global ratio
global im340 
global im380
global bin_340

w = 998;

% img=imread('img0-129.tif');
img=medfilt2(img,[7 7]);
img=adapthisteq(img);
img=adapthisteq(img);


C = semanticseg(img,segnet_rough.net);
bin_340 = C=='cell';

se = strel('disk',3);
bin_340 = imopen(bin_340,se);
bin_340 = imfill(bin_340,'holes');
se = strel('disk',1);
bin_340 = imdilate(bin_340,se);

[m,n]=size(img);
imgrgb=uint8(zeros(m,n,3));
imgrgb(:,:,1)=uint8(img); imgrgb(:,:,2)=imgrgb(:,:,1); imgrgb(:,:,3)=imgrgb(:,:,1)+uint8(bin_340).*255;
axes(handles.axes1)%pongo la imagen en la ventana uno
imshow(imgrgb)

LO=bwconncomp(bin_340);
LBase=labelmatrix(LO);
ratio=(double(im340(:,:,1))./double(im380(:,:,1)));

set(handles.letsgoB,'enable','on');
