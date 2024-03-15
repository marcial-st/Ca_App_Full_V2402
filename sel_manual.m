%% Thesis
%%% Function for automatic ROI selection
function [cell,nroi]=sel_manual(Base,LBase,nroi,fc,A,H,Q,R,ratio_frame)
%%% cell structure is created as follows
%%%    1    cell.user_label  
%%%    2    cell.track_label 
%%%    3    cell.bbox
%%%    4    cell.bbox_wh
%%%    5    cell.x_meas
%%%    6    cell.x_est
%%%    7    cell.ecov_matrix
%%%    8    cell.flag_lost
%%%    9    cell.ca_profile
%%%   10    cell.xy_hist
%%%   11    cell.d_hist
%% Cell Selection
[x,y]= ginput(2);  %% [x,y] cursor position, corners of the ROI 
x=floor(x);
y=floor(y);
dist=sqrt(x.^2+y.^2); % relative distances of chosen points
c_index=find(dist==min(dist)); % find the upper left corner index, [x_i,y_i]
w=abs(x(1)-x(2));  % width
h=abs(y(1)-y(2));  % heigth
cell.bbox=[x(c_index),y(c_index),w,h]; % bounding box [x y x_width y_width]
im_part=imcrop(Base,cell.bbox);  % Mask of ROI in binary Base image
L_part=imcrop(LBase,cell.bbox);  % Mask of ROI in binary labeled Base image
L=bwconncomp(im_part);           % Look for connected objects
LM=labelmatrix(L);         % Labeled Matrix
props=regionprops(LM,'Area','Centroid');   % Area extraction
areas=cell2mat({props.Area});
if isempty(areas)
   warndlg('No cell was found, please select again');
   cell=0;
else
    i_selcel=find(areas==max(areas)); % find main object, cell
    cent_selcel=floor(props(i_selcel).Centroid);
    cell.track_label=L_part(cent_selcel(2),cent_selcel(1));
%% Cell Labeling
    dLaux = strcat('A',num2str(nroi,'%03u'));      %   default label
    def={dLaux};
    dlg_title='Input';
    prompt='Please input label';
    num_lines=1;
    UL= inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(UL)                        %   if user doesn't give one
       cell=0;
       return
    else
         cell.user_label = UL;                           %  user label 
    end
%% Feature Extraction   
    props=regionprops(LBase==cell.track_label,'BoundingBox','Centroid','Perimeter','Area');%bounding box struct
    cent=floor(props.Centroid); % centroids vector, just from selected cells
    
    cell.bbox_wh=[w,h]; % saves the original dimension of bounding box
    cell.x_meas=[cent(1) cent(2) floor(props.Perimeter) floor(props.Area)]; % Features vector, measurements
    [cell.x_est,cell.ecov_matrix]=kalman_v1(A,H,Q,R,cell.x_meas(1:2)',zeros(6,1),eye(6,6),1);% [state vector,error covariance matrix]Kalman estimator
                                       
    cell.xy_est(fc,:)=cell.x_est;
    cell.xy_hist(fc,:)=[cent(1),cent(2)];                             % Initialization of cell's optic flux
    cell.d_hist(fc)=0;                                                % Initialization of cell's classification distance history 
    cell.reliability_counter=0;
    bbox_i=[cell.bbox(1) cell.bbox(2) cell.bbox_wh(1) cell.bbox_wh(2)];
    crop_ratio=imcrop(ratio_frame,bbox_i);
    cell.ca_profile(fc)=mean2(crop_ratio); 
    cell.status='ONLINE';                                                 % Status initialization
    cell.color=[randi([2000,9000])/10000,randi([2000,9000])/10000,randi([2000,9000])/10000];% Assign a color 
    cell.props = props;
    nroi=nroi+1;                                                      % Number of ROIs is increased
end