%% Thesis
%%% Function for automatic ROI selection
function [cell]=sel_auto_import(LBase,nroi,fc,A,H,Q,R,ratio_frame,vector_pos)
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
% cell.track_label=nroi; %% selected cell vector (inner matlab label, 0 <= backgorund)
if vector_pos(nroi,1) ~= -1
    cell.track_label = LBase(vector_pos(nroi,2),vector_pos(nroi,1));
% else
%     cell.track_label = 0;
% end
%%%%%%%%%%%%%%
% if 0==cell.track_label
%    warndlg('No cell was found, please select again');
%    cell=0;
% else
%% Cell Labeling
    dLaux = strcat('A',num2str(nroi,'%03u'));      %   default label
    cell.user_label = dLaux;                        %   assign default label
%% Feature Extraction   
    props=regionprops(LBase==cell.track_label,'BoundingBox','Centroid','Perimeter','Area');%bounding box struct
    cent=floor(props.Centroid); % centroids vector, just from selected cells
    b_box=floor(props.BoundingBox);
    
    cell.bbox=b_box; % bounding box [x y x_width y_width]
    cell.bbox_wh=[b_box(3),b_box(4)]; % saves the original dimension of bounding box
    cell.x_meas=[cent(1) cent(2) floor(props.Perimeter) floor(props.Area)]; % Features vector, measurements
    [cell.x_est,cell.ecov_matrix]=kalman_v1(A,H,Q,R,cell.x_meas(1:2)',zeros(6,1),eye(6,6),1);% [state vector,error covariance matrix]Kalman estimator
    cell.xy_est(fc,:)=cell.x_est;
    cell.reliability_counter=0;
    cell.xy_hist(fc,:)=[cent(1),cent(2)];                             % Initialization of cell's optic flux
    cell.d_hist(fc)=0;                                                % Initialization of cell's classification distance history 
    bbox_i=[cell.bbox(1) cell.bbox(2) cell.bbox_wh(1) cell.bbox_wh(2)];
    cell.mask = imcrop(LBase>0,bbox_i);
    crop_ratio=imcrop(ratio_frame,bbox_i);
    cell.ca_profile(fc)=sum(sum(crop_ratio.*cell.mask))/nnz(cell.mask);
%     cell.ca_profile(fc)=mean2(crop_ratio);  
    cell.status='ONLINE';                                              % Status initialization
    cell.color=[randi([2000,9000])/10000,randi([2000,9000])/10000,randi([2000,9000])/10000];   % Assign a color   
    cell.props = props;
else    
    cell=0;
end