%% This function stores the different data of the current cell in the CELLS structure
function[cell,bbox_i]=store_data(cell,new_tl,fc,LM_fc1_1,A,H,Q,R,ratio_frame)
if new_tl==0 
    cell.track_label=0;
    cell.user_label=cell.user_label;
    cell.bbox=cell.bbox;
    cell.x_meas=cell.x_meas;  
    cell.xy_hist(fc,:)=[floor(cell.x_est(1)),floor(cell.x_est(4))];
    cell.d_hist(fc)=0;
    [cell.x_est,cell.ecov_matrix]=kalman_v1(A,H,Q,R,cell.x_meas',cell.x_est,cell.ecov_matrix,0);% [state vector,error covariance matrix]Kalman estimator
    cell.color=cell.color;
    bbox_i=[cell.bbox(1) cell.bbox(2) cell.bbox_wh(1) cell.bbox_wh(2)];
    crop_ratio=imcrop(ratio_frame,bbox_i);
    cell.ca_profile(fc)=mean2(crop_ratio);
    bbox_i=[cell.bbox(1) cell.bbox(2) cell.bbox_wh(1) cell.bbox_wh(2)]; 
    cell.reliability_counter=cell.reliability_counter+1;    
elseif new_tl== 9999
%     current_cell=logical(LM_fc1_1==new_tl);
%     ccb=bwconncomp(current_cell);              
%     st_aux=regionprops(ccb,'Centroid','BoundingBox','Perimeter','Area');     % centroide de la intersección para buscar en el frame n+1
    cell.track_label=new_tl; % get cell label from frame n+1 "track_label"
%     cell.bbox=st_aux.BoundingBox;
%     cell.x_meas=[floor(st_aux.Centroid(1)) floor(st_aux.Centroid(2)) floor(st_aux.Perimeter) floor(st_aux.Area)]; % Features vector, measurements
    [cell.x_est,cell.ecov_matrix]=kalman_v1(A,H,Q,R,cell.x_meas(1:2)',cell.x_est,cell.ecov_matrix,1);% [state vector,error covariance matrix]Kalman estimator
    cell.xy_est(fc,:)=cell.x_est;    
%     cell.xy_hist(fc,:)=[floor(st_aux.Centroid(1)),floor(st_aux.Centroid(2))];                             % Cell's optic flux
    cell.xy_hist(fc,:)=cell.xy_hist(fc-1,:);
    cell.d_hist(fc)=0;                                           % Distance history
%     bbox_i=[floor(cell.x_est(4)) floor(cell.x_est(1)) cell.bbox_wh(1) cell.bbox_wh(2)];
    bbox_i=[cell.bbox(1) cell.bbox(2) cell.bbox_wh(1) cell.bbox_wh(2)];
    crop_ratio=imcrop(ratio_frame,bbox_i);
    try
%         cell.ca_profile(fc)=sum(sum(crop_ratio.*cell.mask))/nnz(cell.mask);
        cell.ca_profile(fc)=mean2(crop_ratio);
    catch
        disp(strcat("Wrn! Ca profile calculation error, bbox_i=[",num2str(bbox_i),"], fc=[",num2str(fc),"]"))
        cell.mask
        crop_ratio
        cell.ca_profile(fc)=mean2(crop_ratio);
    end
%     bbox_i=[cell.bbox(1) cell.bbox(2) cell.bbox_wh(1) cell.bbox_wh(2)];
    cell.reliability_counter=cell.reliability_counter+1;
else
    current_cell=logical(LM_fc1_1==new_tl);
    ccb=bwconncomp(current_cell); 
    props=regionprops(ccb,'BoundingBox','Centroid','Perimeter','Area','MajoraxisLength','MinoraxisLength','Orientation');% centroide de la intersección para buscar en el frame n+1  
    cell.track_label=new_tl; % get cell label from frame n+1 "track_label"
    cell.bbox=props.BoundingBox;
    cell.x_meas=[floor(props.Centroid(1)) floor(props.Centroid(2)) floor(props.Perimeter) floor(props.Area)]; % Features vector, measurements
    [cell.x_est,cell.ecov_matrix]=kalman_v1(A,H,Q,R,cell.x_meas(1:2)',cell.x_est,cell.ecov_matrix,1);% [state vector,error covariance matrix]Kalman estimator
    cell.xy_est(fc,:)=cell.x_est;    
    cell.xy_hist(fc,:)=[floor(props.Centroid(1)),floor(props.Centroid(2))];                             % Cell's optic flux
    cell.d_hist(fc)=0;                                           % Distance history
    bbox_i=[cell.bbox(1) cell.bbox(2) cell.bbox_wh(1) cell.bbox_wh(2)];
    crop_ratio=imcrop(ratio_frame,bbox_i);    
    try
%         cell.ca_profile(fc)=sum(sum(crop_ratio.*cell.mask))/nnz(cell.mask);
        cell.ca_profile(fc)=mean2(crop_ratio);
    catch
        disp(strcat("Wrn! Ca profile calculation error, bbox_i=[",num2str(bbox_i),"], fc=[",num2str(fc),"]"))
        cell.mask      
        crop_ratio
        disp("size(cell.mask)=")
        size(cell.mask)
        disp("size(crop_ratio)=")
        size(crop_ratio)
        cell.ca_profile(fc)=mean2(crop_ratio);
    end
%     cell.ca_profile(fc)=mean2(crop_ratio);
%     bbox_i=[cell.bbox(1) cell.bbox(2) cell.bbox_wh(1) cell.bbox_wh(2)];
    cell.reliability_counter=0;
end