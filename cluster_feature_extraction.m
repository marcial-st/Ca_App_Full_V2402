function cluster_features = cluster_feature_extraction(IDX,C,K,cells_ci_loc,inj_band,pix2um)

for i_cluster=1:1:K
    cent_clt_i = C(i_cluster,:);
    clt_i_cent_v = ones(length(IDX),2).*cent_clt_i;                                             % cluster i centroid vector
    dist_all_2_clti_v = sqrt(diag((clt_i_cent_v-cells_ci_loc)*(clt_i_cent_v-cells_ci_loc)'));   % distance all cells 2 cluster i centroid
                  
    if ~isempty(inj_band)
        cluster_features(i_cluster).dist2injury = min(abs(inj_band-cent_clt_i(1)));                     % Feat: clt i distance to injury boundary
    else
        disp("[cluster_features] Warning! injury field is empty, so dist2injury is set to 9999")
        cluster_features(i_cluster).dist2injury = 9999;
    end
    cluster_features(i_cluster).dist_cells2centroid = dist_all_2_clti_v(IDX==i_cluster)*pix2um; % Feat: cells clt i distance to clt centroid
    cluster_features(i_cluster).ncells_subclt = sum(IDX==i_cluster);                            % Feat: number of cells per subcluster   
    cluster_cell_idx = find(IDX==i_cluster);                                                    % indices of cells of clt i
    cluster_cell_loc = cells_ci_loc(cluster_cell_idx,:);                                        % location of cells of clt i
                                                                                                % compute mean dist between cells of clt i
    % cell_i_cent_v = ones(length(cluster_cell_idx),2).*cluster_cell_loc(1);                      % centroid of cell_i from clt 1
    dist_between_cells = (sqrt(diag((clt_i_cent_v-cells_ci_loc)*(clt_i_cent_v-cells_ci_loc)')));% mean distance between cell i and other cells of clt i                       
    dist_between_cells=dist_between_cells*pix2um;

    cluster_features(i_cluster).dist_btw_cells_mean = mean(dist_between_cells);                 % Feat: mean distance vector between cells of clt i
    cluster_features(i_cluster).dist_btw_cells_max = max(dist_between_cells);                   % Feat: max distance vector between cells of clt i
    cluster_features(i_cluster).dist_btw_cells_min = min(dist_between_cells);                   % Feat: min distance vector between cells of clt i
    clt_i_cent_v_clt = ones(K,2).*cent_clt_i;
    dist_between_clt = mean(sqrt(diag((clt_i_cent_v_clt-C)*(clt_i_cent_v_clt-C)')))*pix2um;     % mean distance between cell i and other cells of clt i                           
    cluster_features(i_cluster).dist_mean_btw_clt = dist_between_clt;
    cluster_features(i_cluster).location_y = clt_i_cent_v(2)*pix2um;                            % Feat: y coordinate
end

