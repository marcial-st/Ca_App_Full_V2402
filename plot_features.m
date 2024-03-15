function plot_features(class_data,color,exp_title)

nbins=20;

if  ~isempty(class_data)
    feature_names  = fieldnames(class_data.cluster_features);
    n_features = length(feature_names);
    figure;    
    for i_feat = 1:1:n_features
        feature_i = [];
        feature_i = [extractfield(class_data.cluster_features,feature_names{i_feat})];        
        if ~isempty(feature_i)
            feature_i = reshape(feature_i,1,[]);
            subplot(4,2,i_feat);  histogram(feature_i,nbins,"FaceColor",color); ylabel("Frequency"); xlabel("Value");title(feature_names(i_feat),'Interpreter','none'); text(300,2.5,strcat("n=",num2str(length(feature_i))));
        end       
    end
    sgtitle(exp_title)
end
end