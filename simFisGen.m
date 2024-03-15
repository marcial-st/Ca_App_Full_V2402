function [fis]=simFisGen(NCLS)

n_class = max(size(NCLS));

labels = [NCLS(:).label]';
class_dist_mean = [NCLS(:).class_dist_mean]';
tparams=table(labels,class_dist_mean);
tparams = sortrows(tparams, 'class_dist_mean');
dist_x_range = max([NCLS(:).class_dist_max]);

%% MFs setting
fis = mamfis;
%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fis = addInput(fis,[0 dist_x_range],"Name","distance_x");
fis = addMF(fis,"distance_x","zmf",[0.5*tparams.class_dist_mean(1) 1.5*tparams.class_dist_mean(1)],"Name",tparams.labels(1));
for i =2:1:n_class-1
    dstd = abs(tparams.class_dist_mean(i)-tparams.class_dist_mean(i-1))/4;
    fis = addMF(fis,"distance_x","gaussmf",[tparams.class_dist_mean(i) dstd],"Name",tparams.labels(i));
end
fis = addMF(fis,"distance_x","smf",[tparams.class_dist_mean(i) tparams.class_dist_mean(n_class)],"Name",tparams.labels(n_class));

%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fis = addOutput(fis,[0 1],"Name","Output");
step_out = [0:1/(2*n_class):1];
for i =1:1:n_class
    index = 2*i-1;
    fis = addMF(fis,"Output","trimf",[step_out(index) step_out(index+1) step_out(index+2)],"Name",tparams.labels(i),'VariableType',"output");
end
%% Rules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rules = strcat("distance_x==",tparams.labels(1),"=>Output=",tparams.labels(1));
for i=2:1:n_class
    rule1 = strcat("distance_x==",tparams.labels(i),"=>Output=",tparams.labels(i));
    rules = [rules rule1];
end
fis = addRule(fis,rules);
%  fuzzyLogicDesigner(fis)
fis.Rules
figure;
plotmf(fis,'input',1)
figure;
plotmf(fis,'output',1)
