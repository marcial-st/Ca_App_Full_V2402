%%%%%% Thesis
%%%% MST,JMRC,RBR
%%%% Intersection based tracking 
function[new_t_l]=int_tl(cck,LM_fc1_1)
    % LM_fc ---> n
    % LM_fc_1 ---> n+1
ssk=regionprops(cck,'Centroid');     % centroide de la intersección para buscar en el frame n+1
cent=floor(ssk.Centroid);
new_t_l=LM_fc1_1(cent(2),cent(1));