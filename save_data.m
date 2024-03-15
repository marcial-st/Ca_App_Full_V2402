%% Thesis
%% Save data function
% function [all,frames,header,DATA]=save_data(CELLS,k)
function [DATA_table]=save_data(CELLS,k)
% temporary WA to save only 80 ROIs due to xls java limit
[~,sc]=size(CELLS);
if(sc)>227
    CELLS=CELLS(1:227);
end
%
[r,c]=size(CELLS);
DATA_table=cell(k,3*c+1);
frames=(1:1:k-1)';
DATA=zeros(k-1,3*c);
for i=1:1:c
    header(1,3*i-2)=strcat({'Mean ('},CELLS(i).user_label,{')'});
    header(1,3*i-1)=strcat({'X ('},CELLS(i).user_label,{')'});
    header(1,3*i)=strcat({'Y ('},CELLS(i).user_label,{')'});
    if CELLS(i).reliability_counter==0
        DATA(:,3*i-2)=CELLS(i).ca_profile';
        DATA(:,3*i-1:3*i)=CELLS(i).xy_hist;
    else
        [lc_m,lc_n]=size(CELLS(i).ca_profile)
        DATA(1:lc_n,3*i-2)=CELLS(i).ca_profile';
        DATA(1:lc_n,3*i-1:3*i)=CELLS(i).xy_hist;
    end
end
DATA_table(1,2:3*c+1)=header;
DATA_table(2:k,1)=num2cell(frames);
DATA_table(2:k,2:3*c+1)=num2cell(DATA);