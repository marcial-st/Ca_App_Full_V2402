function [NCLS,string_class] = meltNCLS(NCLS,NCLS_2_meld,string_class)

n_classes = length(NCLS);
n_classes_2_meld = length(NCLS_2_meld);

disp('Melting structures ...')

for i_2meld = 1:1:n_classes_2_meld
    NCLS(n_classes+i_2meld).color       = NCLS_2_meld(i_2meld).color;
    % NCLS(n_classes+i_2meld).cPercentage = NCLS_2_meld(i_2meld).cPercentage;
    NCLS(n_classes+i_2meld).cVector     = NCLS_2_meld(i_2meld).cVector;
    NCLS(n_classes+i_2meld).label       = NCLS_2_meld(i_2meld).label;
    NCLS(n_classes+i_2meld).title       = NCLS_2_meld(i_2meld).title;

    string_class(n_classes+i_2meld)       = NCLS_2_meld(i_2meld).label;
end

disp('Structs melt   :)')