function process_list = processListOnlyEnhancement(process_list)

idx = max(size(process_list));    
for i =1:1:idx
selection = string(process_list(i).name);
    switch selection
        case 'MedianFilter'            
            disp(selection)
        case 'Contrast Eq'            
            disp(selection)
        case 'Adaptive Contrast Eq'            
            disp(selection)
        otherwise
            process_list(i).name = [];
            process_list(i).params = [];            
    end
end