%%% Feature extraction of geometric parameters 
clear all; close all;clc

Exp_path = ["G:\Mi unidad\CaExp_2019\100528.1\340",
"G:\Mi unidad\CaExp_2019\100528.4\340",
"G:\Mi unidad\CaExp_2019\100528.5\340",
"G:\Mi unidad\CaExp_2019\100528.6\340",
"G:\Mi unidad\CaExp_2019\100528.8\340",
"G:\Mi unidad\CaExp_2019\100601.1\340"]

load("processList_italy_v200831.mat");
visual_en = 1;

% for i=1:1:length(Files)

local_width       = [];
local_height      = [];
local_orientation = [];
local_axismajor   = [];
local_axisminor   = [];

global_width       = [];
global_height      = [];
global_orientation = [];
global_axismajor   = [];
global_axisminor   = [];
for i_exp = 1:1:length(Exp_path)
    current_exp = Exp_path(i_exp)
    Files = dir(fullfile(current_exp,'*.tif'));
    for i_frame=1:1:1
        img_i = imread(fullfile(current_exp,Files(i_frame).name));
        img_i_bin = processListApply(process_list,img_i);    
        % if visual_en
        %     figure;subplot(1,2,1);imshow(img_i);subplot(1,2,2);imshow(img_i_bin)
        % end
        s = regionprops(img_i_bin,'BoundingBox','MajorAxisLength','MinorAxisLength','Orientation');
        l_nblobs = length(s)
        for k=1:1:length(s)
            s(k).Width = s(k).BoundingBox(3);
            s(k).Heigth = s(k).BoundingBox(4);
        end
    
        local_width       = [local_width;[s.Width]'];
        local_height      = [local_height;[s.Heigth]'];
        local_orientation = [local_orientation;[s.Orientation]'];
        local_axismajor   = [local_axismajor;[s.MajorAxisLength]'];
        local_axisminor   = [local_axisminor;[s.MinorAxisLength]'];
    end
    l_w = max(local_width)
    l_h = max(local_height)
    l_whr_m = mean(local_width./local_height)
    l_whr_s = std(local_width./local_height)


    global_width       = [global_width;      local_width      ];
    global_height      = [global_height;     local_height     ];
    global_orientation = [global_orientation;local_orientation];
    global_axismajor   = [global_axismajor;  local_axismajor  ];
    global_axisminor   = [global_axisminor;  local_axisminor  ];
end

    g_w = max(global_width)
    g_h = max(global_height)
    g_whr_m = mean(global_width./global_height)
    g_whr_s = std(global_width./global_height)


if visual_en
    nbins = 30;
    fsize_xylabel = 30;
    fsize_axis = 20;
    %
    f1=figure; axes1 = axes('Parent',f1);hold(axes1,'on');
    histfit(global_width,nbins,'Normal');
    set(axes1,'FontName','Arial','FontSize',fsize_axis);
    ylabel('Frequency','FontSize',fsize_xylabel,'FontName','Arial');
    xlabel('Width (pixels)','FontSize',fsize_xylabel,'FontName','Arial');
    hold(axes1,'off');
    %
    f2=figure; axes2 = axes('Parent',f2);hold(axes2,'on'); 
    histfit(global_height,nbins,'Normal');
    set(axes2,'FontName','Arial','FontSize',fsize_axis);
    ylabel('Frequency','FontSize',fsize_xylabel,'FontName','Arial');
    xlabel('Height (pixels)','FontSize',fsize_xylabel,'FontName','Arial');
    hold(axes2,'off'); 
    %
    figure; histfit(global_orientation,nbins,'Normal');title("Orientation");xlabel("Degree");ylabel("Frequency")
    figure; histfit(global_axismajor,nbins,'Normal');title("AxisMajor");xlabel("Pixels");ylabel("Frequency")
    figure; histfit(global_axisminor,nbins,'Normal');title("AxisMinor");xlabel("Pixels");ylabel("Frequency")
    %
    f3=figure; axes3 = axes('Parent',f3);hold(axes3,'on'); 
    histogram(global_width./global_height,nbins);
    set(axes3,'FontName','Arial','FontSize',fsize_axis);
    ylabel('Frequency','FontSize',fsize_xylabel,'FontName','Arial');
    xlabel('Width to Height Ratio (a.u.)','FontSize',fsize_xylabel,'FontName','Arial');
    hold(axes3,'off');    
end

pos=floor(reshape([s.BoundingBox],4,[])');
figure;imshow(img_i_bin);hold on
for k=1:1:length(s)
    rectangle('Position',pos(k,:),'EdgeColor','r','LineWidth',2)
end

