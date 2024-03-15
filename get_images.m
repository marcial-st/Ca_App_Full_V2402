% %%% Thesis
% %%% Ratio Sequence V_1
% %%% Load and save 340 nm image sequence
function [D340,D380,numCycles,numFrames340,fpc]=get_images

fpc=100; % frames per cycle

[filename, pathname] = uigetfile('*.tif','Select the 340nm first image');
if filename~=0
    D340 = dir(strcat([pathname, '*.tif']));
    numFrames340=size(D340);
    h1 = waitbar(0,'Loading...');
    for i=1:1:numFrames340(1)
%         nim=strcat(pathname,D340(i).name)
        waitbar(i/numFrames340(1))    
    end
    close(h1);
end

idcs   = strfind(pathname,filesep);
pathname = pathname(1:idcs(end-1));

[filename, pathname] = uigetfile(strcat(pathname,'*.tif'),'Select the 380nm first image');
if filename~=0
    D380 = dir(strcat([pathname, '*.tif']));
    numFrames380=size(D380);
    h1 = waitbar(0,'Loading...');
    for i=1:1:numFrames380(1)
%         nim=strcat(pathname,D380(i).name)
        waitbar(i/numFrames380(1))    
    end
    close(h1);
end

if(numFrames340==numFrames380)
    h1 = waitbar(0,'Series parsing test ...');
        
%     for i=1:1:fr_uf-fr_ui+1
    for i=1:1:numFrames340(1)       

        name340=D340(i).name;
        name380=D380(i).name;
        
        D340(i).number = str2double(regexp(name340,'\d+?\.','once','match'));
        D380(i).number = str2double(regexp(name380,'\d+?\.','once','match'));
        
        if(strcmp(name340(5:length(name340)),name380(5:length(name340)))==0)
            message=strcat('Name missmatch was found, please check the series. ',D340(i).name,'- and -',D380(i).name);
            errordlg(message,'Error')
            break
        end
%         waitbar(i/(fr_uf-fr_ui+1)) 
        waitbar(i/(numFrames340(1)))
    end
    
    [~,idx]=sort([D340.number]);
    D340=D340(idx);
    [~,idx]=sort([D380.number]);
    D380=D380(idx);
    
    close(h1);
    
    fr_l = length(D340);
    prompt = {'Isert range'};
    dlgtitle = 'Number of Frames';
    definput = {['1-',num2str(fr_l)]};
    fr_user = inputdlg(prompt,dlgtitle,[1 40],definput);
    expression = '\-';
    range_split = regexp(fr_user{1},expression,'split');
    fr_ui = str2double(range_split{1})
    fr_uf = str2double(range_split{2})
    
    if fr_ui>fr_uf
        errordlg('Please provide a range in ascendent order')        
    end
    if fr_uf>fr_l
        fr_uf = rf_l;
    end
    
    D340 = D340(fr_ui:fr_uf);
    D380 = D380(fr_ui:fr_uf);
    
%     message='Series matched successfully'
    if (floor(mod((fr_uf-fr_ui+1),fpc))==0)
%         numCycles = 1;
          numCycles = floor((fr_uf-fr_ui+1)/fpc);
    else
        numCycles=floor((fr_uf-fr_ui+1)/fpc+1);
    end
    numFrames340=(fr_uf-fr_ui+1);
else
    errordlg('Number of frames missmatch.','Error')
end