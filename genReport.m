function genReport(spath,exp_name,mtable,data)

import mlreportgen.report.*
import mlreportgen.dom.*;

filePattern = fullfile(spath, '*.png'); % Change to whatever pattern you need.
theFiles = struct2table(dir(filePattern));
theFiles = sortrows(theFiles,"name","ascend");

d = Document(fullfile(spath,strcat("Report_",exp_name)),"pdf");
% open(d);

layout = d.CurrentPageLayout;
layout.PageSize.Height = '8.5in';
layout.PageSize.Width = '11in';
layout.PageSize.Orientation  ='landscape';
layout.PageMargins.Left  = '1cm';
layout.PageMargins.Right = '1cm';
layout.PageMargins.Top = '1cm';
layout.PageMargins.Bottom = '1cm';

pagenumber = PageNumber();
layout.Style = [{pagenumber}];

% Create the footer and add a page number to it
myfooter = PDFPageFooter();
para = Paragraph();
para.HAlign = 'center';
append(para,Page());

% Add the page number to the footer
append(myfooter,para);
layout.PageFooters = myfooter;

% f = waitbar(0,"Generating pdf report ...");
% ch = Chapter(strcat("Chapter Report of Exp: ",exp_name));
% append(d,ch)

p = Paragraph(strcat("Report of Exp: ",exp_name));
p.Bold = true;
append(d,p);
disp(strcat("Generating report, please wait ..."));

for k=1:2
   plot1 = Image(fullfile(theFiles.folder{k},theFiles.name{k}));
   plot1.Style = { HAlign('center') };
   plot1.Height = "3.8in";      
   append(d,plot1); 
end

if ~isempty(mtable)
    table1 = MATLABTable(mtable);
    append(d,table1);
end

if ~isempty(data)
    par=Paragraph(' ');
    par.Style = {OuterMargin('0.5in','0in','0in','12pt')};
    append(d,par);
    fnames = fieldnames(data);
    for kp=1:1:length(fieldnames(data))
        par = Paragraph(strcat(fnames{kp}," = ",num2str(data.(fnames{kp}),'%f04')));
        append(d,par);
    end
end

disp(strcat("     ... 0 %"));
for k=3:length(theFiles.name)
   plot1 = Image(fullfile(theFiles.folder{k},theFiles.name{k}));
   plot1.Style = { HAlign('center') };
   plot1.Height = "3.8in";      
%    plot1.Style = {ScaleToFit};
%    plot1.Width = "4in";
   append(d,plot1); 
%    waitbar(k/length(theFiles.name),f)
   disp(strcat("     ... ",num2str(round(k*100/length(theFiles.name)))," %"));
end
disp(strcat("Opening pdf ..."));
rptview(d);
disp(strcat("Report Finished (:"));
% close(f)
% close(d);
