function cleanSpath(spath)

theFiles = dir(spath);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(spath, baseFileName);
  fprintf(1, 'Now deleting %s\n', fullFileName);
  delete(fullFileName);
end