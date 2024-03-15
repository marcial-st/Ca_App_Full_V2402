function imlegend(colorArr, labelsArr)
% For instance if two legend entries are needed:
% colorArr =
%   Nx3 array of doubles. Each entry should be an RGB percentage value between 0 and 1
%
% labelsArr =
%   1×N cell array
%     {'First name here'}    {'Second name here'}    {'etc'}

hold on;
for ii = 1:length(labelsArr)
  % Make a new legend entry for each label. 'color' contains a 0->255 RGB triplet
  scatter([],[],2, colorArr(ii,:), 'filled', 'DisplayName', labelsArr{ii});
end
hold off;
legend();
end