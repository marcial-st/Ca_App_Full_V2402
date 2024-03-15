function pclass = class_net(cells_online,net)

n_prof = length(cells_online);
pclass = zeros(n_prof,1);

for i=1:1:length(cells_online)
	pmatrix(i,:)=cells_online(i).ca_profile;
end

wb = waitbar(0,'Net classifier progress...');
for i=1:1:n_prof
	cp=cells_online(i).ca_profile;
	ct=1:1:length(cp);

	fn=figure;
    plot(ct,cp,'LineWidth',2.5)
	ylim([min(min(pmatrix)) max(max(pmatrix))])
	xlim([1 length(ct)])
	axis off
	set(fn, 'Color', 'White')	
	pframe = getframe(fn);
	X = rgb2gray(frame2im(pframe));	
	Xr = imresize(X,[150 150],'nearest');
	
    class = classify(net,Xr);
	
    switch class
        case 'c0_unknown'
            pclass(i)=0;
        case 'c1_noresponse'
            pclass(i)=1;
        case 'c2_peak'
            pclass(i)=2;
        case 'c3_peakplateau'
            pclass(i)=3;
        case 'c4_peakoscillations'
            pclass(i)=4;
        case 'c5_peakoscillationsplateau'
            pclass(i)=5;
        case 'c6_oscillations'
            pclass(i)=6;
        case 'c7_oscillationsplateau'
            pclass(i)=7;
        otherwise
            pclass(i)=0;
    end
	waitbar(i/n_prof,wb,'Net classifier progress...');
end
close(wb)