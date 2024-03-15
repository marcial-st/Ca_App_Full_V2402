
function pix2um = get_px2num(im340_1,bar_um,visual_en)

if isempty(visual_en)
    visual_en = 0;
end
if isempty(bar_um)
    bar_um = 40;
end

imcal=rgb2gray(imread("cells_calibration.png"));
if visual_en 
    figure;imshow(imcal);
end
[yimcal,ximcal]=size(imcal);
ycrop = 453;xcrop=382;
imcal_crop = imcrop(imcal,[ycrop xcrop yimcal ximcal])>100;
if visual_en
    figure;imshow(imcal_crop);
end
bar_len = mean(nonzeros(sum(imcal_crop,2)));
[~,x_im340] = size(im340_1);
bar_len_scaled = x_im340*bar_len/ximcal;

pix2um=bar_um/bar_len_scaled;
