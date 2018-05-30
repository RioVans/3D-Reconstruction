function disparityMap = stereo_computation(im1, im2, min_dis, max_dis, win_size, method)
%{  
  Yi Xiao
  3D VISION
  Master in Computer Vision  2017-2018
  Computer Vision Center (Barcelona, Spain)
  ---------------------------------------------
  ---------------------------------------------
  The input parameters are below 5:
  - left image
  - right image
  - minimum disparity
  - maximum disparity
  - window size (e.g. a value of 3 indicates a 3x3 window)
  - method (the method of matching cost, can be chosen between SSD, NCC and BW)

  The output:
  - the disparity map between a pair of rectified images
%}

im1_gray = rgb2gray(im1);
im2_gray = rgb2gray(im2);

[row1,colum1,cn1] = size(im1);
disparityMap = zeros(row1,colum1);
% to get the patch of left image; (i,j) is the index of the centre of
% the patch
if strcmp(method,'SSD')==1
    for i = 1+floor(win_size/2):row1-floor(win_size/2)
        for j = 1+floor(win_size/2)-min_dis:colum1-floor(win_size/2)-max_dis
            im1_window = im1_gray(i-floor(win_size/2):i+floor(win_size/2),j-floor(win_size/2):j+floor(win_size/2));       
            id=1;
            clear ssd;
            for dis = min_dis:max_dis
                im2_window = im2_gray(i-floor(win_size/2):i+floor(win_size/2),j-floor(win_size/2)+dis:j+floor(win_size/2)+dis);
                ssd(id) = sum(sum(abs((double(im1_window)-double(im2_window))).^2));
                id=id+1;
            end
            [m,p] = min(ssd);
            match_p = abs(min_dis+p-1);
            disparityMap(i,j)=match_p;
        end
    end
end

if strcmp(method,'NCC')==1
    for i = 1+floor(win_size/2):row1-floor(win_size/2)
        for j = 1+floor(win_size/2)-min_dis:colum1-floor(win_size/2)-max_dis
            im1_window = im1_gray(i-floor(win_size/2):i+floor(win_size/2),j-floor(win_size/2):j+floor(win_size/2));
            im1_window_mean = sum(sum(im1_window));
            id=1;
            clear ncc;
            for dis = min_dis:max_dis
                im2_window = im2_gray(i-floor(win_size/2):i+floor(win_size/2),j-floor(win_size/2)+dis:j+floor(win_size/2)+dis);
                im2_window_mean= sum(sum(im2_window));
                alpha_left = sqrt(sum(sum((double(im1_window)-im1_window_mean).^2)));
                alpha_right = sqrt(sum(sum((double(im2_window)-im2_window_mean).^2)));
                ncc_up = sum(sum((double(im1_window)-im1_window_mean).*(double(im2_window)-im2_window_mean)));
                ncc_dowm = alpha_left*alpha_right;
                ncc(id)= ncc_up/ncc_dowm;
                id = id+1;
            end
            [m,p] = max(ncc);
            match_p = abs(min_dis+p-1);
            disparityMap(i,j)=match_p;
        end
    end
end

if strcmp(method,'BW')==1
    T=40;
    gamma_c=5;
    gamma_p=win_size/2;
    centre=floor(win_size/2)+1;
    for i = 1+floor(win_size/2):row1-floor(win_size/2)
        for j = 1+floor(win_size/2)-min_dis:colum1-floor(win_size/2)-max_dis
            im1_window = im1_gray(i-floor(win_size/2):i+floor(win_size/2),j-floor(win_size/2):j+floor(win_size/2));
            id=1;
            for dis = min_dis:max_dis
                im2_window = im2_gray(i-floor(win_size/2):i+floor(win_size/2),j-floor(win_size/2)+dis:j+floor(win_size/2)+dis);
                numerator = 0;
                denominator =0;
                for m=1:win_size
                    for n=1:win_size
                        euc_dist= sqrt((m-centre)^2+(n-centre)^2);
                        weight_im1=exp(-((1/3*abs(double(im1_window(m,n))-double(im1_window(centre,centre)))/gamma_c)+(euc_dist/gamma_p)));
                        weight_im2=exp(-((1/3*abs(double(im2_window(m,n))-double(im2_window(centre,centre)))/gamma_c)+(euc_dist/gamma_p)));
                        e=min(abs(double(im1_window(m,n))-double(im2_window(m,n))),T);
                        numerator = numerator+weight_im1*weight_im2*e;
                        denominator = denominator+weight_im1*weight_im2; 
                    end
                end
                cost(id)= numerator/denominator;
                id=id+1;
            end
            [m,p] = min(cost);
            match_p = abs(min_dis+p-1);
            disparityMap(i,j)=match_p;
        end   
    end
end
disparityMap=uint8(disparityMap*255/max(max(disparityMap)));
end