function img = contrast(img)
img = Contrastimg(img,4);


function contrastImg = Contrastimg(img,k)
[imgHei, imgWid] = size(img);
padded_img0 = padarray(img,[k k],'symmetric','pre');
padded_img = padarray(padded_img0,[k k],'symmetric','post');
contrastImg = zeros(imgHei, imgWid);

for row  = (k+1) : (imgHei + k)
    for col = (k+1) : (imgWid + k) 
        temp_img1 = padded_img( (row-1) : (row + 1), (col-1) : (col + 1));
        sum1 = sum(sum(temp_img1));
        mean1 = sum(sum(temp_img1)) / 9;

        temp_img2 = padded_img( (row-2) : (row + 2), (col-2) : (col + 2));
        sum2 = sum(sum(temp_img2));
        mean2 = (sum2-sum1)/(25-9);

        temp_img3 = padded_img( (row-3) : (row + 3), (col-3) : (col + 3));
        sum3 = sum(sum(temp_img3));
        mean3 = (sum3-sum2)/(49-25);
        
        temp_img4 = padded_img( (row-4) : (row + 4), (col-4) : (col + 4));
        sum4 = sum(sum(temp_img4));
        mean4 = (sum4-sum3)/(81-49);

        contrastcoef = [mean1/mean2,mean1/mean3,mean1/mean4];
        contrast = max(contrastcoef);
        if isnan(contrast)
            contrast=2;
        end
        contrastImg(row-k,col-k) =padded_img(row,col).^contrast;

    end
end
end
end