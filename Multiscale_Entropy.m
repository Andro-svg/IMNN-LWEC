function WeightImg_2 = Multiscale_Entropy(img)

WeightImg_2 = Entropy(img,3);

function WeightImg = Entropy(img,k)
[imgHei, imgWid] = size(img);
padded_img0 = padarray(img,[k k],'symmetric','pre');
padded_img = padarray(padded_img0,[k k],'symmetric','post');


WeightImg = zeros(imgHei, imgWid);
sigma = zeros(imgHei, imgWid);
for row  = (k+1) : (imgHei + k)
    for col = (k+1) : (imgWid + k)
        tmp_img = padded_img( (row-k) : (row + k), (col-k) : (col + k));
        G = 1 + double(max(max(img)));
        nk=zeros(G,1);
        %% Weight_Local_Entropy
        for i = 1:(2*k+1)
            for j = 1:(2*k+1)
                Img_level = tmp_img(i,j)+1;  
                nk(Img_level)=nk(Img_level)+1;
                
            end
        end
        meanvalue = mean(mean(tmp_img));
        sigma(row-k,col-k) = sqrt( mean(mean(   (tmp_img - meanvalue).^2 )) );

        idx = find(nk);
        nk(nk==0)=[];
        P = nk/((2*k+1) * (2*k+1));
        for u = 1:length(nk)
            WeightImg(row-k,col-k) = WeightImg(row-k,col-k) - abs((idx(u)-1 - padded_img(row,col)).^2) * P(u) * log(P(u));  
        end
    end
end
WeightImg = (WeightImg./(sigma+1e-2));

end

end