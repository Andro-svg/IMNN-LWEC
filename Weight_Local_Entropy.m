function WeightImg = Weight_Local_Entropy(img, patchSize,slideSize)

[imgHei, imgWid] = size(img);
rowPatchNum = ceil((imgHei - patchSize.M) / slideSize.M) + 1;
colPatchNum = ceil((imgWid - patchSize.N) / slideSize.N) + 1;
rowPosArr = [1 : slideSize.M : (rowPatchNum - 1) * slideSize.M, imgHei - patchSize.M + 1];
colPosArr = [1 : slideSize.N : (colPatchNum - 1) * slideSize.N, imgWid - patchSize.N + 1];
% patchAll = zeros(patchSize, patchSize, rowPatchNum * colPatchNum);
% patchPos = zeros(1, 2, rowPatchNum * colPatchNum);
% k=0;

for row = rowPosArr
    for col = colPosArr
%         k=k+1;
%         patchPos(:, :, k) = [rowPosArr(row) colPosArr(col)];
        tmp_patch = img(row : row + patchSize.M - 1, col : col + patchSize.N - 1);
        G = 256;
        nk=zeros(G,1);
        %% Weight_Local_Entropy
        for i = 1:patchSize.M
            for j = 1:patchSize.N
                Img_level = tmp_patch(i,j)+1;  % 像素种类Img_level
                nk(Img_level)=nk(Img_level)+1;  % 每个种类的个数nk
            end
        end
        WeightPatch = zeros(patchSize.M, patchSize.N);
        idx = find(nk);
        nk(nk==0)=[];
        P = nk/(patchSize.M * patchSize.N); % 每种像素出现的概率
        for i = 1:patchSize.M
            for j = 1:patchSize.N
                for k = 1:length(nk)
                    WeightPatch(i,j) = WeightPatch(i,j) - abs(( idx(k)-1 - tmp_patch(i,j))^2) * P(k) * log(P(k));
                end
            end
        end
        WeightImg(row : row + patchSize.M - 1, col : col + patchSize.N - 1) = WeightPatch;
    end
end
end