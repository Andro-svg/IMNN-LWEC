function [patchTen, patchPosition] = construct_local_patch_ten(img,rowPosArr,colPosArr, patchSize,N,M)

[meshCols, meshRows] = meshgrid(colPosArr, rowPosArr);
idx_fun1 = @(row,col) img(row : row + patchSize - 1, col : col + patchSize - 1);

patchCell = arrayfun(idx_fun1, meshRows, meshCols, 'UniformOutput', false);
patchTen = cat(3, patchCell{:}); 

idx_fun2 = @(row,col) [row,col];
patchPosCell = arrayfun(idx_fun2, meshRows, meshCols, 'UniformOutput', false);
patchPosition = cat(3, patchPosCell{:});


end