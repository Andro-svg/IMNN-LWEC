function [All_Num,time_per_image] = target_detection(readPath, savePath, temporal_step, imgNum, tensor, tuneopts)

if isfield(tuneopts, 'lambdaL');         lambdaL = tuneopts.lambdaL;   end
if isfield(tuneopts, 'omega');           omega = tuneopts.omega;   end


%% Create folders for saving results
if ~exist(savePath)  
    mkdir(savePath);   
end


%% Get all image file names, please make sure that the image file order is correct by this reading way.
filesdir = dir([readPath '/*.jpg']);
if isempty( filesdir )
    filesdir = dir( [readPath '/*.bmp'] );
end
if isempty( filesdir )
    filesdir = dir([readPath '/*.png']);
end
if isempty( filesdir )
    fprintf('\n There is no any image in the folder of %s', readPath);
    return;
end
% get all image file names into a cell array;
files = { filesdir.name };
files = sort_nat(files);

iteration = 0;
time_all=0;
All_Num = length(files);
%% begin to process images using mog based detection method.
t_list =  [1 : temporal_step : length(files)-imgNum + 1, length(files)-imgNum + 1 ];

for t = t_list
    
    iteration = iteration + 1;
    disp('========================================');
    fprintf('%s %d%s\n','Starting', iteration, '-th loop');
    spat_temp_ten = []; 
    priorWeight_ten = [];
    spat_temp_pos = [];
    spat_ten = [];
    tenW = [];
    spat_pos = [];
    T = [];
    BKG = [];
    %% read images and construct the patch image
    t1=clock;
    starttime = clock;
    for tt = 1 : imgNum
        disp([readPath '\' files{tt+t-1}]);
        img = imread([readPath '/' files{tt+t-1}]);
        if size(img, 3) > 1
            img = rgb2gray( img );
        end
        [imgHei, imgWid] = size(img);
        img = double(img);
        imwrite(mat2gray(img), [savePath '/' strtok([files{tt+t-1}],'.') '_ori.jpg']);
        Weight_Entropy = Multiscale_Entropy(img);
        contrastImg = contrast(img);

        lambda1=[];
        lambda2=[];
        lambda11=[];
        lambda22=[];  
        [lambda1, lambda2] = structure_tensor_lambda(Weight_Entropy, 'Gaussian', 3);
        EI = lambda1;
        CI = lambda1.*lambda2./(lambda2+lambda1+0.0001);
        x = 1;
        y = 1;
        prior0 = mat2gray((x+y)*EI.*CI./(x*EI+y*CI+0.001));

        [lambda11, lambda22] = structure_tensor_lambda(contrastImg, 'Gaussian', 3);
        Trace=(lambda11+lambda22);
        Det=(lambda11.*lambda22);

        Trace(Trace==Inf)=[];
        maxT = max(max(Trace));
        minT = min(min(Trace));
        
        Det(Det==Inf)=[];
        maxD = max(max(Det));
        minD = min(min(Det));

        % if not line
        threshT2 =minT +(maxT - minT)*0.13;
        threshD2 = minD +(maxD - minD)*0.01;
        L = zeros([imgHei, imgWid]);
        L = double(not((Trace>= threshT2)  & ( Det <= threshD2)));

        % if point
        threshT1 = minT +(maxT - minT)*0.1;
        threshD1 = minD +(maxD - minD)*0.008;
        P = zeros([imgHei, imgWid]);
        P = double( (Trace>= threshT1) &  ( Det>= threshD1) );
        
        W = 3*P + L;
        prior = prior0.* W;

        %% construct patch tensor
        start0 = clock;
        [tenF, patchNumber, patchPosition] = construct_patch_ten(img, tensor.patchSize, tensor.slideStep);
        spat_temp_ten(:, :, (tt-1) * patchNumber + 1 : tt * patchNumber) = tenF;
        spat_temp_pos(:, :, (tt-1) * patchNumber + 1 : tt * patchNumber) = patchPosition;
        [priorWeight, ~, ~] = construct_patch_ten( prior, tensor.patchSize, tensor.slideStep);
        priorWeight_ten(:, :, (tt-1) * patchNumber + 1 : tt * patchNumber) = priorWeight;
        end0 =clock;
        times=etime(end0,start0);
        disp(['Using ' num2str(tt+t-1) '-th frames to construct and patch tensor:', num2str(times)]);
    end
    endtime = clock;

    %% The LRSD model
    Nway = size(spat_temp_ten);
    Ndim = ndims(spat_temp_ten);

    opts=[];
    opts.max_beta = 1e10*ones(Ndim,Ndim);
    opts.max_rho = 1e4;
    opts.max_iter = 400;
    opts.tol =1e-4;  
    w = 0.001;
    opts.alpha=[0,  w/(2+w),  1/(2+w);
    0,    0,    1/(2+w);
    0,    0,    0]; 
    opts.gamma = 1.5;
    opts.lambda1 = 1e5*lambdaL/(sqrt(max(Nway(1),Nway(2))) * Nway(3));
    opts.lambda2 = 50* opts.lambda1; 
    opts.epsilon = 1e-4;
    opts.omega = omega; % tune
    opts.betaL =[0,    1,    1;
                 0,    0,    1;
                 0,    0,    0]; 

    tenT=[];
    tenB=[];
    tenN=[];
    [tenB, tenT, tenN] = LRSD(spat_temp_ten,priorWeight_ten,opts);

    %% reconstrcut imgNum images
    for kk = 1:imgNum  
        tarImg = zeros([imgHei, imgWid]);
        BKGImg = zeros([imgHei, imgWid]);
        NImg = zeros([imgHei, imgWid]);
        for ss = 1:patchNumber
            position = spat_temp_pos(:,:,(kk-1)*patchNumber+ss);
            row = position(1);
            col = position(2);
            tarImg(row: row + tensor.patchSize-1, col:col + tensor.patchSize-1) = tenT(:,:,(kk-1)*patchNumber+ss);
            BKGImg(row: row + tensor.patchSize-1, col:col + tensor.patchSize-1) = tenB(:,:,(kk-1)*patchNumber+ss);
            NImg(row: row + tensor.patchSize-1, col:col + tensor.patchSize-1) = tenN(:,:,(kk-1)*patchNumber+ss);
        end
        imwrite(mat2gray(tarImg), [savePath '/' strtok([files{kk+t-1}],'.') '_tar.jpg']);
    end
    t2=clock;
    disp(['Programming has running:',num2str(etime(t2,t1))]);
    disp('=====================================================')
    disp(' ')
    time_all = time_all + etime(t2,t1);

    time_per_iter = time_all/(iteration);
    disp(['Each iteration consumes time: ', num2str(time_per_iter)]);
    time_per_image = time_per_iter / imgNum;
    disp(['Each image consumes time: ', num2str(time_per_image)]);
end