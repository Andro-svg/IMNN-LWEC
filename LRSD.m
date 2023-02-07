function [L,S,Noise] = LRSD(X,tenP, opts)

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'max_rho');      max_rho = opts.max_rho;        end
if isfield(opts, 'max_beta');      max_beta = opts.max_beta;        end
if isfield(opts, 'gamma');         gamma = opts.gamma;              end
if isfield(opts, 'lambda1');        lambda1 = opts.lambda1;                end
if isfield(opts, 'lambda2');        lambda2 = opts.lambda2;                end
if isfield(opts, 'epsilon');        epsilon = opts.epsilon;                end
if isfield(opts, 'omega');        omega = opts.omega;                end
if isfield(opts, 'betaL');         betaL = opts.betaL;              end
if isfield(opts, 'alpha');         alpha = opts.alpha;              end


%% Initialization
Nway = size(X);
N = ndims(X);
L = zeros(Nway);
Noise = zeros(Nway);
Noise2 = zeros(Nway);

Z = cell(N,N);
for i=1:N
    for j=1:N
    Z{i,j} = zeros(Nway);
    end
end
P = Z;
S = zeros(Nway);
M = zeros(Nway);

temp = Z;
preNumT = numel(S);
beta = betaL/omega;
rho = 1/omega;
W = ones(size(S))./ tenP;

for iter = 1 : max_iter
    %% update Z
    Lold = L;
    tau = alpha./beta;
    for i=1:N-1
        for j=i+1:N
            tempL   = tensor_permute(L,Nway,i,j);
            tempP   = tensor_permute(P{i,j},Nway,i,j);
            Z{i,j}  = tensor_ipermute(prox_tnn_my( tensor_permute(Z{i,j},Nway,i,j), tempL + tempP/beta(i,j),  tau(i,j),  epsilon),   Nway,i,j);
            temp{i,j} = Z{i,j}-P{i,j}/beta(i,j);
        end
    end

    %% update L
    tempsum = zeros(Nway);
    for i=1:N-1
        for j=1+i:N
        tempsum = tempsum+ beta(i,j)*temp{i,j};
        end
    end
    L = (tempsum+ rho*(X-S-Noise-Noise2+M/rho))/(rho+sum(beta(:)));

    %% update S
    S = prox_l1(X-L-Noise-Noise2+M/rho,W.*lambda1/rho);
    W = 2 ./ ( (abs(S))+ 0.01) ./ tenP ;

    %% update edge
    % L1,1,2 norm
    Q = X-L-S -Noise2 + M/rho;
    for i = 1:Nway(2)
        Noise(:,i,:) = max((1-lambda2/(rho*norm(squeeze(Q(:,i,:))))),0) * Q(:,i,:);
    end

    for m = 1:Nway(1)
        Noise(m,:,:) = max((1-lambda2/(rho*norm(squeeze(Q(m,:,:))))),0) * Q(m,:,:);
    end


    %% update Noise2 
    Noise2 = (rho * (X-L-S-Noise) + M)/(2*lambda2+rho); % Frobenius norm

    %% check the convergence
    dM = X-L-S-Noise-Noise2;
    currNumT = sum(S(:) > 0); 
    chg =norm(X(:)-L(:)-S(:)-Noise(:)-Noise2(:))/norm(X(:));
    fprintf('iter = %d   res=%.10f  \n', iter, chg);    

    if (chg < tol) || (currNumT == preNumT)
        break;
    end
    preNumT = currNumT;  

    %% update M & P
    for i=1:N-1
        for j=i:N 
           P{i,j} = P{i,j}+beta(i,j) * (L-Z{i,j});
        end
    end
    M = M + rho*dM;
    beta = min(gamma*beta,max_beta);  
    rho = min(gamma*rho,max_rho);  
end

