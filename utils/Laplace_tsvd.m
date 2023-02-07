function [U,S,V, Laplace_S] = Laplace_tsvd(A)

% [U,S,V] = tsvd(A,transform,opt) computes the tensor SVD under linear transform, i.e., A=U*S*V^*, where S
% is a f-diagonal tensor, U and V are orthogonal under linear transform.
%
% Input:
%       A       -   n1*n2*n3 tensor
%   transform   -   a structure which defines the linear transform
%       transform.L: the linear transform of two types:
%                  - type I: function handle, i.e., @fft, @dct
%                  - type II: invertible matrix of size n3*n3
%
%       transform.inverseL: the inverse linear transform of transform.L
%                         - type I: function handle, i.e., @ifft, @idct
%                         - type II: inverse matrix of transform.L
%
%       transform.l: a constant which indicates whether the following property holds for the linear transform or not:
%                    L'*L=L*L'=l*I, for some l>0.                           
%                  - transform.l > 0: indicates that the above property holds. Then we set transform.l = l.
%                  - transform.l < 0: indicates that the above property does not hold. Then we can set transform.l = c, for any constant c < 0.
%       If not specified, fft is the default transform, i.e.,
%       transform.L = @fft, transform.l = n3, transform.inverseL = @ifft. 
%
%       opt     -   options for different outputs of U, S and V:
%                   'full': (default) produces full tensor SVD, i.e., A = U*S*V^*, where
%                       U - n1*n1*n3
%                       S - n1*n2*n3
%                       V - n2*n2*n3
%                   'econ': produces the "economy size" decomposition. 
%                       Let m = min(n1,n2). Then, A = U*S*V^*, where
%                       U - n1*m*n3
%                       S - m*m*n3
%                       V - n2*m*n3
%                   'skinny': produces the skinny tensor SVD.
%                       Let r be the tensor tubal rank of A. Then, A = U*S*V^*, where
%                       U - n1*r*n3
%                       S - r*r*n3
%                       V - n2*r*n3
%
% Output: U, S, V

[n1,n2,n3] = size(A);
n12 = min(n1,n2);
A = fft(A,[],3);
U = zeros(n1,n12,n3);
V = zeros(n2,n12,n3);
S = zeros(n12,n12,n3);
for i = 1 : n3
    [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(A(:,:,i),'econ');
end
Laplace_S = Laplace(S, 0.01);
U = ifft(U,[],3);
S = ifft(S,[],3);
V = ifft(V,[],3);

function out = Laplace(x, epsilon)
out = 1 - exp(-x/epsilon);
end

end
