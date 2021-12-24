function [L, S, rL, err] = fused_ls(Z, lambda1, lambda2, mu, tol, max_iter)
%--------------------------------------------------------------------------
% fused_ls implements a linearized alternating direction method of multiplier (ADMM)
% for fused low-rank plus sparse (L+S) decomposition of data matrix Z into
% a low-rank component L, and a sparse component S, with model Z = L + S
%
% Inputs: Z - N x M observed multi-subject connectivity data matrix to be decomposed
%             N - number of connectivity edges, M - number of subjects
%         lambda1  - regularization parameter, controls sparsity of S
%         lambda2  - regularization parameter for penalizing the differences
%                    between successive columns of L
%         mu       - augmented lagrangian parameter
%         tol      - reconstruction error tolerance
%         max_iter - maximum number of iterations
%
% Outputs: L - N x M low-rank matrix
%          S - N x M sparse matrix
%          rL - estimated rank of L
%          err - reconstruction error ||Z - L - S||_F / ||Z||_F
%
% Author: Chee-Ming Ting (2021)
%
% Reference:
% Chee-Ming Ting, Jeremy I Skipper, Steven L Small, Hernando Ombao
%             "Separating Stimulus-Induced and Background Components of 
%             Dynamic Functional Connectivity in Naturalistic fMRI"
%             arXiv preprint arXiv:2102.10331
%--------------------------------------------------------------------------
[N,M] = size(Z);
unobserved = isnan(Z);
Z(unobserved) = 0;
    
d = M-1;
e = eye(M);
D = zeros(d,M);
normZ = norm(Z,'fro');

r = 1;
for i=2:M
     D(r,:) = (e(:,i) - e(:,i-1))';
     r = r + 1;
end 
A = kron(D,speye(N));
    
% default arguments
if nargin < 2
    lambda1 = 1 / sqrt(max(N,M));
end
if nargin < 3
    lambda2 = 1 / sqrt(N*d);
end
if nargin < 4
    mu = 1/norm(Z); % 2-norm or maximum singular value of matrix Z
end
 if nargin < 5
    tol = 1e-6;
end
if nargin < 6
    max_iter = 1000;
end
    
v = 1.5*max(abs(eig(D'*D))); % Parameter v>0 controls proximity to L(k-1)

% initial solution
L     = zeros(N,M);
S     = zeros(N,M);
alpha = zeros(N*d,1);
U     = zeros(N,M);
V     = zeros(N*d,1);
AL    = A*L(:); 
    
% Pre-computed constants
tauL = 1/(mu*(v+1));
tauS = lambda1/mu;
taualpha = lambda2/mu;
    
oldL = L; oldS = S;
    
for iter = 1:max_iter 
    % ADMM step: update L, S and alpha
    C = L(:) - A'*(AL - alpha + V)/v;
    Ctilde = reshape(C,N,M);
    Q = (Z - S + v*Ctilde - U)/(v+1);
    L = svt(double(tauL),double(Q));
    AL = A*L(:);
    S = shrink(tauS,Z-L-U);
    alpha = shrink(taualpha, AL + V);

   % update augmented lagrangian multiplier
    U = U + L + S - Z;
    V = V + AL - alpha;

    % Stopping criterion
    del = norm([L S]- [oldL oldS],'fro') / max([norm([L S],'fro'),1]); % relative changes of parameters
    err = norm(Z - L - S, 'fro') / normZ; % reconstruction error
    if (iter == 1) || (mod(iter, 10) == 0)
        fprintf(1, 'iter: %04d\tdel: %f\terr: %f\trank(L): %d\tcard(S): %d\n', ...
                iter, del, err, rank(L), nnz(S(~unobserved)));
    end        
    if (del < tol && iter>10)
        break; 
    end
    oldL = L;  oldS = S;
    clear C Q
end
rL = rank(L);
end

function X = shrink(tau, Y)
    % soft-thresholding (shrinkage) operator
    X = sign(Y) .* max(abs(Y) - tau, 0);
end

function X = svt(tau, Y)
    % singular value thresholding - shrinkage operator for singular values
    [U, S, V] = svd(Y, 'econ');
    X = U*shrink(tau, S)*V';
end
