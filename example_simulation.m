%==========================================================================
%  Simulation: Fused low-rank plus sparse (L+S) decomposition for separating
%  common/correlated structure and individual-specific background
%  components in multi-subject functional connectivity (FC) networks
%
%  Author: Chee-Ming Ting, 4 Nov 2021
%  
%  Reference: Chee-Ming Ting, Jeremy I Skipper, Steven L Small, Hernando Ombao
%             "Separating Stimulus-Induced and Background Components of 
%             Dynamic Functional Connectivity in Naturalistic fMRI"
%             arXiv preprint arXiv:2102.10331
%  
%  External functions of SBM:
%   - generateSbm.m & bernrnd.m - By Kevin S. Xu
%==========================================================================
clc; clear all; close all;

M = 50;  % Number of subjects
N = 20;  % Number of network nodes
r = 5;   % rank
s = 0.1; % sparsity level

%-------------------------------------------------------------------------%
%             Simulate multi-subject connectivity networks                %
%-------------------------------------------------------------------------%
% Generate basis connectivity matrices (undirected networks with community structure)
K = 2;       % number of communities or blocks
n_k = N/K;   % number of nodes per community
alpha = 0.95; pii = 0.8;
P = pii.*eye(K) + (1-pii).*ones(K); % K x K inter-ccmmunity connection probability matrix
P = alpha.*P;
g = zeros(N,1); % N x 1 community membership vector
for k=1:K
    g((k-1)*n_k+1:(k-1)*n_k+n_k) = k; end

B = zeros(N^2,r); % columns of B: a set of r vectorized basis connectivity matrices
a = -1; b = 1;
for ri = 1:r
    SSigma = generateSbm(g,P); % generate networks using stochastic block model
    SSigma = SSigma.*(a + (b-a)*rand(N)); % edge strenghts from U[-1,1]
    SSigma = (SSigma + SSigma')/2; % Ignore self-connections
    B(:,ri) = SSigma(:); % vectorized connectivity matrix
    clear SSigma
end

% Generate mixing weights for each subject (two groups with different weights)
beta = zeros(r,M);
vars = 0.005;
mu1 = 0.5*ones(r,1);
mu2 = zeros(r,1);
beta(:,1:M/2) = mvnrnd(mu1,vars*eye(r),M/2)'; % represents strong response to stimulus
beta(:,M/2+1:M) =  mvnrnd(mu2,vars*eye(r),M/2)'; % represents weak response to stimulus

% Generate similar/correlated connectivity matrices across subjects
Lt = B*beta; % Lt: matrix of rank r, connectivity profiles in columns of L are similar within each group

% Generate subject-specific sparse background connectivity matrices
St = zeros(N^2,M);
n = N*(N-1)/2;
a = -5; b = 5;
for i =1:M
    BSigma_up = (a + (b-a)*rand(n,1));
    S = sprand(n,1,s); % generate sparse support
    BSigma_up = BSigma_up.*full(spones(S));
    BSigma = tril2mat(BSigma_up,N);
    St(:,i) = BSigma(:); % vectorized connectivity matrix
    clear BSigma_up BSigma
end

% Generate observed connectivity matrix as a sum of low-rank and sparse components
% Columns of Zt: vectorized connectivity of invidividual subjects
Zt = Lt + St;

% Extract upper-triangular part
Zt_up = zeros(N*(N-1)/2,M);
Lt_up = zeros(N*(N-1)/2,M);
St_up = zeros(N*(N-1)/2,M);
for i =1:M
    Zfull = reshape(Zt(:,i),N,N);
    Zt_up(:,i) = Zfull(find(~tril(ones(size(Zfull)))));
    Lfull = reshape(Lt(:,i),N,N);
    Lt_up(:,i) = Lfull(find(~tril(ones(size(Lfull)))));
    Sfull = reshape(St(:,i),N,N);
    St_up(:,i) = Sfull(find(~tril(ones(size(Sfull)))));
    clear Zfull Lfull Sfull
end
Zt = Zt_up; Lt = Lt_up; St = St_up;
clear Zt_up Lt_up St_up

%-------------------------------------------------------------------------%
%            Fused L+S decomposition based on ADMM algorithm              %
%-------------------------------------------------------------------------%
% Parameters
lambda1 = 1 / sqrt(max(N*(N-1)/2,M));
mu = 1/norm(Zt);
lambda2 = 0.001;
max_iter = 5000;
tol = 1e-6;

[Lthat, Sthat, rhat, err] = fused_ls(Zt, lambda1, lambda2, mu, tol, max_iter);
RE_L = relerr(Lt,Lthat);        % relative error of recovered low-rank component L
RE_S = relerr(St,Sthat);        % relative error of recovered sparse component S
shat = nnz(Sthat)/numel(Sthat); % estimated level of sparsity, s
RE_s = abs(shat-s)/s;           % relative error of s

% Plot recovered low-rank and sparse components of functional connectivity (upper-triangular part)
figure('Name','Fused L+S decomposition of multi-subject FC','Color',[1 1 1]);
subplot(3,1,1);
imagesc(Zt); 
caxis([-3 3]); colorbar;
set(gca,'XTick',5:5:M,'XTickLabel',5:5:M,'fontsize',10);
ylabel('Connections', 'fontsize',10);
title('Observed Multi-subject FC','fontsize',10);

subplot(3,1,2);
imagesc(Lthat); 
caxis([-1 1]); colorbar;
set(gca,'XTick',5:5:M,'XTickLabel',5:5:M,'fontsize',10);
ylabel('Connections', 'fontsize',10);
title('Low-rank Component','fontsize',10);

subplot(3,1,3);
imagesc(Sthat); 
caxis([-3 3]); colorbar;
set(gca,'XTick',5:5:M,'XTickLabel',5:5:M,'fontsize',10);
xlabel('Subjects','fontsize',10); ylabel('Connections','fontsize',10);
title('Sparse Component','fontsize',10);
