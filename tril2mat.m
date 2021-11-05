function X = tril2mat(x_up, N)
% X = tril2mat(x_up, N) converts lower triangular part x_up of a symmetric matrix back to
% a N-by-N matrix X
a = triu(ones(N),1);
a(a > 0) = x_up;
X = (a + a')./(eye(N)+1);